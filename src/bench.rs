//! benchmarking force fields

use std::{collections::HashMap, path::Path};

use log::{debug, info};
use openff_toolkit::{
    qcportal::models::Record, qcsubmit::results::ResultCollection,
};

use ligand::{
    forcefield::ForceField,
    molecule::Molecule,
    openmm::{Context, Integrator, Platform},
};
use serde::{Deserialize, Serialize};

use crate::utils::make_rel;

#[allow(unused)]
#[derive(Clone, Deserialize, Serialize)]
struct QMConformerRecord {
    qcarchive_id: String,
    molecule_id: usize,
    mapped_smiles: String,
    energy: f64,
    coordinates: Vec<f64>,
}

impl QMConformerRecord {
    fn from_qcarchive_record(
        molecule_id: usize,
        mapped_smiles: String,
        qc_record: Record,
        coordinates: Vec<f64>,
    ) -> Self {
        const HARTREE2KCALMOL: f64 = 627.5095;
        Self {
            molecule_id,
            mapped_smiles,
            energy: qc_record.get_final_energy() * HARTREE2KCALMOL,
            qcarchive_id: qc_record.id,
            coordinates,
        }
    }
}

#[allow(unused)]
#[derive(Default, Deserialize, Serialize)]
struct MMConformerRecord {
    molecule_id: usize,
    qcarchive_id: String,
    force_field: String,
    mapped_smiles: String,
    coordinates: Vec<f64>,
    energy: f64,
}

#[derive(Deserialize, Serialize)]
struct MoleculeRecord {
    mapped_smiles: String,
    inchi_key: String,
}

impl From<Molecule> for MoleculeRecord {
    fn from(value: Molecule) -> Self {
        Self {
            mapped_smiles: value.to_mapped_smiles(),
            inchi_key: value.to_inchi(),
        }
    }
}

#[derive(Default, Deserialize, Serialize)]
pub struct MoleculeStore {
    qcarchive_records: Vec<QMConformerRecord>,
    molecule_records: Vec<MoleculeRecord>,
    mm_conformers: Vec<MMConformerRecord>,
}

impl MoleculeStore {
    pub fn from_json(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        Ok(serde_json::from_str(&std::fs::read_to_string(path)?)?)
    }

    pub fn optimize_mm(&mut self, forcefield: &str) {
        let inchi_keys = self.get_inchi_keys();
        let mut data = HashMap::new();

        for inchi_key in inchi_keys {
            let molecule_id =
                self.get_molecule_id_by_inchi_key(&inchi_key).unwrap();
            let qm_conformers =
                self.qcarchive_records.clone().into_iter().filter_map(
                    |record| {
                        if record.molecule_id == molecule_id {
                            Some((
                                record.qcarchive_id,
                                record.mapped_smiles,
                                record.coordinates,
                            ))
                        } else {
                            None
                        }
                    },
                );
            for qm_conformer in qm_conformers {
                data.entry(inchi_key.clone())
                    .or_insert(Vec::new())
                    .push(qm_conformer);
            }
        }

        if data.is_empty() {
            // nothing to do, but this is impossible without the DB caching
            // ibstore is doing
            return;
        }

        info!("minimizing blob");

        let minimized_blob = minimize_blob(data, forcefield);

        info!("finished minimizing blob");

        for result in minimized_blob.into_iter() {
            let inchi_key = result.inchi_key;
            let molecule_id =
                self.get_molecule_id_by_inchi_key(&inchi_key).unwrap();
            self.mm_conformers.push(MMConformerRecord {
                molecule_id,
                qcarchive_id: result.qcarchive_id,
                force_field: result.force_field,
                mapped_smiles: result.mapped_smiles,
                coordinates: result.coordinates,
                energy: result.energy,
            });
        }
    }

    pub fn get_dde(&self, _forcefield: &str) -> Vec<(String, f64)> {
        let mut ret = Vec::new();
        let inchi_keys = self.get_inchi_keys();

        debug!("found {} inchi keys", inchi_keys.len());

        for inchi_key in inchi_keys {
            let molecule_id =
                self.get_molecule_id_by_inchi_key(&inchi_key).unwrap();

            let qcarchive_ids =
                self.get_qcarchive_ids_by_molecule_id(molecule_id);
            if qcarchive_ids.len() == 1 {
                // only one conformer, so you can't compute ΔΔE
                debug!("skipping {inchi_key}");
                continue;
            }
            let mut qm_energies =
                self.get_qm_energies_by_molecule_id(molecule_id);
            make_rel(&mut qm_energies);

            let mut mm_energies =
                self.get_mm_energies_by_molecule_id(molecule_id);
            if mm_energies.len() != qm_energies.len() {
                continue;
            }
            make_rel(&mut mm_energies);

            // already check that mm == qm above
            assert_eq!(qcarchive_ids.len(), mm_energies.len());

            let ids = qcarchive_ids.len();

            for i in 0..ids {
                let mm = mm_energies[i];
                let qm = qm_energies[i];
                debug!("{} {} {}", qcarchive_ids[i], mm, qm);
                ret.push((qcarchive_ids[i].clone(), mm - qm));
            }
        }
        ret
    }

    pub fn get_rmsd(&self, _forcefield: &str) -> Vec<(String, f64)> {
        let inchi_keys = self.get_inchi_keys();
        debug!("found {} inchi keys", inchi_keys.len());

        let mut ret = Vec::new();
        for inchi_key in inchi_keys {
            let molecule = Molecule::from_inchi(&inchi_key).unwrap();
            let molecule_id =
                self.get_molecule_id_by_inchi_key(&inchi_key).unwrap();
            let qcarchive_ids =
                self.get_qcarchive_ids_by_molecule_id(molecule_id);
            let qm_conformers =
                self.get_qm_conformers_by_molecule_id(molecule_id);
            let mm_conformers =
                self.get_mm_conformers_by_molecule_id(molecule_id);

            for (i, (mm, qm)) in mm_conformers
                .into_iter()
                .zip(qm_conformers.into_iter())
                .enumerate()
            {
                let rmsd = molecule.get_rmsd(qm, mm);
                ret.push((qcarchive_ids[i].clone(), rmsd));
            }
        }
        ret
    }

    pub fn get_tfd(&self, _forcefield: &str) {
        todo!()
    }

    fn get_inchi_keys(&self) -> Vec<String> {
        self.molecule_records
            .iter()
            .map(|mol| mol.inchi_key.clone())
            .collect()
    }

    fn get_molecule_id_by_inchi_key(&self, inchi_key: &str) -> Option<usize> {
        self.molecule_records
            .iter()
            .position(|m| m.inchi_key == inchi_key)
    }

    pub fn to_json(&self, path: impl AsRef<Path>) {
        std::fs::write(path, serde_json::to_string_pretty(self).unwrap())
            .unwrap()
    }

    fn get_qcarchive_ids_by_molecule_id(&self, id: usize) -> Vec<String> {
        self.qcarchive_records
            .iter()
            .filter_map(|rec| {
                if rec.molecule_id == id {
                    Some(rec.qcarchive_id.clone())
                } else {
                    None
                }
            })
            .collect()
    }

    fn get_qm_energies_by_molecule_id(&self, _molecule_id: usize) -> Vec<f64> {
        self.qcarchive_records
            .iter()
            .flat_map(|rec| {
                if rec.molecule_id == _molecule_id {
                    Some(rec.energy)
                } else {
                    None
                }
            })
            .collect()
    }

    fn get_mm_energies_by_molecule_id(&self, _molecule_id: usize) -> Vec<f64> {
        self.mm_conformers
            .iter()
            .flat_map(|rec| {
                if rec.molecule_id == _molecule_id {
                    Some(rec.energy)
                } else {
                    None
                }
            })
            .collect()
    }

    fn store(&mut self, record: MoleculeRecord) {
        if self.smiles_already_exists(&record.mapped_smiles) {
            return;
        }
        self.molecule_records.push(record);
    }

    fn get_molecule_id_by_smiles(&self, mapped_smiles: String) -> usize {
        self.molecule_records
            .iter()
            .position(|rec| rec.mapped_smiles == mapped_smiles)
            .unwrap()
    }

    fn store_qcarchive(&mut self, record: QMConformerRecord) {
        if !self.qm_conformer_already_exists(&record.qcarchive_id) {
            self.qcarchive_records.push(record);
        }
    }

    /// true if we already have a QCArchive record with this ID
    fn qm_conformer_already_exists(&self, qcarchive_id: &str) -> bool {
        self.qcarchive_records
            .iter()
            .any(|rec| rec.qcarchive_id == qcarchive_id)
    }

    /// true if we already have a molecule with this SMILES string
    fn smiles_already_exists(&self, mapped_smiles: &str) -> bool {
        self.molecule_records
            .iter()
            .any(|rec| rec.mapped_smiles == mapped_smiles)
    }

    fn get_qm_conformers_by_molecule_id(
        &self,
        molecule_id: usize,
    ) -> Vec<Vec<f64>> {
        self.qcarchive_records
            .iter()
            .filter_map(|rec| {
                if rec.molecule_id == molecule_id {
                    Some(rec.coordinates.clone())
                } else {
                    None
                }
            })
            .collect()
    }

    fn get_mm_conformers_by_molecule_id(
        &self,
        molecule_id: usize,
    ) -> Vec<Vec<f64>> {
        self.mm_conformers
            .iter()
            .filter_map(|rec| {
                if rec.molecule_id == molecule_id {
                    Some(rec.coordinates.clone())
                } else {
                    None
                }
            })
            .collect()
    }
}

struct MinimizationInput {
    inchi_key: String,
    qcarchive_id: String,
    force_field: String,
    mapped_smiles: String,
    coordinates: Vec<f64>,
}

impl MinimizationInput {
    fn run_openmm(self) -> anyhow::Result<MinimizationResult> {
        let MinimizationInput {
            inchi_key,
            qcarchive_id,
            force_field,
            mapped_smiles,
            coordinates,
        } = self;
        let molecule = Molecule::from_mapped_smiles(&mapped_smiles)?;
        let forcefield = ForceField::new(&force_field)?;
        let system = forcefield
            .create_interchange(molecule.to_topology())?
            .to_openmm();
        let mut context =
            Context::new(system, Integrator::Verlet(0.1), Platform::Reference);
        context.set_positions(coordinates);
        context.minimize(5e-9, 1500);

        Ok(MinimizationResult {
            inchi_key,
            qcarchive_id,
            force_field,
            mapped_smiles,
            coordinates: context.get_coordinates(),
            energy: context.get_energy(),
        })
    }
}

struct MinimizationResult {
    inchi_key: String,
    qcarchive_id: String,
    force_field: String,
    mapped_smiles: String,
    coordinates: Vec<f64>,
    energy: f64,
}

fn minimize_blob(
    data: HashMap<String, Vec<(String, String, Vec<f64>)>>,
    forcefield: &str,
) -> Vec<MinimizationResult> {
    // define it in this weird way so that we can parallelize later
    let inputs = data.into_iter().flat_map(|(inchi_key, rows)| {
        rows.into_iter().filter_map(
            move |(qcarchive_id, mapped_smiles, coordinates)| {
                let is_isomorphic =
                    Molecule::from_inchi(&inchi_key).unwrap().is_isomorphic(
                        Molecule::from_mapped_smiles(&mapped_smiles).unwrap(),
                    );
                if is_isomorphic {
                    Some(MinimizationInput {
                        inchi_key: inchi_key.clone(),
                        qcarchive_id,
                        force_field: forcefield.to_owned(),
                        mapped_smiles,
                        coordinates,
                    })
                } else {
                    None
                }
            },
        )
    });
    let mut outputs = Vec::new();
    let mut failed = 0;
    for (i, input) in inputs.enumerate() {
        match input.run_openmm() {
            Ok(v) => outputs.push(v),
            Err(_) => failed += 1,
        }
        info!("finished input {i}");
    }
    eprintln!("{failed} minimizations failed, {} succeeded", outputs.len());
    outputs
}

impl From<ResultCollection> for MoleculeStore {
    fn from(collection: ResultCollection) -> Self {
        // let mut molecule_records = Vec::new();
        // let mut qcarchive_records = Vec::new();
        let mut store = Self::default();
        // the parent_ids are supposed to come from collection itself. inside of
        // to_records, and optimization_records, which it calls,
        // collection.into::<CollectionGetResponse>().ids() contains the
        // parent_id values returned by get_qcarchive_ids_by_molecule_id
        for (qcarchive_record, molecule) in collection.to_records() {
            let molecule_record = MoleculeRecord::from(molecule.clone());
            let mapped_smiles = molecule_record.mapped_smiles.clone();
            let smiles = molecule_record.mapped_smiles.clone();
            store.store(molecule_record);
            let molecule_id = store.get_molecule_id_by_smiles(smiles);
            // from qcelemental/checkup_data/physconst.py
            store.store_qcarchive(QMConformerRecord::from_qcarchive_record(
                molecule_id,
                mapped_smiles,
                qcarchive_record,
                molecule.get_conformer(0),
            ));
        }
        store
    }
}
