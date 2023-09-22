//! benchmarking force fields

use std::{collections::HashMap, path::Path};

use openff_toolkit::qcsubmit::results::ResultCollection;

use ligand::{
    forcefield::ForceField,
    molecule::Molecule,
    openmm::{Context, Integrator, Platform},
};
use serde::{Deserialize, Serialize};

#[allow(unused)]
#[derive(Clone, Deserialize, Serialize)]
struct QMConformerRecord {
    qcarchive_id: String,
    molecule_id: usize,
    mapped_smiles: String,
    energy: f64,
    coordinates: Vec<f64>,
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
            inchi_key: value.to_inchikey(),
        }
    }
}

#[derive(Deserialize, Serialize)]
pub(crate) struct MoleculeStore {
    qcarchive_records: Vec<QMConformerRecord>,
    molecule_records: Vec<MoleculeRecord>,
    mm_conformers: Vec<MMConformerRecord>,
}

impl MoleculeStore {
    pub(crate) fn optimize_mm(&mut self, forcefield: &str) {
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

        let minimized_blob = minimize_blob(data, forcefield);
        for result in minimized_blob {
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

    pub(crate) fn get_dde(&self, _forcefield: &str) {
        todo!()
    }

    pub(crate) fn get_rmsd(&self, _forcefield: &str) {
        todo!()
    }

    pub(crate) fn get_tfd(&self, _forcefield: &str) {
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
        rows.into_iter().map(
            move |(qcarchive_id, mapped_smiles, coordinates)| {
                MinimizationInput {
                    inchi_key: inchi_key.clone(),
                    qcarchive_id,
                    force_field: forcefield.to_owned(),
                    mapped_smiles,
                    coordinates,
                }
            },
        )
    });
    let mut outputs = Vec::new();
    let mut failed = 0;
    for input in inputs {
        match input.run_openmm() {
            Ok(v) => outputs.push(v),
            Err(_) => failed += 1,
        }
    }
    eprintln!("{failed} minimizations failed, {} succeeded", outputs.len());
    outputs
}

impl From<ResultCollection> for MoleculeStore {
    fn from(collection: ResultCollection) -> Self {
        let mut molecule_records = Vec::new();
        let mut qcarchive_records = Vec::new();
        for (qcarchive_record, molecule) in collection.to_records() {
            let molecule_record = MoleculeRecord::from(molecule.clone());
            let mapped_smiles = molecule_record.mapped_smiles.clone();
            molecule_records.push(molecule_record);
            // ibstore has some complicated way of retrieving this from the
            // database, but I'll just call the position in the vector its "id"
            let molecule_id = molecule_records.len() - 1;
            // from qcelemental/checkup_data/physconst.py
            const HARTREE2KCALMOL: f64 = 627.5095;
            qcarchive_records.push(QMConformerRecord {
                qcarchive_id: qcarchive_record.id.clone(),
                molecule_id,
                mapped_smiles,
                energy: qcarchive_record.get_final_energy() * HARTREE2KCALMOL,
                coordinates: molecule.get_conformer(0),
            });
        }
        Self {
            qcarchive_records,
            molecule_records,
            mm_conformers: Vec::new(),
        }
    }
}
