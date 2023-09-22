//! benchmarking force fields

use openff_toolkit::{
    qcportal::models::Record, qcsubmit::results::ResultCollection,
};

use ligand::molecule::Molecule;

#[allow(unused)]
struct QMConformerRecord {
    molecule_id: String,
    mapped_smiles: String,
    qc_record: Record,
    coordinates: Vec<f64>,
}

struct MoleculeRecord {
    mapped_smiles: String,
}

impl From<Molecule> for MoleculeRecord {
    fn from(value: Molecule) -> Self {
        Self {
            mapped_smiles: value.to_mapped_smiles(),
        }
    }
}

pub(crate) struct MoleculeStore {}

impl MoleculeStore {
    pub(crate) fn optimize_mm(&self, _forcefield: &str) {
        todo!()
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
}

impl From<ResultCollection> for MoleculeStore {
    fn from(collection: ResultCollection) -> Self {
        let mut molecule_records = Vec::new();
        let mut qcarchive_records = Vec::new();
        for (qcarchive_record, molecule) in collection.to_records() {
            let molecule_record = MoleculeRecord::from(molecule.clone());

            let molecule_id = qcarchive_record.id.clone();
            qcarchive_records.push(QMConformerRecord {
                molecule_id,
                mapped_smiles: molecule_record.mapped_smiles.clone(),
                qc_record: qcarchive_record,
                coordinates: molecule.get_conformer(0),
            });

            molecule_records.push(molecule_record);
        }
        Self {}
    }
}
