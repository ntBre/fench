//! Force field work bENCH

use openff_toolkit::qcsubmit::results::ResultCollection;

use fench::bench::MoleculeStore;

// this is main.py from my benchmarking repo. the goal to emulate
fn main() {
    env_logger::init();

    let forcefield = "openff-2.1.0.offxml";

    let opt = ResultCollection::parse_file("testfiles/filtered-core-opt.json")
        .unwrap();

    // I'm not doing the DB stuff for now. we're doing it live. jokes aside, I
    // don't actually use the DB after a run anyway
    let mut store = MoleculeStore::from(opt);

    log::info!("finished initializing store");

    store.optimize_mm(forcefield);

    store.to_json("testfiles/store.json");

    store.get_dde(forcefield);
    store.get_rmsd(forcefield);
    store.get_tfd(forcefield);
}
