//! Force field work bENCH

use openff_toolkit::qcsubmit::results::ResultCollection;

use crate::bench::MoleculeStore;

mod bench;

// this is main.py from my benchmarking repo. the goal to emulate
fn main() {
    let forcefield = "openff-2.1.0.offxml";

    let opt = ResultCollection::parse_file("testfiles/core-opt.json").unwrap();

    // make sure it's actually loading something
    assert_eq!(400, opt.entries.values().flatten().count());

    // I'm not doing the DB stuff for now. we're doing it live. jokes aside, I
    // don't actually use the DB after a run anyway
    let mut store = MoleculeStore::from(opt);

    store.optimize_mm(forcefield);

    store.get_dde(forcefield);
    store.get_rmsd(forcefield);
    store.get_tfd(forcefield);
}
