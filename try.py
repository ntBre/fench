from ibstore._store import MoleculeStore
from openff.qcsubmit.results import OptimizationResultCollection

data = "try"
store = MoleculeStore.from_qcsubmit_collection(
    OptimizationResultCollection.parse_file(
        "testfiles/filtered-core-opt.json",
    ),
    database_name=f"{data}.sqlite",
)

ff = "openff-2.1.0"
store.optimize_mm(ff)
print(len(store.get_dde(ff)))
