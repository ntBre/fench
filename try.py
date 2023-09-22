import os

from ibstore._store import MoleculeStore
from openff.qcsubmit.results import OptimizationResultCollection

db_name = "try.sqlite"
if os.path.exists(db_name):
    os.remove(db_name)

store = MoleculeStore.from_qcsubmit_collection(
    OptimizationResultCollection.parse_file(
        "testfiles/filtered-core-opt.json",
    ),
    database_name=db_name,
)

ff = "openff-2.1.0"
store.optimize_mm(ff)
print(len(store.get_dde(ff)))
