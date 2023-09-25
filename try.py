import logging
import os

from ibstore._store import MoleculeStore
from openff.qcsubmit.results import OptimizationResultCollection

logging.getLogger("openff").setLevel(logging.ERROR)

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
dde = store.get_dde(ff)
rmsd = store.get_rmsd(ff)
tfd = store.get_tfd(ff)

with open("testfiles/dde.txt", "w") as out:
    for d in sorted(dde, key=lambda x: x.qcarchive_id):
        out.write(f"{d.qcarchive_id} {d.difference}\n")

with open("testfiles/rmsd.txt", "w") as out:
    for d in sorted(rmsd, key=lambda x: x.qcarchive_id):
        out.write(f"{d.qcarchive_id} {d.rmsd}\n")

with open("testfiles/tfd.txt", "w") as out:
    for d in sorted(tfd, key=lambda x: x.qcarchive_id):
        out.write(f"{d.qcarchive_id} {d.tfd}\n")
