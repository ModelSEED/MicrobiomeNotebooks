"""Stream the 1.5 GB KOfam profiles tarball (read-only) and copy only the needed KO HMMs to verify_bio/kofam_subset."""
import sys, tarfile
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from kos import ALL
OUT = HERE.parent / "kofam_subset"
want = {f"profiles/{k}.hmm" for k in ALL}
got = set()
with tarfile.open("/home/freiburger/Documents/EmilyKin/kofam_new/profiles.tar.gz", "r|gz") as tf:
    for m in tf:
        if m.name in want and m.isfile():
            data = tf.extractfile(m).read()
            (OUT / Path(m.name).name).write_bytes(data)
            got.add(m.name)
            if got == want:
                break
print("wanted", len(want), "extracted", len(got), "missing", sorted(w.split('/')[-1][:-4] for w in want - got))
