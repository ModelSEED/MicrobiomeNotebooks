"""Build the DRAM 1.5.0 databases from the verified raw files fetched by fetch_dbs.sh.

Usage: python setup_dbs.py RAW_DIR OUT_DIR THREADS PHASE      PHASE = nopfam | pfam

Two phases through DRAM's own select_db mechanism, so everything except Pfam is built (and any
failure seen) while the 40.5 GB Pfam-A.full.gz is still downloading:
  nopfam: clear_config=True, select_db = every non-UniRef, non-Pfam database + the 5 sheets
  pfam:   select_db = [pfam, pfam_hmm] added to the same config and description db

Same databases as `DRAM-setup.py prepare_databases --skip_uniref` (KOfam, Pfam, dbCAN V11, MEROPS,
RefSeq viral, VOGDB, the five DRAM sheets), with two deviations forced by dead URLs:

* DRAM downloads the dbCAN family-activity and subfamily-EC tables from bcb.unl.edu, which now
  redirects every path to the dbCAN2 home page (HTML, HTTP 200). prepare_databases() only honours
  arguments ending in `_loc`, so `dbcan_fam_activities` / `dbcan_subfam_ec` cannot be supplied by
  the user (they are always downloaded). The two download functions are replaced below by copies of
  the same files from pro.unl.edu (same names, V11, 08062022), fetched and checked by fetch_dbs.sh.
* vog_annotations is handled the same way (DRAM's URL redirects correctly, but this pins it to the
  md5-checked copy).
* VOGDB release 237's vog.hmm.tar.gz nests the models under hmm/ (hmm/VOG00001.hmm ...), but DRAM's
  process_vogdb() globs only <dir>/VOG*.hmm, which would merge nothing and fail hmmpress after the
  Pfam build. process_vogdb() is replaced by the same steps with a recursive glob.

Everything else goes through prepare_databases() unchanged. Raw inputs are not modified; the dbCAN
HMM file is copied first because process_dbcan() moves its input.
"""
import os
import shutil
import sys
import tarfile
from glob import glob

from mag_annotator import database_processing as dp

raw, out, threads, phase = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
assert phase in ("nopfam", "pfam"), phase
R = lambda f: os.path.join(raw, f)


def _copy_into(src):
    def f(output_dir, logger=dp.LOGGER, verbose=True, **kw):
        dst = os.path.join(output_dir, os.path.basename(src))
        shutil.copyfile(src, dst)
        logger.info(f"(salterns patch) copied verified {src} -> {dst} instead of downloading")
        return dst
    return f


dp.download_dbcan_fam_activities = _copy_into(R("CAZyDB.08062022.fam-activities.txt"))
dp.download_dbcan_subfam_ec = _copy_into(R("CAZyDB.08062022.fam.subfam.ec.txt"))
dp.download_vog_annotations = _copy_into(R("vog.annotations.tsv.gz"))

def process_vogdb(vog_hmm_targz, output_dir=".", logger=dp.LOGGER, version=dp.DEFAULT_VOGDB_VERSION,
                  threads=1, verbose=True):
    hmm_dir = os.path.join(output_dir, "vogdb_hmms")
    os.mkdir(hmm_dir)
    with tarfile.open(vog_hmm_targz) as t:
        t.extractall(hmm_dir, filter="data")
    vog_hmms = os.path.join(output_dir, f"vog_{version}_hmms.txt")
    files = sorted(glob(os.path.join(hmm_dir, "**", "VOG*.hmm"), recursive=True))
    assert len(files) > 1000, f"only {len(files)} VOG hmms found"
    logger.info(f"(salterns patch) merging {len(files)} VOG hmms (recursive glob)")
    dp.merge_files(files, vog_hmms)
    dp.run_process(["hmmpress", "-f", vog_hmms], logger, verbose=verbose)
    dp.LOGGER.info("VOGdb database processed")
    return {"vogdb": vog_hmms}


dp.process_vogdb = process_vogdb

os.makedirs(out, exist_ok=True)
common = dict(output_dir=out, skip_uniref=True, keep_database_files=False, threads=threads, verbose=True)

if phase == "nopfam":
    stage = os.path.join(raw, "for_setup")
    os.makedirs(stage, exist_ok=True)
    dbcan_copy = os.path.join(stage, "dbCAN-HMMdb-V11.txt")
    shutil.copyfile(R("dbCAN-HMMdb-V11.txt"), dbcan_copy)
    for f in ("kofam_profiles.tar.gz", "kofam_ko_list.gz", "pepunit.lib", "viral.1.protein.faa.gz", "vog.hmm.tar.gz"):
        assert os.path.getsize(R(f)) > 0, f
    dp.prepare_databases(
        kofam_hmm_loc=R("kofam_profiles.tar.gz"),
        kofam_ko_list_loc=R("kofam_ko_list.gz"),
        dbcan_loc=dbcan_copy,
        dbcan_version="11",
        dbcan_date="08062022",
        viral_loc=R("viral.1.protein.faa.gz"),
        peptidase_loc=R("pepunit.lib"),
        vogdb_loc=R("vog.hmm.tar.gz"),
        select_db=["kofam_hmm", "kofam_ko_list", "dbcan", "dbcan_fam_activities", "dbcan_subfam_ec", "vogdb",
                   "vog_annotations", "viral", "peptidase", "genome_summary_form", "module_step_form",
                   "function_heatmap_form", "amg_database", "etc_module_database"],
        clear_config=True,
        **common,
    )
else:
    for f in ("Pfam-A.full.gz", "Pfam-A.hmm.dat.gz"):
        assert os.path.getsize(R(f)) > 0, f
    dp.prepare_databases(
        pfam_loc=R("Pfam-A.full.gz"),
        pfam_hmm_loc=R("Pfam-A.hmm.dat.gz"),
        select_db=["pfam", "pfam_hmm"],
        clear_config=False,
        **common,
    )
print("SETUP_FINISHED", phase)
