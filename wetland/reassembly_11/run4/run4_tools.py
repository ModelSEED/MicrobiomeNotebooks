#!/usr/bin/env python3
"""Run 4 (per Andrew, 2026-11-17: MAGs of the 11 re-assembled libraries; previous results = run-3 merged). Derived from run3_tools.py.
Run 4 (per Andrew, 2026-11-12): classify the MAGs extracted from the 7 DAS Tool runs repeated on 2026-11-03
(248 expected) with the same local GTDB-Tk 2.8.0 / GTDB R232 as runs 1-2, plus control MAGs from run 2, then
merge them into the run-2 tables so writeback_r232.py sees one result set for every MAG of narrative 272312.
Local files only (no KBase calls).

  select  all_mags_run4_full.tsv (make_all_list.py) -> all_mags_run4.tsv = MAGs with no run-2 result + N_CONTROLS
          run-2 MAGs (fixed seed, archaea included) to check reproducibility.
  gate    every run-3 MAG has a summary row (or GTDB-Tk logged why not); every control gets exactly its run-2
          lineage. Exit 4 otherwise.
  merge   merged_run4/gtdbtk_run4/ (name = the run-3 out_dir, which writeback_r232 checks) (run-2 summaries + new rows, filtered/failed lists of both runs, run-3 gtdbtk.json and
          trees) + mag_gtdb_lineages_run4_merged.tsv + all_manifest_run4_merged.tsv for writeback_r232.py (WB_* vars).
"""
import csv, glob, json, os, random, shutil, sys

H = "/home/freiburger/gtdbtk_kbase/"
N_CONTROLS, N_ARCHAEA, SEED = 60, 10, 20261112
read = lambda p: list(csv.DictReader(open(H + p), delimiter="\t"))


def write(path, rows, cols):
    with open(H + path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def select():
    full = read("all_mags_run4_full.tsv")
    run2 = {r["name"]: r for r in read("mag_gtdb_lineages_run3_merged.tsv")}
    new = [r for r in full if r["name"] not in run2]
    old = sorted(r["name"] for r in full if r["name"] in run2)
    rng = random.Random(SEED)
    arc = [n for n in old if run2[n]["classification"].startswith("d__Archaea")]
    ctl = set(rng.sample(arc, N_ARCHAEA))
    ctl |= set(rng.sample([n for n in old if n not in ctl], N_CONTROLS - N_ARCHAEA))
    rows = new + [r for r in full if r["name"] in ctl]
    write("all_mags_run4.tsv", rows, list(full[0].keys()))
    with open(H + "run4_controls.txt", "w") as fh:
        fh.write("\n".join(sorted(ctl)) + "\n")
    print(f"select: {len(full)} MAGs in KBase, {len(old)} with a run-2 result, {len(new)} new; "
          f"+ {len(ctl)} controls ({N_ARCHAEA} archaea) -> all_mags_run4.tsv ({len(rows)})")


def summary_rows(out):
    rows = {}
    for f in sorted(set(os.path.realpath(p) for p in glob.glob(f"{H}{out}/**/gtdbtk.*.summary.tsv", recursive=True))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            rows[r["user_genome"]] = (os.path.basename(f), r)
    return rows


def reasons(out):
    rs = {}
    for f in glob.glob(f"{H}{out}/align/*filtered.tsv") + glob.glob(f"{H}{out}/identify/*failed*.tsv"):
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if p and p[0]:
                rs[p[0]] = os.path.basename(f)
    return rs


def gate():
    man = {r["name"] for r in read("all_manifest_run4.tsv")}
    run4 = {r["name"]: r for r in read("mag_gtdb_lineages_run4.tsv")}
    run2 = {r["name"]: r for r in read("mag_gtdb_lineages_run3_merged.tsv")}
    ctl = [l.strip() for l in open(H + "run4_controls.txt") if l.strip()]
    rs = reasons("gtdbtk_run4")
    missing = sorted(man - set(run4) - set(rs))
    skipped = sorted((man - set(run4)) & set(rs))
    diff = [n for n in ctl if n not in run4 or run4[n]["classification"] != run2[n]["classification"]]
    new = [n for n in run4 if n not in run2]
    ok = [n for n in new if run4[n]["classification"] and not run4[n]["classification"].startswith("Unclassified")]
    print(f"run 4: {len(run4)}/{len(man)} MAGs in the summaries; new MAGs {len(new)}, {len(ok)} with a lineage; "
          f"skipped by GTDB-Tk: {skipped}")
    print(f"controls: {len(ctl) - len(diff)}/{len(ctl)} identical to run 2" + (f"; DIFFER: {diff}" if diff else ""))
    if missing or diff:
        print(f"GATE FAILED: missing {missing[:5]}, control differences {diff[:5]}")
        sys.exit(4)
    print("GATE PASSED")


def merge():
    out, src2, src3 = H + "merged_run4/gtdbtk_run4", "merged_run3/gtdbtk_run3", "gtdbtk_run4"
    if os.path.exists(out):
        shutil.rmtree(out)
    for d in ("classify", "align", "identify"):
        os.makedirs(f"{out}/{d}", exist_ok=True)
    run2_names = {r["name"] for r in read("mag_gtdb_lineages_run3_merged.tsv")}
    s2, s3 = summary_rows(src2), summary_rows(src3)
    new = {n: v for n, v in s3.items() if n not in run2_names}           # controls keep their run-2 rows
    for fname in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
        rows = [r for f, r in s2.values() if f == fname] + [r for f, r in new.values() if f == fname]
        cols = next(iter(r.keys() for f, r in list(s2.values()) + list(s3.values()) if f == fname))
        with open(f"{out}/classify/{fname}", "w", newline="") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    for src in (src2, src3):
        for f in glob.glob(f"{H}{src}/align/*filtered.tsv") + glob.glob(f"{H}{src}/identify/*failed*.tsv"):
            sub = os.path.basename(os.path.dirname(f))
            with open(f"{out}/{sub}/{os.path.basename(f)}", "a") as fh:
                fh.writelines(l for l in open(f) if src == src2 or l.split("\t")[0] not in run2_names)
    for t in glob.glob(f"{H}{src3}/classify/*.tree"):
        shutil.copy(t, f"{out}/classify/run4_{os.path.basename(t)}")
    shutil.copy(f"{H}{src3}/gtdbtk.json", f"{out}/gtdbtk.json")
    open(f"{out}/README.txt", "w").write(
        "Merged GTDB-Tk 2.8.0 / GTDB R232 results for narrative 272312: the run-3 merged set (12,859 MAGs,\n"
        "merged_run3/gtdbtk_run3) + the MAGs first classified in run 4 (gtdbtk_run4: MAGs of the 11 libraries re-assembled\n"
        "2026-11-15/16; control MAGs keep their identical earlier rows). Trees run4_* are run 4's placements.\n")
    lin2 = read("mag_gtdb_lineages_run3_merged.tsv")
    lin3 = [r for r in read("mag_gtdb_lineages_run4.tsv") if r["name"] not in run2_names]
    write("mag_gtdb_lineages_run4_merged.tsv", lin2 + lin3, list(lin2[0].keys()))
    m2 = read("all_manifest_run3_merged.tsv")
    m3 = [r for r in read("all_manifest_run4.tsv") if r["name"] not in run2_names]
    write("all_manifest_run4_merged.tsv", m2 + m3, list(m2[0].keys()))
    print(f"merge: {len(lin2)} run-2 + {len(lin3)} run-3 lineages; manifest {len(m2)} + {len(m3)}; "
          f"summaries {len(s2)} + {len(new)} -> merged_run4/gtdbtk_run4")


if __name__ == "__main__":
    {"select": select, "gate": gate, "merge": merge}[sys.argv[1]]()
