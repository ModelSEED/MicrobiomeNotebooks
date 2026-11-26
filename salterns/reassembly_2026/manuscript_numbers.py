#!/usr/bin/env python
"""Snapshot of the 2026 re-assembly numbers for the manuscript (re-runnable while the poplar rollout progresses).

Collects, for every rollout sample that is complete on poplar:
  NEW  = rollout kept MAGs (merge + Tadpole + default-mode SPAdes; site-coverage binning; DAS Tool;
         keep = CheckM2 >= 50 % complete, < 10 % contaminated, GUNC pass)        <sample>/keep/bins.tsv
  2026 = the KBase 2026 pipeline's DAS Tool bins of the same unit, scored on poplar with the same CheckM2 + GUNC
         and filtered by the same rule                                           qc/old/<sample>/
  2020 = the original MAGs (310 dereplicated + 140 non-dereplicated, Zhou et al.) of the same sample, scored with the
         same CheckM2 + GUNC (vs2020/checkm2*, gunc*) and filtered by the same rule
plus skani recovery (ANI >= 95 %, larger aligned fraction >= 50 %) of the 2020 MAGs by NEW kept MAGs of the same
sample (vs2020/skani_completed.tsv; run vs2020/skani_completed.sh on poplar first).
Two "lost" sets of medium/high-quality original MAGs (both from comparison/mags_old_vs_new.tsv + skani_edges.tsv):
  orig2020_lost_in_first_pass / _recovered : the 37 whose species (skani 95 % ANI / 50 % AF cluster) has no first-pass
                                             2026 MAG (the manuscript's "whose species the first pass had not recovered")
  orig2020_lost30 / orig2020_lost30_recovered : the project's "30 lost" (tadpole/RESULTS.md; make_mag_list.py): the
                                             subset with no first-pass MAG even above 80 % ANI
Writes manuscript_numbers_<date>.json and table_s7_partial_<date>.tsv next to this script.
Run: ~/Documents/py_venv/bin/python manuscript_numbers.py
"""
import csv, io, json, os, statistics as st, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.dirname(HERE)
P = "/scratch1/afreiburger/salterns_2026"
SSH = ["ssh", "-o", "BatchMode=yes", "poplar"]

REMOTE = r'''
import csv, json, os, glob, statistics as st
R = "/scratch1/afreiburger/salterns_2026/rollout"
def n50(p):
    L = []; n = 0
    for l in open(p):
        if l[0] == ">":
            if n: L.append(n)
            n = 0
        else: n += len(l.strip())
    if n: L.append(n)
    L.sort(reverse=True); t = sum(L); a = 0
    for x in L:
        a += x
        if a >= t / 2: return x
    return 0
state = json.load(open(f"{R}/STATE.json"))["samples"]
out = {"complete": sorted(k for k, v in state.items() if v["state"] == "complete"), "n_samples": len(state),
       "meta_fallback": [], "new": {}, "old": {}}
for s in out["complete"]:
    am = os.path.join(R, s, "asm", "asm_meta.json")
    if os.path.exists(am) and "meta_fallback" in open(am).read(): out["meta_fallback"].append(s)
    rows = list(csv.DictReader(open(f"{R}/{s}/keep/bins.tsv"), delimiter="\t"))
    out["new"][s] = [{k: r[k] for k in ("bin", "keep", "checkm2_completeness", "checkm2_contamination", "gunc_pass", "n50", "size_bp")} for r in rows]
    c2 = f"{R}/qc/old/{s}/checkm2/quality_report.tsv"; g = f"{R}/qc/old/{s}/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv"
    if os.path.exists(c2) and os.path.exists(g):
        C = {r["Name"]: r for r in csv.DictReader(open(c2), delimiter="\t")}
        G = {r["genome"]: r for r in csv.DictReader(open(g), delimiter="\t")}
        fs = {os.path.basename(f).rsplit(".", 1)[0]: f for f in glob.glob(f"{R}/inputs/old/{s}/*")}
        out["old"][s] = [{"bin": n, "checkm2_completeness": r["Completeness"], "checkm2_contamination": r["Contamination"],
                          "gunc_pass": G.get(n, {}).get("pass.GUNC", ""), "n50": n50(fs[n]) if n in fs else ""}
                         for n, r in C.items()]
print(json.dumps(out))
'''


def keep(c, x, g):
    return float(c) >= 50 and float(x) < 10 and str(g) == "True"


def hq(c, x):
    return float(c) >= 90 and float(x) < 5


def med(v):
    v = [int(float(x)) for x in v if str(x) not in ("", "nan")]
    return int(st.median(v)) if v else None


def n50(p):
    L = []; n = 0
    for l in open(p):
        if l[0] == ">":
            if n: L.append(n)
            n = 0
        else:
            n += len(l.strip())
    if n: L.append(n)
    L.sort(reverse=True); t = sum(L); a = 0
    for x in L:
        a += x
        if a >= t / 2: return x
    return 0


def main():
    remote = json.loads(subprocess.run(SSH + ["python3", "-"], input=REMOTE, capture_output=True, text=True,
                                       check=True).stdout)
    tmp = os.path.join(HERE, "manuscript_numbers_tmp"); os.makedirs(tmp, exist_ok=True)
    for f in ("checkm2/quality_report.tsv", "checkm2_rest/quality_report.tsv", "gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",
              "gunc_rest/GUNC.progenomes_2.1.maxCSS_level.tsv", "skani_completed.tsv"):
        subprocess.run(["rsync", "-a", "-e", "ssh -o BatchMode=yes", f"poplar:{P}/vs2020/{f}", os.path.join(tmp, f.replace("/", "__"))],
                       stdin=subprocess.DEVNULL, capture_output=True)
    # 2020 originals: CheckM2/GUNC on poplar, N50 locally
    C, G = {}, {}
    for f in ("checkm2__quality_report.tsv", "checkm2_rest__quality_report.tsv"):
        p = os.path.join(tmp, f)
        if os.path.exists(p):
            C.update({r["Name"]: r for r in csv.DictReader(open(p), delimiter="\t")})
    for f in ("gunc__GUNC.progenomes_2.1.maxCSS_level.tsv", "gunc_rest__GUNC.progenomes_2.1.maxCSS_level.tsv"):
        p = os.path.join(tmp, f)
        if os.path.exists(p):
            G.update({r["genome"]: r for r in csv.DictReader(open(p), delimiter="\t")})
    idx = list(csv.DictReader(open(os.path.join(S, "MAG_contigs", "index.tsv")), delimiter="\t"))
    import re
    orig = {}
    for r in idx:
        s = re.search(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool", r["assembly_name"]).group(1).replace("R2_restored_DShore", "R2_restored_D_shore")
        g = r["file"].rsplit(".", 1)[0]
        orig.setdefault(s, []).append({"genome": g, "population": r["population"], "assembly_name": r["assembly_name"],
                                       "n50": n50(os.path.join(S, "MAG_contigs", r["file"])), "qc": C.get(g),
                                       # GUNC 1.1.1 names a genome by its file name cut at the first '.fa'
                                       # ('...contigs.fa_assembly.fa' -> '...contigs'); the files on poplar are the local .fna linked as .fa
                                       "gunc": G.get(r["file"][:-4].__add__(".fa").split(".fa")[0] if r["file"].endswith(".fna") else r["file"].split(".fa")[0], {}).get("pass.GUNC", "")})
    # skani recovery of 2020 MAGs by NEW kept MAGs of the same sample
    rec = set()
    sk = os.path.join(tmp, "skani_completed.tsv")
    if os.path.exists(sk):
        for r in csv.DictReader(open(sk), delimiter="\t"):
            a, q = r["Ref_file"], r["Query_file"]
            if float(r["ANI"]) < 95 or max(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"])) < 50:
                continue
            ao, qo = "/originals_all/" in a, "/originals_all/" in q
            if ao == qo:
                continue
            o_, m_ = (a, q) if ao else (q, a)
            rec.add((os.path.basename(o_).rsplit(".", 1)[0], os.path.basename(m_).rsplit(".", 1)[0]))
    # original MAGs whose species the first (KBase) pass did not recover (comparison/README.md: 'lost')
    lost = set()
    cmpf = os.path.join(HERE, "comparison", "mags_old_vs_new.tsv")
    if os.path.exists(cmpf):
        cm = list(csv.DictReader(open(cmpf), delimiter="\t"))
        newc = {r["species_cluster"] for r in cm if r["genome"].startswith("new|")}
        lost = {r["genome"][4:] for r in cm if r["genome"].startswith("old|") and r["species_cluster"] not in newc
                and r["mimag"] != "below MQ"}
    # the project's "30 lost" originals (no first-pass MAG above 80 % ANI), as listed by tadpole/scripts/make_mag_list.py
    gl = os.path.join(HERE, "annotation_runs", "tadpole", "scripts", "genome_check_mags.tsv")
    lost30 = {r["mag"] for r in csv.DictReader(open(gl), delimiter="\t") if r["population"] == "lost30"}
    assert len(lost30) == 30 and lost30 <= lost, (len(lost30), len(lost30 - lost))
    rows, tot = [], {}
    for s in remote["complete"]:
        new = remote["new"][s]; nk = [b for b in new if b["keep"] == "yes"]
        old = remote["old"].get(s, []); ok_ = [b for b in old if keep(b["checkm2_completeness"], b["checkm2_contamination"], b["gunc_pass"])]
        og = orig.get(s, [])
        oq = [o for o in og if o["qc"] is not None]
        ok20 = [o for o in oq if keep(o["qc"]["Completeness"], o["qc"]["Contamination"], o["gunc"])]
        recovered = [o for o in og if any(g == o["genome"] and m.startswith(s + "__") for g, m in rec)]
        row = {"sample": s, "site": s.split("_")[0] if not s.startswith("R2_restored") else "R2_restored",
               "new_bins": len(new), "new_kept": len(nk), "new_hq": sum(hq(b["checkm2_completeness"], b["checkm2_contamination"]) for b in nk),
               "new_med_n50": med(b["n50"] for b in nk),
               "kbase2026_bins": len(old), "kbase2026_kept": len(ok_), "kbase2026_hq": sum(hq(b["checkm2_completeness"], b["checkm2_contamination"]) for b in ok_),
               "kbase2026_med_n50": med(b["n50"] for b in ok_),
               "orig2020_mags": len(og), "orig2020_scored": len(oq), "orig2020_kept": len(ok20) if len(oq) == len(og) else None,
               "orig2020_hq": sum(hq(o["qc"]["Completeness"], o["qc"]["Contamination"]) for o in ok20) if len(oq) == len(og) else None,
               "orig2020_med_n50": med(o["n50"] for o in (ok20 if len(oq) == len(og) else og)),
               "orig2020_recovered": len(recovered) if os.path.exists(sk) else None,
               "orig2020_lost_in_first_pass": sum(1 for o in og if o["assembly_name"] in lost),
               "orig2020_lost_recovered": sum(1 for o in recovered if o["assembly_name"] in lost) if os.path.exists(sk) else None,
               "orig2020_lost30": sum(1 for o in og if o["assembly_name"] in lost30),
               "orig2020_lost30_recovered": sum(1 for o in recovered if o["assembly_name"] in lost30) if os.path.exists(sk) else None}
        rows.append(row)
    keys = ("new_kept", "new_hq", "kbase2026_kept", "kbase2026_hq", "orig2020_mags", "orig2020_kept", "orig2020_hq", "orig2020_recovered", "orig2020_lost_in_first_pass", "orig2020_lost_recovered",
            "orig2020_lost30", "orig2020_lost30_recovered")
    for k in keys:
        vals = [r[k] for r in rows]
        tot[k] = sum(vals) if all(v is not None for v in vals) else None
    date = time.strftime("%Y-%m-%d", time.gmtime())
    res = {"snapshot_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "n_complete": len(rows),
           "n_samples": remote["n_samples"], "meta_fallback": remote["meta_fallback"], "rows": rows, "totals": tot,
           "pooled_new_med_n50": med(b["n50"] for s in remote["complete"] for b in remote["new"][s] if b["keep"] == "yes"),
           "pooled_kbase2026_med_n50": med(b["n50"] for s in remote["complete"] for b in remote["old"].get(s, [])
                                           if keep(b["checkm2_completeness"], b["checkm2_contamination"], b["gunc_pass"])),
           "pooled_orig2020_med_n50": med(o["n50"] for s in remote["complete"] for o in orig.get(s, [])
                                          if o["qc"] and keep(o["qc"]["Completeness"], o["qc"]["Contamination"], o["gunc"]))}
    json.dump(res, open(os.path.join(HERE, f"manuscript_numbers_{date}.json"), "w"), indent=1)
    with open(os.path.join(HERE, f"table_s7_partial_{date}.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(json.dumps({k: res[k] for k in ("snapshot_utc", "n_complete", "meta_fallback", "totals", "pooled_new_med_n50",
                                         "pooled_kbase2026_med_n50", "pooled_orig2020_med_n50")}, indent=1))


if __name__ == "__main__":
    main()
