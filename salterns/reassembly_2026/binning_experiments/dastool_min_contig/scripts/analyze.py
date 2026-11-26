#!/usr/bin/env python
"""Metrics, recovery (skani), DAS Tool-vs-CheckM calibration and 43-unit extrapolation for the min-contig experiment.

Inputs: runs/<unit>/<variant>/out/ (run_unit.py), results/originals.tsv and results/current_mags.tsv (prep_refs.py),
the production DAS Tool outputs and KBase CheckM tables in ../../work/ (read only) and ../../ledger.json (read only).
Outputs in results/: bins.tsv, variant_metrics.tsv, variant_totals.tsv, recovery_originals.tsv, recovery_current.tsv,
lost_originals.tsv, organism_turnover.tsv, calibration_bins.tsv, calibration.tsv, maxbin2_posthoc.tsv,
extrapolation_43.tsv, costs.tsv, report_tables.md, skani_all.tsv.
Match rule (as in comparison/): skani ANI >= 95 and the larger aligned fraction >= 50.
Usage: ~/Documents/py_venv/bin/python -I analyze.py
"""
import glob
import json
import os
import statistics
import subprocess
import sys
from datetime import datetime

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
RA = os.path.dirname(os.path.dirname(EXP))
RES = os.path.join(EXP, "results")
sys.path.insert(0, HERE)
from units import UNITS, unit_samples  # noqa: E402

SKANI = os.path.expanduser("~/Documents/gtdbtk_venv/bin/skani")
VARS = ["V0", "V1", "V2", "V3", "V4", "V5", "V6"]
VLAB = {"V0": "V0 as run", "V1": "V1 MaxBin2 >=2500", "V2": "V2 all >=3000", "V3": "V3 all >=5000",
        "V4": "V4 no MaxBin2", "V5": "V5 as run + filtered extras", "V6": "V6 as run + >=5000 extras"}
CO = lambda u: u.endswith("_sediment")
BANDS, BLAB = [49.99, 69.99, 89.99, 101], ["50-70", "70-90", "ge90"]   # DAS completeness bands of passing bins


def n50(lens):
    lens = sorted(lens, reverse=True)
    tot, acc = sum(lens), 0
    for L in lens:
        acc += L
        if acc > tot / 2:          # DAS Tool's definition (cumsum > half)
            return L
    return 0


def production_attempt(unit):
    for a in sorted(glob.glob(os.path.join(RA, "work/dastool", unit, "attempt*")), reverse=True):
        if json.load(open(os.path.join(a, "status.json"))).get("state") == "done":
            return a


def seqlen(unit):
    a = production_attempt(unit)
    t = pd.read_csv(os.path.join(a, "out", f"{unit}.seqlength"), sep="\t", header=None, names=["c", "l"],
                    dtype={"c": str})
    t["c"] = t.c.str.split().str[0]
    return t.set_index("c").l


# ---------------------------------------------------------------- bins of every variant
def load_bins():
    rows = []
    for u in UNITS:
        L = seqlen(u)
        for v in VARS:
            od = os.path.join(EXP, "runs", u, v, "out")
            s = pd.read_csv(os.path.join(od, f"{u}_DASTool_summary.tsv"), sep="\t")
            c2b = pd.read_csv(os.path.join(od, f"{u}_DASTool_contig2bin.tsv"), sep="\t", header=None, names=["c", "b"],
                              dtype=str)
            c2b["l"] = c2b.c.map(L)
            lt = c2b[c2b.l < 2500].groupby("b").l.sum() / c2b.groupby("b").l.sum()
            for r in s.itertuples():
                rows.append({"unit": u, "variant": v, "bin": r.bin, "bin_set": r.bin_set, "SCG_set": r.SCG_set,
                             "das_completeness": r.SCG_completeness, "das_redundancy": r.SCG_redundancy,
                             "bin_score": r.bin_score, "size_bp": r.size, "n_contigs": r.contigs, "n50": r.N50,
                             "frac_bp_lt2500": round(float(lt.get(r.bin, 0.0)), 4),
                             "fasta": os.path.join(od, f"{u}_DASTool_bins", f"{r.bin}.fa")})
    b = pd.DataFrame(rows)
    b["pass_mq"] = (b.das_completeness >= 50) & (b.das_redundancy < 10)
    b["pass_hq"] = (b.das_completeness >= 90) & (b.das_redundancy < 5)
    assert all(os.path.exists(f) for f in b.fasta), "missing bin FASTA"
    b["gid"] = "bin|" + b.unit + "|" + b.variant + "|" + b.bin
    return b


# ---------------------------------------------------------------- skani over everything (one triangle)
def run_skani(b, orig, cur):
    out = os.path.join(RES, "skani_all.tsv")
    paths = pd.concat([b.fasta, orig.fasta, cur.fasta]).drop_duplicates()
    lst = os.path.join(RES, "skani_inputs.txt")
    if not os.path.exists(out) or open(lst).read().split() != list(paths):
        paths.to_csv(lst, index=False, header=False)
        subprocess.run(["nice", "-n", "10", SKANI, "triangle", "-l", lst, "-t", "24", "-E", "--min-af", "15", "-o", out],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    e = pd.read_csv(out, sep="\t")
    f2g = dict(zip(b.fasta, b.gid)) | dict(zip(orig.fasta, orig.genome)) | dict(zip(cur.fasta, cur.genome))
    e["a"], e["b"] = e.Ref_file.map(f2g), e.Query_file.map(f2g)
    assert e.a.notna().all() and e.b.notna().all()
    e["af_a"], e["af_b"] = e.Align_fraction_ref, e.Align_fraction_query
    e = e[["a", "b", "ANI", "af_a", "af_b"]]
    both = pd.concat([e, e.rename(columns={"a": "b", "b": "a", "af_a": "af_b", "af_b": "af_a"})], ignore_index=True)
    both["match"] = (both.ANI >= 95) & (both[["af_a", "af_b"]].max(axis=1) >= 50)
    return both           # directed: a -> b, af_a = fraction of a aligned


# ---------------------------------------------------------------- calibration: DAS Tool SCG vs KBase CheckM (43 units)
def calibration():
    rows = []
    for d in sorted(glob.glob(os.path.join(RA, "work/checkm/*/"))):
        u = os.path.basename(d.rstrip("/"))
        a = production_attempt(u)
        if a is None or not os.path.exists(os.path.join(d, f"checkm_{u}.tsv")):
            continue
        s = pd.read_csv(os.path.join(a, "out", f"{u}_DASTool_summary.tsv"), sep="\t")
        c = pd.read_csv(os.path.join(d, f"checkm_{u}.tsv"), sep="\t")
        c["bin"] = c.assembly_name.str.split(".fasta__").str[0]
        m = s.merge(c[["bin", "completeness", "contamination", "n50", "size_bp"]], on="bin", how="outer", indicator=True)
        m["unit"] = u
        rows.append(m)
    t = pd.concat(rows, ignore_index=True)
    return t


def ledger_costs():
    L = json.load(open(os.path.join(RA, "ledger.json")))["tasks"]
    rows = []
    for k, t in L.items():
        for a in t["attempts"]:
            if a.get("status") == "completed" and a.get("finished_utc") and a.get("submitted_utc"):
                h = (datetime.fromisoformat(a["finished_utc"]) - datetime.fromisoformat(a["submitted_utc"])).total_seconds() / 3600
                rows.append({"stage": t["stage"], "unit": t["sample"], "hours_submit_to_finish": h})
    return pd.DataFrame(rows)


def main():
    orig = pd.read_csv(os.path.join(RES, "originals.tsv"), sep="\t", keep_default_na=False)
    cur = pd.read_csv(os.path.join(RES, "current_mags.tsv"), sep="\t", keep_default_na=False)
    b = load_bins()
    e = run_skani(b, orig, cur)
    b.drop(columns=["gid"]).to_csv(os.path.join(RES, "bins.tsv"), sep="\t", index=False)
    passing = set(b[b.pass_mq].gid)
    M = e[e.match]
    # best match per (a, unit-variant) among passing bins
    eb = e[e.b.str.startswith("bin|")].copy()
    eb["unit"], eb["variant"] = eb.b.str.split("|").str[1], eb.b.str.split("|").str[2]
    eb["b_pass"] = eb.b.isin(passing)

    # ---------- recovery (a): same-sample originals; (b): current MAGs of the unit
    rec_o, rec_c = [], []
    for u in UNITS:
        o = orig[orig.unit_sample.isin(unit_samples(u))]
        c = cur[cur["sample"] == u]
        for v in VARS:
            x = eb[(eb.unit == u) & (eb.variant == v)]
            for ref, out, tag in ((o, rec_o, "original"), (c, rec_c, "current")):
                for g in ref.genome:
                    y = x[(x.a == g)]
                    yp = y[y.b_pass & y.match]
                    best = (yp if len(yp) else y[y.b_pass]).sort_values(["ANI", "af_a"], ascending=False).head(1)
                    out.append({"unit": u, "variant": v, "genome": g, "recovered": len(yp) > 0,
                                "n_matching_passing_bins": len(yp),
                                "best_bin": best.b.iloc[0].split("|")[-1] if len(best) else "",
                                "best_ani": best.ANI.iloc[0] if len(best) else np.nan,
                                "best_af_of_ref": best.af_a.iloc[0] if len(best) else np.nan,
                                "best_af_of_bin": best.af_b.iloc[0] if len(best) else np.nan})
    rec_o = pd.DataFrame(rec_o).merge(orig[["genome", "sample", "mimag", "genus", "species_cluster", "n50", "lost58",
                                             "lost30"]], on="genome")
    rec_c = pd.DataFrame(rec_c).merge(cur[["genome", "mimag", "genus", "species_cluster", "n50"]], on="genome")
    bn50 = dict(zip(b.unit + "|" + b.variant + "|" + b.bin, b.n50))
    for r in (rec_o, rec_c):
        r["best_bin_n50"] = [bn50.get(f"{u}|{v}|{x}", np.nan) for u, v, x in zip(r.unit, r.variant, r.best_bin)]
    rec_o.to_csv(os.path.join(RES, "recovery_originals.tsv"), sep="\t", index=False)
    rec_c.to_csv(os.path.join(RES, "recovery_current.tsv"), sep="\t", index=False)

    # ---------- organism turnover vs V0 (passing bins matched between variants)
    tv = []
    for u in UNITS:
        p0 = b[(b.unit == u) & (b.variant == "V0") & b.pass_mq]
        for v in VARS[1:]:
            pk = b[(b.unit == u) & (b.variant == v) & b.pass_mq]
            mk = M[M.a.isin(p0.gid) & M.b.isin(pk.gid)]
            kept = set(mk.a)
            new = set(pk.gid) - set(mk.b)
            for r in p0.itertuples():
                tv.append({"unit": u, "variant": v, "v0_bin": r.bin, "v0_bin_set": r.bin_set, "v0_compl": r.das_completeness,
                           "v0_red": r.das_redundancy, "v0_n50": r.n50, "v0_size": r.size_bp, "status": "kept" if r.gid in kept else "dropped",
                           "match_bin": ";".join(sorted(x.split("|")[-1] for x in mk[mk.a == r.gid].b)),
                           "match_n50": mk[mk.a == r.gid].b.map(lambda g: bn50[g.split("|", 1)[1]]).max() if r.gid in kept else np.nan,
                           "match_compl": b.set_index("gid").loc[list(mk[mk.a == r.gid].b), "das_completeness"].max() if r.gid in kept else np.nan})
            c0 = pd.read_csv(os.path.join(EXP, "runs", u, "V0", "out", f"{u}_DASTool_contig2bin.tsv"), sep="\t",
                             header=None, names=["c", "b"], dtype=str)
            c0p = set(c0[c0.b.isin(p0.bin)].c)
            ck = pd.read_csv(os.path.join(EXP, "runs", u, v, "out", f"{u}_DASTool_contig2bin.tsv"), sep="\t", header=None,
                             names=["c", "b"], dtype=str)
            curu = set(cur[cur["sample"] == u].genome)
            for g in new:
                r = b.set_index("gid").loc[g]
                nbr = set(M[M.a == g].b)
                tv.append({"unit": u, "variant": v, "v0_bin": "", "v0_bin_set": "", "status": "new", "match_bin": r.bin,
                           "match_n50": r.n50, "match_compl": r.das_completeness, "new_red": r.das_redundancy,
                           "new_size": r.size_bp, "new_bin_set": r.bin_set,
                           "new_matches_any_V0_bin": bool(nbr & set(b[(b.unit == u) & (b.variant == "V0")].gid)),
                           "new_matches_current_MAG": bool(nbr & curu),
                           "new_matches_original": any(x.startswith("old|") for x in nbr),
                           "new_bp_in_V0_passing_bins": round(float(ck[ck.b == r.bin].c.isin(c0p).mean()), 4)})
    tv = pd.DataFrame(tv)
    tv.to_csv(os.path.join(RES, "organism_turnover.tsv"), sep="\t", index=False)

    # ---------- lost originals: any bin (passing or not) of any variant; and presence in the assembly
    lo = orig[orig.lost58].copy()
    lrows = []
    for g in lo.genome:
        r = lo.set_index("genome").loc[g]
        row = {"genome": g, "sample": r["sample"], "genus": r.genus, "mimag": r.mimag, "completeness": r.completeness,
               "contamination": r.contamination, "lost30": r.lost30}
        own = [u for u in UNITS if r.unit_sample in unit_samples(u)]
        row["experiment_unit"] = ";".join(own)
        for u in own:            # several units when a core's original is also in scope of its site co-assembly
            sk = os.path.join(EXP, "runs", u, "skani_originals_vs_assembly.tsv")
            s = pd.read_csv(sk, sep="\t")
            s = s[s.Query_file == r.fasta]
            row[f"in_assembly_{u}"] = (f"ANI {s.ANI.max():.1f} AF {s.Align_fraction_query.max():.0f}" if len(s) else "none")
            row["in_assembly_af"] = max(row.get("in_assembly_af", 0.0), s.Align_fraction_query.max() if len(s) else 0.0)
            row["in_assembly_ani"] = np.nanmax([row.get("in_assembly_ani", np.nan), s.ANI.max() if len(s) else np.nan])
        for v in VARS:           # any unit of the experiment (a co-assembly may recover another core's genome)
            x = eb[(eb.a == g) & (eb.variant == v)]
            xm = x[x.match]
            row[f"{v}_any_bin_match"] = len(xm) > 0
            row[f"{v}_passing_bin_match"] = bool((xm.b_pass).any())
            bb = x.sort_values(["ANI", "af_a"], ascending=False).head(1)
            row[f"{v}_best"] = (f"{bb.b.iloc[0].split('|', 1)[1]} ANI {bb.ANI.iloc[0]:.1f} AF {bb.af_a.iloc[0]:.0f}"
                                if len(bb) else "")
        lrows.append(row)
    lost = pd.DataFrame(lrows)
    lost.to_csv(os.path.join(RES, "lost_originals.tsv"), sep="\t", index=False)

    # ---------- per unit x variant metrics
    rows = []
    for u in UNITS:
        for v in VARS:
            x = b[(b.unit == u) & (b.variant == v)]
            p = x[x.pass_mq]
            ro = rec_o[(rec_o.unit == u) & (rec_o.variant == v)]
            ro0 = rec_o[(rec_o.unit == u) & (rec_o.variant == "V0")]
            rc = rec_c[(rec_c.unit == u) & (rec_c.variant == v)]
            rc0 = rec_c[(rec_c.unit == u) & (rec_c.variant == "V0")]
            g_o = set(ro[ro.recovered].genome) - set(ro0[ro0.recovered].genome)
            l_o = set(ro0[ro0.recovered].genome) - set(ro[ro.recovered].genome)
            g_c = set(rc[rc.recovered].genome) - set(rc0[rc0.recovered].genome)
            l_c = set(rc0[rc0.recovered].genome) - set(rc[rc.recovered].genome)
            t = tv[(tv.unit == u) & (tv.variant == v)] if v != "V0" else None
            rows.append({
                "unit": u, "variant": v, "das_bins": len(x), "pass_50_10": len(p), "pass_90_5": int(x.pass_hq.sum()),
                "pass_from_maxbin2": int(p.bin_set.str.startswith("maxbin2").sum()),
                "pass_from_filtered_set": int(p.bin_set.str.contains("_ge").sum()),
                "median_n50_pass": p.n50.median(), "median_contigs_pass": p.n_contigs.median(),
                "median_size_mbp_pass": round(p.size_bp.median() / 1e6, 3), "total_mbp_pass": round(p.size_bp.sum() / 1e6, 1),
                "median_das_compl_pass": p.das_completeness.median(), "median_frac_bp_lt2500_pass": p.frac_bp_lt2500.median(),
                "orig_n": len(ro), "orig_mq_n": int((ro.mimag != "below MQ").sum()),
                "orig_recovered": int(ro.recovered.sum()), "orig_mq_recovered": int((ro.recovered & (ro.mimag != "below MQ")).sum()),
                "orig_species_recovered": ro[ro.recovered].species_cluster.nunique(),
                "orig_gained_vs_V0": len(g_o), "orig_lost_vs_V0": len(l_o),
                "orig_median_af_recovered": ro[ro.recovered].best_af_of_ref.median(),
                "orig_median_n50_ratio_recovered": (ro[ro.recovered].best_bin_n50 / ro[ro.recovered].n50).median(),
                "cur_n": len(rc), "cur_recovered": int(rc.recovered.sum()), "cur_gained_vs_V0": len(g_c),
                "cur_lost_vs_V0": len(l_c),
                "v0_pass_kept": int((t.status == "kept").sum()) if t is not None else len(p),
                "v0_pass_dropped": int((t.status == "dropped").sum()) if t is not None else 0,
                "new_pass_vs_V0": int((t.status == "new").sum()) if t is not None else 0,
                "lost58_any_bin": int(sum(((eb.a == g) & (eb.unit == u) & (eb.variant == v) & eb.match).any()
                                          for g in lost.genome[[u in s.split(";") for s in lost.experiment_unit]])),
                "lost58_in_scope": int(sum(u in s.split(";") for s in lost.experiment_unit)),
            })
    vm = pd.DataFrame(rows)
    vm.to_csv(os.path.join(RES, "variant_metrics.tsv"), sep="\t", index=False)

    # ---------- calibration (all 43 units, production V0 bins)
    cal = calibration()
    cal.to_csv(os.path.join(RES, "calibration_bins.tsv"), sep="\t", index=False)
    j = cal[cal._merge == "both"].copy()
    j["das_mq"] = (j.SCG_completeness >= 50) & (j.SCG_redundancy < 10)
    j["ck_mq"] = (j.completeness >= 50) & (j.contamination < 10)
    j["das_hq"] = (j.SCG_completeness >= 90) & (j.SCG_redundancy < 5)
    j["ck_hq"] = (j.completeness >= 90) & (j.contamination < 5)
    from scipy.stats import pearsonr, spearmanr
    calrows = {
        "units": j.unit.nunique(), "bins_joined": len(j), "bins_das_only": int((cal._merge == "left_only").sum()),
        "bins_checkm_only": int((cal._merge == "right_only").sum()),
        "pearson_completeness": pearsonr(j.SCG_completeness, j.completeness)[0],
        "spearman_completeness": spearmanr(j.SCG_completeness, j.completeness)[0],
        "pearson_redundancy_vs_contamination": pearsonr(j.SCG_redundancy, j.contamination)[0],
        "spearman_redundancy_vs_contamination": spearmanr(j.SCG_redundancy, j.contamination)[0],
        "median_checkm_minus_das_completeness": (j.completeness - j.SCG_completeness).median(),
        "median_checkm_minus_das_redundancy": (j.contamination - j.SCG_redundancy).median(),
        "slope_checkm_on_das_completeness": np.polyfit(j.SCG_completeness, j.completeness, 1)[0],
        "intercept_checkm_on_das_completeness": np.polyfit(j.SCG_completeness, j.completeness, 1)[1],
        "das_mq_and_ck_mq": int((j.das_mq & j.ck_mq).sum()), "das_mq_not_ck_mq": int((j.das_mq & ~j.ck_mq).sum()),
        "ck_mq_not_das_mq": int((~j.das_mq & j.ck_mq).sum()), "neither_mq": int((~j.das_mq & ~j.ck_mq).sum()),
        "das_hq_and_ck_hq": int((j.das_hq & j.ck_hq).sum()), "das_hq_not_ck_hq": int((j.das_hq & ~j.ck_hq).sum()),
        "ck_hq_not_das_hq": int((~j.das_hq & j.ck_hq).sum()),
        "p_ck_mq_given_das_mq": (j.das_mq & j.ck_mq).sum() / j.das_mq.sum(),
        "p_ck_mq_given_not_das_mq": (~j.das_mq & j.ck_mq).sum() / (~j.das_mq).sum(),
        "p_ck_hq_given_das_hq": (j.das_hq & j.ck_hq).sum() / j.das_hq.sum(),
        "p_ck_hq_given_not_das_hq": (~j.das_hq & j.ck_hq).sum() / (~j.das_hq).sum(),
        **{f"p_ck_mq_given_das_mq_compl_{k}": v for k, v in
           j[j.das_mq].groupby(pd.cut(j[j.das_mq].SCG_completeness, BANDS, labels=BLAB), observed=True).ck_mq.mean().items()},
        **{f"n_das_mq_compl_{k}": int(v) for k, v in
           j[j.das_mq].groupby(pd.cut(j[j.das_mq].SCG_completeness, BANDS, labels=BLAB), observed=True).size().items()},
        "n50_identical_frac": float((j.N50 == j.n50).mean()), "size_identical_frac": float((j.size == j.size_bp).mean()),
        "total_checkm_mq_43_units": int(j.ck_mq.sum()), "total_checkm_hq_43_units": int(j.ck_hq.sum()),
    }
    pd.Series(calrows).to_csv(os.path.join(RES, "calibration.tsv"), sep="\t", header=False)

    # ---------- expected KBase CheckM passes per variant (calibrated) and totals over the 8 units
    p1, p0 = calrows["p_ck_mq_given_das_mq"], calrows["p_ck_mq_given_not_das_mq"]
    h1, h0 = calrows["p_ck_hq_given_das_hq"], calrows["p_ck_hq_given_not_das_hq"]
    # MQ: a DAS-passing bin's CheckM pass rate depends on its DAS completeness (50-70: ~55 %, 70-90: ~84 %, >=90: ~94 %),
    # so the expected count is summed per bin by completeness band; a DAS-failing bin gets the overall failing rate.
    band_rate = {k: calrows[f"p_ck_mq_given_das_mq_compl_{k}"] for k in BLAB if f"p_ck_mq_given_das_mq_compl_{k}" in calrows}
    b["exp_ck_mq"] = np.where(b.pass_mq, pd.cut(b.das_completeness, BANDS, labels=BLAB).map(band_rate).astype(float), p0)
    vm["exp_checkm_mq_two_rate"] = vm.pass_50_10 * p1 + (vm.das_bins - vm.pass_50_10) * p0
    vm["exp_checkm_mq"] = [b[(b.unit == u) & (b.variant == v)].exp_ck_mq.sum() for u, v in zip(vm.unit, vm.variant)]
    vm["exp_checkm_hq"] = vm.pass_90_5 * h1 + (vm.das_bins - vm.pass_90_5) * h0
    ck = j.groupby("unit").agg(ck_mq=("ck_mq", "sum"), ck_hq=("ck_hq", "sum"))
    vm["actual_checkm_mq_V0"] = vm.unit.map(ck.ck_mq)
    vm["actual_checkm_hq_V0"] = vm.unit.map(ck.ck_hq)
    vm.to_csv(os.path.join(RES, "variant_metrics.tsv"), sep="\t", index=False)
    tot = []
    for v in VARS:
        x = vm[vm.variant == v]
        p = b[(b.variant == v) & b.pass_mq]
        ro = rec_o[rec_o.variant == v]
        tot.append({"variant": VLAB[v], "das_bins": x.das_bins.sum(), "pass_50_10": x.pass_50_10.sum(),
                    "pass_90_5": x.pass_90_5.sum(), "pass_from_maxbin2": x.pass_from_maxbin2.sum(),
                    "pass_from_filtered_set": x.pass_from_filtered_set.sum(),
                    "exp_checkm_mq": round(x.exp_checkm_mq.sum(), 1),
                    "exp_checkm_mq_two_rate": round(x.exp_checkm_mq_two_rate.sum(), 1),
                    "exp_checkm_hq": round(x.exp_checkm_hq.sum(), 1),
                    "median_n50_pass": p.n50.median(), "median_contigs_pass": p.n_contigs.median(),
                    "median_size_mbp_pass": round(p.size_bp.median() / 1e6, 2), "total_mbp_pass": round(p.size_bp.sum() / 1e6, 1),
                    "median_das_compl_pass": p.das_completeness.median(),
                    "orig_recovered": x.orig_recovered.sum(), "orig_mq_recovered": x.orig_mq_recovered.sum(),
                    "orig_n": x.orig_n.sum(), "orig_mq_n": x.orig_mq_n.sum(),
                    "orig_gained_vs_V0": x.orig_gained_vs_V0.sum(), "orig_lost_vs_V0": x.orig_lost_vs_V0.sum(),
                    "orig_median_af_recovered": ro[ro.recovered].best_af_of_ref.median(),
                    "orig_median_n50_ratio_recovered": (ro[ro.recovered].best_bin_n50 / ro[ro.recovered].n50).median(),
                    "cur_recovered": x.cur_recovered.sum(), "cur_n": x.cur_n.sum(),
                    "v0_pass_kept": x.v0_pass_kept.sum(), "v0_pass_dropped": x.v0_pass_dropped.sum(),
                    "new_pass_vs_V0": x.new_pass_vs_V0.sum(), "lost58_any_bin": x.lost58_any_bin.sum(),
                    "lost58_in_scope": x.lost58_in_scope.sum()})
    tot = pd.DataFrame(tot)
    tot.to_csv(os.path.join(RES, "variant_totals.tsv"), sep="\t", index=False)

    # ---------- MaxBin2 post-hoc filter: how lossy is V1 as a stand-in for MaxBin2 at 2500?
    mp = []
    for u in UNITS:
        a = production_attempt(u)
        L = seqlen(u)
        mb = pd.read_csv(os.path.join(a, "c2b", "maxbin2.contigs2bin.tsv"), sep="\t", header=None, names=["c", "b"], dtype=str)
        mb["c"] = mb.c.str.split().str[0]
        mb["l"] = mb.c.map(L)
        scg = pd.concat([pd.read_csv(os.path.join(a, "out", f"{u}_proteins.faa.{k}.scg"), sep="\t", header=None,
                                     names=["p", "g"], dtype=str).assign(set=k) for k in ("bacteria", "archaea")])
        scg["c"] = scg.p.str.replace(r"_[0-9]+$", "", regex=True)
        scg["l"] = scg.c.map(L)
        bac = scg[scg.set == "bacteria"]
        # MaxBin2 sets its number of bins from marker-gene copy numbers on contigs >= min_contig_length (median count
        # of its marker set); proxy: the median copy number of DAS Tool's 51 bacterial SCGs at each length floor
        med = {m: statistics.median(bac[bac.l >= m].g.value_counts().reindex(bac.g.unique(), fill_value=0)) for m in (1000, 2500)}
        per_bin = mb.groupby("b").apply(lambda t: pd.Series({"bp": t.l.sum(), "bp_ge2500": t.l[t.l >= 2500].sum()}),
                                        include_groups=False)
        scg_b = scg[scg.c.isin(set(mb.c))]
        mp.append({"unit": u, "maxbin2_bins": mb.b.nunique(), "maxbin2_contigs": len(mb), "maxbin2_mbp": round(mb.l.sum() / 1e6, 1),
                   "frac_contigs_lt2500": round((mb.l < 2500).mean(), 3),
                   "frac_bp_lt2500": round(mb.l[mb.l < 2500].sum() / mb.l.sum(), 3),
                   "frac_binned_scg_hits_on_lt2500": round((scg_b.l < 2500).mean(), 3),
                   "median_per_bin_frac_bp_removed": round(float((1 - per_bin.bp_ge2500 / per_bin.bp).median()), 3),
                   "bins_left_empty_by_filter": int((per_bin.bp_ge2500 == 0).sum()),
                   "bins_after_filter": int((per_bin.bp_ge2500 > 0).sum()),
                   "median_scg_copies_ge1000": med[1000], "median_scg_copies_ge2500": med[2500],
                   "assembly_frac_bp_lt2500": round(L[L < 2500].sum() / L.sum(), 3)})
    mp = pd.DataFrame(mp)
    mp.to_csv(os.path.join(RES, "maxbin2_posthoc.tsv"), sep="\t", index=False)

    # ---------- 43-unit extrapolation
    allv0 = cal[cal._merge != "right_only"].copy()
    allv0["das_mq"] = (allv0.SCG_completeness >= 50) & (allv0.SCG_redundancy < 10)
    allv0["ck_mq"] = (allv0.completeness >= 50) & (allv0.contamination < 10)
    allv0["ck_hq"] = (allv0.completeness >= 90) & (allv0.contamination < 5)
    allv0["co"] = allv0.unit.map(CO)
    base = allv0.groupby("co").agg(units=("unit", "nunique"), das_bins=("bin", "size"), ck_mq=("ck_mq", "sum"),
                                   ck_hq=("ck_hq", "sum"), maxbin2_selected=("bin_set", lambda s: (s == "maxbin2").sum()))
    base_mb_mq = allv0[allv0.ck_mq].groupby("co").bin_set.apply(lambda s: (s == "maxbin2").sum())
    ex = []
    for v in VARS:
        for co in (False, True):
            x = vm[(vm.variant == v) & (vm.unit.map(CO) == co)]
            x0 = vm[(vm.variant == "V0") & (vm.unit.map(CO) == co)]
            r_mq = x.exp_checkm_mq.sum() / x0.exp_checkm_mq.sum()
            r_hq = x.exp_checkm_hq.sum() / x0.exp_checkm_hq.sum()
            pk = b[(b.variant == v) & b.pass_mq & (b.unit.map(CO) == co)]
            pk0 = b[(b.variant == "V0") & b.pass_mq & (b.unit.map(CO) == co)]
            ex.append({"variant": VLAB[v], "units": "co-assemblies (4)" if co else "per-sample (39)",
                       "ratio_exp_checkm_mq_vs_V0_8units": round(r_mq, 3), "ratio_exp_checkm_hq_vs_V0_8units": round(r_hq, 3),
                       "actual_checkm_mq_V0_43": int(base.loc[co, "ck_mq"]), "actual_checkm_hq_V0_43": int(base.loc[co, "ck_hq"]),
                       "projected_checkm_mq_43": round(base.loc[co, "ck_mq"] * r_mq), "projected_checkm_hq_43": round(base.loc[co, "ck_hq"] * r_hq),
                       "median_n50_ratio_vs_V0_8units": round(pk.n50.median() / pk0.n50.median(), 2),
                       "basis_units": ",".join(sorted(x.unit))})
    ex = pd.DataFrame(ex)
    tot43 = ex.groupby("variant", sort=False)[["actual_checkm_mq_V0_43", "actual_checkm_hq_V0_43", "projected_checkm_mq_43",
                                               "projected_checkm_hq_43"]].sum().reset_index()
    tot43["units"] = "all 43"
    ex = pd.concat([ex, tot43], ignore_index=True)
    ex.to_csv(os.path.join(RES, "extrapolation_43.tsv"), sep="\t", index=False)

    # replacement rate of MaxBin2-derived V0 passing bins in V4 (8 units), applied to all 43 units
    t4 = tv[(tv.variant == "V4") & (tv.status != "new")]
    mbt = t4[t4.v0_bin_set == "maxbin2"]   # V0 bin sets are the three binners as run
    repl = {"v0_pass_from_maxbin2_8units": len(mbt), "replaced_in_V4": int((mbt.status == "kept").sum()),
            "replacement_rate": (mbt.status == "kept").mean(),
            "v0_ckmq_from_maxbin2_43units": int(base_mb_mq.sum()),
            "das_selected_from_maxbin2_43units": int(base.maxbin2_selected.sum()), "das_bins_43units": int(base.das_bins.sum())}
    repl["expected_net_loss_43units_V4"] = round(repl["v0_ckmq_from_maxbin2_43units"] * (1 - repl["replacement_rate"]))

    # ---------- costs
    lc = ledger_costs()
    cost = lc.groupby("stage").hours_submit_to_finish.describe()[["count", "50%", "max"]].rename(
        columns={"50%": "median_h", "max": "max_h"})
    cost.to_csv(os.path.join(RES, "costs.tsv"), sep="\t")
    chk = {u: json.load(open(os.path.join(EXP, "runs", u, "run_checks.json"))) for u in UNITS}
    rt = pd.DataFrame({u: {"fetch_s": c["fetch_seconds"], "selection_only_s (5 variants concurrently)": c["V0"]["seconds"],
                           "fresh_prodigal_diamond_V0_s (-t 24)": c.get("V0fresh", {}).get("seconds")} for u, c in chk.items()}).T

    # ---------- report tables
    md = []
    w = md.append
    w("# Report tables (generated by scripts/analyze.py)\n")
    w("## Reproduction of the production DAS Tool run (V0)\n")
    rep = pd.DataFrame({u: {"stored bins": c["V0"]["vs_stored"]["stored_bins"],
                            "V0 (stored SCG, --resume) identical": c["V0"]["vs_stored"]["identical"],
                            "V0fresh (prodigal+diamond -t 24) identical": c.get("V0fresh", {}).get("vs_stored", {}).get("identical"),
                            "fresh SCG lines differing (bac/arc)": "{}/{}".format(
                                *[c["V0fresh"]["scg_vs_stored"][k]["only_stored"] + c["V0fresh"]["scg_vs_stored"][k]["only_fresh"]
                                  for k in ("bacteria", "archaea")]) if "V0fresh" in c else ""} for u, c in chk.items()}).T
    w(rep.to_markdown()); w("")
    w("## Totals over the 8 units (DAS Tool SCG estimates; passing = DAS completeness >= 50 and redundancy < 10)\n")
    w(tot.set_index("variant").T.to_markdown()); w("")
    w("## Per unit\n")
    w(vm[["unit", "variant", "das_bins", "pass_50_10", "pass_90_5", "pass_from_maxbin2", "median_n50_pass",
          "median_contigs_pass", "median_size_mbp_pass", "total_mbp_pass", "orig_n", "orig_recovered", "orig_gained_vs_V0",
          "orig_lost_vs_V0", "cur_n", "cur_recovered", "v0_pass_dropped", "new_pass_vs_V0", "lost58_in_scope",
          "lost58_any_bin"]].to_markdown(index=False)); w("")
    w("## Calibration: DAS Tool SCG estimates vs KBase CheckM (production V0 bins, all units)\n")
    w(pd.Series(calrows).to_frame("value").to_markdown()); w("")
    w("## MaxBin2 post-hoc filter (V1) vs a MaxBin2 re-run at 2500 bp\n")
    w(mp.to_markdown(index=False)); w("")
    w("## Extrapolation to all 43 units (KBase CheckM >= 50/<10 and >= 90/<5)\n")
    w(ex.drop(columns=["basis_units"]).to_markdown(index=False)); w("")
    w("MaxBin2 replacement (V4): " + json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in repl.items()}) + "\n")
    w("## Lost original MAGs in scope (species cluster without any re-assembly MAG)\n")
    cols = ["genome", "sample", "genus", "mimag", "completeness", "contamination", "lost30", "in_assembly_ani", "in_assembly_af"] + \
           [f"{v}_any_bin_match" for v in VARS] + ["V0_best", "V3_best"]
    cols += [c for c in lost.columns if c.startswith("in_assembly_") and c not in cols]
    w(lost[lost.experiment_unit != ""][cols].assign(genome=lambda d: d.genome.str.replace("old|Salt_Pond_MetaG", "", regex=False).str[:45])
      .to_markdown(index=False)); w("")
    w("## Run times and KBase stage costs\n")
    w(rt.to_markdown()); w("")
    w(cost.round(2).to_markdown()); w("")
    open(os.path.join(RES, "report_tables.md"), "w").write("\n".join(md))
    json.dump({"calibration": calrows, "replacement": repl}, open(os.path.join(RES, "summary.json"), "w"), indent=1,
              default=float)
    print("\n".join(md))


if __name__ == "__main__":
    main()
