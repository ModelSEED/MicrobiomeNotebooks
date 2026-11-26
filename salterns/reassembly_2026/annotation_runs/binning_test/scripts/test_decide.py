#!/usr/bin/env python3
"""Synthetic tests of decide.py against an independent exact-arithmetic oracle of the pre-registered rule.

Usage: python -I test_decide.py /path/to/decide.py [workdir]
Builds one fake <root> per scenario (FASTA bins, CheckM2/CheckM1/GUNC/skani tables, originals.tsv, STATE.json),
runs decide.py, and compares every criterion's pass flag and the overall decision with the oracle.
"""
import json
import os
import shutil
import statistics
import subprocess
import sys
from fractions import Fraction as F

DECIDE = sys.argv[1]
WORK = sys.argv[2] if len(sys.argv) > 2 else "/tmp/test_decide_work"
PRIMARY, SECONDARY = "R1_A_D1", "SF2_C_D1"


def G(name, comp, cont, n50=10000, gunc=True, c1=None, checkm2=True):
    """A genome: CheckM2 completeness/contamination as printed strings, contigs giving the wanted N50."""
    return {"name": name, "comp": str(comp), "cont": str(cont), "contigs": [n50] * 3, "gunc": gunc,
            "c1": c1 if c1 is not None else (str(comp), str(cont)), "checkm2": checkm2}


def gunc_name(g):
    return (g + ".fa").split(".fa")[0]          # GUNC 1.1.1: basename.split(file_suffix)[0]


def write_set(root, gset, s, genomes, pops=None):
    d = (os.path.join(root, s, "dastool", "out", f"{s}_DASTool_bins") if gset == "new"
         else os.path.join(root, "inputs", gset, s))
    os.makedirs(d, exist_ok=True)
    for g in genomes:
        with open(os.path.join(d, g["name"] + ".fa"), "w") as fh:
            for i, L in enumerate(g["contigs"]):
                fh.write(f">{g['name']}_c{i}\n{'A' * L}\n")
    if gset == "orig":
        with open(os.path.join(d, "originals.tsv"), "w") as fh:
            fh.write("mag\tpopulation\tmd5\n")
            for g in genomes:
                fh.write(f"{g['name']}\t{pops[g['name']]}\tx\n")
    qc = os.path.join(root, "qc", gset, s)
    for k in ("checkm2", "checkm1", "gunc"):
        os.makedirs(os.path.join(qc, k), exist_ok=True)
    if not genomes:
        for k in ("checkm2", "checkm1", "gunc"):
            open(os.path.join(qc, k, "EMPTY"), "w").close()
        return d
    with open(os.path.join(qc, "checkm2", "quality_report.tsv"), "w") as fh:
        fh.write("Name\tCompleteness\tContamination\tCompleteness_Model_Used\n")
        for g in genomes:
            if g["checkm2"]:
                fh.write(f"{g['name']}\t{g['comp']}\t{g['cont']}\tNeural Network\n")
    with open(os.path.join(qc, "checkm1", "qa.tsv"), "w") as fh:
        fh.write("Bin Id\tMarker lineage\tCompleteness\tContamination\n")
        for g in genomes:
            fh.write(f"{g['name']}\tk__Bacteria (UID203)\t{g['c1'][0]}\t{g['c1'][1]}\n")
    with open(os.path.join(qc, "gunc", "GUNC.progenomes_2.1.maxCSS_level.tsv"), "w") as fh:
        fh.write("genome\tn_genes_called\tclade_separation_score\tcontamination_portion\tpass.GUNC\n")
        for g in genomes:
            if g["gunc"] is not None:
                fh.write(f"{gunc_name(g['name'])}\t100\t0.1\t0.0\t{g['gunc']}\n")
    return d


def build(root, spec, states):
    shutil.rmtree(root, ignore_errors=True)
    os.makedirs(root)
    for s, sp in spec.items():
        dirs = {}
        dirs["new"] = write_set(root, "new", s, sp["new"])
        dirs["old"] = write_set(root, "old", s, sp["old"])
        dirs["orig"] = write_set(root, "orig", s, sp["orig"], sp["pops"])
        sk = os.path.join(root, "qc", "skani", s)
        os.makedirs(sk, exist_ok=True)
        with open(os.path.join(sk, "skani_triangle.tsv"), "w") as fh:
            fh.write("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n")
            for (sa, ga, sb, gb, ani, afr, afq) in sp["skani"]:
                fh.write(f"{os.path.join(dirs[sa], ga + '.fa')}\t{os.path.join(dirs[sb], gb + '.fa')}\t{ani}\t{afr}\t{afq}\tx\ty\n")
    cfg = {"root": root, "sf2_deadline_utc": "2026-11-14T06:18:00Z",
           "samples": {PRIMARY: {"role": "primary"}, SECONDARY: {"role": "secondary"}}}
    json.dump(cfg, open(os.path.join(root, "config.json"), "w"))
    json.dump({"samples": {s: {"state": st, "info": None} for s, st in states.items()}},
              open(os.path.join(root, "STATE.json"), "w"))


# ------------------------------------------------------------------ oracle (exact rational arithmetic)
def oracle_sample(sp):
    derep = {g for g, p in sp["pops"].items() if p == "derep310"}
    out = {}
    for gs in ("new", "old"):
        G_ = [g for g in sp[gs] if g["checkm2"]]
        mq = [g for g in G_ if F(g["comp"]) >= 50 and F(g["cont"]) < 10]
        hq = [g for g in G_ if F(g["comp"]) >= 90 and F(g["cont"]) < 5]
        mqn = {g["name"] for g in mq}
        rec = set()
        for (sa, ga, sb, gb, ani, afr, afq) in sp["skani"]:
            if F(str(ani)) >= 95 and max(F(str(afr)), F(str(afq))) >= 50:
                for (x, gx, y, gy) in ((sa, ga, sb, gb), (sb, gb, sa, ga)):
                    if x == gs and gx in mqn and y == "orig" and gy in derep:
                        rec.add(gy)
        n50s = sorted(F(sorted(g["contigs"], reverse=True)[0]) for g in mq)    # equal contigs -> N50 = length
        out[gs] = {"mq": len(mq), "hq": len(hq),
                   "n50": statistics.median(n50s) if mq else None,
                   "cont": statistics.median([F(g["cont"]) for g in mq]) if mq else None,
                   "gunc": F(sum(1 for g in mq if g["gunc"] is True), len(mq)) if mq else None,
                   "rec": len(rec)}
    N, O = out["new"], out["old"]

    def und(n, o, t):
        if n is None:
            return False
        if o is None:
            return True
        return t(n, o)
    c = {"a": F(N["mq"]) >= F(6, 5) * O["mq"],
         "b": F(N["hq"]) >= F(6, 5) * O["hq"] or N["hq"] >= O["hq"] + 5,
         "c": und(N["n50"], O["n50"], lambda n, o: n >= F(3, 2) * o),
         "d1": und(N["cont"], O["cont"], lambda n, o: n <= o + 1),
         "d2": und(N["gunc"], O["gunc"], lambda n, o: n >= o - F(1, 20)),
         "e": N["rec"] >= O["rec"]}
    c["d"] = c["d1"] and c["d2"]
    c["pass"] = all(c[k] for k in "abcde")
    return c


def oracle_decision(spec, states):
    p, q = states[PRIMARY], states[SECONDARY]
    ok = {s: (oracle_sample(spec[s])["pass"] if states[s] == "complete" else None) for s in spec}
    if q in ("asm_failed", "asm_timeout"):
        if p != "complete":
            return "INCOMPLETE"
        return "R1_ONLY_PASS" if ok[PRIMARY] else "NO_GO"
    if p != "complete" or q != "complete":
        return "INCOMPLETE"
    return "GO" if ok[PRIMARY] and ok[SECONDARY] else "NO_GO"


# ------------------------------------------------------------------ scenario construction
def base_old():
    """OLD: 10 MQ+ (4 HQ), contamination 0.36 x10, N50 10 kb, GUNC 9/10, plus 3 low-quality bins."""
    gs = [G(f"metabat2.bin.{i:03d}", 95 if i < 4 else 70, "0.36", 10000, gunc=(i != 9)) for i in range(10)]
    gs += [G(f"concoct.bin.{i:03d}_sub", 40, "2.0", 3000) for i in range(3)]
    return gs


def base_new():
    """NEW: 14 MQ+ (10 HQ), contamination 0.5, N50 20 kb, GUNC all pass."""
    gs = [G(f"metabat2.bin.{i:03d}", 96 if i < 10 else 60, "0.5", 20000) for i in range(14)]
    gs += [G(f"maxbin2.Bin.{i:03d}", 30, "1.0", 3000) for i in range(2)]
    return gs


def base_orig():
    gs = [G(f"Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.{i}.contigs__", 98, "1.0", 50000) for i in range(4)]
    gs += [G("Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.concoct_out.3.contigs.fa_assembly", 99, "0.8", 50000)]
    pops = {g["name"]: ("derep310" if "contigs__" in g["name"] else "nonderep140") for g in gs}
    return gs, pops


def sample(new=None, old=None, rec_old=2, rec_new=3, extra_skani=()):
    orig, pops = base_orig()
    old = old if old is not None else base_old()
    new = new if new is not None else base_new()
    sk = []
    for i in range(rec_old):
        sk.append(("old", old[i]["name"], "orig", orig[i]["name"], 98.5, 80.0, 75.0))
    for i in range(rec_new):
        sk.append(("orig", orig[i]["name"], "new", new[i]["name"], 99.1, 60.0, 90.0))
    sk += list(extra_skani)
    return {"new": new, "old": old, "orig": orig, "pops": pops, "skani": sk}


def mod(gs, idx, **kw):
    gs = [dict(g) for g in gs]
    for i in idx:
        g = gs[i]
        for k, v in kw.items():
            if k == "n50":
                g["contigs"] = [v] * 3
            else:
                g[k] = str(v) if k in ("comp", "cont") else v
    return gs


def scenarios():
    S = {}
    good = sample()
    S["baseline: both pass -> GO"] = ({PRIMARY: good, SECONDARY: sample()}, {PRIMARY: "complete", SECONDARY: "complete"})
    # a) OLD mq 10: NEW 12 passes, 11 fails (NEW 14 MQ+, drop to 12 / 11 by making bins non-MQ)
    S["a tie: NEW 12 MQ+ vs OLD 10 -> pass"] = ({PRIMARY: sample(new=mod(base_new(), [12, 13], comp=40)), SECONDARY: sample()}, None)
    S["a fail: NEW 11 MQ+ vs OLD 10"] = ({PRIMARY: sample(new=mod(base_new(), [11, 12, 13], comp=40)), SECONDARY: sample()}, None)
    S["a: contamination 10.0 is not MQ"] = ({PRIMARY: sample(new=mod(base_new(), [10, 11, 12], cont="10.00")), SECONDARY: sample()}, None)
    # b) OLD hq 4 -> 1.2x = 4.8, +5 = 9
    S["b: NEW 5 HQ vs OLD 4 -> pass (1.2x)"] = ({PRIMARY: sample(new=mod(base_new(), range(5, 10), comp=89)), SECONDARY: sample()}, None)
    S["b: NEW 4 HQ vs OLD 4 -> fail"] = ({PRIMARY: sample(new=mod(base_new(), range(4, 10), comp=89)), SECONDARY: sample()}, None)
    S["b: HQ contamination 5.0 is not HQ"] = ({PRIMARY: sample(new=mod(base_new(), range(0, 10), cont="5.0")), SECONDARY: sample()}, None)
    big_old = [G(f"metabat2.bin.{i:03d}", 95, "0.36", 10000) for i in range(30)]
    new35 = [G(f"metabat2.bin.{i:03d}", 95, "0.5", 20000) for i in range(36)]
    S["b: NEW 35 HQ vs OLD 30 -> pass (+5)"] = ({PRIMARY: sample(old=big_old, new=mod(new35, [35], comp=60)), SECONDARY: sample()}, None)
    S["b: NEW 34 HQ vs OLD 30 -> fail"] = ({PRIMARY: sample(old=big_old, new=mod(new35, [34, 35], comp=60)), SECONDARY: sample()}, None)
    # c) OLD median N50 10000 -> 15000 passes, 14999 fails; even-count median with .5
    S["c tie: N50 15000 vs 10000 -> pass"] = ({PRIMARY: sample(new=mod(base_new(), range(14), n50=15000)), SECONDARY: sample()}, None)
    S["c fail: N50 14999 vs 10000"] = ({PRIMARY: sample(new=mod(base_new(), range(14), n50=14999)), SECONDARY: sample()}, None)
    old_half = mod(base_old(), range(5), n50=10001)          # OLD median (10000+10001)/2 = 10000.5 -> 15000.75
    S["c: N50 15001 vs 10000.5 -> pass"] = ({PRIMARY: sample(old=old_half, new=mod(base_new(), range(14), n50=15001)), SECONDARY: sample()}, None)
    S["c: N50 15000 vs 10000.5 -> fail"] = ({PRIMARY: sample(old=old_half, new=mod(base_new(), range(14), n50=15000)), SECONDARY: sample()}, None)
    # d1) OLD median contamination 0.36 -> NEW 1.36 passes (float tie), 1.37 fails
    S["d1 tie: cont 1.36 vs 0.36 -> pass"] = ({PRIMARY: sample(new=mod(base_new(), range(14), cont="1.36")), SECONDARY: sample()}, None)
    S["d1 fail: cont 1.37 vs 0.36"] = ({PRIMARY: sample(new=mod(base_new(), range(14), cont="1.37")), SECONDARY: sample()}, None)
    old57 = mod(base_old(), range(10), cont="0.57")
    S["d1 tie: cont 1.57 vs 0.57 -> pass"] = ({PRIMARY: sample(old=old57, new=mod(base_new(), range(14), cont="1.57")), SECONDARY: sample()}, None)
    old_mix = mod(mod(base_old(), range(5), cont="1.35"), range(5, 10), cont="1.37")   # median 1.36 (two-value)
    S["d1 tie: cont 2.36 vs median(1.35,1.37) -> pass"] = ({PRIMARY: sample(old=old_mix, new=mod(base_new(), range(14), cont="2.36")), SECONDARY: sample()}, None)
    # d2) OLD GUNC 1/5 = 0.2 -> NEW 3/20 = 0.15 passes (float tie), 2/20 fails. 5 OLD MQ+ (4 HQ ... )
    old5 = [G(f"metabat2.bin.{i:03d}", 95, "0.36", 10000, gunc=(i == 0)) for i in range(5)]
    new20 = [G(f"metabat2.bin.{i:03d}", 95, "0.5", 20000, gunc=(i < 3)) for i in range(20)]
    S["d2 tie: GUNC 3/20 vs 1/5 -> pass"] = ({PRIMARY: sample(old=old5, new=new20), SECONDARY: sample()}, None)
    S["d2 fail: GUNC 2/20 vs 1/5"] = ({PRIMARY: sample(old=old5, new=mod(new20, [2], gunc=False)), SECONDARY: sample()}, None)
    S["d2: GUNC row missing counts as fail"] = ({PRIMARY: sample(new=mod(base_new(), range(2), gunc=None)), SECONDARY: sample()}, None)
    # e) originals
    S["e tie: NEW 2 vs OLD 2 -> pass"] = ({PRIMARY: sample(rec_new=2), SECONDARY: sample()}, None)
    S["e fail: NEW 1 vs OLD 2"] = ({PRIMARY: sample(rec_new=1), SECONDARY: sample()}, None)
    orig, _ = base_orig()
    S["e: match by non-MQ NEW bin does not count"] = ({PRIMARY: sample(rec_new=1, extra_skani=[
        ("new", "maxbin2.Bin.000", "orig", orig[2]["name"], 99, 90, 90)]), SECONDARY: sample()}, None)
    S["e: ANI 94.99 does not count"] = ({PRIMARY: sample(rec_new=1, extra_skani=[
        ("new", "metabat2.bin.005", "orig", orig[2]["name"], 94.99, 90, 90)]), SECONDARY: sample()}, None)
    S["e: max AF 50 (query side) counts"] = ({PRIMARY: sample(rec_new=1, extra_skani=[
        ("new", "metabat2.bin.005", "orig", orig[2]["name"], 95.0, 20, 50.0)]), SECONDARY: sample()}, None)
    S["e: nonderep140 original does not count"] = ({PRIMARY: sample(rec_new=1, extra_skani=[
        ("new", "metabat2.bin.005", "orig", orig[4]["name"], 99, 90, 90)]), SECONDARY: sample()}, None)
    S["e: two NEW MAGs on one original count once"] = ({PRIMARY: sample(rec_new=1, extra_skani=[
        ("new", "metabat2.bin.005", "orig", orig[0]["name"], 99, 90, 90)]), SECONDARY: sample()}, None)
    # undefined medians / empty sets
    no_mq_new = mod(base_new(), range(14), comp=40)
    S["empty: NEW no MQ+ -> c,d fail"] = ({PRIMARY: sample(new=no_mq_new, rec_new=0), SECONDARY: sample()}, None)
    S["empty: OLD no MQ+ -> c,d pass"] = ({PRIMARY: sample(old=mod(base_old(), range(10), comp=40), rec_old=0), SECONDARY: sample()}, None)
    S["empty: NEW no bins at all"] = ({PRIMARY: sample(new=[], rec_new=0), SECONDARY: sample()}, None)
    S["missing CheckM2 row -> not MQ"] = ({PRIMARY: sample(new=mod(base_new(), [11, 12, 13], checkm2=False)), SECONDARY: sample()}, None)
    # decision branches
    fail = sample(rec_new=1)
    S["R1 pass, SF2 fail -> NO_GO"] = ({PRIMARY: good, SECONDARY: fail}, {PRIMARY: "complete", SECONDARY: "complete"})
    S["R1 fail, SF2 pass -> NO_GO"] = ({PRIMARY: fail, SECONDARY: good}, {PRIMARY: "complete", SECONDARY: "complete"})
    S["SF2 asm_timeout, R1 pass -> R1_ONLY_PASS"] = ({PRIMARY: good, SECONDARY: good}, {PRIMARY: "complete", SECONDARY: "asm_timeout"})
    S["SF2 asm_failed, R1 fail -> NO_GO"] = ({PRIMARY: fail, SECONDARY: good}, {PRIMARY: "complete", SECONDARY: "asm_failed"})
    S["SF2 asm_timeout, R1 incomplete -> INCOMPLETE"] = ({PRIMARY: good, SECONDARY: good}, {PRIMARY: "incomplete", SECONDARY: "asm_timeout"})
    S["SF2 incomplete (not asm) -> INCOMPLETE"] = ({PRIMARY: good, SECONDARY: good}, {PRIMARY: "complete", SECONDARY: "incomplete"})
    S["R1 incomplete -> INCOMPLETE"] = ({PRIMARY: good, SECONDARY: good}, {PRIMARY: "incomplete", SECONDARY: "complete"})
    return S


def main():
    os.makedirs(WORK, exist_ok=True)
    rows, bad = [], 0
    for i, (title, (spec, states)) in enumerate(scenarios().items()):
        states = states or {PRIMARY: "complete", SECONDARY: "complete"}
        root = os.path.join(WORK, f"s{i:02d}")
        build(root, spec, states)
        r = subprocess.run([sys.executable, "-I", DECIDE, "--config", os.path.join(root, "config.json"),
                            "--state", os.path.join(root, "STATE.json"), "--out", root], capture_output=True, text=True)
        dec = json.load(open(os.path.join(root, "DECISION.json")))
        want = oracle_decision(spec, states)
        diffs = []
        if dec["decision"] != want:
            diffs.append(f"decision {dec['decision']} != oracle {want}")
        if dec["decision"] == "ERROR":
            diffs.append(dec.get("error", "") + " " + r.stderr[-500:])
        for s in (PRIMARY, SECONDARY):
            e = dec.get("samples", {}).get(s, {})
            if "criteria" not in e:
                continue
            o = oracle_sample(spec[s])
            got = {k: e["criteria"][k]["pass"] for k in "abce"}
            got["d1"] = e["criteria"]["d"]["d1_contamination"]["pass"]
            got["d2"] = e["criteria"]["d"]["d2_gunc_pass_fraction"]["pass"]
            got["d"] = e["criteria"]["d"]["pass"]
            got["pass"] = e["pass"]
            for k, v in got.items():
                if v != o[k]:
                    diffs.append(f"{s} {k}: decide={v} oracle={o[k]}")
        p_oracle = oracle_sample(spec[PRIMARY])
        flags = "".join(k if p_oracle[k] else "-" for k in ("a", "b", "c", "d1", "d2", "e"))
        bad += bool(diffs)
        rows.append((title, flags, dec["decision"], want, "; ".join(diffs) or "ok"))
    w = max(len(r[0]) for r in rows)
    print(f"{'scenario':{w}s}  R1 oracle  decision      oracle        check")
    for t, f, d, o, c in rows:
        print(f"{t:{w}s}  {f:9s}  {d:12s}  {o:12s}  {c}")
    print(f"\n{len(rows)} scenarios, {bad} with a mismatch")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
