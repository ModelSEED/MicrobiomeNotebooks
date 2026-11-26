#!/usr/bin/env python3
"""Pre-registered decision for the salterns binning test (stdlib only). Writes DECISION.json, REPORT.md and mags.tsv.

PRE-REGISTERED RULE (fixed before any result was seen; implemented literally):
Per sample S in {R1_A_D1, SF2_C_D1}, compare NEW (Arm C assembly binned like 2026) vs OLD (the 2026 pipeline's DAS Tool
output for the same unit: all its DAS Tool bins), both scored here with the same tools: CheckM2 (primary), CheckM1
lineage_wf --reduced_tree (secondary), GUNC.
 a) MQ+ MAGs (CheckM2 completeness >= 50, contamination < 10): NEW >= 1.2 x OLD
 b) HQ MAGs (>= 90, < 5): NEW >= 1.2 x OLD or NEW >= OLD + 5
 c) median N50 of MQ+ MAGs: NEW >= 1.5 x OLD
 d) median CheckM2 contamination of MQ+ MAGs: NEW <= OLD + 1.0 point; and GUNC pass fraction of MQ+ MAGs:
    NEW >= OLD - 0.05
 e) original species recovered (skani ANI >= 95 and max AF >= 50 between an MQ+ MAG and an original MAG of sample S;
    originals = derep310 MAGs whose sample is S): NEW >= OLD
 S passes if a, b, c, d and e all pass. CheckM1-based a/b are reported (informational).
 GO for the full salterns rollout only if BOTH samples pass. If only R1_A_D1 passes: NO_GO (the SF2-type decision
 needs Andrew). If SF2_C_D1's Arm C assembly fails or never finishes within 24 h: decide on R1_A_D1 alone and mark
 the decision R1_ONLY_PASS (if R1_A_D1 passes) rather than GO.

Conventions fixed here (not part of the rule text, stated in DECISION.json):
 - A median or fraction over an empty MQ+ set is undefined. Undefined OLD with defined NEW -> the criterion passes;
   undefined NEW -> it fails.
 - GUNC pass = pass.GUNC == True in GUNC.<db>.maxCSS_level.tsv; a genome missing from GUNC's output counts as not
   passing.
 - A genome missing from the CheckM2 report is not MQ+.
 - A sample whose processing failed (a required step failed after retries) is INCOMPLETE; the overall decision is then
   INCOMPLETE (needs Andrew), except that a failed/timed-out SF2 assembly gives R1_ONLY_PASS / NO_GO as the rule says.
 - (review, 2026-11-13) Every threshold is evaluated in exact rational arithmetic on the values as printed (CheckM2
   completeness/contamination as decimal strings, GUNC pass fractions as ratios, N50 medians), so a value exactly at a
   threshold gets the inclusive inequality the rule states (floats made e.g. 1.36 <= 0.36 + 1.0 and 0.15 >= 0.2 - 0.05
   false). Reported numbers are floats of the exact values.
 - (review, 2026-11-13) GUNC 1.1.1 names a genome basename.split(file_suffix)[0], i.e. it cuts the file name at its
   first '.fa'; GUNC rows are joined to genomes on that key (it changes only originals named '*.fa_assembly'). Two
   genomes of one set with the same key would get no GUNC value (not passing) and are listed in qc_status.

Usage: python -I decide.py --config scripts/config.json --state STATE.json --out <root>
"""
import argparse
import csv
import datetime
import json
import os
import statistics
import sys
import traceback
from fractions import Fraction

RULE = __doc__.split("Conventions fixed here")[0].split("PRE-REGISTERED RULE")[1].strip()
CONVENTIONS = __doc__.split("Conventions fixed here")[1].split("Usage:")[0].strip()


def ts():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fasta_stats(path):
    ls, n = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if n is not None:
                    ls.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n is not None:
        ls.append(n)
    srt, tot, acc, n50 = sorted(ls, reverse=True), sum(ls), 0, 0
    for x in srt:
        acc += x
        if 2 * acc >= tot:
            n50 = x
            break
    return {"size_bp": tot, "n_contigs": len(ls), "n50": n50}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def exact(s):
    """Exact rational value of a number as printed in a QC table ('1.36' -> 34/25)."""
    s = str(s).strip()
    try:
        return Fraction(s)
    except ValueError:
        return Fraction(float(s))          # raises for nan/inf -> decision ERROR, never a silent pass


def gunc_key(genome):
    """GUNC 1.1.1 (run_from_fnas): genome name = os.path.basename(fna).split(file_suffix)[0], file_suffix '.fa'."""
    return (genome + ".fa").split(".fa")[0]


def gdir(cfg, gset, s):
    if gset == "new":
        return os.path.join(cfg["root"], s, "dastool", "out", f"{s}_DASTool_bins")
    return os.path.join(cfg["root"], "inputs", gset, s)


def load_set(cfg, gset, s):
    """{genome: {...}} with stats and QC for one genome set of one sample."""
    d = gdir(cfg, gset, s)
    out = {}
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            if fn.endswith(".fa"):
                g = fn[:-3]
                out[g] = {"genome": g, "set": gset, "sample": s, "path": os.path.realpath(os.path.join(d, fn)),
                          **fasta_stats(os.path.join(d, fn))}
    qc = os.path.join(cfg["root"], "qc", gset, s)
    status = {}
    # CheckM2
    f = os.path.join(qc, "checkm2", "quality_report.tsv")
    status["checkm2"] = "ok" if os.path.exists(f) else ("empty" if os.path.exists(os.path.join(qc, "checkm2", "EMPTY")) else "missing")
    if os.path.exists(f):
        for r in read_tsv(f):
            if r["Name"] in out:
                out[r["Name"]].update(checkm2_completeness=float(r["Completeness"]),
                                      checkm2_contamination=float(r["Contamination"]),
                                      checkm2_model=r.get("Completeness_Model_Used", ""),
                                      x_c2_comp=exact(r["Completeness"]), x_c2_cont=exact(r["Contamination"]))
    # CheckM1
    f = os.path.join(qc, "checkm1", "qa.tsv")
    status["checkm1"] = "ok" if os.path.exists(f) else ("empty" if os.path.exists(os.path.join(qc, "checkm1", "EMPTY")) else "missing")
    if os.path.exists(f):
        for r in read_tsv(f):
            g = r.get("Bin Id")
            if g in out:
                out[g].update(checkm1_completeness=float(r["Completeness"]), checkm1_contamination=float(r["Contamination"]),
                              checkm1_marker_lineage=r.get("Marker lineage", ""))
    # GUNC
    gd = os.path.join(qc, "gunc")
    fs = [x for x in os.listdir(gd) if x.startswith("GUNC.") and x.endswith(".maxCSS_level.tsv")] if os.path.isdir(gd) else []
    status["gunc"] = "ok" if fs else ("empty" if os.path.exists(os.path.join(gd, "EMPTY")) else "missing")
    if fs:
        keymap = {}
        for g in out:
            keymap.setdefault(gunc_key(g), []).append(g)
        status["gunc_name_collisions"] = sorted(k for k, v in keymap.items() if len(v) > 1)
        for r in read_tsv(os.path.join(gd, fs[0])):
            hit = keymap.get(r.get("genome"), [])
            if len(hit) == 1:
                out[hit[0]].update(gunc_pass=str(r.get("pass.GUNC", "")).strip().lower() == "true",
                                   gunc_css=r.get("clade_separation_score", ""),
                                   gunc_contamination_portion=r.get("contamination_portion", ""))
    for g in out.values():
        c, k = g.get("x_c2_comp"), g.get("x_c2_cont")
        g["mq"] = c is not None and c >= 50 and k < 10
        g["hq"] = c is not None and c >= 90 and k < 5
        c1, k1 = g.get("checkm1_completeness"), g.get("checkm1_contamination")
        g["mq_checkm1"] = c1 is not None and c1 >= 50 and k1 < 10
        g["hq_checkm1"] = c1 is not None and c1 >= 90 and k1 < 5
    return out, status


def originals_meta(cfg, s):
    f = os.path.join(cfg["root"], "inputs", "orig", s, "originals.tsv")
    return {r["mag"]: r["population"] for r in read_tsv(f)} if os.path.exists(f) else {}


def skani_matches(cfg, s, sets):
    """{(set, genome): set(original derep310 names matched at ANI>=95 and max AF>=50)}"""
    f = os.path.join(cfg["root"], "qc", "skani", s, "skani_triangle.tsv")
    path2 = {}
    for gs, genomes in sets.items():
        for g in genomes.values():
            path2[g["path"]] = (gs, g["genome"])
    res, n_edges = {}, 0
    if not os.path.exists(f):
        return res, "missing", 0
    for r in read_tsv(f):
        a = path2.get(os.path.realpath(r["Ref_file"]))
        b = path2.get(os.path.realpath(r["Query_file"]))
        if not a or not b:
            continue
        ani, af = float(r["ANI"]), max(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"]))
        n_edges += 1
        if ani < 95 or af < 50:
            continue
        for x, y in ((a, b), (b, a)):
            if y[0] == "orig" and x[0] in ("new", "old"):
                res.setdefault(x, set()).add(y[1])
    return res, "ok", n_edges


def med(xs):
    return statistics.median(xs) if xs else None


def tofloat(x):
    return None if x is None else float(x)


def crit_ratio(new, old, factor):
    return new >= factor * old


def evaluate(cfg, s):
    sets, status = {}, {}
    for gs in ("new", "old", "orig"):
        sets[gs], status[gs] = load_set(cfg, gs, s)
    orig_pop = originals_meta(cfg, s)
    derep = sorted(g for g, p in orig_pop.items() if p == "derep310")
    matches, skani_status, n_edges = skani_matches(cfg, s, sets)
    m, X = {}, {}
    for gs in ("new", "old"):
        G = list(sets[gs].values())
        mq = [g for g in G if g["mq"]]
        rec = set()
        for g in mq:
            rec |= matches.get((gs, g["genome"]), set()) & set(derep)
        # exact values used by the criteria (Fractions); m[] holds their float copies for the report
        X[gs] = {"n50": med([Fraction(g["n50"]) for g in mq]),
                 "cont": med([g["x_c2_cont"] for g in mq]),
                 "gunc": Fraction(sum(bool(g.get("gunc_pass")) for g in mq), len(mq)) if mq else None}
        m[gs] = {
            "n_bins": len(G),
            "mq": len(mq),
            "hq": sum(g["hq"] for g in G),
            "median_n50_mq": tofloat(X[gs]["n50"]),
            "median_contamination_mq": tofloat(X[gs]["cont"]),
            "gunc_pass_fraction_mq": tofloat(X[gs]["gunc"]),
            "gunc_missing_mq": sum("gunc_pass" not in g for g in mq),
            "originals_recovered": len(rec),
            "originals_recovered_names": sorted(rec),
            "mq_checkm1": sum(g["mq_checkm1"] for g in G),
            "hq_checkm1": sum(g["hq_checkm1"] for g in G),
            "mbp_in_mq": round(sum(g["size_bp"] for g in mq) / 1e6, 2),
            "median_completeness_mq": med([g["checkm2_completeness"] for g in mq]),
            "median_size_mb_mq": med([g["size_bp"] / 1e6 for g in mq]),
            "median_contigs_mq": med([g["n_contigs"] for g in mq]),
        }
    N, O = m["new"], m["old"]

    def undefined_rule(newv, oldv, test):
        if newv is None:
            return False, "NEW undefined (no MQ+ MAG) -> fail"
        if oldv is None:
            return True, "OLD undefined (no MQ+ MAG) -> pass"
        return test(newv, oldv), ""

    XN, XO = X["new"], X["old"]
    crit = {}
    # a, b on integers: NEW >= 1.2 x OLD  <=>  5 NEW >= 6 OLD (exact)
    crit["a"] = {"metric": "MQ+ MAGs (CheckM2 completeness >= 50, contamination < 10)", "new": N["mq"], "old": O["mq"],
                 "threshold": f"NEW >= 1.2 x OLD = {1.2 * O['mq']:.2f}", "pass": 5 * N["mq"] >= 6 * O["mq"]}
    crit["b"] = {"metric": "HQ MAGs (CheckM2 >= 90, < 5)", "new": N["hq"], "old": O["hq"],
                 "threshold": f"NEW >= 1.2 x OLD = {1.2 * O['hq']:.2f} or NEW >= OLD + 5 = {O['hq'] + 5}",
                 "pass": 5 * N["hq"] >= 6 * O["hq"] or N["hq"] >= O["hq"] + 5}
    ok, note = undefined_rule(XN["n50"], XO["n50"], lambda n, o: n >= Fraction(3, 2) * o)
    crit["c"] = {"metric": "median N50 of MQ+ MAGs (bp)", "new": N["median_n50_mq"], "old": O["median_n50_mq"],
                 "threshold": f"NEW >= 1.5 x OLD = {float(Fraction(3, 2) * XO['n50']):.1f}" if XO["n50"] is not None else "NEW >= 1.5 x OLD",
                 "pass": ok, "note": note}
    ok1, note1 = undefined_rule(XN["cont"], XO["cont"], lambda n, o: n <= o + 1)
    ok2, note2 = undefined_rule(XN["gunc"], XO["gunc"], lambda n, o: n >= o - Fraction(1, 20))
    crit["d"] = {"metric": "median CheckM2 contamination of MQ+ MAGs (%) AND GUNC pass fraction of MQ+ MAGs",
                 "d1_contamination": {"new": N["median_contamination_mq"], "old": O["median_contamination_mq"],
                                      "threshold": (f"NEW <= OLD + 1.0 = {O['median_contamination_mq'] + 1.0:.2f}"
                                                    if O["median_contamination_mq"] is not None else "NEW <= OLD + 1.0"),
                                      "pass": ok1, "note": note1},
                 "d2_gunc_pass_fraction": {"new": N["gunc_pass_fraction_mq"], "old": O["gunc_pass_fraction_mq"],
                                           "threshold": (f"NEW >= OLD - 0.05 = {O['gunc_pass_fraction_mq'] - 0.05:.3f}"
                                                         if O["gunc_pass_fraction_mq"] is not None else "NEW >= OLD - 0.05"),
                                           "pass": ok2, "note": note2},
                 "pass": ok1 and ok2}
    crit["e"] = {"metric": f"originals (derep310 of {s}, n={len(derep)}) recovered by an MQ+ MAG at skani ANI >= 95 and max AF >= 50",
                 "new": N["originals_recovered"], "old": O["originals_recovered"], "threshold": "NEW >= OLD",
                 "pass": N["originals_recovered"] >= O["originals_recovered"]}
    passed = all(crit[k]["pass"] for k in "abcde")
    info = {
        "checkm1_a": {"new": N["mq_checkm1"], "old": O["mq_checkm1"], "threshold": f"NEW >= 1.2 x OLD = {1.2 * O['mq_checkm1']:.2f}",
                      "pass": 5 * N["mq_checkm1"] >= 6 * O["mq_checkm1"]},
        "checkm1_b": {"new": N["hq_checkm1"], "old": O["hq_checkm1"],
                      "threshold": f"NEW >= 1.2 x OLD = {1.2 * O['hq_checkm1']:.2f} or NEW >= OLD + 5",
                      "pass": 5 * N["hq_checkm1"] >= 6 * O["hq_checkm1"] or N["hq_checkm1"] >= O["hq_checkm1"] + 5},
        "kbase_checkm_old": kbase_link(cfg, s, sets["old"]),
        "originals_qc": {g["genome"]: {k: g.get(k) for k in ("checkm2_completeness", "checkm2_contamination",
                                                             "checkm1_completeness", "checkm1_contamination", "gunc_pass",
                                                             "n50", "n_contigs", "size_bp")} | {"population": orig_pop.get(g["genome"])}
                         for g in sets["orig"].values()},
        "qc_status": status, "skani_status": skani_status, "skani_edges_among_sets": n_edges,
        "derep310_originals": derep,
    }
    rows = [g for gs in ("new", "old", "orig") for g in sets[gs].values()]
    return {"pass": passed, "criteria": crit, "metrics": m, "informational": info}, rows


def kbase_link(cfg, s, old):
    """KBase CheckM 1.1.2 (kb_Msuite, reduced tree) on the same OLD bins vs CheckM1/CheckM2 here."""
    f = os.path.join(cfg["root"], "inputs", "kbase_checkm", f"checkm_{s}.tsv")
    if not os.path.exists(f):
        return {"status": "missing"}
    kb = {}
    for r in read_tsv(f):
        g = r["assembly_name"].split(".fasta__")[0]
        kb[g] = (float(r["completeness"]), float(r["contamination"]))
    pairs = [(kb[g], old[g]) for g in kb if g in old]
    d_c1 = [abs(k[0] - o["checkm1_completeness"]) for k, o in pairs if "checkm1_completeness" in o]
    d_k1 = [abs(k[1] - o["checkm1_contamination"]) for k, o in pairs if "checkm1_contamination" in o]
    return {"status": "ok", "bins_matched": len(pairs), "bins_in_kbase_table": len(kb),
            "kbase_mq": sum(1 for k in kb.values() if k[0] >= 50 and k[1] < 10),
            "kbase_hq": sum(1 for k in kb.values() if k[0] >= 90 and k[1] < 5),
            "checkm1_here_mq": sum(1 for _, o in pairs if o["mq_checkm1"]),
            "checkm2_here_mq": sum(1 for _, o in pairs if o["mq"]),
            "median_abs_diff_completeness_kbase_vs_checkm1": med(d_c1),
            "median_abs_diff_contamination_kbase_vs_checkm1": med(d_k1)}


def decide(cfg, state):
    samples = cfg["samples"]
    primary = next(s for s, S in samples.items() if S["role"] == "primary")
    secondary = next(s for s, S in samples.items() if S["role"] == "secondary")
    res, all_rows = {}, []
    for s in samples:
        st = state["samples"].get(s, {}).get("state", "unknown")
        entry = {"state": st, "state_info": state["samples"].get(s, {}).get("info")}
        if st == "complete":
            ev, rows = evaluate(cfg, s)
            entry.update(ev)
            all_rows += rows
        else:
            entry["pass"] = None
        res[s] = entry
    p, q = res[primary], res[secondary]
    notes = []
    if q["state"] in ("asm_failed", "asm_timeout"):
        if p["state"] != "complete":
            decision = "INCOMPLETE"
            notes.append(f"{secondary} assembly {q['state']}; {primary} processing incomplete")
        else:
            decision = "R1_ONLY_PASS" if p["pass"] else "NO_GO"
            notes.append(f"{secondary} Arm C assembly {q['state']}: decided on {primary} alone (pre-registered).")
    elif p["state"] != "complete" or q["state"] != "complete":
        decision = "INCOMPLETE"
        notes.append("a required processing step failed; see samples.*.state_info; needs Andrew")
    elif p["pass"] and q["pass"]:
        decision = "GO"
        notes.append("both samples pass: GO for the full salterns rollout")
    else:
        decision = "NO_GO"
        if p["pass"] and not q["pass"]:
            notes.append(f"only {primary} passes: NO_GO; the SF2-type decision needs Andrew")
        elif q["pass"] and not p["pass"]:
            notes.append(f"only {secondary} passes")
        else:
            notes.append("neither sample passes")
    return {"decision": decision, "generated_utc": ts(), "sf2_deadline_utc": cfg["sf2_deadline_utc"],
            "rule": RULE, "conventions": CONVENTIONS, "notes": notes, "samples": res}, all_rows


def fmt(v):
    if v is None:
        return "–"
    if isinstance(v, float):
        return f"{v:,.2f}" if abs(v) < 1000 else f"{v:,.0f}"
    return str(v)


def report(dec):
    L = [f"# Salterns binning test: decision {dec['decision']}", "", f"Generated {dec['generated_utc']}. "
         f"SF2 assembly deadline {dec['sf2_deadline_utc']}.", ""]
    L += [f"- {n}" for n in dec["notes"]] + [""]
    for s, e in dec["samples"].items():
        L.append(f"## {s}: {e['state']}" + ("" if e.get("pass") is None else (" — PASS" if e["pass"] else " — FAIL")))
        L.append("")
        if "criteria" not in e:
            L.append(f"Not evaluated: `{json.dumps(e.get('state_info'))}`")
            L.append("")
            continue
        L.append("| criterion | metric | NEW | OLD | threshold | pass |")
        L.append("|---|---|---|---|---|---|")
        for k in "abce":
            c = e["criteria"][k]
            L.append(f"| {k} | {c['metric']} | {fmt(c['new'])} | {fmt(c['old'])} | {c['threshold']} | {c['pass']} |")
        d = e["criteria"]["d"]
        for kk, lab in (("d1_contamination", "median CheckM2 contamination of MQ+ (%)"),
                        ("d2_gunc_pass_fraction", "GUNC pass fraction of MQ+")):
            c = d[kk]
            L.append(f"| d | {lab} | {fmt(c['new'])} | {fmt(c['old'])} | {c['threshold']} | {c['pass']} |")
        L.append("")
        inf = e["informational"]
        L.append(f"CheckM1 (informational): MQ+ NEW {inf['checkm1_a']['new']} vs OLD {inf['checkm1_a']['old']} "
                 f"(a {inf['checkm1_a']['pass']}); HQ NEW {inf['checkm1_b']['new']} vs OLD {inf['checkm1_b']['old']} "
                 f"(b {inf['checkm1_b']['pass']}).")
        L.append(f"KBase CheckM link (OLD bins): `{json.dumps(inf['kbase_checkm_old'])}`")
        L.append("")
        L.append("| set | bins | MQ+ | HQ | Mbp in MQ+ | median completeness | median contamination | median N50 | median contigs | GUNC pass | originals recovered |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for gs in ("new", "old"):
            m = e["metrics"][gs]
            L.append(f"| {gs.upper()} | {m['n_bins']} | {m['mq']} | {m['hq']} | {m['mbp_in_mq']} | {fmt(m['median_completeness_mq'])} | "
                     f"{fmt(m['median_contamination_mq'])} | {fmt(m['median_n50_mq'])} | {fmt(m['median_contigs_mq'])} | "
                     f"{fmt(m['gunc_pass_fraction_mq'])} | {m['originals_recovered']} {m['originals_recovered_names']} |")
        L.append("")
    L.append("## Rule (pre-registered)\n\n```\n" + dec["rule"] + "\n```\n\n## Conventions\n\n```\n" + dec["conventions"] + "\n```")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--state", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out_json = os.path.join(a.out, "DECISION.json")
    try:
        cfg = json.load(open(a.config))
        state = json.load(open(a.state))
        dec, rows = decide(cfg, state)
        cols = ["sample", "set", "genome", "size_bp", "n_contigs", "n50", "checkm2_completeness", "checkm2_contamination",
                "checkm2_model", "mq", "hq", "checkm1_completeness", "checkm1_contamination", "checkm1_marker_lineage",
                "mq_checkm1", "hq_checkm1", "gunc_pass", "gunc_css", "gunc_contamination_portion", "path"]
        with open(os.path.join(a.out, "mags.tsv"), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(fmt_cell(r.get(c)) for c in cols) + "\n")
        with open(os.path.join(a.out, "REPORT.md"), "w") as fh:
            fh.write(report(dec))
        tmp = out_json + ".part"
        with open(tmp, "w") as fh:
            json.dump(dec, fh, indent=1, default=list)
        os.replace(tmp, out_json)
        print(json.dumps({"decision": dec["decision"], "notes": dec["notes"]}, indent=1))
    except Exception as e:
        err = {"decision": "ERROR", "generated_utc": ts(), "error": f"{type(e).__name__}: {e}",
               "traceback": traceback.format_exc(), "rule": RULE}
        with open(out_json, "w") as fh:
            json.dump(err, fh, indent=1)
        print(err["traceback"], file=sys.stderr)
        sys.exit(1)


def fmt_cell(v):
    return "" if v is None else str(v)


if __name__ == "__main__":
    main()
