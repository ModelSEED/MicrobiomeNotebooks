#!/usr/bin/env python3
"""Stdlib helpers for the salterns Arm C rollout (driver.py task scripts call these with python -I).
btlib.py (copied unchanged from the reviewed binning test) still does the binning-specific steps
(merge-depth, abund, c2b-*, das-inputs); this file adds what the rollout needs on top.

Subcommands:
  asm-stats SCAFFOLDS MINLEN OUT_FA OUT_JSON   keep sequences >= MINLEN (length with N's, as btlib filter-fasta) and
                                               write size/contiguity stats of the input and of the kept set
  spades-meta --dir D --mode M --attempt K --mem GB --threads T --out JSON
                                               k list, SPAdes version, peak memory and time from params.txt/spades.log
  oom-check FILE...                            exit 0 (and print the lines) if any file shows an out-of-memory failure
  verify-import --bins DIR --checkm2 F --gunc F [--checkm1 F]
                                               every bin has exactly one CheckM2 / GUNC (/ CheckM1) row
  keep --sample S --bins DIR --checkm2 F --gunc F [--checkm1 F] --outdir D --kept-dir K
                                               apply the KEEP rule, write D/bins.tsv + D/kept.txt, copy kept FASTAs to
                                               K/<S>__<bin>.fa
  report --config CFG [--interim]              pooled tables + per-sample summary + REPORT.md under <root>/results/

KEEP rule (Andrew, 2026-11-13): CheckM2 completeness >= 50 AND CheckM2 contamination < 10 AND GUNC pass.
Comparisons use exact decimals on the printed CheckM2 values. GUNC pass = pass.GUNC == True in
GUNC.progenomes_2.1.maxCSS_level.tsv; a bin with no GUNC row, or a non-True value, does not pass (binning-test rule).
HQ = kept and CheckM2 completeness >= 90 and contamination < 5.
"""
import argparse
import csv
import datetime
import glob
import json
import os
import re
import shutil
import statistics
import sys
from decimal import Decimal, InvalidOperation

THRESH = (1000, 2500, 5000, 10000, 20000)
OOM_PATTERNS = [r"Cannot allocate memory", r"unable to allocate", r"std::bad_alloc", r"OutOfMemoryError",
                r"OS return value: -9\b", r"\bKilled\b", r"memory limit (exceeded|reached)", r"Out of memory"]


# ------------------------------------------------------------------------------------------------ FASTA helpers
def fasta_lengths(path):
    out, cur, n = [], None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out.append((cur, n))
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if cur is not None:
        out.append((cur, n))
    return out


def n50(lengths):
    ls = sorted(lengths, reverse=True)
    tot, acc = sum(ls), 0
    for x in ls:
        acc += x
        if 2 * acc >= tot:
            return x
    return 0


def len_stats(ls):
    d = {"n": len(ls), "bp": sum(ls), "n50": n50(ls), "max": max(ls) if ls else 0}
    for t in THRESH:
        sel = [x for x in ls if x >= t]
        d[f"n_ge{t}"] = len(sel)
        d[f"Mbp_ge{t}"] = round(sum(sel) / 1e6, 3)
    d["n50_ge1000"] = n50([x for x in ls if x >= 1000])
    d["n50_ge2500"] = n50([x for x in ls if x >= 2500])
    return d


def asm_stats(inp, minlen, out_fa, out_json):
    minlen = int(minlen)
    all_l, kept_l = [], []
    tmp = out_fa + ".tmp"
    with open(inp) as fh, open(tmp, "w") as fo:
        hdr, lines, n = None, [], 0

        def flush():
            if hdr is None:
                return
            all_l.append(n)
            if n >= minlen:
                kept_l.append(n)
                fo.write(hdr)
                fo.writelines(lines)

        for line in fh:
            if line.startswith(">"):
                flush()
                hdr, lines, n = line, [], 0
            else:
                lines.append(line)
                n += len(line.strip())
        flush()
    os.replace(tmp, out_fa)
    res = {"input": os.path.realpath(inp), "min_len": minlen, "all": len_stats(all_l), "kept": len_stats(kept_l),
           "utc": utc()}
    with open(out_json, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ------------------------------------------------------------------------------------------------ SPAdes
def _mem_gb(tok):
    m = re.match(r"^(\d+(?:\.\d+)?)([KMGT]?)$", tok)
    if not m:
        return None
    v = float(m.group(1))
    return v * {"K": 1e-6, "M": 1e-3, "G": 1.0, "T": 1e3, "": 1e-9}[m.group(2)]


def spades_meta(a):
    d = a.dir
    res = {"mode": a.mode, "attempt": a.attempt, "mem_limit_gb": a.mem, "threads": a.threads, "spades_dir": d,
           "utc": utc()}
    p = os.path.join(d, "params.txt")
    if os.path.exists(p):
        txt = open(p, errors="replace").read()
        m = re.search(r"^\s*k: \[([0-9, ]+)\]", txt, re.M)
        res["k_list"] = m.group(1).replace(" ", "") if m else None
        m = re.search(r"SPAdes version: (\S+)", txt)
        res["spades_version"] = m.group(1) if m else None
        m = re.search(r"^Command line: (.*)$", txt, re.M)
        res["command_line"] = m.group(1).replace("\t", " ").strip() if m else None
        m = re.search(r"Memory limit \(in Gb\): (\d+)", txt)
        res["mem_limit_gb_params"] = int(m.group(1)) if m else None
    peak, last = 0.0, None
    lg = os.path.join(d, "spades.log")
    if os.path.exists(lg):
        for line in open(lg, errors="replace"):
            m = re.match(r"^\s*(\d+:\d+:\d+\.\d+)\s+(\S+)\s*/\s*(\S+)\s+(INFO|WARN|ERROR|DEBUG)", line)
            if m:
                g = _mem_gb(m.group(3))
                if g:
                    peak = max(peak, g)
            if "Assembling time" in line or "SPAdes pipeline finished" in line:
                last = line.strip()
    res["peak_rss_gb_log"] = round(peak, 1)
    res["last_line"] = last
    w = os.path.join(d, "warnings.log")
    res["n_warnings"] = sum(1 for _ in open(w, errors="replace")) if os.path.exists(w) else 0
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


def oom_check(files):
    hits = []
    for f in files:
        if not os.path.exists(f):
            continue
        for line in open(f, errors="replace"):
            for pat in OOM_PATTERNS:
                if re.search(pat, line, re.I):
                    hits.append(f"{f}: {line.strip()[:200]}")
                    break
    for h in hits[:10]:
        print(h)
    sys.exit(0 if hits else 1)


# ------------------------------------------------------------------------------------------------ QC tables
def dec(x):
    try:
        return Decimal(str(x).strip())
    except (InvalidOperation, AttributeError):
        return None


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bin_files(d):
    return sorted(f for f in glob.glob(os.path.join(d, "*.fa")))


def bin_name(path):
    return os.path.basename(path)[:-3]


def gunc_key(name):
    return os.path.basename(name).split(".fa")[0]


def load_qc(checkm2, gunc, checkm1):
    c2 = {r["Name"]: r for r in read_tsv(checkm2)}
    gu = {}
    for r in read_tsv(gunc):
        gu[gunc_key(r["genome"])] = r
    c1 = {}
    for r in read_tsv(checkm1):
        if "Bin Id" in r:
            c1[r["Bin Id"]] = r
    return c2, gu, c1


def classify(name, c2, gu):
    """(keep, quality, reasons, fields) for one bin under the KEEP rule."""
    reasons = []
    r2, rg = c2.get(name), gu.get(gunc_key(name))
    comp = dec(r2["Completeness"]) if r2 else None
    cont = dec(r2["Contamination"]) if r2 else None
    if r2 is None or comp is None or cont is None:
        reasons.append("checkm2_missing")
    else:
        if not comp >= 50:
            reasons.append("completeness<50")
        if not cont < 10:
            reasons.append("contamination>=10")
    if rg is None:
        reasons.append("gunc_missing")
    elif str(rg.get("pass.GUNC", "")).strip() != "True":
        reasons.append("gunc_fail" if str(rg.get("pass.GUNC", "")).strip() == "False" else "gunc_na")
    keep = not reasons
    q = ""
    if keep:
        q = "HQ" if (comp >= 90 and cont < 5) else "MQ"
    return keep, q, reasons, comp, cont


def bins_table(sample, bins_dir, checkm2, gunc, checkm1):
    c2, gu, c1 = load_qc(checkm2, gunc, checkm1)
    rows = []
    for f in bin_files(bins_dir):
        b = bin_name(f)
        ls = [x for _, x in fasta_lengths(f)]
        keep, q, reasons, comp, cont = classify(b, c2, gu)
        r2, rg, r1 = c2.get(b, {}), gu.get(gunc_key(b), {}), c1.get(b, {})
        rows.append({
            "sample": sample, "bin": b, "mag_id": f"{sample}__{b}", "size_bp": sum(ls), "n_contigs": len(ls),
            "n50": n50(ls), "max_contig": max(ls) if ls else 0,
            "checkm2_completeness": r2.get("Completeness", ""), "checkm2_contamination": r2.get("Contamination", ""),
            "checkm2_model": r2.get("Completeness_Model_Used", ""),
            "checkm1_completeness": r1.get("Completeness", ""), "checkm1_contamination": r1.get("Contamination", ""),
            "checkm1_lineage": r1.get("Marker lineage", ""),
            "gunc_pass": rg.get("pass.GUNC", ""), "gunc_css": rg.get("clade_separation_score", ""),
            "gunc_contamination_portion": rg.get("contamination_portion", ""),
            "gunc_rrs": rg.get("reference_representation_score", ""), "gunc_level": rg.get("taxonomic_level", ""),
            "keep": "yes" if keep else "no", "quality": q, "drop_reasons": ";".join(reasons),
            "fasta": f,
        })
    return rows


BIN_COLS = ["sample", "bin", "mag_id", "keep", "quality", "drop_reasons", "checkm2_completeness",
            "checkm2_contamination", "checkm2_model", "gunc_pass", "gunc_css", "gunc_contamination_portion",
            "gunc_rrs", "gunc_level", "checkm1_completeness", "checkm1_contamination", "checkm1_lineage", "size_bp",
            "n_contigs", "n50", "max_contig", "fasta"]


def write_tsv(path, rows, cols):
    tmp = path + ".part"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def verify_import(a):
    bins = [bin_name(f) for f in bin_files(a.bins)]
    assert bins, f"no bins in {a.bins}"
    c2, gu, c1 = load_qc(a.checkm2, a.gunc, a.checkm1)
    n2 = len(read_tsv(a.checkm2))
    ng = len(read_tsv(a.gunc))
    msg = {"bins": len(bins), "checkm2_rows": n2, "gunc_rows": ng}
    assert n2 == len(bins) and set(c2) == set(bins), f"CheckM2 rows do not match the bins: {msg}"
    assert ng == len(bins) and set(gu) == {gunc_key(b) for b in bins}, f"GUNC rows do not match the bins: {msg}"
    if a.checkm1:
        n1 = len(read_tsv(a.checkm1))
        msg["checkm1_rows"] = n1
        assert n1 == len(bins) and set(c1) == set(bins), f"CheckM1 rows do not match the bins: {msg}"
    print(json.dumps(msg))


def keep_cmd(a):
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(a.kept_dir, exist_ok=True)
    rows = bins_table(a.sample, a.bins, a.checkm2, a.gunc, a.checkm1)
    # idempotent: remove this sample's earlier copies first
    for f in glob.glob(os.path.join(a.kept_dir, f"{a.sample}__*.fa")):
        os.remove(f)
    kept = []
    for r in rows:
        if r["keep"] == "yes":
            dst = os.path.join(a.kept_dir, f"{r['mag_id']}.fa")
            shutil.copyfile(r["fasta"], dst + ".part")
            os.replace(dst + ".part", dst)
            kept.append(r["mag_id"])
    write_tsv(os.path.join(a.outdir, "bins.tsv"), rows, BIN_COLS)
    with open(os.path.join(a.outdir, "kept.txt"), "w") as fh:
        fh.writelines(k + "\n" for k in kept)
    s = summarize_rows(rows)
    with open(os.path.join(a.outdir, "summary.json"), "w") as fh:
        json.dump(s, fh, indent=1)
    print(json.dumps(s, indent=1))


def med(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def summarize_rows(rows):
    kept = [r for r in rows if r["keep"] == "yes"]
    reasons = {}
    for r in rows:
        for x in filter(None, r["drop_reasons"].split(";")):
            reasons[x] = reasons.get(x, 0) + 1
    fl = lambda v: float(v) if v not in ("", None) else None
    return {"bins": len(rows), "kept": len(kept), "kept_HQ": sum(r["quality"] == "HQ" for r in kept),
            "kept_MQ": sum(r["quality"] == "MQ" for r in kept), "dropped": len(rows) - len(kept),
            "drop_reasons": reasons,
            "dropped_only_gunc": sum(1 for r in rows if r["drop_reasons"] in ("gunc_fail", "gunc_na", "gunc_missing")),
            "kept_median_n50": med([r["n50"] for r in kept]),
            "kept_Mbp": round(sum(r["size_bp"] for r in kept) / 1e6, 2),
            "kept_median_completeness": med([fl(r["checkm2_completeness"]) for r in kept]),
            "kept_median_contamination": med([fl(r["checkm2_contamination"]) for r in kept])}


# ------------------------------------------------------------------------------------------------ report
def jload(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def report(a):
    cfg = json.load(open(a.config))
    root = cfg["root"]
    res = os.path.join(root, "results")
    os.makedirs(res, exist_ok=True)
    done = lambda n: os.path.exists(os.path.join(root, ".done", n))
    failed = lambda n: os.path.exists(os.path.join(root, ".failed", n))
    # GTDB-Tk
    gt = {}
    for f in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
        for r in read_tsv(os.path.join(root, "gtdbtk", f)):
            gt[r["user_genome"]] = r
    # 2026 KBase-selected MAGs (mags_old_vs_new.tsv; KBase CheckM 1.1.2 values)
    old_tab = read_tsv(cfg.get("old_table"))
    kb = {}
    for r in old_tab:
        if r.get("set") == "re-assembly (2026)":
            kb.setdefault(r["sample"], []).append(r)
    state = (jload(os.path.join(root, "STATE.json")) or {}).get("samples", {})
    all_bins, old_bins, summ = [], [], []
    for s in sorted(cfg["samples"], key=lambda x: cfg["samples"][x]["order"]):
        S = cfg["samples"][s]
        unit = S.get("unit", s)
        I = S.get("import_from") or {}
        row = {"sample": s, "unit_2026": unit, "site": S["group_name"], "n_coverage_libs": len(S["group"]),
               "imported_from_binning_test": "yes" if (I and not I.get("staged")) else "no",
               "processed_on": I.get("source", "") if I.get("staged") else cfg.get("host", "")}
        row["status"] = "complete" if done(f"keep:{s}") else state.get(s, {}).get("state", "running")
        am = jload(os.path.join(root, s, "asm", "asm_meta.json")) or {}
        row.update({"assembly_mode": am.get("mode", ""), "spades_attempt": am.get("attempt", ""),
                    "k_list": am.get("k_list", ""), "spades_mem_limit_gb": am.get("mem_limit_gb", ""),
                    "spades_threads": am.get("threads", ""), "spades_peak_rss_gb": am.get("peak_rss_gb_log", "")})
        ast = (jload(os.path.join(root, s, "asm", "asm_stats.json")) or {}).get("all", {})
        for t in (1000, 2500, 10000):
            row[f"asm_Mbp_ge{t}"] = ast.get(f"Mbp_ge{t}", "")
        row["asm_n50_ge1000"] = ast.get("n50_ge1000", "")
        row["asm_n50_ge2500"] = ast.get("n50_ge2500", "")
        rows = read_tsv(os.path.join(root, s, "keep", "bins.tsv"))
        for r in rows:
            r["size_bp"] = int(r["size_bp"])
            r["n50"] = int(r["n50"])
            g = gt.get(r["mag_id"], {})
            r["gtdb_classification"] = g.get("classification", "")
            r["gtdb_closest_reference"] = g.get("closest_genome_reference", "")
            r["gtdb_closest_ani"] = g.get("closest_genome_ani", "")
            r["gtdb_note"] = g.get("note", "")
            r["gtdb_warnings"] = g.get("warnings", "")
        all_bins += rows
        sm = summarize_rows(rows) if rows else {}
        row.update({"new_bins": sm.get("bins", ""), "new_kept": sm.get("kept", ""), "new_kept_HQ": sm.get("kept_HQ", ""),
                    "new_kept_MQ": sm.get("kept_MQ", ""), "new_dropped": sm.get("dropped", ""),
                    "drop_completeness<50": sm.get("drop_reasons", {}).get("completeness<50", "") if sm else "",
                    "drop_contamination>=10": sm.get("drop_reasons", {}).get("contamination>=10", "") if sm else "",
                    "drop_gunc": (sm.get("drop_reasons", {}).get("gunc_fail", 0) + sm.get("drop_reasons", {}).get(
                        "gunc_na", 0) + sm.get("drop_reasons", {}).get("gunc_missing", 0)) if sm else "",
                    "dropped_only_by_gunc": sm.get("dropped_only_gunc", ""),
                    "new_kept_median_n50": sm.get("kept_median_n50", ""), "new_kept_Mbp": sm.get("kept_Mbp", ""),
                    "new_kept_median_completeness": sm.get("kept_median_completeness", ""),
                    "new_kept_median_contamination": sm.get("kept_median_contamination", "")})
        # 2026 side, (i) KBase-selected MAGs with KBase CheckM 1.1.2 values
        k = kb.get(unit, [])
        row["old2026_kbase_mags"] = len(k)
        row["old2026_kbase_HQ"] = sum(1 for r in k if r.get("quality", "").startswith("high"))
        row["old2026_kbase_median_n50"] = med([int(r["n50"]) for r in k if r.get("n50")])
        # (ii) all 2026 DAS Tool bins of the unit, scored here with CheckM2 + GUNC (same KEEP rule)
        od = os.path.join(root, "inputs", "old", s)
        c2f = os.path.join(root, "qc", "old", s, "checkm2", "quality_report.tsv")
        guf = glob.glob(os.path.join(root, "qc", "old", s, "gunc", "GUNC.*.maxCSS_level.tsv"))
        if os.path.isdir(od) and os.path.exists(c2f) and guf:
            orows = bins_table(unit, od, c2f, guf[0], None)
            for r in orows:
                r["rollout_sample"] = s
            old_bins += orows
            osm = summarize_rows(orows)
            row.update({"old2026_checkm2_source": "binning_test" if (I and not I.get("staged")) else "rollout",
                        "old2026_dastool_bins": osm["bins"], "old2026_checkm2_kept": osm["kept"],
                        "old2026_checkm2_kept_HQ": osm["kept_HQ"], "old2026_checkm2_kept_median_n50": osm["kept_median_n50"],
                        "old2026_checkm2_kept_Mbp": osm["kept_Mbp"]})
            if sm:
                row["delta_kept_new_minus_old"] = sm["kept"] - osm["kept"]
                row["delta_HQ_new_minus_old"] = sm["kept_HQ"] - osm["kept_HQ"]
                if osm["kept_median_n50"] and sm.get("kept_median_n50"):
                    row["n50_ratio_new_over_old"] = round(sm["kept_median_n50"] / osm["kept_median_n50"], 2)
        else:
            row["old2026_checkm2_source"] = "not available"
        summ.append(row)
    scols = ["sample", "unit_2026", "site", "n_coverage_libs", "status", "imported_from_binning_test", "processed_on",
             "assembly_mode",
             "spades_attempt", "k_list", "spades_mem_limit_gb", "spades_threads", "spades_peak_rss_gb",
             "asm_Mbp_ge1000", "asm_Mbp_ge2500", "asm_Mbp_ge10000", "asm_n50_ge1000", "asm_n50_ge2500",
             "new_bins", "new_kept", "new_kept_HQ", "new_kept_MQ", "new_dropped", "drop_completeness<50",
             "drop_contamination>=10", "drop_gunc", "dropped_only_by_gunc", "new_kept_median_n50", "new_kept_Mbp",
             "new_kept_median_completeness", "new_kept_median_contamination", "old2026_kbase_mags", "old2026_kbase_HQ",
             "old2026_kbase_median_n50", "old2026_checkm2_source", "old2026_dastool_bins", "old2026_checkm2_kept",
             "old2026_checkm2_kept_HQ", "old2026_checkm2_kept_median_n50", "old2026_checkm2_kept_Mbp",
             "delta_kept_new_minus_old", "delta_HQ_new_minus_old", "n50_ratio_new_over_old"]
    write_tsv(os.path.join(res, "samples_summary.tsv"), summ, scols)
    bcols = BIN_COLS[:-1] + ["gtdb_classification", "gtdb_closest_reference", "gtdb_closest_ani", "gtdb_note",
                             "gtdb_warnings"]
    write_tsv(os.path.join(res, "bins_all.tsv"), all_bins, bcols)
    write_tsv(os.path.join(res, "mags_kept.tsv"), [r for r in all_bins if r["keep"] == "yes"], bcols)
    write_tsv(os.path.join(res, "bins_dropped.tsv"), [r for r in all_bins if r["keep"] != "yes"], bcols)
    write_tsv(os.path.join(res, "old2026_bins_checkm2_gunc.tsv"), old_bins, ["rollout_sample"] + BIN_COLS[:-1])
    write_md(res, cfg, summ, all_bins, old_bins, gt, interim=a.interim)
    print(f"report: {len(summ)} samples, {len(all_bins)} bins, {sum(r['keep'] == 'yes' for r in all_bins)} kept -> {res}")


def _fmt(v):
    if v in (None, ""):
        return "–"
    if isinstance(v, float):
        return f"{v:,.1f}" if abs(v) < 1e5 else f"{v:,.0f}"
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)


def write_md(res, cfg, summ, all_bins, old_bins, gt, interim=False):
    kept = [r for r in all_bins if r["keep"] == "yes"]
    L = [f"# Salterns Arm C rollout: summary{' (INTERIM)' if interim else ''}", "",
         f"Generated {utc()} by `rolib.py report` from `{cfg['root']}`.", "",
         "KEEP rule: CheckM2 completeness >= 50 AND CheckM2 contamination < 10 AND GUNC pass "
         "(`pass.GUNC == True`; a missing GUNC row does not pass). HQ = kept with CheckM2 >= 90 / < 5. "
         "CheckM1 `--reduced_tree` values are informational only.", "",
         "Per-sample assemblies only. The 2026 site co-assemblies (R1/R2/R2A/SF2 sediment) are NOT redone; their "
         "2026 MAGs remain as they are.", ""]
    n_s = len(summ)
    comp = sum(1 for r in summ if r["status"] == "complete")
    L += [f"- Samples: {n_s}; complete: {comp}; failed: {sum(1 for r in summ if r['status'] == 'failed')}; "
          f"running: {sum(1 for r in summ if r['status'] == 'running')}.",
          f"- NEW DAS Tool bins: {len(all_bins):,}; kept: {len(kept):,} (HQ {sum(r['quality'] == 'HQ' for r in kept):,}, "
          f"MQ {sum(r['quality'] == 'MQ' for r in kept):,}); dropped: {len(all_bins) - len(kept):,}.",
          f"- Kept MAGs with a GTDB-Tk (R232) classification: {sum(1 for r in kept if r.get('gtdb_classification')):,}.",
          f"- Assembly mode: " + ", ".join(f"{m or 'none yet'} {sum(1 for r in summ if r['assembly_mode'] == m)}"
                                          for m in sorted({r['assembly_mode'] for r in summ})) + ".", ""]
    elsewhere = sorted({r["processed_on"] for r in summ if r.get("processed_on") and r["processed_on"] != cfg.get("host")})
    for h in elsewhere:
        mine = [r for r in summ if r.get("processed_on") == h]
        L += [f"- Assembled, binned and QC'd on {h} with the same driver and recipe, then imported here (md5-verified; "
              f"the KEEP rule re-applied here gave the same kept set): "
              + (", ".join(r["sample"] for r in mine if r["status"] == "complete") or "none yet") + "."
              + ("" if all(r["status"] == "complete" for r in mine) else
                 " Not imported (yet): " + ", ".join(f"{r['sample']} ({r['status']})" for r in mine
                                                       if r["status"] != "complete") + "."), ""]
    okb = [r for r in old_bins]
    if okb:
        osum = summarize_rows(okb)
        L += [f"- 2026 DAS Tool bins of the same units scored here with CheckM2 + GUNC: {osum['bins']:,} bins, "
              f"{osum['kept']:,} pass the same KEEP rule (HQ {osum['kept_HQ']:,}).", ""]
    L += ["## Per sample", "",
          "| sample | status | mode | asm Mbp >=2.5 kb | asm N50 >=1 kb | NEW bins | kept | HQ | kept median N50 | "
          "2026 KBase MAGs (HQ) | 2026 bins CheckM2-kept (HQ) | 2026 kept median N50 |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in summ:
        L.append(f"| {r['sample']} | {r['status']} | {r['assembly_mode'] or '–'} | {_fmt(r.get('asm_Mbp_ge2500'))} | "
                 f"{_fmt(r.get('asm_n50_ge1000'))} | {_fmt(r.get('new_bins'))} | {_fmt(r.get('new_kept'))} | "
                 f"{_fmt(r.get('new_kept_HQ'))} | {_fmt(r.get('new_kept_median_n50'))} | "
                 f"{_fmt(r.get('old2026_kbase_mags'))} ({_fmt(r.get('old2026_kbase_HQ'))}) | "
                 f"{_fmt(r.get('old2026_checkm2_kept'))} ({_fmt(r.get('old2026_checkm2_kept_HQ'))}) | "
                 f"{_fmt(r.get('old2026_checkm2_kept_median_n50'))} |")
    L += ["", "Columns: '2026 KBase MAGs' = the KBase-CheckM-selected 2026 MAGs of the same unit "
          "(`comparison/mags_old_vs_new.tsv`, KBase CheckM 1.1.2). '2026 bins CheckM2-kept' = all 2026 DAS Tool bins of "
          "the unit re-scored here with CheckM2 + GUNC under the same KEEP rule (CheckM2 on both sides).", "",
          "Files: `samples_summary.tsv`, `bins_all.tsv`, `mags_kept.tsv`, `bins_dropped.tsv`, "
          "`old2026_bins_checkm2_gunc.tsv`; kept FASTAs in `../mags_kept/`; GTDB-Tk output in `../gtdbtk/`."]
    with open(os.path.join(res, "REPORT.md.part"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    os.replace(os.path.join(res, "REPORT.md.part"), os.path.join(res, "REPORT.md"))


# ------------------------------------------------------------------------------------------------ main
def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    cmd = sys.argv[1]
    if cmd == "asm-stats":
        asm_stats(*sys.argv[2:6])
        return
    if cmd == "oom-check":
        oom_check(sys.argv[2:])
        return
    ap = argparse.ArgumentParser(prog=f"rolib.py {cmd}")
    if cmd == "spades-meta":
        for x in ("--dir", "--mode", "--attempt", "--mem", "--threads", "--out"):
            ap.add_argument(x, required=True)
        spades_meta(ap.parse_args(sys.argv[2:]))
    elif cmd in ("verify-import", "keep"):
        ap.add_argument("--bins", required=True)
        ap.add_argument("--checkm2", required=True)
        ap.add_argument("--gunc", required=True)
        ap.add_argument("--checkm1")
        if cmd == "keep":
            ap.add_argument("--sample", required=True)
            ap.add_argument("--outdir", required=True)
            ap.add_argument("--kept-dir", required=True)
            keep_cmd(ap.parse_args(sys.argv[2:]))
        else:
            verify_import(ap.parse_args(sys.argv[2:]))
    elif cmd == "report":
        ap.add_argument("--config", required=True)
        ap.add_argument("--interim", action="store_true")
        report(ap.parse_args(sys.argv[2:]))
    else:
        raise SystemExit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
