#!/usr/bin/env python3
"""Evaluate the tadpole-pilot assemblies against the (iii) MAG-source reconstruction and the 2026 KBase assembly.

Stdlib only (runs with /scratch1/afreiburger/salterns_2026/envs/spades3153/bin/python, Python 3.9).

For every <root>/<sample>/arm{A,B,C}/spades/{scaffolds,contigs}.fasta that exists:
  n and Mbp in sequences >= 1, 2.5, 5, 10, 20 kb; N50 of sequences >= 1 kb and >= 2.5 kb; max length;
  length-weighted median SPAdes k-mer coverage of sequences >= 2.5 kb (from NODE_..._cov_X names);
  ratios to the (iii) estimates (orig_*_est) and to the 2026 KBase values (new_*) from targets.tsv.
'go' = Mbp >= 2.5 kb AND Mbp >= 10 kb both reach GO_FRACTION (0.85) of the (iii) estimate (README decision rule).
Primary comparison: scaffolds.fasta (the 2026 KBase assemblies are scaffolds.fasta >= 1 kb); contigs.fasta is
also reported because whether (iii) used contigs or scaffolds is unknown. Lengths include N's, as in KBase dna_size.

Also collects read-preparation statistics from <sample>/logs/*.log (BBDuk, BBMerge, Tadpole) and the SPAdes
insert-size estimate from spades.log.

Usage: python eval_assemblies.py --root /scratch1/afreiburger/salterns_2026/tadpole_pilot \
           --targets scripts/targets.tsv --out-prefix eval/pilot
Writes <out-prefix>_assemblies.tsv and <out-prefix>_prep.tsv and prints a summary.
"""
import argparse
import csv
import glob
import os
import re
import sys

THRESH = [1000, 2500, 5000, 10000, 20000]
GO_FRACTION = 0.85


def read_lengths(path):
    lens, covs = [], []
    name, n = None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens.append(n)
                    covs.append(_cov(name))
                name, n = line[1:].strip(), 0
            else:
                n += len(line.strip())
        if name is not None:
            lens.append(n)
            covs.append(_cov(name))
    return lens, covs


def _cov(name):
    m = re.search(r"_cov_([0-9.]+)", name)
    return float(m.group(1)) if m else None


def n50(lens):
    s = sorted(lens, reverse=True)
    half, acc = sum(s) / 2.0, 0
    for x in s:
        acc += x
        if acc >= half:
            return x
    return 0


def weighted_median(pairs):
    pairs = sorted((c, l) for c, l in pairs if c is not None)
    if not pairs:
        return None
    half, acc = sum(l for _, l in pairs) / 2.0, 0
    for c, l in pairs:
        acc += l
        if acc >= half:
            return c
    return None


def stats(lens, covs):
    out = {}
    for t in THRESH:
        sel = [x for x in lens if x >= t]
        out["n_ge%d" % t] = len(sel)
        out["Mbp_ge%d" % t] = round(sum(sel) / 1e6, 2)
    out["N50_ge1000"] = n50([x for x in lens if x >= 1000])
    out["N50_ge2500"] = n50([x for x in lens if x >= 2500])
    out["max_len"] = max(lens) if lens else 0
    wm = weighted_median([(c, l) for l, c in zip(lens, covs) if l >= 2500])
    out["kmer_cov_wmedian_ge2500"] = round(wm, 2) if wm is not None else ""
    return out


def ratio(a, b):
    try:
        a, b = float(a), float(b)
        return round(a / b, 3) if b else ""
    except (TypeError, ValueError):
        return ""


def grep_last(path, pattern):
    """Last match of pattern (regex with groups) in a text file, or None."""
    if not os.path.exists(path):
        return None
    hit = None
    rx = re.compile(pattern, re.IGNORECASE)
    with open(path, errors="replace") as fh:
        for line in fh:
            m = rx.search(line)
            if m:
                hit = m
    return hit


def prep_stats(sdir):
    rows = []
    logs = os.path.join(sdir, "logs")
    for path in sorted(glob.glob(os.path.join(logs, "*.log"))):
        step = os.path.basename(path)[:-4]
        rec = {"step": step}
        m = grep_last(path, r"^Input:\s+(\d+) reads\s+(\d+) bases")
        if m:
            rec["input_reads"], rec["input_bases"] = m.group(1), m.group(2)
        m = grep_last(path, r"^Result:\s+(\d+) reads \(([\d.]+)%\)\s+(\d+) bases \(([\d.]+)%\)")
        if m:
            rec.update(result_reads=m.group(1), result_pct_reads=m.group(2),
                       result_bases=m.group(3), result_pct_bases=m.group(4))
        m = grep_last(path, r"^(KTrimmed|Contaminants|Total Removed):\s+(\d+) reads \(([\d.]+)%\)\s+(\d+) bases \(([\d.]+)%\)")
        if m:
            rec["removed_or_trimmed"] = "%s %s reads (%s%%) %s bases (%s%%)" % m.groups()
        m = grep_last(path, r"^Pairs:\s+(\d+)")
        if m:
            rec["bbmerge_pairs"] = m.group(1)
        m = grep_last(path, r"^Joined:\s+(\d+)\s+([\d.]+)%")
        if m:
            rec["bbmerge_joined"], rec["bbmerge_joined_pct"] = m.group(1), m.group(2)
        m = grep_last(path, r"^Ambiguous:\s+(\d+)\s+([\d.]+)%")
        if m:
            rec["bbmerge_ambiguous_pct"] = m.group(2)
        m = grep_last(path, r"^No Solution:\s+(\d+)\s+([\d.]+)%")
        if m:
            rec["bbmerge_nosolution_pct"] = m.group(2)
        m = grep_last(path, r"^Avg Insert:\s+([\d.]+)")
        if m:
            rec["bbmerge_avg_insert"] = m.group(1)
        m = grep_last(path, r"^Errors detected:\s+(\d+)")
        if m:
            rec["errors_detected"] = m.group(1)
        m = grep_last(path, r"^Errors corrected:\s+(\d+)")
        if m:
            rec["errors_corrected"] = m.group(1)
        m = grep_last(path, r"^Reads with errors detected:\s+(\d+)\s+\(([\d.]+)%")
        if m:
            rec["reads_with_errors_detected"], rec["reads_with_errors_detected_pct"] = m.group(1), m.group(2)
        m = grep_last(path, r"^Reads fully corrected:\s+(\d+)\s+\(([\d.]+)% of detected")
        if m:
            rec["reads_fully_corrected"], rec["reads_fully_corrected_pct_of_detected"] = m.group(1), m.group(2)
        m = grep_last(path, r"^Unique Kmers:\s+(\d+)")
        if m:
            rec["unique_kmers"] = m.group(1)
        m = grep_last(path, r"^(?:Total time|Time):\s+([\d.]+) seconds")
        if m:
            rec["time_s"] = m.group(1)
        m = grep_last(path, r"^Output:\s+(\d+) reads\s+(\d+) bases")
        if m:
            rec["output_reads"], rec["output_bases"] = m.group(1), m.group(2)
        rows.append(rec)
    return rows


def _mem_gb(tok):
    """SPAdes memory token ('512M', '24G') -> GB."""
    v, u = float(tok[:-1]), tok[-1].upper()
    return v / 1024.0 if u == "M" else (v * 1024.0 if u == "T" else v)


def resources(sdir):
    """Wall time per step from the driver logs (driver*.log, logs/armA_chain.out: '[ts] start X' / '[ts] done  X'),
    peak memory per step: BBTools = largest 'used=' Java heap figure in the step log (the heap itself is fixed by
    -Xmx/-Xms, so this is what the k-mer tables actually needed); SPAdes = largest peak-RSS field ('cur / PEAK') in
    spades.log (per-k processes, so the max over the log is the run's peak)."""
    import datetime
    starts, rows = {}, {}
    logs = sorted(glob.glob(os.path.join(sdir, "driver*.log"))) + [os.path.join(sdir, "logs", "armA_chain.out")]
    rx = re.compile(r"^\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\] (start|done ) (\S+)")
    for path in logs:
        if not os.path.exists(path):
            continue
        for line in open(path, errors="replace"):
            m = rx.match(line)
            if not m:
                continue
            t = datetime.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
            if m.group(2) == "start":
                starts[m.group(3)] = t
            elif m.group(3) in starts:
                st = starts[m.group(3)]
                rows[m.group(3)] = dict(step=m.group(3), start=st.strftime("%F %T"), end=t.strftime("%F %T"),
                                        wall_h=round((t - st).total_seconds() / 3600.0, 2))
    # Fallback for A4_metaspades: relaunching run_sample.sh (e.g. to add Arm C later) truncates logs/armA_chain.out,
    # which held its start/done lines. Start = the driver's 'Arm A chain ... PID' line, end = mtime of .done/A4_metaspades.
    marker = os.path.join(sdir, ".done", "A4_metaspades")
    if "A4_metaspades" not in rows and os.path.exists(marker):
        chain = None
        for path in logs:
            if os.path.exists(path):
                for line in open(path, errors="replace"):
                    m = re.match(r"^\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\] Arm A chain", line)
                    if m and chain is None:
                        chain = datetime.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
        end = datetime.datetime.fromtimestamp(os.path.getmtime(marker))
        rows["A4_metaspades"] = dict(step="A4_metaspades", start=chain.strftime("%F %T") if chain else "",
                                     end=end.strftime("%F %T") + " (.done mtime)",
                                     wall_h=round((end - chain).total_seconds() / 3600.0, 2) if chain else "")
        starts.pop("A4_metaspades", None)
    for name, st in starts.items():                  # started, not (yet) done
        rows.setdefault(name, dict(step=name, start=st.strftime("%F %T"), end="running", wall_h=""))
    spades_dirs = {"A4_metaspades": "armA", "B3_metaspades": "armB", "C1_spades_default": "armC"}
    for name, r in rows.items():
        peak = None
        if name in spades_dirs:
            log = os.path.join(sdir, spades_dirs[name], "spades", "spades.log")
            if os.path.exists(log):
                for line in open(log, errors="replace"):
                    m = re.match(r"^\s*\d+:\d+:\d+\.\d+\s+(\d+[MGT]) / (\d+[MGT])\s", line)
                    if m:
                        peak = max(peak or 0, _mem_gb(m.group(2)))
            r["peak_mem_kind"] = "SPAdes peak RSS"
        else:
            log = os.path.join(sdir, "logs", name + ".log")
            if os.path.exists(log):
                for m in re.finditer(r"used=(\d+)m", open(log, errors="replace").read()):
                    peak = max(peak or 0, int(m.group(1)) / 1024.0)
            r["peak_mem_kind"] = "Java heap used (max reported)"
            m = grep_last(log, r"-Xmx(\d+[gGmM])")
            r["java_xmx"] = m.group(1) if m else ""
        r["peak_mem_GB"] = round(peak, 1) if peak is not None else ""
    return sorted(rows.values(), key=lambda r: r["start"])


def spades_log_info(spdir):
    log = os.path.join(spdir, "spades.log")
    info = {}
    m = grep_last(log, r"Insert size = ([\d.]+), deviation = ([\d.]+)")
    if m:
        info["insert_size"], info["insert_sd"] = m.group(1), m.group(2)
    m = grep_last(log, r"SPAdes version: (\S+)|SPAdes genome assembler v(\S+)")
    if m:
        info["spades_version"] = m.group(1) or m.group(2)
    m = grep_last(log, r"(\d+:\d+:\d+\.\d+)\s+\S+\s+\S+\s+INFO.*")
    if m:
        info["last_log_time"] = m.group(1)
    info["finished"] = "yes" if grep_last(log, r"SPAdes pipeline finished") else "no"
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--samples", nargs="*", help="default: every <root>/<sample> with an arm dir")
    a = ap.parse_args()

    targets = {r["sample"]: r for r in csv.DictReader(open(a.targets), delimiter="\t")}
    samples = a.samples or sorted(
        os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(a.root, "*", "armA")))

    asm_rows, prep_rows, res_rows = [], [], []
    for s in samples:
        sdir = os.path.join(a.root, s)
        t = targets.get(s, {})
        for rec in prep_stats(sdir):
            prep_rows.append(dict(sample=s, **rec))
        for rec in resources(sdir):
            res_rows.append(dict(sample=s, **rec))
        for arm in ("armA", "armB", "armC"):
            spdir = os.path.join(sdir, arm, "spades")
            if not os.path.isdir(spdir):
                continue
            info = spades_log_info(spdir)
            for kind in ("scaffolds", "contigs"):
                fa = os.path.join(spdir, kind + ".fasta")
                row = dict(sample=s, arm=arm, kind=kind, spades_finished=info.get("finished", "no"),
                           spades_version=info.get("spades_version", ""),
                           insert_size=info.get("insert_size", ""), insert_sd=info.get("insert_sd", ""))
                if not os.path.exists(fa):
                    row["status"] = "not finished"
                    asm_rows.append(row)
                    continue
                lens, covs = read_lengths(fa)
                row["status"] = "ok"
                row.update(stats(lens, covs))
                for key, okey, nkey in (("Mbp_ge2500", "orig_Mbp_ge2500_est", "new_Mbp_ge2500"),
                                        ("Mbp_ge5000", "orig_Mbp_ge5000_est", "new_Mbp_ge5000"),
                                        ("Mbp_ge10000", "orig_Mbp_ge10000_est", "new_Mbp_ge10000"),
                                        ("Mbp_ge20000", "orig_Mbp_ge20000_est", "new_Mbp_ge20000"),
                                        ("N50_ge2500", "orig_N50_ge2500_est", "new_N50_ge2500")):
                    row["iii_" + key] = t.get(okey, "")
                    row["kbase2026_" + key] = t.get(nkey, "")
                    row["over_iii_" + key] = ratio(row[key], t.get(okey))
                    row["over_kbase2026_" + key] = ratio(row[key], t.get(nkey))
                row["kbase2026_Mbp_ge1000"] = t.get("new_Mbp_ge1000", "")
                row["over_kbase2026_Mbp_ge1000"] = ratio(row["Mbp_ge1000"], t.get("new_Mbp_ge1000"))
                r1, r2 = row["over_iii_Mbp_ge2500"], row["over_iii_Mbp_ge10000"]
                row["go_85pct_of_iii"] = ("yes" if (r1 != "" and r2 != "" and r1 >= GO_FRACTION and r2 >= GO_FRACTION)
                                          else ("no" if r1 != "" else ""))
                asm_rows.append(row)

    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    for rows, suffix in ((asm_rows, "_assemblies.tsv"), (prep_rows, "_prep.tsv"), (res_rows, "_resources.tsv")):
        if not rows:
            continue
        cols = []
        for r in rows:
            cols += [k for k in r if k not in cols]
        with open(a.out_prefix + suffix, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    print("Read preparation")
    for r in prep_rows:
        bits = [r["sample"], r["step"]]
        for k in ("result_pct_bases", "bbmerge_joined_pct", "bbmerge_avg_insert", "reads_with_errors_detected_pct",
                  "errors_corrected", "time_s"):
            if r.get(k):
                bits.append("%s=%s" % (k, r[k]))
        print("  " + "  ".join(bits))
    print("\nAssemblies (Mbp >= 2.5 kb | >= 10 kb | >= 20 kb | N50 >= 2.5 kb) vs (iii) and 2026 KBase")
    for r in asm_rows:
        if r.get("status") != "ok":
            print("  %-9s %s %-9s %s (SPAdes finished: %s)" % (r["sample"], r["arm"], r["kind"], r["status"],
                                                              r["spades_finished"]))
            continue
        print("  %-9s %s %-9s %7.1f | %6.1f | %6.1f | %6d   (iii) %s | %s | %s | %s   2026 %s | %s | %s | %s   "
              "x(iii) %s / %s   go=%s" % (
                  r["sample"], r["arm"], r["kind"], r["Mbp_ge2500"], r["Mbp_ge10000"], r["Mbp_ge20000"],
                  r["N50_ge2500"], r["iii_Mbp_ge2500"], r["iii_Mbp_ge10000"], r["iii_Mbp_ge20000"],
                  r["iii_N50_ge2500"], r["kbase2026_Mbp_ge2500"], r["kbase2026_Mbp_ge10000"],
                  r["kbase2026_Mbp_ge20000"], r["kbase2026_N50_ge2500"], r["over_iii_Mbp_ge2500"],
                  r["over_iii_Mbp_ge10000"], r["go_85pct_of_iii"]))
    print("\nWall time and peak memory per step")
    for r in res_rows:
        print("  %-9s %-20s %s -> %-19s %6s h  %6s GB (%s)" % (r["sample"], r["step"], r["start"], r["end"],
                                                              r["wall_h"], r.get("peak_mem_GB", ""),
                                                              r.get("peak_mem_kind", "")))
    print("\nWrote %s_assemblies.tsv, %s_prep.tsv and %s_resources.tsv" % (a.out_prefix, a.out_prefix, a.out_prefix))
    return 0


if __name__ == "__main__":
    sys.exit(main())
