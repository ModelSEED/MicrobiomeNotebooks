#!/usr/bin/env python
"""READ-ONLY: in every tracked csv/tsv (and tsv-like .txt) that has a column of bin names, find
  (a) MAG-ID columns: cells equal to (or embedding) the current (R220) MAG iterativeID of the row's bin,
      with the number of rows that are current / stale (another ID string) / blank, and how many of the
      51 changing MAGs are in the file -> the extension list for scripts/apply_mag_ids.py;
  (b) taxonomy columns: cells equal to the row MAG's R220 name at some rank (or its full lineage), and how
      many of those cells change in R232.
Writes scan/bin_keyed_columns.tsv and scan/bin_keyed_taxonomy.tsv.
"""
import collections, csv, io, json, os, re, subprocess
csv.field_size_limit(1 << 30)
EK = os.path.expanduser("~/Documents/EmilyKin")
OUT = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/data"
SCAN = os.path.join(OUT, "scan")
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PFX = dict(zip(RANKS, ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]))
idmap = {r["MAG"]: r for r in csv.DictReader(open(os.path.join(OUT, "mag_id_map_r220_r232.tsv")), delimiter="\t")}
lin = {r["MAG"]: r for r in csv.DictReader(open(os.path.join(OUT, "mag_lineage_r220_r232.tsv")), delimiter="\t")}
bins = set(idmap)
changed_bins = {b for b, r in idmap.items() if r["changed"] == "True"}
BINRE = re.compile(r"(?<![A-Za-z0-9_])((?:CAN_[1-5]|coasm)_bin\.\d+)(?![0-9])")

files = subprocess.run(["git", "-C", EK, "grep", "-I", "-l", "-E", r"(CAN_[1-5]|coasm)_bin\.[0-9]", "--", "*.csv", "*.tsv", "*.txt"],
                       capture_output=True).stdout.decode().split("\n")
files = [f for f in files if f]
idrows, taxrows = [], []
for f in files:
    p = os.path.join(EK, f)
    if os.path.getsize(p) > 300_000_000: continue
    text = open(p, encoding="utf8", errors="replace").read()
    first = text.split("\n", 1)[0]
    delim = "\t" if first.count("\t") >= max(1, first.count(",")) else ","
    rows = list(csv.reader(io.StringIO(text), delimiter=delim))
    if len(rows) < 2: continue
    hdr = rows[0]; ncol = max(len(r) for r in rows)
    # bin column: column whose cells are mostly exact bin names
    bc = collections.Counter()
    for r in rows[1:]:
        for i, c in enumerate(r):
            if c.strip() in bins: bc[i] += 1
    if not bc: continue
    bcol, nb = bc.most_common(1)[0]
    if nb < 1: continue
    stats = collections.defaultdict(lambda: collections.Counter())
    tax = collections.defaultdict(lambda: [0, 0, set()])
    rows_with_bin = 0; changed_in_file = set(); nonblank = collections.Counter()
    for r in rows[1:]:
        if bcol >= len(r): continue
        b = r[bcol].strip()
        if b not in bins: continue
        rows_with_bin += 1
        if b in changed_bins: changed_in_file.add(b)
        cur = idmap[b]["old_id"]
        L = lin[b]
        for i, c in enumerate(r):
            if i == bcol: continue
            v = c.strip()
            if not v: stats[i]["blank"] += 1; continue
            nonblank[i] += 1
            if v == cur: stats[i]["current_exact"] += 1
            elif re.search(r"(?<![A-Za-z0-9.])" + re.escape(cur) + r"(?![0-9])" + ("" if cur.endswith("_m") else r"(?!_m)"), v):
                stats[i]["current_embedded"] += 1
            elif re.fullmatch(r"[A-Za-z][A-Za-z0-9_\-]*\.\d+(?:_ma?)?", v) and v != b:
                stats[i]["other_id_like"] += 1
            else:
                stats[i]["other"] += 1
            # taxonomy
            if v.startswith("d__") and ";" in v:
                if v == L["r220_classification"]:
                    tax[(i, "lineage")][0] += 1
                    if L["r220_classification"] != L["r232_classification"]: tax[(i, "lineage")][1] += 1; tax[(i, "lineage")][2].add(b)
                continue
            for R in RANKS[1:]:
                a = L[f"r220_{R}"]
                if a and (v == a or v == PFX[R] + a):
                    tax[(i, R)][0] += 1
                    if a != L[f"r232_{R}"]: tax[(i, R)][1] += 1; tax[(i, R)][2].add(b)
                    break
    for i, s in stats.items():
        n_id = s["current_exact"] + s["current_embedded"]
        if n_id == 0 and s["other_id_like"] < max(3, 0.5 * rows_with_bin): continue
        if n_id == 0 and s["other_id_like"] and not re.search(r"iter|id|name|label|tip|mag", hdr[i] if i < len(hdr) else "", re.I):
            continue
        idrows.append(dict(file=f, delim="tab" if delim == "\t" else "comma", bin_col=hdr[bcol] if bcol < len(hdr) else bcol,
                           id_col=hdr[i] if i < len(hdr) else f"col{i}", rows_with_bin=rows_with_bin,
                           current_exact=s["current_exact"], current_embedded=s["current_embedded"],
                           stale_or_other_id=s["other_id_like"], blank=s["blank"], other=s["other"],
                           changed_mags_in_file=len(changed_in_file)))
    for (i, R), (n, nch, bs) in tax.items():
        if n < 1: continue
        taxrows.append(dict(file=f, bin_col=hdr[bcol] if bcol < len(hdr) else bcol, tax_col=hdr[i] if i < len(hdr) else f"col{i}",
                            rank=R, n_cells_r220=n, n_cells_change_r232=nch, n_mags_change=len(bs), rows_with_bin=rows_with_bin,
                            nonblank_cells=nonblank[i], gtdb_match_fraction=round(n / max(1, nonblank[i]), 3),
                            mags_change=";".join(sorted(bs))[:400]))
for name, rows in (("bin_keyed_columns.tsv", idrows), ("bin_keyed_taxonomy.tsv", taxrows)):
    with open(os.path.join(SCAN, name), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
print(len(files), "bin files;", len(idrows), "id columns;", len(taxrows), "taxonomy columns")
