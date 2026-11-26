#!/usr/bin/env python
"""READ-ONLY scan of every tracked file of ~/Documents/EmilyKin for the R220 -> R232 MAG rename.

For each tracked text file (csv/tsv/json/ipynb/md/txt/py/R/sh/log/html/svg/fa/... and .gz of text) it finds
ID-like tokens and classifies each token against the sets built by build_maps.py:
  bin            one of the 291 bin names (CAN_1_bin.107, coasm_bin.185, ...)
  old_minted     old ID of a changed MAG ending in _m (30): MAG-only string, always MAG context
  old_inherit    old ID of a changed MAG that is also an ASV iterativeID (21): context must be decided
  old_display    'Ca. Competibacter.32' / 'Candidatus Falkowbacteria.7' / '*Ca.* ...' display forms of the two Ca_ old IDs
  new_asv        new ID that already exists as an ASV iterativeID (Simplicispira.1, Holophagaceae.9, Taibaiella.5,
                 Hoeflea.1): today these occurrences are ASV context
  new_other      any other new ID already present (should be none; Chitinophagaceae.1_m is counted as old_minted)
  mag_only       unchanged MAG ID that is not an ASV ID (column signal: MAG)
  shared         unchanged MAG ID that is also an ASV ID (no signal)
  asv_only       ASV iterativeID that is no MAG ID (column signal: ASV)
Each occurrence of an old_inherit token gets a context from its column (tabular / JSON path) or line (text):
MAG / ASV / ambiguous, using the other tokens of that column / line (see decide()).
Also: per-MAG taxonomy cells (R220 names at the MAG's own rank, keyed by the row's bin or MAG ID) and R220-only
taxon-name occurrences in text. Writes rename_plan/data/scan/*.tsv. Never writes to the repo.
"""
import collections, csv, gzip, io, json, os, re, subprocess, sys, zipfile

csv.field_size_limit(1 << 30)
EK = os.path.expanduser("~/Documents/EmilyKin")
OUT = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/data"
SCAN = os.path.join(OUT, "scan"); os.makedirs(SCAN, exist_ok=True)
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PFX = dict(zip(RANKS, ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]))

# ---------------- sets -----------------------------------------------------------------------------------------
idmap = list(csv.DictReader(open(os.path.join(OUT, "mag_id_map_r220_r232.tsv")), delimiter="\t"))
lin = {r["MAG"]: r for r in csv.DictReader(open(os.path.join(OUT, "mag_lineage_r220_r232.tsv")), delimiter="\t")}
bins = {r["MAG"] for r in idmap}
old_of = {r["MAG"]: r["old_id"] for r in idmap}
new_of = {r["MAG"]: r["new_id"] for r in idmap}
asv_ids = set(json.load(open(os.path.join(EK, "iterativeIDs.json"))))
cur_mag_ids = set(old_of.values())
mag_only_ids = cur_mag_ids - asv_ids
shared_ids = cur_mag_ids & asv_ids
asv_only_ids = asv_ids - cur_mag_ids
changed = {r["old_id"]: r["MAG"] for r in idmap if r["changed"] == "True"}
old_minted = {o for o in changed if o.endswith("_m")}
old_inherit = set(changed) - old_minted
new_ids = {r["new_id"]: r["MAG"] for r in idmap if r["changed"] == "True"}
new_asv = {n for n in new_ids if n in asv_ids}
id2bin = {v: k for k, v in old_of.items()}          # MAG interpretation of a current ID
display_roots = {o[3:]: o for o in changed if o.startswith("Ca_")}   # Competibacter.32 -> Ca_Competibacter.32
assert not (old_inherit - asv_ids)

def kind(t):
    if t in bins: return "bin"
    if t in old_minted: return "old_minted"
    if t in old_inherit: return "old_inherit"
    if t in new_asv: return "new_asv"
    if t in new_ids: return "new_other"
    if t in mag_only_ids: return "mag_only"
    if t in shared_ids: return "shared"
    if t in asv_only_ids: return "asv_only"
    return None

TOK = re.compile(r"(?<![A-Za-z0-9_\-.])([A-Za-z][A-Za-z0-9_\-]*(?:_sp\.)?\.\d+(?:_ma?(?![A-Za-z0-9]))?)(?![0-9])")
DISP = re.compile(r"(?:Ca\.?[ _]?|\*Ca\.\*\s?|_Ca\._\s?|Candidatus\s|Ca\s)$")

def tokens(s):
    """yield (token, kind) for every ID-like token in s; resolves embedded names (emilykin_..._b1_Ca_X.7)."""
    for mt in TOK.finditer(s):
        t = mt.group(1)
        k = kind(t)
        if k is None and "_" in t:
            parts = t.split("_")
            for i in range(1, len(parts)):
                tt = "_".join(parts[i:])
                k2 = kind(tt)
                if k2:
                    t, k = tt, k2; break
        if k is None and t in display_roots and DISP.search(s[max(0, mt.start() - 12):mt.start()]):
            yield display_roots[t], "old_display"; continue
        if k:
            yield t, k

MAGHDR = re.compile(r"mag|bin|genome|user_genome|tip|best_mag", re.I)
ASVHDR = re.compile(r"asv|member|otu|feature|hash|amplicon", re.I)

def decide(col_counts, header=""):
    """context for old_inherit occurrences of one column / path from the other tokens in it."""
    mag_sig = col_counts["mag_only"] + col_counts["old_minted"] + col_counts["new_other"]
    asv_sig = col_counts["asv_only"] + col_counts["new_asv"]
    if mag_sig and not asv_sig: return "MAG", "column holds MAG-only IDs, no ASV-only IDs"
    if asv_sig and not mag_sig: return "ASV", "column holds ASV-only IDs, no MAG-only IDs"
    if mag_sig and asv_sig: return "ambiguous", f"mixed column (MAG-only {mag_sig}, ASV-only {asv_sig})"
    h = str(header)
    if MAGHDR.search(h) and not ASVHDR.search(h): return "MAG", f"only shared IDs; header '{h[:40]}' says MAG"
    if ASVHDR.search(h) and not MAGHDR.search(h): return "ASV", f"only shared IDs; header '{h[:40]}' says ASV"
    return "ambiguous", "only shared IDs, header uninformative"

# R220 taxon names that change (for text / code scan)
tx = list(csv.DictReader(open(os.path.join(OUT, "taxon_name_changes_r220_r232.tsv")), delimiter="\t"))
r220_changed = {}
for r in tx:
    if r["r220_name"]:
        r220_changed.setdefault(r["r220_name"], set()).add((r["rank"], r["r232_name"], r["r220_name_is_asv_taxon_name"] == "True",
                                                           r["r220_name_still_used_in_r232_mag_taxonomy"] == "True"))
names_sorted = sorted(r220_changed, key=len, reverse=True)
NAMERE = re.compile(r"(?<![A-Za-z0-9])(" + "|".join(re.escape(n) for n in names_sorted) + r")(?![A-Za-z0-9])")

# ---------------- per-file scanning -----------------------------------------------------------------------------
col_rows = []      # file, column, kinds..., context counts
tax_rows = []      # file, column, rank, n_mag_rows, n_r220_match, n_need_change
name_rows = []     # file, name, rank, n, n_with_prefix, context
file_rows = {}

def newcounts(): return collections.Counter()

def add_col(path, col, cnt, ctx, why, example, binjoin=0, header=""):
    col_rows.append(dict(file=path, column=col, header=str(header)[:60],
                         old_minted=cnt["old_minted"], old_inherit=cnt["old_inherit"], old_display=cnt["old_display"],
                         new_asv_preexisting=cnt["new_asv"], new_other_preexisting=cnt["new_other"],
                         bins=cnt["bin"], mag_only=cnt["mag_only"], shared=cnt["shared"], asv_only=cnt["asv_only"],
                         inherit_context=ctx, inherit_context_basis=why, inherit_bin_joined=binjoin,
                         old_ids_found=";".join(sorted(cnt["_ids"])) if cnt.get("_ids") else "",
                         example=example[:160]))

def mag_of_row(cells_tokens, mag_cols):
    bs = [t for tok in cells_tokens for (t, k) in tok if k == "bin"]
    if bs: return bs[0]
    for ci in mag_cols:
        if ci < len(cells_tokens):
            for t, k in cells_tokens[ci]:
                if k in ("mag_only", "old_minted", "shared", "old_inherit") and t in id2bin:
                    return id2bin[t]
    return None

def tax_check(mag, value, acc, col):
    """compare a cell to the MAG's R220 taxonomy; acc[col] counts"""
    v = value.strip()
    if not v or mag not in lin: return
    L = lin[mag]
    if v.startswith("d__") and ";" in v:
        if v == L["r220_classification"]:
            acc[(col, "lineage")][0] += 1; acc[(col, "lineage")][1] += L["r220_classification"] != L["r232_classification"]
        return
    for R in RANKS[1:]:
        a = L[f"r220_{R}"]
        if a and (v == a or v == PFX[R] + a):
            acc[(col, R)][0] += 1; acc[(col, R)][1] += a != L[f"r232_{R}"]
            return

def scan_table(path, text, delim):
    rows = list(csv.reader(io.StringIO(text), delimiter=delim))
    if not rows: return
    hdr = rows[0]
    ncol = max(len(r) for r in rows)
    colcnt = [newcounts() for _ in range(ncol)]
    hdrcnt = newcounts()
    ex = [""] * ncol
    rowtoks = []
    for r in rows[1:]:
        rt = [list(tokens(c)) if c else [] for c in r]
        rowtoks.append(rt)
        for ci, tl in enumerate(rt):
            for t, k in tl:
                colcnt[ci][k] += 1
                if k.startswith("old"):
                    colcnt[ci].setdefault("_ids", set()).add(t)
                    if not ex[ci]: ex[ci] = r[ci]
    ht = [list(tokens(c)) for c in hdr]
    for tl in ht:
        for t, k in tl:
            hdrcnt[k] += 1
            if k.startswith("old"): hdrcnt.setdefault("_ids", set()).add(t)
    # MAG-key columns: columns whose decided context is MAG or which hold bins
    ctxs = []
    for ci in range(ncol):
        h = hdr[ci] if ci < len(hdr) else f"col{ci}"
        ctx, why = decide(colcnt[ci], h)
        ctxs.append((ctx, why))
    mag_cols = [ci for ci in range(ncol) if ctxs[ci][0] == "MAG"]
    # bin-joined inherit occurrences (row's bin owns that old ID) -- supporting evidence only
    binjoin = collections.Counter()
    acc = collections.defaultdict(lambda: [0, 0])
    for r, rt in zip(rows[1:], rowtoks):
        rowbins = {t for tl in rt for (t, k) in tl if k == "bin"}
        for ci, tl in enumerate(rt):
            for t, k in tl:
                if k == "old_inherit" and changed[t] in rowbins: binjoin[ci] += 1
        m = mag_of_row(rt, mag_cols)
        if m:
            for ci, c in enumerate(r):
                tax_check(m, c, acc, hdr[ci] if ci < len(hdr) else f"col{ci}")
    for ci in range(ncol):
        c = colcnt[ci]
        if any(c[k] for k in ("old_minted", "old_inherit", "old_display", "new_asv", "new_other")):
            h = hdr[ci] if ci < len(hdr) else f"col{ci}"
            add_col(path, f"col[{ci}]={h}", c, ctxs[ci][0] if c["old_inherit"] or c["old_display"] else "",
                    ctxs[ci][1] if c["old_inherit"] or c["old_display"] else "", ex[ci], binjoin[ci], h)
    if any(hdrcnt[k] for k in ("old_minted", "old_inherit", "old_display", "new_asv", "new_other")):
        ctx, why = decide(hdrcnt, "<header row>")
        add_col(path, "<header row>", hdrcnt, ctx if hdrcnt["old_inherit"] else "", why if hdrcnt["old_inherit"] else "",
                ",".join(sorted(hdrcnt.get("_ids", [])))[:160])
    for (col, R), (n, nch) in acc.items():
        tax_rows.append(dict(file=path, column=col, rank=R, n_mag_cells=n, n_need_change=nch))
    file_rows[path]["n_bin_tokens"] = sum(c["bin"] for c in colcnt) + hdrcnt["bin"]
    file_rows[path]["format"] = "table"

def walk_json(obj, path, out, mag=None):
    """flatten JSON to (normpath, string, owning_mag) for keys and string values."""
    if isinstance(obj, dict):
        bm = None
        for k, v in obj.items():
            if isinstance(v, str) and v in bins and re.search(r"bin|mag|genome", k, re.I): bm = v
        m2 = bm or mag
        for k, v in obj.items():
            ks = str(k)
            kk = kind(ks) if TOK.fullmatch(ks) else None
            norm = "<BIN>" if kk == "bin" else ("<ID>" if kk else ks)
            out.append((path + "/{keys}", ks, m2, True))
            km = ks if kk == "bin" else (id2bin.get(ks) if kk in ("mag_only", "old_minted") else None)
            walk_json(v, path + "/" + norm, out, km or m2)
    elif isinstance(obj, list):
        for v in obj:
            walk_json(v, path + "[]", out, mag)
    elif isinstance(obj, str):
        out.append((path, obj, mag, False))
    elif obj is not None:
        out.append((path, str(obj), mag, False))

def scan_json(path, obj):
    flat = []
    walk_json(obj, "", flat)
    colcnt = collections.defaultdict(newcounts)
    ex = {}
    acc = collections.defaultdict(lambda: [0, 0])
    nb = 0
    for p, s, m, iskey in flat:
        for t, k in tokens(s):
            colcnt[p][k] += 1
            nb += k == "bin"
            if k.startswith("old"):
                colcnt[p].setdefault("_ids", set()).add(t)
                ex.setdefault(p, s)
        if m and not iskey:
            tax_check(m, s, acc, p)
    for p, c in colcnt.items():
        if any(c[k] for k in ("old_minted", "old_inherit", "old_display", "new_asv", "new_other")):
            ctx, why = decide(c, p) if (c["old_inherit"] or c["old_display"]) else ("", "")
            add_col(path, p, c, ctx, why, ex.get(p, ""), 0, p)
    for (col, R), (n, nch) in acc.items():
        tax_rows.append(dict(file=path, column=col, rank=R, n_mag_cells=n, n_need_change=nch))
    file_rows[path]["n_bin_tokens"] = nb
    file_rows[path]["format"] = "json"

LINE_MAG = re.compile(r"\bMAGs?\b|\bbins?\b|GTDB|mag_iterativeID", re.I)
LINE_ASV = re.compile(r"\bASVs?\b|\bmembers?\b|16S|amplicon|MiDAS", re.I)

def scan_text(path, text, colname="<text>"):
    c = newcounts(); ctxc = collections.Counter(); ex = ""
    nb = 0
    for line in text.splitlines():
        tl = list(tokens(line))
        if not tl: continue
        ks = collections.Counter(k for _, k in tl)
        nb += ks["bin"]
        for t, k in tl:
            c[k] += 1
            if k.startswith("old"): c.setdefault("_ids", set()).add(t)
        if ks["old_inherit"] or ks["old_display"]:
            mag_sig = ks["bin"] + ks["mag_only"] + ks["old_minted"] + bool(LINE_MAG.search(line))
            asv_sig = ks["asv_only"] + ks["new_asv"] + bool(LINE_ASV.search(line))
            ctx = "MAG" if mag_sig and not asv_sig else "ASV" if asv_sig and not mag_sig else "ambiguous"
            ctxc[ctx] += ks["old_inherit"] + ks["old_display"]
        if (ks["old_minted"] or ks["old_inherit"] or ks["old_display"]) and not ex:
            ex = line.strip()
    if any(c[k] for k in ("old_minted", "old_inherit", "old_display", "new_asv", "new_other")):
        ctx = ";".join(f"{k}:{v}" for k, v in sorted(ctxc.items()))
        add_col(path, colname, c, ctx, "per-line: bins / MAG-only IDs / 'MAG' vs ASV-only IDs / 'ASV','member'" if ctx else "", ex)
    # taxon names
    if colname in ("<text>", "source"):
        nc = collections.Counter(); npfx = collections.Counter()
        for mt in NAMERE.finditer(text):
            n = mt.group(1); nc[n] += 1
            if text[max(0, mt.start() - 3):mt.start()] in ("p__", "c__", "o__", "f__", "g__", "s__"): npfx[n] += 1
        for n, v in nc.items():
            name_rows.append(dict(file=path, part=colname, r220_name=n, n=v, n_gtdb_prefixed=npfx[n],
                                  ranks=";".join(sorted({f"{a}->{b or '(none)'}" for a, b, _, _ in r220_changed[n]})),
                                  also_asv_taxon=any(x[2] for x in r220_changed[n]),
                                  still_used_r232=any(x[3] for x in r220_changed[n])))
    return nb

TEXT_EXT = {".csv", ".tsv", ".json", ".ipynb", ".md", ".txt", ".py", ".r", ".sh", ".log", ".html", ".svg", ".yaml", ".yml",
            ".out", ".err", ".text", ".fa", ".faa", ".fasta", ".ffn", ".nwk", ".tree", ".treefile", ".contree", ".nex",
            ".support", ".besttree", ".bootstraps", ".mltrees", ".starttree", ".bestmodel", ".iqtree", ".bionj", ".mldist",
            ".phy", ".bed", ".tex", ".bib", ".js", ".css", ".xml", ".sbml", ".patch", ".rmd", ".qmd", ".toml", ".cfg", ".ini",
            ".env", ".done", ".gitignore", ".lock", ".sha256", ".hmm", ".rba", ".dat"}
ZIP_EXT = {".xlsx", ".docx", ".pptx"}
BIN_EXT = {".png", ".pdf", ".jpg", ".jpeg", ".gif", ".tif", ".tiff", ".qza", ".qzv", ".npy", ".npz", ".pkl", ".h5", ".rds"}
bin_rows = []

def classify_ext(p):
    b = os.path.basename(p).lower()
    e = os.path.splitext(b)[1]
    if e == ".gz": return "gz", os.path.splitext(b[:-3])[1]
    if "raxml." in b: return "text", ".tree"
    if e in TEXT_EXT or e == "": return "text", e
    if e in ZIP_EXT: return "zip", e
    return "binary", e

def scan_file(p):
    full = os.path.join(EK, p)
    kindf, e = classify_ext(p)
    try: size = os.path.getsize(full)
    except OSError: return
    file_rows[p] = dict(file=p, ext=e, size=size, format="", n_bin_tokens=0)
    if kindf == "zip":
        try:
            z = zipfile.ZipFile(full)
            txt = "\n".join(re.sub(r"<[^>]+>", " ", z.read(n).decode("utf8", "ignore")) for n in z.namelist()
                            if n.endswith(".xml") and ("sheet" in n or "sharedStrings" in n or "word/document" in n
                                                       or "slides/slide" in n or "footnotes" in n or "comments" in n
                                                       or "embeddings" in n or "chart" in n))
        except Exception as ex:
            txt = ""; file_rows[p]["format"] = f"zip-error {ex}"
        c = collections.Counter(k for _, k in tokens(txt))
        nm = collections.Counter(m.group(1) for m in NAMERE.finditer(txt))
        bin_rows.append(dict(file=p, ext=e, size=size, old_minted=c["old_minted"], old_inherit=c["old_inherit"],
                             old_display=c["old_display"], bins=c["bin"], new_asv_preexisting=c["new_asv"],
                             r220_names=";".join(f"{k}:{v}" for k, v in nm.most_common(12)), note="text extracted from OOXML"))
        file_rows[p]["format"] = "ooxml"
        return
    if kindf == "binary":
        note = ""
        c = collections.Counter(); nm = collections.Counter()
        if e == ".pdf":
            try:
                txt = subprocess.run(["pdftotext", "-q", full, "-"], capture_output=True, text=True, timeout=120).stdout
                c = collections.Counter(k for _, k in tokens(txt)); nm = collections.Counter(m.group(1) for m in NAMERE.finditer(txt))
                note = "pdftotext"
            except Exception as ex:
                note = f"pdftotext failed {ex}"
        bin_rows.append(dict(file=p, ext=e, size=size, old_minted=c["old_minted"], old_inherit=c["old_inherit"],
                             old_display=c["old_display"], bins=c["bin"], new_asv_preexisting=c["new_asv"],
                             r220_names=";".join(f"{k}:{v}" for k, v in nm.most_common(12)), note=note or "image/binary: not text-searchable"))
        file_rows[p]["format"] = "binary"
        return
    if size > 400_000_000: return
    raw = open(full, "rb").read()
    if kindf == "gz":
        try: raw = gzip.decompress(raw)
        except Exception: return
    if b"\0" in raw[:8192]:
        file_rows[p]["format"] = "binary?"; return
    text = raw.decode("utf8", "replace")
    if e in (".csv", ".tsv"):
        delim = "\t" if e == ".tsv" else ","
        first = text[:5000].split("\n", 1)[0]
        if e == ".csv" and first.count("\t") > first.count(","): delim = "\t"
        scan_table(p, text, delim)
    elif e == ".json":
        try:
            scan_json(p, json.loads(text))
        except Exception:
            file_rows[p]["format"] = "text(json-fail)"
            file_rows[p]["n_bin_tokens"] = scan_text(p, text)
    elif e == ".ipynb":
        try:
            nb = json.loads(text)
            src = "\n".join("".join(c.get("source", [])) for c in nb.get("cells", []))
            outs = []
            for c in nb.get("cells", []):
                for o in c.get("outputs", []):
                    outs.append("".join(o.get("text", [])) if "text" in o else
                                "\n".join("".join(v) if isinstance(v, list) else str(v)
                                          for kk, v in o.get("data", {}).items() if kk.startswith("text/")))
            n1 = scan_text(p, src, "source"); n2 = scan_text(p, "\n".join(outs), "outputs")
            file_rows[p]["n_bin_tokens"] = n1 + n2; file_rows[p]["format"] = "ipynb"
        except Exception:
            file_rows[p]["n_bin_tokens"] = scan_text(p, text); file_rows[p]["format"] = "text(ipynb-fail)"
    else:
        file_rows[p]["n_bin_tokens"] = scan_text(p, text); file_rows[p]["format"] = "text"

def main():
    files = subprocess.run(["git", "-C", EK, "ls-files", "-z"], capture_output=True).stdout.decode().split("\0")
    files = [f for f in files if f]
    # candidate text files: any ID token kind of interest / bin / R220 name (prefilter with git grep for speed)
    pat = os.path.join(SCAN, "_prefilter_patterns.txt")
    with open(pat, "w") as fh:
        for s in sorted(set(changed) | set(new_ids) | set(display_roots) | set(r220_changed)):
            fh.write(s + "\n")
        for s in ("coasm_bin.", "CAN_1_bin.", "CAN_2_bin.", "CAN_3_bin.", "CAN_4_bin.", "CAN_5_bin."):
            fh.write(s + "\n")
    cand = set(subprocess.run(["git", "-C", EK, "grep", "-I", "-l", "-F", "-f", pat], capture_output=True).stdout.decode().split("\n"))
    cand.discard("")
    todo = []
    for f in files:
        k, e = classify_ext(f)
        if k in ("zip", "binary", "gz") or f in cand:
            todo.append(f)
    print(f"tracked {len(files)}, candidates {len(todo)}", file=sys.stderr)
    for i, f in enumerate(todo):
        try:
            scan_file(f)
        except Exception as ex:
            print("ERR", f, ex, file=sys.stderr)
        if i % 500 == 0: print(i, f, file=sys.stderr)
    def w(name, rows, fields=None):
        if not rows: rows = [{}]
        fields = fields or list(rows[0])
        with open(os.path.join(SCAN, name), "w") as fh:
            ww = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            ww.writeheader(); ww.writerows(r for r in rows if r)
    w("columns_raw.tsv", col_rows)
    w("taxonomy_cells_raw.tsv", tax_rows)
    w("r220_names_raw.tsv", name_rows)
    w("binary_raw.tsv", bin_rows)
    w("files_raw.tsv", list(file_rows.values()), ["file", "ext", "size", "format", "n_bin_tokens"])
    print(f"columns {len(col_rows)} taxonomy {len(tax_rows)} names {len(name_rows)} binary {len(bin_rows)}", file=sys.stderr)

if __name__ == "__main__":
    main()
