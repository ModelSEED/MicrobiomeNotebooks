"""Independent bin-keyed verification of the three docs-stage workbooks against their .pre_r232.bak.xlsx backups.

Expected values are re-derived from primary sources, not from the stage's change list:
  ID columns           -> mag_iterativeID_old_to_new.json (committed R232 map)
  rank / lineage cols  -> repo gtdbtk.bac120.summary.tsv (R232); the cell format is learned from the backup against
                          HEAD's gtdbtk.bac120.summary.tsv (R220), and must reproduce the backup for 291/291 first
  S6 B..T              -> builder MAG_target_gene_filled.xlsx 'Table S6_ANI' by bin (current); HEAD parity checked
  S6 U..AA             -> Table_S6_skani_novelty.tsv by bin (current); HEAD parity checked
Every changed cell must be explained; every expected change must be present; non-bin rows and other sheets unchanged.
"""
import json, csv, io, subprocess, re, sys, collections, zipfile
import openpyxl
from openpyxl.utils import get_column_letter as L

REPO = "/home/freiburger/Documents/EmilyKin"
OUT = sys.argv[1]
m = json.load(open(f"{REPO}/mag_iterativeID_old_to_new.json"))

def summ(text):
    return {r["user_genome"]: r for r in csv.DictReader(io.StringIO(text), delimiter="\t")}
S232 = summ(open(f"{REPO}/gtdbtk.bac120.summary.tsv", newline="").read().replace("\r\n", "\n"))
S220 = summ(subprocess.run(["git", "-C", REPO, "show", "HEAD:gtdbtk.bac120.summary.tsv"], capture_output=True,
                           text=True).stdout.replace("\r\n", "\n"))
assert len(S232) == len(S220) == 291

RANKS = {"Phylum": "p", "Class": "c", "Order": "o", "Family": "f", "Genus": "g", "genus": "g", "Species": "s"}
def rank(S, b, k):
    for x in S[b]["classification"].split(";"):
        if x.startswith(k + "__"): return x[3:]
    return ""

def norm(v):
    if v is None: return ""
    if isinstance(v, float):
        return repr(round(v, 9))
    if isinstance(v, int): return repr(float(v))
    s = str(v)
    return s

def tsv_rows(text):
    return {r["MAG"]: r for r in csv.DictReader(io.StringIO(text.replace("\r\n", "\n")), delimiter="\t")}
T6_cur = tsv_rows(open(f"{REPO}/Table_S6_skani_novelty.tsv", newline="").read())
T6_head = tsv_rows(subprocess.run(["git", "-C", REPO, "show", "HEAD:Table_S6_skani_novelty.tsv"], capture_output=True,
                                  text=True).stdout)

def builder_s6(path_or_bytes):
    wb = openpyxl.load_workbook(path_or_bytes, data_only=False)
    ws = wb["Table S6_ANI"]
    hdr = [c.value for c in ws[1]]
    res = {}
    for row in ws.iter_rows(min_row=2):
        b = row[0].value
        if b in m: res[b] = {hdr[i]: row[i].value for i in range(len(hdr))}
    return res
B6_cur = builder_s6(f"{REPO}/MAG_target_gene_filled.xlsx")
B6_head = builder_s6(io.BytesIO(subprocess.run(["git", "-C", REPO, "show", "HEAD:MAG_target_gene_filled.xlsx"],
                                               capture_output=True).stdout))

def num_eq(a, b):
    """value equality allowing str<->number typing differences ('97.58' vs 97.58)."""
    if norm(a) == norm(b): return True
    try:
        return abs(float(a) - float(b)) < 1e-9
    except (TypeError, ValueError):
        return False

report = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s); report.append(s)

def sheet_hdr(ws):
    """header row = the row above the first bin row; ID column found in that row or the rows above."""
    first = None
    for row in ws.iter_rows():
        if isinstance(row[0].value, str) and row[0].value.strip() in m:
            first = row[0].row; break
    if first is None: return None, {}
    hdr = {}
    for r in range(1, first):
        for c in range(1, ws.max_column + 1):
            v = ws.cell(r, c).value
            if isinstance(v, str) and v.strip(): hdr[c] = v  # lowest header row wins
    return first, hdr

def lineage_fmt_learn(bak_vals, col):
    """find f with f(R220 classification) == backup cell for all bins"""
    cands = {
        "full": lambda S, b: S[b]["classification"],
        "full_space": lambda S, b: S[b]["classification"].replace(";", "; "),
    }
    for name, f in cands.items():
        if all(norm(bak_vals[b]) == norm(f(S220, b)) for b in bak_vals):
            return name, f
    bad = [(b, bak_vals[b], S220[b]["classification"]) for b in bak_vals if bak_vals[b] != S220[b]["classification"]][:3]
    return None, bad

FILES = [("MAG_target_gene_filled_UPDATED (1).xlsx", "MAG_target_gene_filled_UPDATED (1).pre_r232.bak.xlsx"),
         ("cladeexcel.xlsx", "cladeexcel.pre_r232.bak.xlsx"),
         ("abundance and coverage (1).xlsx", "abundance and coverage (1).pre_r232.bak.xlsx")]
summary = collections.Counter()
all_rows = []
for cur_f, bak_f in FILES:
    wc = openpyxl.load_workbook(f"{REPO}/{cur_f}", data_only=False)
    wbk = openpyxl.load_workbook(f"{REPO}/{bak_f}", data_only=False)
    log(f"##### {cur_f}")
    log(f"sheets: {wc.sheetnames == wbk.sheetnames} {wc.sheetnames}")
    log(f"defined names equal: {sorted(wc.defined_names.keys()) == sorted(wbk.defined_names.keys())}")
    n_changed_file = 0
    for name in wbk.sheetnames:
        a, c = wbk[name], wc[name]
        # structure
        st = {
            "dims": (a.dimensions, c.dimensions),
            "merged": (sorted(map(str, a.merged_cells.ranges)), sorted(map(str, c.merged_cells.ranges))),
            "freeze": (a.freeze_panes, c.freeze_panes),
            "autofilter": (a.auto_filter.ref, c.auto_filter.ref),
            "condfmt": (len(a.conditional_formatting), len(c.conditional_formatting)),
            "dv": (len(a.data_validations.dataValidation), len(c.data_validations.dataValidation)),
            "colwidths": ({k: v.width for k, v in a.column_dimensions.items()},
                          {k: v.width for k, v in c.column_dimensions.items()}),
            "rowheights": ({k: v.height for k, v in a.row_dimensions.items() if v.height},
                           {k: v.height for k, v in c.row_dimensions.items() if v.height}),
        }
        diffs = [k for k, (x, y) in st.items() if x != y]
        # cells: values + styles + formulas
        maxr, maxc = max(a.max_row, c.max_row), max(a.max_column, c.max_column)
        first, hdr = sheet_hdr(a)
        changed = []
        styled = 0
        fa = fc = 0
        for r in range(1, maxr + 1):
            for col in range(1, maxc + 1):
                x, y = a.cell(r, col), c.cell(r, col)
                if isinstance(x.value, str) and x.value.startswith("="): fa += 1
                if isinstance(y.value, str) and y.value.startswith("="): fc += 1
                if x.value != y.value or type(x.value) != type(y.value):
                    changed.append((r, col, x.value, y.value))
                if x.has_style != y.has_style or (x.has_style and (x._style != y._style)):
                    styled += 1
        log(f"[{name}] changed cells {len(changed)}; style diffs {styled}; formulas {fa}->{fc}; "
            f"structure diffs {diffs or 'none'}")
        n_changed_file += len(changed)
        if not changed and first is None: continue
        # expected values per bin row
        binrow = {}
        if first:
            for r in range(first, a.max_row + 1):
                v = a.cell(r, 1).value
                if isinstance(v, str) and v.strip() in m: binrow[r] = v.strip()
        exp = {}  # (r,col) -> expected value
        why = {}
        for col, h in hdr.items():
            hs = h.strip()
            if hs in ("MAG iteration ID", "iterative_ID", "iterativeID", "iterative ID"):
                for r, b in binrow.items(): exp[(r, col)] = m[b]; why[(r, col)] = "ID"
            elif hs in RANKS and name != "Table S6_ANI":
                k = RANKS[hs]
                # learn blank convention from backup vs R220
                okfmt = all(norm(a.cell(r, col).value) == norm(rank(S220, b, k)) for r, b in binrow.items())
                if not okfmt:
                    bad = [(b, a.cell(r, col).value, rank(S220, b, k)) for r, b in binrow.items()
                           if norm(a.cell(r, col).value) != norm(rank(S220, b, k))][:4]
                    log(f"   !! {name} col {L(col)} {hs}: backup != f(R220) for some bins: {bad}")
                for r, b in binrow.items():
                    exp[(r, col)] = rank(S232, b, k) or None; why[(r, col)] = hs
            elif hs in ("Taxonomy", "Taxonomy ", "GTDB taxonomy") and name != "Table S6_ANI":
                vals = {b: a.cell(r, col).value for r, b in binrow.items()}
                fname, f = lineage_fmt_learn(vals, col)
                if fname is None: log(f"   !! {name} col {L(col)} lineage format not learned: {f}")
                else:
                    for r, b in binrow.items(): exp[(r, col)] = f(S232, b); why[(r, col)] = "lineage"
        if name == "Table S6_ANI":
            # HEAD parity: backup == HEAD builder (A..T) + HEAD TSV (U..AA)
            par_bad = 0
            for r, b in binrow.items():
                for col, h in hdr.items():
                    if col == 1: continue
                    if col <= 20:
                        e_head, e_cur = B6_head[b].get(h), B6_cur[b].get(h)
                    else:
                        e_head, e_cur = T6_head[b].get(h), T6_cur[b].get(h)
                    if not num_eq(a.cell(r, col).value, e_head if e_head not in ("",) else None) and \
                       not (norm(a.cell(r, col).value) == "" and norm(e_head) == ""):
                        par_bad += 1
                        if par_bad <= 3: log(f"   S6 HEAD-parity miss {b} {h}: bak={a.cell(r,col).value!r} head={e_head!r}")
                    exp[(r, col)] = e_cur; why[(r, col)] = "S6:" + ("builder" if col <= 20 else "tsv")
            log(f"   S6 HEAD parity (backup == HEAD builder B..T + HEAD TSV U..AA): {par_bad} mismatching cells")
        # classify changes
        unexplained, wrong = [], []
        for r, col, x, y in changed:
            key = (r, col)
            if key not in exp:
                unexplained.append((L(col) + str(r), a.cell(r, 1).value, hdr.get(col), x, y)); continue
            e = exp[key]
            if not (num_eq(y, e) or (norm(y) == "" and norm(e) == "")):
                wrong.append((L(col) + str(r), binrow.get(r), hdr.get(col), x, y, e))
            summary[(cur_f, name, why[key])] += 1
            all_rows.append((cur_f, name, L(col) + str(r), binrow.get(r), hdr.get(col).strip(), x, y))
        missing = []
        chset = {(r, col) for r, col, _, _ in changed}
        for key, e in exp.items():
            if key in chset: continue
            v = a.cell(*key).value
            if not (num_eq(v, e) or (norm(v) == "" and norm(e) == "")):
                missing.append((L(key[1]) + str(key[0]), binrow.get(key[0]), hdr.get(key[1]), v, e))
        log(f"   explained {len(changed) - len(unexplained) - len(wrong)}, wrong value {len(wrong)}, "
            f"unexplained {len(unexplained)}, missing expected changes {len(missing)}")
        for w in wrong[:8]: log("   WRONG", w)
        for u in unexplained[:8]: log("   UNEXPLAINED", u)
        for x in missing[:8]: log("   MISSING", x)
        # D1: ID column unique and equal to map
        for col, h in hdr.items():
            if h.strip() in ("MAG iteration ID", "iterative_ID"):
                ids = [c.cell(r, col).value for r in binrow]
                okmap = sum(1 for r, b in binrow.items() if c.cell(r, col).value == m[b])
                log(f"   D1 {L(col)}: {okmap}/{len(binrow)} equal map, duplicates {len(ids) - len(set(ids))}")
    log(f"total changed cells in {cur_f}: {n_changed_file}")
    # raw zip check: <f> count and <v> presence per sheet part
    for fname in (cur_f, bak_f):
        z = zipfile.ZipFile(f"{REPO}/{fname}")
        parts = sorted(n for n in z.namelist() if n.startswith("xl/worksheets/sheet"))
        fc = {p: len(re.findall(rb"<f[ >]", z.read(p))) for p in parts}
        fv = {p: len(re.findall(rb"<f[^>]*>[^<]*</f><v>", z.read(p))) for p in parts}
        log(f"   raw {fname[:28]}: <f> per part {[(p.split('/')[-1], n) for p, n in fc.items() if n]}, "
            f"<f>+cached<v> {sum(fv.values())}")
    za, zb = zipfile.ZipFile(f"{REPO}/{cur_f}"), zipfile.ZipFile(f"{REPO}/{bak_f}")
    log(f"   zip parts equal: {sorted(za.namelist()) == sorted(zb.namelist())}; "
        f"differing parts: {[n for n in zb.namelist() if n in za.namelist() and za.read(n) != zb.read(n)]}")

log("SUMMARY of changed cells by (file, sheet, reason):")
for k, v in sorted(summary.items()): log(f"   {v:5d} {k}")
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["file", "sheet", "cell", "bin", "column", "before", "after"])
    w.writerows(all_rows)
