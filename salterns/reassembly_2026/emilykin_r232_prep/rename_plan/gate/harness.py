#!/usr/bin/env python3
"""Baseline gate harness (critic C15 / correction 17).

Runs one table builder or figure renderer inside the HEAD worktree (9e087c3), detects every file it wrote,
compares each against the committed version (HEAD blob) or, for untracked outputs such as the gitignored PDFs,
against the main repo's untracked copy, and then restores the worktree to HEAD so the next step starts clean.

Usage:
  harness.py run STEP [--cwd REL] [--env K=V ...] [--strace] [--keep] -- CMD ...
  harness.py restore              # restore the worktree to HEAD (tracked files) and drop new files

Verdict tiers per written file:
  A  byte-identical
  B  equivalent: pixel-identical PNG, text-identical PDF, cell-identical XLSX, text/media-identical DOCX,
     parsed-equal JSON, numerically identical CSV/TSV (|rel diff| <= 1e-9), or EOL-only difference
  C  differs (how is reported)
  N  new file with no committed counterpart
Nothing here writes to ~/Documents/EmilyKin.
"""
import argparse, csv, difflib, hashlib, io, json, math, os, re, shutil, subprocess, sys, time, zipfile
from pathlib import Path

RP = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan")
WT = RP / "head_worktree"
MAIN = Path.home() / "Documents" / "EmilyKin"
SCR = RP / "gate_scratch"
HEAD = "9e087c3"
PY = str(Path.home() / "Documents" / "py_venv" / "bin" / "python")
SKIP_DIRS = {".git", "__pycache__", ".ipynb_checkpoints"}
INPUT_LIST = SCR / "copied_inputs.tsv"       # relpaths copied from the main repo into the worktree
csv.field_size_limit(sys.maxsize)


def tracked_set():
    out = subprocess.run(["git", "-C", str(WT), "ls-files", "-z"], capture_output=True, check=True).stdout
    return set(out.decode().split("\0")) - {""}


def snapshot():
    snap = {}
    for root, dirs, files in os.walk(WT):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for f in files:
            p = os.path.join(root, f)
            try:
                st = os.lstat(p)
            except FileNotFoundError:
                continue
            snap[os.path.relpath(p, WT)] = (st.st_mtime_ns, st.st_size)
    return snap


def head_blob(rel):
    r = subprocess.run(["git", "-C", str(WT), "show", f"{HEAD}:{rel}"], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def copied_inputs():
    if not INPUT_LIST.exists():
        return set()
    return {l.split("\t")[0] for l in INPUT_LIST.read_text().splitlines() if l and not l.startswith("#")}


# ---------------------------------------------------------------- comparisons
def _num(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def _close(a, b):
    if a == b:
        return True
    fa, fb = _num(a), _num(b)
    if fa is None or fb is None:
        return False
    return abs(fa - fb) <= 1e-9 * max(1.0, abs(fa), abs(fb))


def cmp_table(old, new, delim):
    t_old = old.decode("utf-8", "replace").replace("\r\n", "\n")
    t_new = new.decode("utf-8", "replace").replace("\r\n", "\n")
    if t_old == t_new:
        return "B", "EOL-only difference (CRLF vs LF)"
    ro = list(csv.reader(io.StringIO(t_old), delimiter=delim))
    rn = list(csv.reader(io.StringIO(t_new), delimiter=delim))
    if not ro or not rn:
        return "C", f"empty table (old {len(ro)} rows, new {len(rn)} rows)"
    ho, hn = ro[0], rn[0]
    notes = []
    if ho != hn:
        notes.append(f"header differs: -{[c for c in ho if c not in hn][:8]} +{[c for c in hn if c not in ho][:8]}"
                     + ("" if set(ho) != set(hn) else " (column order only)"))
    if len(ro) != len(rn):
        notes.append(f"rows {len(ro)-1} -> {len(rn)-1}")
    # align by column name, rows by order (and also try by first column key)
    common = [c for c in ho if c in hn]
    io_, in_ = {c: i for i, c in enumerate(ho)}, {c: i for i, c in enumerate(hn)}
    def cell(r, i):
        return r[i] if i < len(r) else ""
    ndiff = nnum = 0; cols = {}; ex = []
    key_old = {cell(r, 0): r for r in ro[1:]}
    key_new = {cell(r, 0): r for r in rn[1:]}
    by_key = len(key_old) == len(ro) - 1 and len(key_new) == len(rn) - 1 and set(key_old) == set(key_new)
    pairs = ([(key_old[k], key_new[k]) for k in key_old] if by_key else list(zip(ro[1:], rn[1:])))
    order_same = [cell(r, 0) for r in ro[1:]] == [cell(r, 0) for r in rn[1:]]
    for a, b in pairs:
        for c in common:
            va, vb = cell(a, io_[c]), cell(b, in_[c])
            if va != vb:
                if _close(va, vb):
                    nnum += 1
                else:
                    ndiff += 1; cols[c] = cols.get(c, 0) + 1
                    if len(ex) < 6:
                        ex.append(f"{cell(a,0)}|{c}: {va[:60]!r} -> {vb[:60]!r}")
    if by_key and not order_same:
        notes.append("row order differs (compared by first-column key)")
    if ndiff == 0 and not notes:
        return "B", f"numerically identical ({nnum} cells differ only in float repr <=1e-9)" if nnum else "parsed-identical"
    if ndiff == 0 and notes == ["row order differs (compared by first-column key)"]:
        return "C", "same cells, different row order" + (f"; {nnum} float-repr cells" if nnum else "")
    top = ", ".join(f"{c}:{n}" for c, n in sorted(cols.items(), key=lambda x: -x[1])[:10])
    return "C", "; ".join(notes + [f"{ndiff} cells differ in {len(cols)} columns ({top})"] + (["e.g. " + " ; ".join(ex)] if ex else []))


def _jdiff(a, b, path, out, lim=12):
    if len(out) >= lim:
        return
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a.keys() | b.keys():
            if k not in a: out.append(f"+{path}/{k}")
            elif k not in b: out.append(f"-{path}/{k}")
            else: _jdiff(a[k], b[k], f"{path}/{k}", out, lim)
            if len(out) >= lim: return
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: list len {len(a)} -> {len(b)}"); return
        for i, (x, y) in enumerate(zip(a, b)):
            _jdiff(x, y, f"{path}[{i}]", out, lim)
            if len(out) >= lim: return
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not _close(a, b): out.append(f"{path}: {a} -> {b}")
    elif a != b:
        out.append(f"{path}: {str(a)[:60]!r} -> {str(b)[:60]!r}")


def cmp_json(old, new):
    try:
        a, b = json.loads(old), json.loads(new)
    except Exception as e:
        return cmp_text(old, new)
    d = []
    _jdiff(a, b, "", d)
    if not d:
        return "B", "parsed-identical (formatting/key order/float repr only)"
    return "C", f"{len(d)}+ differing paths: " + " ; ".join(d[:8])


def cmp_text(old, new):
    a = old.decode("utf-8", "replace").replace("\r\n", "\n")
    b = new.decode("utf-8", "replace").replace("\r\n", "\n")
    if a == b:
        return "B", "EOL-only difference (CRLF vs LF)"
    al, bl = a.splitlines(), b.splitlines()
    dl = [l for l in difflib.unified_diff(al, bl, lineterm="", n=0) if l[:1] in "+-" and l[:3] not in ("+++", "---")]
    return "C", f"{sum(l[0]=='-' for l in dl)} lines removed / {sum(l[0]=='+' for l in dl)} added; first: " + " | ".join(x[:110] for x in dl[:6])


def cmp_png(old, new):
    from PIL import Image
    import numpy as np
    try:
        A = np.asarray(Image.open(io.BytesIO(old)).convert("RGBA")).astype(int)
        B = np.asarray(Image.open(io.BytesIO(new)).convert("RGBA")).astype(int)
    except Exception as e:
        return "C", f"image read error {e}"
    if A.shape != B.shape:
        return "C", f"size {A.shape[1]}x{A.shape[0]} -> {B.shape[1]}x{B.shape[0]}"
    if (A == B).all():
        return "B", "pixel-identical (PNG metadata/compression differs)"
    diff = (A != B).any(axis=2)
    return "C", f"same size {A.shape[1]}x{A.shape[0]}; {diff.mean()*100:.3f}% pixels differ; mean |d| {abs(A-B).mean():.3f}"


def _pdftext(b):
    p = SCR / "_tmp_cmp.pdf"
    p.write_bytes(b)
    r = subprocess.run(["pdftotext", "-layout", str(p), "-"], capture_output=True)
    n = subprocess.run(["pdfinfo", str(p)], capture_output=True, text=True).stdout
    pages = re.search(r"Pages:\s+(\d+)", n)
    return r.stdout.decode("utf-8", "replace"), (pages.group(1) if pages else "?")


def cmp_pdf(old, new):
    ta, pa = _pdftext(old); tb, pb = _pdftext(new)
    if ta == tb and pa == pb:
        return "B", f"text-identical ({pa} pages; pdftotext -layout)"
    if " ".join(ta.split()) == " ".join(tb.split()) and pa == pb:
        return "B", f"text-identical up to whitespace/layout ({pa} pages)"
    v, how = cmp_text(ta.encode(), tb.encode())
    return "C", f"pages {pa}->{pb}; text: {how}"


def cmp_xlsx(old, new):
    import openpyxl, warnings
    warnings.filterwarnings("ignore")
    try:
        wa = openpyxl.load_workbook(io.BytesIO(old), data_only=False)
        wb = openpyxl.load_workbook(io.BytesIO(new), data_only=False)
    except Exception as e:
        return "C", f"xlsx read error {e}"
    notes = []
    if wa.sheetnames != wb.sheetnames:
        notes.append(f"sheets {wa.sheetnames} -> {wb.sheetnames}")
    tot = 0
    for sn in wa.sheetnames:
        if sn not in wb.sheetnames:
            continue
        a, b = wa[sn], wb[sn]
        mr, mc = max(a.max_row, b.max_row), max(a.max_column, b.max_column)
        nd = 0; ex = []; cols = {}
        for r in range(1, mr + 1):
            for c in range(1, mc + 1):
                va, vb = a.cell(r, c).value, b.cell(r, c).value
                if va == vb or (va in (None, "") and vb in (None, "")):
                    continue
                if isinstance(va, (int, float)) and isinstance(vb, (int, float)) and _close(va, vb):
                    continue
                nd += 1; cols[c] = cols.get(c, 0) + 1
                if len(ex) < 3:
                    ex.append(f"{openpyxl.utils.get_column_letter(c)}{r}: {str(va)[:40]!r}->{str(vb)[:40]!r}")
        if nd or a.max_row != b.max_row:
            hdrs = []
            for c, n in sorted(cols.items(), key=lambda x: -x[1])[:6]:
                h = next((a.cell(rr, c).value for rr in (3, 2, 1) if a.cell(rr, c).value), None)
                hdrs.append(f"{openpyxl.utils.get_column_letter(c)}[{str(h)[:25]}]:{n}")
            notes.append(f"'{sn}': rows {a.max_row}->{b.max_row}, {nd} cells differ ({', '.join(hdrs)}) e.g. {'; '.join(ex)}")
            tot += nd
    # drawings/images/charts survive?
    za, zb = zipfile.ZipFile(io.BytesIO(old)), zipfile.ZipFile(io.BytesIO(new))
    ma = sorted(n for n in za.namelist() if re.search(r"(media|drawings|charts)/", n))
    mb = sorted(n for n in zb.namelist() if re.search(r"(media|drawings|charts)/", n))
    if ma != mb:
        notes.append(f"drawing/media parts {len(ma)} -> {len(mb)}")
    if not notes:
        return "B", "cell-identical (all sheets, formulas/values; zip container differs)"
    return "C", " || ".join(notes)


def cmp_docx(old, new):
    import docx
    def dump(b):
        d = docx.Document(io.BytesIO(b))
        paras = [p.text for p in d.paragraphs]
        tabs = [[[c.text for c in row.cells] for row in t.rows] for t in d.tables]
        z = zipfile.ZipFile(io.BytesIO(b))
        media = {n: hashlib.md5(z.read(n)).hexdigest() for n in z.namelist() if n.startswith("word/media/")}
        return paras, tabs, media
    pa, ta, ma = dump(old); pb, tb, mb = dump(new)
    notes = []
    if pa != pb:
        v, how = cmp_text("\n".join(pa).encode(), "\n".join(pb).encode())
        notes.append(f"paragraphs {len(pa)}->{len(pb)}: {how}")
    if ta != tb:
        nd = sum(1 for x, y in zip(ta, tb) if x != y) + abs(len(ta) - len(tb))
        notes.append(f"tables {len(ta)}->{len(tb)}, {nd} differ")
    if ma != mb:
        same = sum(1 for k in ma if mb.get(k) == ma[k])
        notes.append(f"media {len(ma)}->{len(mb)} ({same} identical)")
    if not notes:
        return "B", "text/tables/media identical (docx container differs)"
    return "C", " || ".join(notes)


def compare(rel, old, new):
    if old == new:
        return "A", "byte-identical"
    ext = Path(rel).suffix.lower()
    if rel.endswith(".gz"):
        import gzip
        try:
            o, n = gzip.decompress(old), gzip.decompress(new)
            if o == n:
                return "B", "gzip payload byte-identical (gzip header/mtime differs)"
            return compare(rel[:-3], o, n)
        except Exception:
            pass
    try:
        if ext in (".csv",):
            return cmp_table(old, new, ",")
        if ext in (".tsv", ".bed"):
            return cmp_table(old, new, "\t")
        if ext == ".json":
            return cmp_json(old, new)
        if ext in (".png", ".jpg", ".jpeg"):
            return cmp_png(old, new)
        if ext == ".pdf":
            return cmp_pdf(old, new)
        if ext in (".xlsx", ".xlsm"):
            return cmp_xlsx(old, new)
        if ext == ".docx":
            return cmp_docx(old, new)
        if ext in (".md", ".txt", ".log", ".py", ".html", ".svg", ".nwk", ".tree", ".support", ".R", ".sh", ".out", ""):
            return cmp_text(old, new)
    except Exception as e:
        return "C", f"compare error ({type(e).__name__}: {e})"
    return "C", f"binary differs ({len(old)} -> {len(new)} bytes)"


# ---------------------------------------------------------------- run / restore
def restore(written, tracked, inputs, save_dir=None):
    to_checkout = []
    for rel, kind in written:
        p = WT / rel
        if kind == "deleted":
            if rel in tracked:
                to_checkout.append(rel)
            elif rel in inputs and (MAIN / rel).exists():
                p.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(MAIN / rel, p)
            continue
        if rel in tracked:
            to_checkout.append(rel)
        elif rel in inputs:
            shutil.copy2(MAIN / rel, p)              # a script overwrote a copied input: put it back
        else:
            if p.is_symlink() or p.is_file():
                p.unlink()
    for i in range(0, len(to_checkout), 200):
        subprocess.run(["git", "-C", str(WT), "checkout", HEAD, "--"] + to_checkout[i:i + 200], check=True,
                       capture_output=True)


def clean_check():
    st = subprocess.run(["git", "-C", str(WT), "status", "--porcelain"], capture_output=True, text=True).stdout
    tracked_mod = [l for l in st.splitlines() if not l.startswith("??")]
    return tracked_mod


def run(args):
    step = args.step
    out = SCR / step
    if out.exists():
        shutil.rmtree(out)
    (out / "out").mkdir(parents=True)
    pre_dirty = clean_check()
    if pre_dirty:
        sys.exit(f"worktree not clean before {step}: {pre_dirty[:5]}")
    tracked = tracked_set(); inputs = copied_inputs()
    patched = []
    for pf in args.patch or []:
        pf = Path(pf).resolve()
        r = subprocess.run(["git", "-C", str(WT), "apply", "--whitespace=nowarn", str(pf)], capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"patch {pf} failed: {r.stderr}")
        shutil.copy2(pf, out / pf.name)
        patched += [l[6:] for l in pf.read_text().splitlines() if l.startswith("+++ b/")]
    before = snapshot()
    env = dict(os.environ)
    env.update({"OMP_NUM_THREADS": "24", "OPENBLAS_NUM_THREADS": "24", "MKL_NUM_THREADS": "24",
                "NUMEXPR_MAX_THREADS": "24", "NUMBA_NUM_THREADS": "24", "MPLBACKEND": "Agg",
                "PYTHONDONTWRITEBYTECODE": "1"})
    for kv in args.env or []:
        k, v = kv.split("=", 1); env[k] = v
    cwd = WT / (args.cwd or ".")
    cmd = list(args.cmd)
    if cmd and cmd[0] == "--":
        cmd = cmd[1:]
    if cmd and cmd[0] in ("python", "py"):
        cmd[0] = PY
    full = ["nice", "-n", "10", "/usr/bin/time", "-v", "-o", str(out / "time.txt")]
    if args.strace:
        full += ["strace", "-f", "-qq", "-e", "trace=openat,stat,newfstatat,access", "-e", "status=failed",
                 "-o", str(out / "strace.txt")]
    full += cmd
    t0 = time.time()
    with open(out / "run.log", "w") as log:
        log.write(f"# cwd={cwd}\n# cmd={' '.join(cmd)}\n# env+={args.env}\n")
        log.flush()
        rc = subprocess.run(full, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    after = snapshot()
    written = []
    for rel, sig in after.items():
        if rel.startswith(("gate_tmp/",)):
            continue
        if rel not in before:
            written.append((rel, "new"))
        elif before[rel] != sig:
            written.append((rel, "modified"))
    for rel in before:
        if rel not in after:
            written.append((rel, "deleted"))
    res = []
    for rel, kind in sorted(written):
        p = WT / rel
        if kind == "deleted":
            res.append(dict(file=rel, kind=kind, tier="C", how="file deleted by the script", ref="")); continue
        new = p.read_bytes()
        if rel in tracked:
            old = head_blob(rel); ref = "HEAD blob"
        elif rel in inputs:
            old = (MAIN / rel).read_bytes(); ref = "copied input (main repo)"
        elif (MAIN / rel).exists() and (MAIN / rel).is_file():
            old = (MAIN / rel).read_bytes(); ref = "main-repo untracked copy"
        else:
            old = None; ref = "none"
        if old is None:
            tier, how = "N", f"new file, no committed counterpart ({len(new)} bytes)"
        else:
            tier, how = compare(rel, old, new)
        if tier != "A":
            dst = out / "out" / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dst)
        res.append(dict(file=rel, kind=kind, tier=tier, how=how, ref=ref))
    tm = (out / "time.txt").read_text() if (out / "time.txt").exists() else ""
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", tm)
    missing = []
    if args.strace and (out / "strace.txt").exists():
        seen = set()
        for line in (out / "strace.txt").read_text(errors="replace").splitlines():
            m = re.search(r'"(/[^"]+)"', line)
            if not m or "ENOENT" not in line:
                continue
            p = m.group(1)
            if not p.startswith(str(WT) + "/") or "__pycache__" in p or p.endswith((".pyc", ".so")):
                continue
            rel = os.path.relpath(p, WT)
            if rel in seen:
                continue
            seen.add(rel)
            if (MAIN / rel).exists():
                missing.append(rel)
    rec = dict(step=step, cwd=str(cwd.relative_to(WT)), cmd=" ".join(cmd), env=args.env or [], rc=rc, patched=patched,
               wall_s=round(wall, 1), max_rss_mb=round(int(rss.group(1)) / 1024, 1) if rss else None,
               files=res, missing_untracked_inputs=missing)
    (out / "result.json").write_text(json.dumps(rec, indent=1))
    if not args.keep:
        restore(written + [(p_, "modified") for p_ in patched], tracked, inputs)
        left = clean_check()
        if left:
            print("WARNING: worktree still dirty after restore:", left[:5])
    # print summary
    print(f"== {step}: rc={rc} wall={wall:.1f}s rss={rec['max_rss_mb']}MB files={len(res)}")
    for r in res:
        print(f"  [{r['tier']}] {r['file']} ({r['kind']}; vs {r['ref']}): {r['how'][:400]}")
    if missing:
        print("  MISSING untracked inputs present in main repo:", missing[:30])
    return rc


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd_")
    r = sub.add_parser("run")
    r.add_argument("step")
    r.add_argument("--cwd")
    r.add_argument("--env", action="append")
    r.add_argument("--strace", action="store_true")
    r.add_argument("--keep", action="store_true", help="do not restore the worktree afterwards")
    r.add_argument("--patch", action="append", help="git-apply this diff (worktree only) before the run; reverted after")
    sub.add_parser("restore")
    argv = sys.argv[1:]
    cmd = []
    if "--" in argv:
        i = argv.index("--"); cmd = argv[i + 1:]; argv = argv[:i]
    a = ap.parse_args(argv)
    a.cmd = cmd
    if a.cmd_ == "run":
        sys.exit(0 if run(a) == 0 else 1)
    elif a.cmd_ == "restore":
        st = subprocess.run(["git", "-C", str(WT), "status", "--porcelain", "-uall"], capture_output=True, text=True).stdout
        print(st[:3000])


if __name__ == "__main__":
    main()
