#!/usr/bin/env python
"""READ-ONLY: for each file listed in scan/files_raw.tsv, find the tracked scripts (py/R/ipynb/sh/Rmd/js) that
mention its basename, and guess which of them WRITE it (heuristic: the line, or the 3 lines after it, contain a
write call / an output-path variable). Writes scan/generators.tsv (file, writers, readers, n_refs).
"""
import collections, csv, os, re, subprocess

EK = os.path.expanduser("~/Documents/EmilyKin")
SCAN = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/data/scan"
files = [r["file"] for r in csv.DictReader(open(os.path.join(SCAN, "files_raw.tsv")), delimiter="\t")]
# basenames to search; too-generic names (summary.json, README.md ...) are still searched but flagged
base = collections.defaultdict(list)
for f in files:
    base[os.path.basename(f)].append(f)
pat = os.path.join(SCAN, "_basenames.txt")
open(pat, "w").write("\n".join(sorted(base)) + "\n")
CODE = ["*.py", "*.R", "*.r", "*.ipynb", "*.sh", "*.Rmd", "*.js"]
res = subprocess.run(["git", "-C", EK, "grep", "-n", "-I", "-F", "-f", pat, "--"] + CODE, capture_output=True).stdout.decode("utf8", "replace")
WRITE = re.compile(r"to_csv|to_json|to_excel|to_parquet|savefig|\.write\(|writerows?|json\.dump|\.dump\(|write_text|write_bytes|"
                   r"open\([^)]*[\"'](?:w|a|wb|w\+)[\"']|ggsave|write\.(?:csv|table|tsv)|write_tsv|fwrite|saveRDS|png\(|pdf\(|"
                   r"\bOUT\w*\s*=|\bout\w*\s*=|\boutput\w*\s*=|--out|\bdst\b|\bsave\(|>\s*\"?\$?\w*/?[\w./-]*$|tee ", re.I)
odd = [b for b in base if not re.fullmatch(r"[\w.\-%()]+", b)]
lines_by_file = {}
hits = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))   # basename -> script -> [refs, writes]
for ln in res.splitlines():
    try:
        script, lno, text = ln.split(":", 2)
    except ValueError:
        continue
    if script not in lines_by_file:
        try:
            lines_by_file[script] = open(os.path.join(EK, script), encoding="utf8", errors="replace").read().splitlines()
        except Exception:
            lines_by_file[script] = []
    L = lines_by_file[script]; i = int(lno) - 1
    ctx = " ".join(L[i:i + 4]) if L else text
    cand = {t for t in re.findall(r"[\w.\-%()]+", text) if t in base}
    cand |= {t.split("/")[-1] for t in re.findall(r"[\w./\-%()]+", text) if t.split("/")[-1] in base}
    cand |= {b for b in odd if b in text}
    for b in cand:
        if b in text:
            h = hits[b][script]; h[0] += 1
            if WRITE.search(ctx): h[1] += 1
rows = []
for f in files:
    b = os.path.basename(f)
    hs = {s: v for s, v in hits.get(b, {}).items() if s != f}
    writers = sorted(s for s, v in hs.items() if v[1])
    readers = sorted(s for s, v in hs.items() if not v[1])
    rows.append(dict(file=f, basename_shared_by=len(base[b]), writers=";".join(writers[:12]) + (" ..." if len(writers) > 12 else ""),
                     readers=";".join(readers[:12]) + (" ..." if len(readers) > 12 else ""), n_scripts=len(hs)))
with open(os.path.join(SCAN, "generators.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
print(len(rows), "files;", sum(1 for r in rows if r["writers"]), "with a writer candidate;", sum(1 for r in rows if r["n_scripts"]), "with any script ref")
