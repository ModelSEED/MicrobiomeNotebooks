#!/usr/bin/env python
"""READ-ONLY: hard-coded MAG IDs and R220 taxon names in tracked code (py / R / sh / js / ipynb source cells).
One row per (file, line, string). Uses the token classifier of scan_repo.py.
Writes scan/code_hardcoded_raw.tsv."""
import csv, importlib.util, json, os, re, subprocess, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("scan", os.path.join(HERE, "scan_repo.py"))
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
EK = S.EK
files = subprocess.run(["git", "-C", EK, "ls-files", "-z", "--", "*.py", "*.R", "*.r", "*.sh", "*.js", "*.ipynb", "*.Rmd"],
                       capture_output=True).stdout.decode().split("\0")
rows = []
# extra literal checks (GTDB names written into code in either form)
EXTRA = re.compile(r"Patescibacteri\w*|Myxococcota\w*|Verrucomicrobi\w*|Paceibacteria|Minisyncoccia|Ignavibacteri\w*|"
                   r"Burkholderiaceae\w*|Flavobacteriales\w*")
for f in filter(None, files):
    p = os.path.join(EK, f)
    try: text = open(p, encoding="utf8", errors="replace").read()
    except Exception: continue
    if f.endswith(".ipynb"):
        try:
            nb = json.loads(text)
            lines = []
            for ci, c in enumerate(nb.get("cells", [])):
                for li, l in enumerate("".join(c.get("source", [])).splitlines()):
                    lines.append((f"cell{ci}:{li + 1}", l))
        except Exception:
            lines = [(str(i + 1), l) for i, l in enumerate(text.splitlines())]
    else:
        lines = [(str(i + 1), l) for i, l in enumerate(text.splitlines())]
    for ln, l in lines:
        if len(l) > 3000: l = l[:3000]
        toks = list(S.tokens(l))
        ks = {k for _, k in toks}
        for t, k in toks:
            if k in ("old_minted", "old_inherit", "old_display", "new_asv", "new_other"):
                ctx = "MAG" if k == "old_minted" else ""
                if k in ("old_inherit", "old_display"):
                    mag = ("bin" in ks) or ("mag_only" in ks) or bool(S.LINE_MAG.search(l))
                    asv = ("asv_only" in ks) or bool(S.LINE_ASV.search(l))
                    ctx = "MAG" if mag and not asv else "ASV" if asv and not mag else "ambiguous"
                if k == "new_asv": ctx = "ASV (today)"
                rows.append(dict(file=f, line=ln, kind=k, string=t, context=ctx, snippet=l.strip()[:200]))
        for m in S.NAMERE.finditer(l):
            n = m.group(1)
            pre = l[max(0, m.start() - 3):m.start()]
            info = S.r220_changed[n]
            rows.append(dict(file=f, line=ln, kind="r220_taxon_name" + ("_gtdb_prefixed" if pre.endswith("__") else ""), string=n,
                             context=("GTDB" if pre.endswith("__") else "ASV-name-collision" if any(x[2] for x in info) else "GTDB-only name"),
                             snippet=l.strip()[:200]))
with open(os.path.join(S.SCAN, "code_hardcoded_raw.tsv"), "w") as fh:
    w = csv.DictWriter(fh, fieldnames=["file", "line", "kind", "string", "context", "snippet"], delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
print(len(rows), "rows in", len({r['file'] for r in rows}), "code files")
