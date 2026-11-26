"""Collect gate_scratch/*/result.json into apply_log/gate_results.tsv (one row per written file) and print a step table."""
import json, csv
from pathlib import Path
RP = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan")
rows = []; steps = []
for p in sorted((RP / "gate_scratch").glob("*/result.json")):
    r = json.loads(p.read_text())
    if not isinstance(r, dict):
        continue
    tiers = {}
    for f in r["files"]:
        tiers[f["tier"]] = tiers.get(f["tier"], 0) + 1
        rows.append([r["step"], r["cmd"].replace(str(Path.home()) + "/Documents/py_venv/bin/", ""), r["cwd"], ";".join(r["env"]),
                     ";".join(r.get("patched", [])), r["rc"], r["wall_s"], r["max_rss_mb"], f["file"], f["kind"], f["ref"], f["tier"], f["how"]])
    if not r["files"]:
        rows.append([r["step"], r["cmd"].replace(str(Path.home()) + "/Documents/py_venv/bin/", ""), r["cwd"], ";".join(r["env"]),
                     ";".join(r.get("patched", [])), r["rc"], r["wall_s"], r["max_rss_mb"], "(no file written)", "", "", "", ""])
    steps.append((r["step"], r["rc"], r["wall_s"], r["max_rss_mb"], " ".join(f"{k}{v}" for k, v in sorted(tiers.items()))))
with open(RP / "apply_log" / "gate_results.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["step", "command", "cwd", "env", "worktree_patches", "rc", "wall_s", "max_rss_MB", "file", "kind",
                "compared_against", "tier", "how"])
    w.writerows(rows)
for s in steps:
    print("\t".join(map(str, s)))
print(len(rows), "rows")
