"""Record md5 of ASV-keyed / 'leave' data files of ~/Documents/EmilyKin (HEAD 9e087c3) before the R232 run.
Rows: every files_plan.tsv row whose action starts with 'leave', plus core ASV files with zero old-ID hits
(not listed in files_plan.tsv because it only lists files that contain an old ID)."""
import csv, glob, hashlib, os, subprocess, sys
REPO = os.path.expanduser("~/Documents/EmilyKin")
RP = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
OUT = f"{RP}/apply_log/asv_md5_before.tsv"
tracked = set(subprocess.run(["git", "-C", REPO, "ls-files", "-z"], capture_output=True, check=True).stdout.decode().split("\0"))
plan = list(csv.DictReader(open(f"{RP}/data/files_plan.tsv"), delimiter="\t"))
plan_action = {r["file"]: r["action"] for r in plan}
rows = {}
ASV_ACTIONS = ("leave (ASV context)", "leave (frozen run provenance; members are ASVs)")
for r in plan:
    if r["action"].startswith("leave"):
        rows[r["file"]] = dict(source="files_plan.tsv leave", action=r["action"], area=r["area"],
                               asv_keyed="yes" if r["action"] in ASV_ACTIONS else "no (leave, not ASV-keyed)")
core = ["abundances.csv", "abundances (1).csv", "taxonomy.json", "taxonIDs.json", "iterativeIDs.json",
        "iterativeID_taxonomy.json", "iterativeID_phylums.json", "iterativeID_color_map.json", "iterativeID_levels.json",
        "relative_abundance.json", "asv_to_mag_mapping.json", "asvs.fasta", "asv_iterativeID_old_to_new.csv",
        "asv_iterativeID_old_to_new.json", "Genus_color_map.json", "Species_color_map.json", "Phylum_color_map.json",
        "phylum_colors.json", "phylum_base_overrides.json", "proteo_class_color.json", "nonzero_per_day.json",
        "correlation_matrix_root.csv", "pvalue_matrix_root.csv", "FDR_passing_pairs.npy",
        "significantly_connected_organisms.json", "table_rel_export.csv", "sample_days.json", "phase_day_ranges.json",
        "can_sample_phase.json", "network/cooccurrence_edges_p_value_FDR.csv"]
core += sorted(glob.glob("network/network_module_membership_p_value_FDR*.json", root_dir=REPO))
core += sorted(glob.glob("network/network_summary_p_value_FDR*.json", root_dir=REPO))
for d in ["correlations", "dbRDA", "asv_mag_mapping", "network_method_comparison", "abundace_correlations", "rel_ab_heatmaps"]:
    core += sorted(f for f in tracked if f.startswith(d + "/"))
for f in core:
    if not f or f in rows or f not in tracked:
        continue
    act = plan_action.get(f)
    if act and not act.startswith("leave"):
        # in the plan with a non-leave action: not an unchanged ASV file; record but flag
        rows[f] = dict(source="core list, but files_plan.tsv action is NOT leave", action=act, area="", asv_keyed="no (plan changes it)")
    else:
        rows[f] = dict(source="core ASV file added (zero old-ID hits, not in files_plan.tsv)", action="leave (implied)", area="", asv_keyed="yes")
# C13: the all-phases network is MAG-only, not ASV context
C13 = {"network/network_module_membership_p_value_FDR_allphases.json", "network/network_summary_p_value_FDR_allphases.json"}
def eol(b):
    if b"\0" in b[:8192]:
        return "binary"
    crlf = b.count(b"\r\n"); lf = b.count(b"\n") - crlf
    return "CRLF" if crlf and not lf else ("LF" if lf and not crlf else ("mixed" if crlf and lf else "none"))
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["file", "md5", "bytes", "eol", "asv_keyed", "action", "source", "note"])
    for f in sorted(rows):
        p = os.path.join(REPO, f)
        if not os.path.exists(p):
            w.writerow([f, "MISSING", "", "", rows[f]["asv_keyed"], rows[f]["action"], rows[f]["source"], ""]); continue
        b = open(p, "rb").read()
        note = "C13: MAG-only network (140 nodes); D2-frozen keeps its results, but the apply may add bin/R232-ID fields" if f in C13 else ""
        w.writerow([f, hashlib.md5(b).hexdigest(), len(b), eol(b), rows[f]["asv_keyed"], rows[f]["action"], rows[f]["source"], note])
n = len(rows); y = sum(r["asv_keyed"] == "yes" for r in rows.values())
print(f"wrote {OUT}: {n} files ({y} ASV-keyed)")
from collections import Counter
print(Counter(r["source"] for r in rows.values()))
