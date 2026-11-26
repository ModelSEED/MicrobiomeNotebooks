"""Validation (2026-11-24): local CheckM 1.1.2 (hw.run_checkm, FM or TF) on bins that KBase's kb_Msuite already
CheckM'd, compared bin by bin with the KBase values (work/checkm/<unit>/checkm_<unit>.tsv). Saves nothing to KBase.
Usage: python checkm_vs_kbase.py <unit> <AssemblySet ref> <out dir> [TF]"""
import csv, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import hw
unit, target, out = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
os.makedirs(out, exist_ok=True)
spec = {"workdir": out, "unit": unit, "threads": 8, "remote": sys.argv[4] if len(sys.argv) > 4 else None}
items = [it["ref"] for it in hw.get_obj(target)["data"]["items"]]
infos = hw.rpc("ws", "Workspace.get_object_info3", [{"objects": [{"ref": r} for r in items]}])[0]["infos"]
bdir = os.path.join(out, "bins"); os.makedirs(bdir, exist_ok=True)
for r, i in zip(items, infos):
    hw.fetch_assembly_fasta(r, os.path.join(bdir, i[1] + ".fna"))
ext, qa, _ = hw.run_checkm(spec, out, bdir, {i[1] for i in infos})
kb = {r["assembly_name"]: r for r in csv.DictReader(open(
    f"/home/freiburger/Documents/MicrobiomeNotebooks/grow/assembly_binning_2026/work/checkm/{unit}/checkm_{unit}.tsv"), delimiter="\t")}
d_comp = [abs(float(ext[n]["Completeness"]) - float(kb[n]["completeness"])) for n in kb]
d_cont = [abs(float(ext[n]["Contamination"]) - float(kb[n]["contamination"])) for n in kb]
lin = sum(1 for n in kb if str(ext[n].get("marker lineage")) == kb[n]["marker_lineage"])
pas = sum(1 for n in kb if (float(ext[n]["Completeness"]) >= 50 and float(ext[n]["Contamination"]) < 10) == (kb[n]["pass_50_10"] == "1"))
res = {"unit": unit, "bins": len(kb), "host": spec["remote"] or "FM", "max_abs_diff_completeness": max(d_comp),
       "max_abs_diff_contamination": max(d_cont), "identical_values": sum(1 for a, b in zip(d_comp, d_cont) if a < 1e-9 and b < 1e-9),
       "same_marker_lineage": lin, "same_50_10_call": pas}
print(json.dumps(res, indent=1)); json.dump(res, open(os.path.join(out, "comparison.json"), "w"), indent=1)
