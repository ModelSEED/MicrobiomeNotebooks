"""Smoke test (2026-11-24) of hw.py's three binner stages on a small real dataset (1,500 deadrun_0006 contigs, 1M
read pairs), with KBase fetch/save stubbed out. Two 'libraries' (the same reads twice) exercise the multi-library
paths that sitecov uses. Saves nothing to KBase. Usage: python smoke_binners.py <stage> [nlibs]"""
import json, os, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import hw
stage, nlibs = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1
D = os.path.join(HERE, "smoke")
wd = os.path.join(D, f"{stage}_{nlibs}lib"); shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
hw.fetch_assembly_fasta = lambda ref, path: (shutil.copy(os.path.join(D, "asm.fa"), path), {})
hw.fetch_reads = lambda ref, cache: (os.path.join(D, "reads.inter.fastq.gz"), {})
hw.check_against_assembly = lambda ci, ref: None
def fake_save(spec, st, bin_dir, p, ci):
    data = hw.binned_contigs_from_dir(bin_dir, p["assembly_ref"], ci)
    if not data["bins"]:
        raise hw.Empty("no bins")
    json.dump(data, open(os.path.join(spec["workdir"], "binned_contigs.json"), "w"))
    return {"n_bins": len(data["bins"]), "binned_bp": sum(c["len"] for b in data["bins"] for c in b["contigs"].values()),
            "total_contig_len": data["total_contig_len"], "example_bin": {k: v for k, v in data["bins"][0].items() if k != "contigs"}}
hw.save_binned = fake_save
base = {"assembly_ref": "0/0/0", "reads_list": [f"9/{i}/1" for i in range(nlibs)], "workspace_name": "x"}
P = {"metabat": {**base, "min_contig_length": 2500},
     "maxbin": {**base, "thread": 32, "prob_threshold": 0.8, "markerset": "40", "min_contig_length": 1000, "plotmarker": 0},
     "concoct": {**base, "read_mapping_tool": "bowtie2_default", "min_contig_length": 2500, "contig_split_size": 10000,
                 "contig_split_overlap": 0, "kmer_size": 4, "max_clusters_for_vgmm": 400, "max_iterations_for_vgmm": 500,
                 "total_percentage_pca": 90, "no_cov_normalization": 0, "no_total_coverage": 0}}[stage]
spec = {"stage": stage, "unit": "smoke", "workdir": wd, "threads": 4, "mem_gb": 8, "params": P, "reads_cache": D}
tmp = os.path.join(wd, "tmp"); os.makedirs(tmp)
try:
    out = hw.STAGES[stage](spec, tmp)
    print(json.dumps({"stage": stage, "nlibs": nlibs, "status": "done", **out}, indent=1))
except hw.Empty as e:
    print(json.dumps({"stage": stage, "nlibs": nlibs, "status": "empty", "error": str(e)}))
