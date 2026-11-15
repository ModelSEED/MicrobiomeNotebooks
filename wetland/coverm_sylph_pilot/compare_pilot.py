#!/usr/bin/env python3
"""Compare the CoverM and sylph pilot outputs (poplar run_pilot.sh) for 2 wetland metagenomes.

Inputs (in the results directory, argv[1]): coverm_<s>.tsv (CoverM genome, defaults: min covered fraction 10%),
coverm_<s>.raw.tsv (same BAM, min covered fraction 0: counts, reads_per_base), sylph_derep.tsv (sylph profile -u vs
the same 188 dereplicated MAGs), sylph_all.tsv (vs all wetland MAGs), time_*.txt (/usr/bin/time -v), and the GTDB R232
lineage table for labels. Writes comparison tables + summary.json + summary.md next to the inputs.
"""
import glob, json, os, re, sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

D = sys.argv[1]
LIN = "/home/freiburger/gtdbtk_kbase/mag_gtdb_lineages_run3_merged.tsv"
SAMPLES = ["SRR21120499", "SRR21120498"]
PROJECT = {f"SRR21120{n}" for n in range(495, 504)}
lin = pd.read_csv(LIN, sep="\t", usecols=["name", "classification"]).set_index("name")["classification"]


def strip(g):
    return re.sub(r"\.(fna|fa|fasta)$", "", os.path.basename(str(g)))


def times():
    out = {}
    for f in glob.glob(os.path.join(D, "time_*.txt")):
        t = open(f).read()
        wall = re.search(r"Elapsed \(wall clock\) time.*: (\S+)", t).group(1)
        parts = [float(x) for x in wall.split(":")]
        sec = sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
        rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1)) / 1048576
        cpu = re.search(r"Percent of CPU this job got: (\S+)", t).group(1)
        out[os.path.basename(f)[5:-4]] = {"wall_min": round(sec / 60, 1), "max_rss_GB": round(rss, 1), "cpu": cpu}
    return out


def coverm(s):
    d = pd.read_csv(os.path.join(D, f"coverm_{s}.tsv"), sep="\t")
    d.columns = ["genome"] + [re.sub(r"^.*?\s(?=[A-Z])", "", c, count=1) for c in d.columns[1:]]
    raw = pd.read_csv(os.path.join(D, f"coverm_{s}.raw.tsv"), sep="\t")
    raw.columns = ["genome"] + [re.sub(r"^.*?\s(?=[A-Z])", "", c, count=1) for c in raw.columns[1:]]
    return d, raw


def sylph(path):
    if not os.path.exists(path):
        return None
    lines = [l for l in open(path) if l.strip()]
    hdr = next(i for i, l in enumerate(lines) if l.startswith("Sample_file"))
    rows = [l.rstrip("\n").split("\t") for l in lines[hdr:]]
    df = pd.DataFrame(rows[1:], columns=rows[0])
    for c in ("Taxonomic_abundance", "Sequence_abundance", "Adjusted_ANI", "True_cov", "Naive_ANI"):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["genome"] = df["Genome_file"].map(strip)
    df["sample"] = df["Sample_file"].map(lambda x: os.path.basename(x).split("_")[0])
    return df


def bray_curtis(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.abs(a - b).sum() / (a + b).sum()) if (a + b).sum() else float("nan")


summary = {"timing": times(), "samples": {}}
sd, sa = sylph(os.path.join(D, "sylph_derep.tsv")), sylph(os.path.join(D, "sylph_all.tsv"))
md = []
for s in SAMPLES:
    cv, raw = coverm(s)
    unm = float(cv.loc[cv.genome == "unmapped", "Relative Abundance (%)"].iloc[0])
    cvg = cv[cv.genome != "unmapped"].copy()
    rawg = raw[raw.genome != "unmapped"].copy()
    rawg["Read Count"] = pd.to_numeric(rawg["Read Count"])
    mapped_any = float(100 - raw.loc[raw.genome == "unmapped", "Relative Abundance (%)"].iloc[0])
    det_c = set(cvg.loc[cvg["Relative Abundance (%)"] > 0, "genome"])
    ys = sd[sd["sample"] == s].set_index("genome") if sd is not None else pd.DataFrame()
    det_s = set(ys.index)
    m = cvg.set_index("genome")[["Relative Abundance (%)", "Trimmed Mean", "Mean", "Covered Fraction"]].join(
        ys[["Taxonomic_abundance", "Sequence_abundance", "Adjusted_ANI", "True_cov"]], how="outer").fillna(0)
    m["coverm_rel_norm"] = m["Relative Abundance (%)"] / max(m["Relative Abundance (%)"].sum(), 1e-9) * 100
    m["lineage"] = [lin.get(g, "") for g in m.index]
    m.to_csv(os.path.join(D, f"compare_{s}.tsv"), sep="\t")
    both = m[(m["Relative Abundance (%)"] > 0) & (m["Taxonomic_abundance"] > 0)]
    r = {
        "coverm_reads_mapped_pct_default": round(100 - unm, 2),
        "coverm_reads_mapped_pct_any_coverage": round(mapped_any, 2),
        "sylph_sequence_abundance_sum_pct": round(float(ys["Sequence_abundance"].sum()), 2) if len(ys) else None,
        "detected_coverm": len(det_c), "detected_sylph": len(det_s), "detected_both": len(det_c & det_s),
        "only_coverm": sorted(det_c - det_s)[:20], "n_only_coverm": len(det_c - det_s),
        "only_sylph": sorted(det_s - det_c)[:20], "n_only_sylph": len(det_s - det_c),
        "spearman_relabund_both": round(spearmanr(both["coverm_rel_norm"], both["Taxonomic_abundance"])[0], 3) if len(both) > 2 else None,
        "pearson_relabund_both": round(pearsonr(both["coverm_rel_norm"], both["Taxonomic_abundance"])[0], 3) if len(both) > 2 else None,
        "spearman_trimmedmean_vs_truecov": round(spearmanr(both["Trimmed Mean"], both["True_cov"])[0], 3) if len(both) > 2 else None,
        "bray_curtis_relabund_all": round(bray_curtis(m["coverm_rel_norm"], m["Taxonomic_abundance"]), 3),
        "sylph_adjusted_ani_quantiles": ys["Adjusted_ANI"].quantile([0, .25, .5, .75, 1]).round(2).tolist() if len(ys) else None,
        "sylph_detected_ani_lt_99": int((ys["Adjusted_ANI"] < 99).sum()) if len(ys) else None,
        "own_sample_mags_detected_sylph": int(sum(1 for g in det_s if g.endswith(s))),
    }
    if sa is not None:
        ya = sa[sa["sample"] == s]
        own = ya["genome"].str.endswith("__" + s)
        inproj = ya["genome"].map(lambda g: g.split("__")[-1] in PROJECT)
        r.update({
            "sylph_all_detected": int(len(ya)),
            "sylph_all_detected_own_sample": int(own.sum()),
            "sylph_all_detected_same_project_other_samples": int((inproj & ~own).sum()),
            "sylph_all_detected_other_projects": int((~inproj).sum()),
            "sylph_all_sequence_abundance_sum_pct": round(float(ya["Sequence_abundance"].sum()), 2),
            "sylph_all_seqabund_from_other_projects_pct": round(float(ya.loc[~inproj, "Sequence_abundance"].sum()), 2),
            "sylph_all_top_other_project": ya.loc[~inproj].nlargest(5, "Sequence_abundance")[
                ["genome", "Sequence_abundance", "Adjusted_ANI"]].assign(
                lineage=lambda x: [lin.get(g, "") for g in x.genome]).to_dict("records"),
        })
    summary["samples"][s] = r
json.dump(summary, open(os.path.join(D, "summary.json"), "w"), indent=1, default=str)
print(json.dumps(summary, indent=1, default=str))
