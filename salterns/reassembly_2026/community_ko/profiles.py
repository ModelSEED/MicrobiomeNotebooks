"""Community gene profiles of the methylated-compound panel (ko_panel_methylated.tsv) for Figure 1, Figure S4 and Table S1.

read_profile()      read-based survey (the version named in ../read_survey_current.json, or an explicit results dir):
                    per KO × sample: reads, copies per genome equivalent (GE), and, where the survey gives them, the
                    archaeal-only reads and copies per GE.
assembly_profile()  assembly-based KOfam profile (assembly_based/results/ko_by_sample.tsv[.gz]) when present.
panel()             the panel with the figure annotations (compound class, process type).
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R26 = os.path.dirname(HERE)
SED = [f"{s}_{c}_{d}" for s in ("R1", "R2", "SF2", "R2A") for c in "ABC" for d in ("D1", "D2")]

# compound classes (figure annotation), in row order; keys are the panel's `group` values
CLASS = [
    ("TMA and TMAO", ["TMA", "TMAO", "TMA/DMA"]),
    ("DMA", ["DMA"]),
    ("MMA", ["MMA", "methylamines"]),
    ("Choline and carnitine", ["choline", "carnitine", "betaine aldehyde"]),
    ("GB, DMG and sarcosine", ["glycine betaine", "dimethylglycine", "sarcosine", "glycine", "serine"]),
    ("Methanol", ["methanol"]),
    ("Formaldehyde and C1", ["formaldehyde", "formate", "C1 (H4F)"]),
    ("DMS and sulfonates", ["DMS", "DMS/methanethiol", "DMSO", "DMSO2", "methanethiol", "methanesulfonate"]),
    ("DMSP and MMPA", ["DMSP", "MMPA"]),
    ("Phosphonates", ["MPn", "phosphonates"]),
    ("Methane", ["methane"]),
    ("Acetate", ["acetate", "acetyl-CoA"]),
    ("H2/CO2", ["CO2/H2", "methyl-H4MPT"]),
]
METHANOGENESIS = {"mttB", "mttC", "mtbB", "mtbC", "mtmB", "mtmC", "mtbA", "mtaB", "mtaC", "mtaA", "mtgB", "mtsA",
                  "mcrA", "mcrB", "mcrG"}
TMA_PRODUCTION = {"cutC", "cutD", "cntA", "cntB", "grdI", "grdH", "torA", "torZ"}
SYNTHESIS = {"betA", "betB", "gbsB", "CMO", "gsdmt", "sdmt", "mpnS", "dmsA", "mddA"}
TRANSPORT = {"amtB", "phnD"}


def panel():
    p = pd.read_csv(os.path.join(HERE, "ko_panel_methylated.tsv"), sep="\t")
    cmap = {g: c for c, gs in CLASS for g in gs}
    p["compound_class"] = p["group"].map(cmap)
    assert p.compound_class.notna().all(), p.loc[p.compound_class.isna(), ["ko", "group"]]

    def ptype(r):
        if r.status == "FigS4":
            return "methanogenesis (C1 core)" if r.compound_class == "H2/CO2" else "acetate / CODH"
        g = r.gene
        if g in METHANOGENESIS:
            return "methanogenesis"
        if g in TMA_PRODUCTION:
            return "TMA production"
        if g in SYNTHESIS:
            return "biosynthesis / production"
        if g in TRANSPORT:
            return "transport"
        return "degradation"
    p["ptype"] = p.apply(ptype, axis=1)
    p["label"] = p.gene
    dup = p.gene.duplicated(keep="first")
    p.loc[dup, "label"] = p.loc[dup, "gene"].str.split(" ").str[0] + "_" + p.loc[dup, "ko"]
    assert p.label.is_unique
    p["order"] = np.arange(len(p))
    return p


def survey_dir(results=None):
    if results:
        return results
    cur = json.load(open(os.path.join(R26, "read_survey_current.json")))
    return os.path.join(R26, cur["results"])


def read_profile(results=None):
    """Long table: ko, gene, sample, reads, copies_per_GE (+ reads_archaeal, copies_per_GE_archaeal when given)."""
    d = survey_dir(results)
    t = pd.read_csv(os.path.join(d, "per_sample_KO.tsv"), sep="\t")
    p = panel()
    if "ko" in t.columns and t["ko"].notna().any():
        t = t[t["ko"].isin(p.ko)].copy()
    else:                                                     # v1: gene-symbol labels only
        sym = dict(zip(p.label, p.ko))
        sym.update({g: k for g, k in zip(p.gene, p.ko) if g not in sym})
        t = t[t.label.isin(sym)].copy()
        t["ko"] = t.label.map(sym)
    keep = ["ko", "sample", "reads", "copies_per_GE"] + [c for c in ("reads_archaeal", "copies_per_GE_archaeal")
                                                           if c in t.columns]
    t = t[keep].groupby(["ko", "sample"], as_index=False).sum()
    return t


def wide(t, value="copies_per_GE", samples=SED, kos=None):
    w = t.pivot(index="ko", columns="sample", values=value).reindex(columns=samples)
    if kos is not None:
        w = w.reindex(kos)
    return w


def assembly_profile():
    base = os.path.join(HERE, "assembly_based", "results", "ko_by_sample.tsv")
    for f in (base, base + ".gz"):
        if os.path.exists(f):
            return pd.read_csv(f, sep="\t")
    return None
