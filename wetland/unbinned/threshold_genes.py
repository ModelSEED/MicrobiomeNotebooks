#!/usr/bin/env python3
"""Threshold sensitivity, part 2b: what genes and functions does each contig-length band add? (per Andrew 2026-11-25)

Per sample, classes as in threshold_sensitivity.py:
  binned   = contigs of the sample's admitted MAGs (all >= 2.5 kb)
  unb_2.5  = unbinned >= 2,500 bp            (the 2.5 kb paradigm: in the assembly, not in a MAG)
  unb_2.0  = 2,000-2,499 bp                  (what a 2 kb floor adds; the KBase assemblies stop here)
  unb_1.0  = 1,000-1,999 bp, from the 1 kb re-assembly (reassembly/), minus contigs already contained in a KBase
             contig (compare_reassembly.py)
Annotation: Bakta 1.11.4, DB 6.0-full, --meta, as for the MAGs. The KBase-assembly part and the re-assembly's 1-2 kb
part are separate Bakta runs with identical settings; Bakta annotates each contig independently (pyrodigal meta mode
picks its model per sequence), so the split does not change the calls.

Caveat measured here, not assumed: Bakta calls genes on linear contigs with pyrodigal closed=True (bakta/features/
cds.py), so a gene running off a contig edge is NOT called. Short contigs lose relatively more genes this way. Each
band is therefore also gene-called with pyrodigal (meta mode) closed=False: genes, edge-partial genes, and the share
of open-ended genes Bakta's closed mode leaves out.

Per class (per sample, then pooled over samples):
  - structure: contigs, bp, Bakta CDS per Mbp; pyrodigal open-ended genes per Mbp and their partial share
  - annotation yield of Bakta CDS: hypothetical, UniRef90, UniRef50-only, KO, EC, COG, Pfam (%)
  - non-coding: tRNA, tmRNA, rRNA 5S/16S/23S, ncRNA, ncRNA-region, CRISPR (per Mbp, and the share of the sample's
    16S/23S genes in each class)
  - functions found only in this class: KO, EC, COG and named products present in the class but absent from every
    >= 2.5 kb contig of the sample (binned or not), for unb_2.0 and unb_1.0. Null: a random set of >= 2.5 kb contigs
    of the same total bp, removed from the >= 2.5 kb pool, scored the same way (|F(S) - F(rest)|, mean of NULL_REPS
    draws). For unb_2.5 the comparison is against the MAGs only (what the unbinned >= 2.5 kb contigs add), no null.
Inputs per run: ~/gtdbtk_kbase/bakta_thr/<run>/kbase.tsv and reasm12.tsv (Bakta feature tables), the KBase assembly,
~/gtdbtk_kbase/reasm/<run>.1to2kb.fa and <run>.classes.tsv.gz.
Outputs: unbinned/threshold/genes_per_class.tsv, genes_only_in.tsv, genes_summary.json.
Usage: ~/Documents/py_venv/bin/python -I threshold_genes.py RUN [RUN ...] [--threads 8]
"""
import argparse, concurrent.futures as cf, glob, json, os, random
from collections import defaultdict
import pandas as pd
import pyrodigal

ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
ASM = os.path.expanduser("~/gtdbtk_kbase/assemblies"); MAGS = os.path.expanduser("~/gtdbtk_kbase/genomes")
RE = os.path.expanduser("~/gtdbtk_kbase/reasm"); BK = os.path.expanduser("~/gtdbtk_kbase/bakta_thr")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "threshold")
CL = ("binned", "unb_2.5", "unb_2.0", "unb_1.0"); NULL_REPS = 20
COLS = ["contig", "type", "start", "stop", "strand", "locus_tag", "gene", "product", "dbxrefs"]
FUNCS = ("KO", "EC", "COG", "product")


def fasta(path):
    name = None; seq = []
    for l in open(path):
        if l[0] == ">":
            if name is not None: yield name.split()[0], "".join(seq)
            name = l[1:]; seq = []
        else: seq.append(l.strip())
    if name is not None: yield name.split()[0], "".join(seq)


def bakta(path):
    d = pd.read_csv(path, sep="\t", comment="#", header=None, names=COLS, dtype=str, keep_default_na=False)
    x = d.dbxrefs.str.split(r",\s*")
    pick = lambda pre: x.map(lambda v: {t[len(pre):] for t in v if t.startswith(pre)})
    d["KO"], d["EC"] = pick("KEGG:"), pick("EC:")
    d["COG"] = x.map(lambda v: {t[4:] for t in v if t.startswith("COG:COG")})
    d["PFAM"] = x.map(lambda v: {t for t in v if t.startswith("PFAM:")})
    d["U90"] = d.dbxrefs.str.contains("UniRef90_"); d["U50"] = d.dbxrefs.str.contains("UniRef50_") & ~d.U90
    d["hyp"] = d["product"].str.contains("hypothetical protein", case=False)
    d["product_set"] = [set() if h else {p} for h, p in zip(d.hyp, d["product"])]
    return d


finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=True)


def open_genes(seq):
    g = finder.find_genes(seq.encode())
    return len(g), sum(1 for x in g if x.partial_begin or x.partial_end)


per_rows, only_rows = [], []
for run in a.runs:
    binned = set()
    for f in glob.glob(f"{MAGS}/*__{run}.fna"):
        binned |= {l[1:].split()[0] for l in open(f) if l[0] == ">"}
    cls, seqs = {}, {}
    for n, s in fasta(glob.glob(f"{ASM}/{run}_*.fna")[0]):
        L = len(s); cls[n] = "binned" if n in binned else "unb_2.5" if L >= 2500 else "unb_2.0"; seqs[n] = s
    rc = pd.read_csv(f"{RE}/{run}.classes.tsv.gz", sep="\t", keep_default_na=False)
    new12 = set(rc[(rc.length < 2000) & (rc.match == "")].contig)          # not already inside a KBase contig
    for n, s in fasta(f"{RE}/{run}.1to2kb.fa"):
        if n in new12: cls["R:" + n] = "unb_1.0"; seqs["R:" + n] = s
    f1 = bakta(f"{BK}/{run}/kbase.tsv"); f2 = bakta(f"{BK}/{run}/reasm12.tsv"); f2["contig"] = "R:" + f2.contig
    f2 = f2[f2.contig.isin(cls)]
    ft = pd.concat([f1, f2], ignore_index=True); ft["cls"] = ft.contig.map(cls)
    assert ft.cls.notna().all(), ft[ft.cls.isna()].contig.head()
    names = list(seqs)
    with cf.ThreadPoolExecutor(a.threads) as ex:
        og = dict(zip(names, ex.map(lambda n: open_genes(seqs[n]), names)))
    bp = defaultdict(int); nct = defaultdict(int); ogn = defaultdict(int); ogp = defaultdict(int)
    for n, c in cls.items():
        bp[c] += len(seqs[n]); nct[c] += 1; ogn[c] += og[n][0]; ogp[c] += og[n][1]
    cds = ft[ft.type == "cds"]; rr = ft[ft.type == "rRNA"]
    tot16 = max(1, rr["product"].str.startswith("16S").sum()); tot23 = max(1, rr["product"].str.startswith("23S").sum())
    for c in CL:
        x = cds[cds.cls == c]; n = max(1, len(x)); mb = max(1e-9, bp[c] / 1e6); y = ft[ft.cls == c]; r = rr[rr.cls == c]
        per_rows.append({
            "run": run, "class": c, "contigs": nct[c], "bp": bp[c], "bakta_CDS": len(x), "bakta_CDS_per_Mbp": round(len(x) / mb, 1),
            "open_genes": ogn[c], "open_genes_per_Mbp": round(ogn[c] / mb, 1), "open_partial_pct": round(ogp[c] / max(1, ogn[c]) * 100, 2),
            "bakta_vs_open_pct": round(len(x) / max(1, ogn[c]) * 100, 2),
            "hypothetical_pct": round(x.hyp.sum() / n * 100, 2), "UniRef90_pct": round(x.U90.sum() / n * 100, 2),
            "UniRef50_only_pct": round(x.U50.sum() / n * 100, 2), "KO_pct": round((x.KO.str.len() > 0).sum() / n * 100, 2),
            "EC_pct": round((x.EC.str.len() > 0).sum() / n * 100, 2), "COG_pct": round((x.COG.str.len() > 0).sum() / n * 100, 2),
            "Pfam_pct": round((x.PFAM.str.len() > 0).sum() / n * 100, 2),
            **{f"{t}_per_Mbp": round((y.type == t).sum() / mb, 3) for t in ("tRNA", "tmRNA", "ncRNA", "ncRNA-region", "crispr", "sorf")},
            "rRNA_5S": int(r["product"].str.startswith("5S").sum()), "rRNA_16S": int(r["product"].str.startswith("16S").sum()),
            "rRNA_23S": int(r["product"].str.startswith("23S").sum()),
            "share_of_sample_16S_pct": round(r["product"].str.startswith("16S").sum() / tot16 * 100, 2),
            "share_of_sample_23S_pct": round(r["product"].str.startswith("23S").sum() / tot23 * 100, 2)})
    # functions only in a class, vs every >= 2.5 kb contig, with a bp-matched null
    fset = lambda df, k: set().union(*(df["product_set"] if k == "product" else df[k])) if len(df) else set()
    long_ctg = [n for n, c in cls.items() if c in ("binned", "unb_2.5")]
    long_f = {k: fset(cds[cds.cls.isin(("binned", "unb_2.5"))], k) for k in FUNCS}
    rnd = random.Random(f"{run}-null")
    for c in ("unb_2.0", "unb_1.0", "unb_2.5"):
        ref = {k: fset(cds[cds.cls == "binned"], k) for k in FUNCS} if c == "unb_2.5" else long_f
        obs = {k: len(fset(cds[cds.cls == c], k) - ref[k]) for k in FUNCS}
        pool = list(long_ctg); null = defaultdict(list)
        for _ in range(NULL_REPS if c != "unb_2.5" else 0):   # unb_2.5 vs binned: no bp-matched null (it outweighs binned)
            rnd.shuffle(pool); S, acc = [], 0
            for n in pool:
                if acc >= bp[c]: break
                S.append(n); acc += len(seqs[n])
            Sset = set(S); sdf = cds[cds.contig.isin(Sset)]
            rest = cds[cds.cls.isin(("binned", "unb_2.5")) & ~cds.contig.isin(Sset)]
            for k in FUNCS: null[k].append(len(fset(sdf, k) - fset(rest, k)))
        only_rows.append({"run": run, "class": c, "vs": "binned" if c == "unb_2.5" else "all_ge2.5kb", "class_bp": bp[c],
                          **{f"{k}_only": obs[k] for k in FUNCS},
                          **{f"{k}_only_null_mean": (round(sum(null[k]) / len(null[k]), 1) if null[k] else None) for k in FUNCS},
                          **{f"{k}_in_class": len(fset(cds[cds.cls == c], k)) for k in FUNCS}})
    print(run, "done", flush=True)
per = pd.DataFrame(per_rows); only = pd.DataFrame(only_rows)
per.to_csv(f"{OUT}/genes_per_class.tsv", sep="\t", index=False); only.to_csv(f"{OUT}/genes_only_in.tsv", sep="\t", index=False)
pool = {}
for c, g in per.groupby("class"):
    mb = g.bp.sum() / 1e6; cds_n = g.bakta_CDS.sum()
    w = lambda col: round(float((g[col] * g.bakta_CDS).sum() / max(1, cds_n)), 2)
    pool[c] = {"samples": len(g), "contigs": int(g.contigs.sum()), "bp": int(g.bp.sum()),
               "bp_pct_of_all": round(g.bp.sum() / per.bp.sum() * 100, 2),
               "bakta_CDS_per_Mbp": round(cds_n / mb, 1), "open_genes_per_Mbp": round(g.open_genes.sum() / mb, 1),
               "open_partial_pct": round(float((g.open_partial_pct * g.open_genes).sum() / g.open_genes.sum()), 2),
               "bakta_vs_open_pct": round(cds_n / g.open_genes.sum() * 100, 2),
               **{k: w(k) for k in ("hypothetical_pct", "UniRef90_pct", "UniRef50_only_pct", "KO_pct", "EC_pct", "COG_pct", "Pfam_pct")},
               "rRNA_16S": int(g.rRNA_16S.sum()), "rRNA_23S": int(g.rRNA_23S.sum()),
               "tRNA_per_Mbp": round(float((g["tRNA_per_Mbp"] * g.bp).sum() / g.bp.sum()), 3)}
summ = {"runs": a.runs, "pooled": pool, "only_in": only.groupby("class")[[c for c in only.columns if c.endswith("_only") or c.endswith("_null_mean")]].sum().round(1).to_dict("index")}
json.dump(summ, open(f"{OUT}/genes_summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
