"""Per-genome gene table and annotation evidence for the four arms of the pre-registered comparison.

Re-runnable; reads only. For each of the 57 modelling representatives writes
  ../annotation_comparison/evidence/<mag_id>.json.gz  with
    genes:  model gene id -> [protein_id, contig, start, end, strand]   (Prodigal genes = RAST features)
    rast:   gene -> [RAST roles]   (rast_service/functions, split like fliu's builder: '; ', ' / ', ' @ ')
    sources: {kofam|eggnog|bakta: gene -> [[evidence type, term, [ModelSEED rxn ids]], ...]}
    bakta_products: gene -> [product strings]   (matched to template role names at build time)
    bakta_extra_genes: Bakta CDS with no Prodigal partner (enter as extra features)
    stats: counts, Bakta overlap matching
Model gene ids are the Prodigal protein ids without the '<kofam_name>|' prefix (cobra parses '|' as OR).

Evidence rules (ANNOTATION_CHOICE_CRITERIA.md):
  kofam : DRAM-rule best KO (kofam_id) -> ko2rxn; accepted Sec read-through rescues replace the KO of their pieces.
  eggnog: KEGG_ko (KO), EC, BiGG_Reaction (BiGG) -> maps (runs on the same Prodigal proteins: identity join).
  bakta : EC -> ec2rxn; product -> template roles (at build time). Bakta calls its own genes: a Bakta CDS is
          attached to a Prodigal gene on the same contig and strand with >= 80 % reciprocal coordinate overlap
          (best overlap if several); unmatched CDS become extra features 'bakta_<locus_tag>'. Pseudogenes skipped.
usage: python evidence.py [mag_id ...]
"""
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
M26 = HERE.parent
R26 = M26.parent
ANN = M26 / "annotations"
OUT = M26 / "annotation_comparison" / "evidence"
MAPS = json.load(open(HERE / "maps/term_maps.json"))
KO2, EC2, BIGG2 = MAPS["ko2rxn"], MAPS["ec2rxn"], MAPS["bigg2rxn"]
EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
MIN_RECIP = 0.80


def split_functions(fn):
    """fliu's split (cliff_mags/mac/salterns/build_models.ipynb cell 25), as in the scouting common.py."""
    return [x for x in re.split("; | / | @", fn) if x]


def gid(protein_id):
    return protein_id.split("|", 1)[1]


def ec_list(s):
    if not isinstance(s, str) or not s:
        return []
    return [e.strip() for e in re.split(r"[;,]", s) if EC_RE.match(e.strip())]


def manifest():
    m = pd.read_csv(ANN / "scripts/mags57.tsv", sep="\t")
    return m.set_index("mag_id")


def overlap_match(prod, bak):
    """prod/bak: DataFrames with contig,start,end,strand (+ gene / locus_tag). Returns {locus_tag: gene}."""
    out = {}
    pg = {k: g.sort_values("start") for k, g in prod.groupby(["contig", "strand"])}
    for _, b in bak.iterrows():
        g = pg.get((b.contig, b.strand))
        if g is None:
            continue
        cand = g[(g.start <= b.end) & (g.end >= b.start)]
        best, bestv = None, 0.0
        for _, p in cand.iterrows():
            ov = min(p.end, b.end) - max(p.start, b.start) + 1
            f1, f2 = ov / (p.end - p.start + 1), ov / (b.end - b.start + 1)
            if f1 >= MIN_RECIP and f2 >= MIN_RECIP and min(f1, f2) > bestv:
                best, bestv = p.gene, min(f1, f2)
        if best is not None:
            out[b.locus_tag] = best
    return out


def build(mag_id, man):
    row = man.loc[mag_id]
    kof = pd.read_csv(R26 / row.local_kofam_per_mag, sep="\t", dtype={"kofam_id": str})
    rast = pd.read_csv(M26 / "rast_service/functions" / f"{mag_id}.tsv", sep="\t", dtype=str).fillna("")
    assert set(rast.protein_id) == set(kof.gene), mag_id
    kof["g"] = kof.gene.map(gid)
    assert kof.g.is_unique, mag_id
    genes = {r.g: [r.gene, r.contig, int(r.start), int(r.end), r.strand] for r in kof.itertuples()}
    rast_roles = {}
    for r in rast.itertuples():
        roles = split_functions(r.function) if r.function else []
        if roles:
            rast_roles[gid(r.protein_id)] = roles
    src = {"kofam": defaultdict(list), "eggnog": defaultdict(list), "bakta": defaultdict(list)}
    stats = {"prodigal_genes": len(genes), "rast_genes_with_function": len(rast_roles)}

    # --- KOfam (+ Sec rescue) ---
    ko_of = {r.g: r.kofam_id for r in kof.itertuples() if isinstance(r.kofam_id, str) and r.kofam_id}
    sec = pd.read_csv(ANN / "kofam/sec_readthrough_57.tsv", sep="\t", dtype=str)
    sec = sec[(sec.mag_id == mag_id) & (sec.accepted == "True")]
    n_sec = 0
    for r in sec.itertuples():
        for pc in str(r.pieces).split(";"):
            if pc and pc != "nan":
                ko_of[gid(pc)] = r.ko
                n_sec += 1
    for g, ko in ko_of.items():
        src["kofam"][g].append(["KO", ko, KO2.get(ko, [])])
    stats.update(kofam_genes_with_ko=len(ko_of), kofam_sec_pieces_reassigned=n_sec)

    # --- eggNOG (same Prodigal proteins) ---
    egg = pd.read_csv(ANN / f"eggnog/{mag_id}/{mag_id}.eggnog_genes.tsv", sep="\t", dtype=str).fillna("")
    assert set(egg.gene) == set(kof.gene), mag_id
    for r in egg.itertuples():
        g = gid(r.gene)
        for ko in [k for k in r.ko.split(",") if k]:
            src["eggnog"][g].append(["KO", ko, KO2.get(ko, [])])
        for ec in ec_list(r.EC):
            src["eggnog"][g].append(["EC", ec, EC2.get(ec, [])])
        for b in [x for x in r.BiGG_Reaction.split(",") if x]:
            src["eggnog"][g].append(["BiGG", b, BIGG2.get(b, [])])
    stats["eggnog_genes_with_term"] = len(src["eggnog"])

    # --- Bakta (own gene calls -> Prodigal by reciprocal overlap) ---
    bk = pd.read_csv(ANN / f"bakta/{mag_id}/{mag_id}.bakta_genes.tsv", sep="\t", dtype=str).fillna("")
    bk = bk[(bk.type == "cds") & (bk.pseudogene != "True")].copy()
    bk["start"], bk["end"] = bk.start.astype(int), bk.end.astype(int)
    prod = kof[["g", "contig", "start", "end", "strand"]].rename(columns={"g": "gene"})
    match = overlap_match(prod, bk)
    products, extra = defaultdict(list), {}
    for r in bk.itertuples():
        g = match.get(r.locus_tag)
        if g is None:
            g = f"bakta_{r.locus_tag}"
            extra[g] = [r.locus_tag, r.contig, int(r.start), int(r.end), r.strand]
        ecs = ec_list(r.ec)
        for ec in ecs:
            src["bakta"][g].append(["EC", ec, EC2.get(ec, [])])
        if r.product and r.product != "hypothetical protein":
            products[g].append(r.product)
            for ec in ecs:                       # RAST-style role spelling "product (EC x)"
                products[g].append(f"{r.product} (EC {ec})")
    stats.update(bakta_cds=len(bk), bakta_matched=len(match), bakta_unmatched=len(extra),
                 bakta_genes_with_ec=sum(1 for v in src["bakta"].values() if v),
                 bakta_genes_with_product=len(products))
    out = {"mag_id": mag_id, "genes": genes, "rast": rast_roles,
           "sources": {k: dict(v) for k, v in src.items()}, "bakta_products": dict(products),
           "bakta_extra_genes": extra, "stats": stats,
           "inputs": {"kofam": str(row.local_kofam_per_mag), "rast": f"modelling_2026/rast_service/functions/{mag_id}.tsv",
                      "eggnog": f"modelling_2026/annotations/eggnog/{mag_id}/{mag_id}.eggnog_genes.tsv",
                      "bakta": f"modelling_2026/annotations/bakta/{mag_id}/{mag_id}.bakta_genes.tsv",
                      "sec": "modelling_2026/annotations/kofam/sec_readthrough_57.tsv",
                      "maps": "modelling_2026/pipeline/maps/term_maps.json"}}
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / f"{mag_id}.json.gz", "wt") as fh:
        json.dump(out, fh)
    return stats


def main(argv):
    man = manifest()
    mags = argv or list(man.index)
    rows = []
    for m in mags:
        s = build(m, man)
        rows.append({"mag_id": m, **s})
        print(m, s, flush=True)
    if not argv:
        pd.DataFrame(rows).to_csv(OUT.parent / "evidence_stats.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
