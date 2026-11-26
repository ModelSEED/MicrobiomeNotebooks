#!/usr/bin/env python
"""Split / tabulate the annot57 outputs (collected from poplar) into per-MAG coordinate tables and coverage counts.

Inputs (annotations dir A):
  scripts/mags57.tsv                         mag_id, set, kofam_name, ... (kofam_name = gene-ID prefix of the Prodigal genes)
  prodigal/<mag_id>.genes.tsv                KOfam-track Prodigal genes (-p meta -g 11), the eggNOG input
  eggnog/_pooled/g{1,2}/annot57_g*.emapper.{annotations,hits,seed_orthologs}
  bakta/<mag_id>/<mag_id>.json               Bakta 1.12.0 (--keep-contig-headers)
  kofam/<mag_id>.tsv.gz, kofam/sec_readthrough_57.tsv
  logs/bakta_<mag_id>.time, logs/emapper_g*.time
Outputs:
  eggnog/<mag_id>/<mag_id>.emapper.{annotations,hits,seed_orthologs}   (rows of that MAG, original format)
  eggnog/<mag_id>/<mag_id>.eggnog_genes.tsv  every Prodigal gene + coordinates + eggNOG fields ('-' -> empty)
  bakta/<mag_id>/<mag_id>.bakta_genes.tsv    every Bakta feature + coordinates (original contig IDs) + KO/EC/evidence
  coverage_per_genome.tsv, runtimes.tsv
usage: summarize.py <annotations dir> <modelling_species.tsv>
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

A = Path(sys.argv[1])
species = pd.read_csv(sys.argv[2], sep="\t")
species["mag_id"] = species.rep_genome.str.split("|").str[1]
mags = pd.read_csv(A / "scripts/mags57.tsv", sep="\t")
k2m = dict(zip(mags.kofam_name, mags.mag_id))
assert len(mags) == 57 and mags.mag_id.isin(species.mag_id).all()


def blank(x):
    return "" if pd.isna(x) or x == "-" else x


# ---------------- eggNOG ----------------
EGG_COLS = None
ann, hits, seeds = [], [], []
for g in (1, 2):
    d = A / f"eggnog/_pooled/g{g}"
    p = d / f"annot57_g{g}.emapper.annotations"
    with open(p) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                assert EGG_COLS in (None, cols)
                EGG_COLS = cols
                continue
            ann.append(line.rstrip("\n").split("\t"))
    for kind, store in (("hits", hits), ("seed_orthologs", seeds)):
        with open(d / f"annot57_g{g}.emapper.{kind}") as fh:
            for line in fh:
                if line.startswith("#qseqid"):
                    SEED_HDR = line
                if line.startswith("#"):
                    continue
                store.append(line)
# .emapper.hits has no header: DIAMOND --outfmt 6 as run by emapper (see the '##' diamond command in seed_orthologs)
HITS_HDR = "#qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovhsp\tscovhsp\n"
ann = pd.DataFrame(ann, columns=EGG_COLS)
ann["mag_id"] = ann["query"].str.split("|").str[0].map(k2m)
assert ann.mag_id.notna().all() and not ann["query"].duplicated().any()
hits_by = {}
for line in hits:
    hits_by.setdefault(k2m[line.split("|", 1)[0]], []).append(line)
seeds_by = {}
for line in seeds:
    seeds_by.setdefault(k2m[line.split("|", 1)[0]], []).append(line)

# ---------------- per MAG ----------------
rows = []
for m in mags.itertuples():
    mid = m.mag_id
    genes = pd.read_csv(A / f"prodigal/{mid}.genes.tsv", sep="\t", dtype={"partial": str, "rbs_motif": str})
    # eggNOG split + coordinate table
    od = A / f"eggnog/{mid}"
    od.mkdir(parents=True, exist_ok=True)
    a = ann[ann.mag_id == mid].drop(columns="mag_id")
    assert a["query"].isin(genes.gene).all()
    with open(od / f"{mid}.emapper.annotations", "w") as fo:
        fo.write("#" + "\t".join(EGG_COLS) + "\n")
        a.to_csv(fo, sep="\t", index=False, header=False)
    for kind, src, hdr in (("hits", hits_by, HITS_HDR), ("seed_orthologs", seeds_by, SEED_HDR)):
        with open(od / f"{mid}.emapper.{kind}", "w") as fo:
            fo.write(hdr + "".join(src.get(mid, [])))
    e = a.rename(columns={"query": "gene"}).copy()
    for c in e.columns:
        if c != "gene":
            e[c] = e[c].map(blank)
    e["ko"] = e.KEGG_ko.str.replace("ko:", "", regex=False)
    eg = genes[["gene", "contig", "contig_len", "start", "end", "strand", "partial", "aa_len"]].merge(e, on="gene", how="left")
    eg.to_csv(od / f"{mid}.eggnog_genes.tsv", sep="\t", index=False)
    has = lambda c: eg[c].fillna("").ne("")
    egg_prod = has("Description") | has("Preferred_name")

    # Bakta coordinate table from the JSON
    bj = A / f"bakta/{mid}/{mid}.json"
    b = None
    if bj.exists():
        J = json.load(open(bj))
        seq_ids = {s["id"] for s in J["sequences"]}
        fasta_ids = {l[1:].split()[0] for l in open(A.parent.parent / ("work/rollout_mags_kept" if m.set == "rollout"
                                                                         else "work/coasm_mags_kept") / f"{mid}.fa")
                     if l.startswith(">")}
        assert seq_ids <= fasta_ids, (mid, list(seq_ids - fasta_ids)[:3])  # --keep-contig-headers took effect
        br = []
        for f in J["features"]:
            xr = f.get("db_xrefs") or []
            if f.get("ups") or f.get("ips"):
                ev = "UPS/IPS"
            elif f.get("psc") and f["psc"].get("valid"):
                ev = "PSC"
            elif f.get("pscc"):
                ev = "PSCC"
            elif f["type"] in ("cds", "sorf"):
                ev = "none"
            else:
                ev = ""
            exc = f.get("exception")
            br.append(dict(
                locus_tag=f.get("locus", ""), type=f["type"], contig=f["sequence"], start=f["start"], end=f["stop"],
                strand=f["strand"], aa_len=len(f["aa"]) if f.get("aa") else "", gene=f.get("gene") or "",
                product=f.get("product") or "", hypothetical=bool(f.get("hypothetical")),
                pseudogene=bool(f.get("pseudogene")), truncated=f.get("truncated") or "", evidence=ev,
                ko=";".join(x[5:] for x in xr if x.startswith("KEGG:")),
                ec=";".join(x[3:] for x in xr if x.startswith("EC:")),
                go=";".join(x for x in xr if x.startswith("GO:")),
                cog=";".join(x[4:] for x in xr if x.startswith("COG:")),
                uniref=";".join(x[7:] for x in xr if x.startswith("UniRef:")),
                pfam=";".join(x[5:] for x in xr if x.startswith("PFAM:")),
                exception=f"{exc['type']}@aa{exc.get('codon_position', '')}" if exc else "",
                db_xrefs=";".join(xr)))
        b = pd.DataFrame(br)
        b.to_csv(A / f"bakta/{mid}/{mid}.bakta_genes.tsv", sep="\t", index=False)

    # KOfam
    k = pd.read_csv(A / f"kofam/{mid}.tsv.gz", sep="\t", dtype=str).fillna("")
    assert len(k) == len(genes) and (k.gene.values == genes.gene.values).all()
    sec = pd.read_csv(A / "kofam/sec_readthrough_57.tsv", sep="\t")
    sec = sec[(sec.mag_id == mid) & sec.accepted]
    kofam_ko = k.kofam_id.ne("")

    r = dict(mag_id=mid, set=m.set, genus=species.set_index("mag_id").genus[mid],
             prodigal_cds=len(genes),
             kofam_ko=int(kofam_ko.sum()),
             kofam_ec=int(k.kofam_definition.str.contains(r"\[EC:", regex=True).sum()),
             kofam_product=int(k.kofam_definition.ne("").sum()),
             kofam_sec_accepted=len(sec), kofam_sec_new=int((sec.status == "sec_readthrough").sum()),
             eggnog_hit=int(has("seed_ortholog").sum()),
             eggnog_ko=int(has("ko").sum()), eggnog_ec=int(has("EC").sum()),
             eggnog_bigg=int(has("BiGG_Reaction").sum()), eggnog_kegg_reaction=int(has("KEGG_Reaction").sum()),
             eggnog_product=int(egg_prod.sum()),
             kofam_or_eggnog_ko=int((kofam_ko.values | has("ko").values).sum()))
    if b is not None:
        c = b[b.type == "cds"]
        r.update(bakta_cds=len(c), bakta_sorf=int((b.type == "sorf").sum()), bakta_pseudo=int(c.pseudogene.sum()),
                 bakta_ko=int(c.ko.ne("").sum()), bakta_ec=int(c.ec.ne("").sum()),
                 bakta_product=int((~c.hypothetical & c["product"].ne("") &
                                    c["product"].ne("hypothetical protein")).sum()),
                 bakta_ups_ips=int((c.evidence == "UPS/IPS").sum()), bakta_psc=int((c.evidence == "PSC").sum()),
                 bakta_pscc=int((c.evidence == "PSCC").sum()), bakta_nohit=int((c.evidence == "none").sum()),
                 bakta_sec_cds=int(c.exception.str.startswith("selenocysteine").sum()),
                 bakta_trna=int((b.type == "tRNA").sum()), bakta_rrna=int((b.type == "rRNA").sum()))
    rows.append(r)

cov = pd.DataFrame(rows)


# ---------------- run times (GNU time -v) ----------------
def gtime(p):
    t = {}
    if not p.exists():
        return t
    for line in open(p):
        line = line.strip()
        if line.startswith("Elapsed (wall clock)"):
            v = line.split(": ", 1)[1].split(":")
            v = [float(x) for x in v]
            t["wall_min"] = round((v[0] * 3600 + v[1] * 60 + v[2]) / 60 if len(v) == 3 else (v[0] * 60 + v[1]) / 60, 1)
        elif line.startswith("User time"):
            t["user_s"] = float(line.split(": ")[1])
        elif line.startswith("System time"):
            t["sys_s"] = float(line.split(": ")[1])
        elif line.startswith("Maximum resident"):
            t["max_rss_gb"] = round(int(line.split(": ")[1]) / 1024 ** 2, 1)
        elif line.startswith("Exit status"):
            t["exit"] = int(line.split(": ")[1])
    if "user_s" in t:
        t["cpu_min"] = round((t.pop("user_s") + t.pop("sys_s")) / 60, 1)
    return t


rt = []
for mid in mags.mag_id:
    rt.append(dict(job=f"bakta {mid}", **gtime(A / f"logs/bakta_{mid}.time")))
for g in (1, 2):
    rt.append(dict(job=f"emapper g{g}", **gtime(A / f"logs/emapper_g{g}.time")))
rt = pd.DataFrame(rt)
rt.to_csv(A / "runtimes.tsv", sep="\t", index=False)
bt = rt[rt.job.str.startswith("bakta")].assign(mag_id=lambda x: x.job.str[6:])
cov = cov.merge(bt[["mag_id", "wall_min", "cpu_min", "max_rss_gb"]].rename(
    columns={"wall_min": "bakta_wall_min", "cpu_min": "bakta_cpu_min", "max_rss_gb": "bakta_max_rss_gb"}), on="mag_id", how="left")
cov.to_csv(A / "coverage_per_genome.tsv", sep="\t", index=False)

# gene-caller concordance, Bakta CDS vs KOfam-track Prodigal genes (same contig IDs): shared stop codon / identical coords
conc = []
for mid in mags.mag_id:
    bp = A / f"bakta/{mid}/{mid}.bakta_genes.tsv"
    if not bp.exists():
        continue
    b = pd.read_csv(bp, sep="\t", usecols=["type", "contig", "start", "end", "strand"]).query("type == 'cds'")
    p = pd.read_csv(A / f"prodigal/{mid}.genes.tsv", sep="\t", usecols=["contig", "start", "end", "strand"])
    stop = lambda x: set(zip(x.contig, x.strand, x.end.where(x.strand == "+", x.start)))
    bs, ps = stop(b), stop(p)
    same = set(zip(b.contig, b.strand, b.start, b.end)) & set(zip(p.contig, p.strand, p.start, p.end))
    conc.append(dict(mag_id=mid, bakta_cds=len(b), prodigal_genes=len(p), shared_stop=len(bs & ps),
                     frac_bakta_with_prodigal_stop=round(len(bs & ps) / len(bs), 4),
                     frac_prodigal_with_bakta_stop=round(len(bs & ps) / len(ps), 4),
                     identical_coords=len(same)))
if conc:
    pd.DataFrame(conc).to_csv(A / "gene_caller_concordance.tsv", sep="\t", index=False)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
num = cov.select_dtypes("number")
print(num.describe().T[["min", "50%", "max"]].round(1).to_string())
frac = pd.DataFrame({c: cov[c] / cov.prodigal_cds for c in ["kofam_ko", "kofam_ec", "eggnog_hit", "eggnog_ko", "eggnog_ec",
                                                             "eggnog_bigg", "eggnog_product", "kofam_or_eggnog_ko"]})
if "bakta_cds" in cov:
    for c in ["bakta_ko", "bakta_ec", "bakta_product", "bakta_ups_ips", "bakta_psc", "bakta_pscc", "bakta_nohit"]:
        frac[c] = cov[c] / cov.bakta_cds
print((100 * frac).describe().T[["min", "50%", "max"]].round(1).to_string())
print(rt.to_string())
print("totals:", cov[["prodigal_cds", "kofam_ko", "eggnog_ko", "eggnog_bigg"] + (["bakta_cds", "bakta_ko"] if "bakta_cds" in cov else [])].sum().to_dict())
