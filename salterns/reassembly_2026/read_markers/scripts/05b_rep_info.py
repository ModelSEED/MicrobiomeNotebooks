#!/usr/bin/env python
"""Step 5b: one row per DIAMOND-DB representative (99 % / 95 % cluster):
label, role, members by source, salterns MAG members (MAG ids, species clusters, GTDB R232 genus), external genera,
and a genus for the representative:
  genus_source = 'salterns_MAG'  : majority GTDB genus of the binned salterns members (rollout kept MAGs / co-assembly MAGs)
               = 'external'      : majority genus of the KEGG / UniProt members
               = 'transferred'   : salterns-only, unbinned cluster -> genus of the best DIAMOND blastp hit among all
                                   external proteins of target / reported labels (pident recorded; genus kept at any pident,
                                   the read-level tables bin by it)
Also the median complete-gene length (aa) per label (KEGG members if >= 5; else complete salterns genes if >= 5; else
the longest half of all members, because UniProt sets contain PCR-amplicon fragments), used for RPKM and
copies per genome equivalent. Writes refs/rep_info.tsv.gz and refs/label_lengths.tsv.
"""
import os, subprocess, gzip, re
from collections import Counter
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers"
DM = "/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond"
os.chdir(W)
meta = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
c99 = pd.read_csv("refs/clu/clusters99.tsv", sep="\t", header=None, names=["label", "rep", "member"], dtype=str)
m = c99.merge(meta.drop(columns=["label"]), left_on="member", right_on="pid", how="left")
assert m.pid.notna().all()


def majority(xs):
    xs = [x for x in xs if x]
    return Counter(xs).most_common(1)[0][0] if xs else ""


rows = []
for rep, g in m.groupby("rep"):
    sal = g[g.source == "salterns"]
    binned = sal[sal.mag != ""]
    ext = g[g.source != "salterns"]
    rows.append(dict(rep=rep, label=g.label.iloc[0], role=g.role.iloc[0], n_members=len(g), n_salterns=len(sal),
                     n_salterns_binned=len(binned), n_kegg=(g.source == "kegg").sum(), n_uniprot=(g.source == "uniprot").sum(),
                     samples=",".join(sorted(set(sal["sample"]))), mags=",".join(sorted(set(binned.mag))),
                     species_clusters=",".join(sorted(set(binned.species_cluster) - {""})),
                     genus_mag=majority(binned.genus), family_mag=majority(binned.family), phylum_mag=majority(binned.phylum),
                     domain_mag=majority(binned.domain),
                     genus_ext=majority(ext.genus), family_ext=majority(ext.family), phylum_ext=majority(ext.phylum),
                     domain_ext=majority(ext.domain), organism_ext=majority(ext.organism)))
ri = pd.DataFrame(rows)

# transfer taxonomy to salterns-only unbinned reps of counted classes
counted = ri.role.isin(["target", "decoy_reported"])
need = ri[counted & (ri.n_salterns_binned == 0) & (ri.n_kegg + ri.n_uniprot == 0)]
ext_ids = meta[(meta.source != "salterns") & meta.role.isin(["target", "decoy_reported"])].pid
os.makedirs("refs/tx", exist_ok=True)
open("refs/tx/q.ids", "w").write("\n".join(need.rep) + "\n")
open("refs/tx/r.ids", "w").write("\n".join(ext_ids) + "\n")
A = "/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin"
subprocess.run(f"{A}/seqkit grep -f refs/tx/q.ids db/survey.faa > refs/tx/q.faa && "
               f"{A}/seqkit grep -f refs/tx/r.ids refs/labeled.faa > refs/tx/r.faa && "
               f"nice -n 10 {DM} makedb --in refs/tx/r.faa -d refs/tx/r --quiet && "
               f"nice -n 10 {DM} blastp -q refs/tx/q.faa -d refs/tx/r -o refs/tx/hits.tsv -k 1 --sensitive --threads 8 --quiet "
               f"--outfmt 6 qseqid sseqid pident length evalue bitscore", shell=True, check=True)
h = pd.read_csv("refs/tx/hits.tsv", sep="\t", header=None, names="q s pid len e bits".split(), dtype={"q": str, "s": str})
mi = meta.set_index("pid")
tr = {q: (mi.at[s, "genus"], mi.at[s, "family"], mi.at[s, "phylum"], mi.at[s, "domain"], p, mi.at[s, "label"])
      for q, s, p in zip(h.q, h.s, h.pid)}
ri["genus_transferred"] = [tr.get(r, ("",) * 6)[0] for r in ri.rep]
ri["family_transferred"] = [tr.get(r, ("",) * 6)[1] for r in ri.rep]
ri["phylum_transferred"] = [tr.get(r, ("",) * 6)[2] for r in ri.rep]
ri["domain_transferred"] = [tr.get(r, ("",) * 6)[3] for r in ri.rep]
ri["transfer_pident"] = [tr.get(r, ("",) * 6)[4] for r in ri.rep]
ri["transfer_hit_label"] = [tr.get(r, ("",) * 6)[5] for r in ri.rep]


def pick(r):
    if r.n_salterns_binned > 0 and r.genus_mag:
        return r.genus_mag, r.family_mag, r.phylum_mag, r.domain_mag, "salterns_MAG"
    if r.n_salterns_binned > 0:
        return "", r.family_mag, r.phylum_mag, r.domain_mag, "salterns_MAG"
    if r.n_kegg + r.n_uniprot > 0:
        return r.genus_ext, r.family_ext, r.phylum_ext, r.domain_ext, "external"
    if r.genus_transferred or r.family_transferred:
        return r.genus_transferred, r.family_transferred, r.phylum_transferred, r.domain_transferred, "transferred"
    return "", "", "", "", "none"


p = ri.apply(pick, axis=1, result_type="expand")
ri[["genus", "family", "phylum", "domain", "genus_source"]] = p
ri.to_csv("refs/rep_info.tsv.gz", sep="\t", index=False)

# label lengths
ln = {}
for line in open("refs/labeled.faa"):
    if line.startswith(">"):
        cur = line[1:].strip()
    else:
        ln[cur] = len(line.strip())
meta["aa_len"] = meta.pid.map(ln)
partial = {}
for base in sorted(set(meta[meta.source == "salterns"].set)):
    for line in open(f"genes/{base}/{base}.faa"):
        if line.startswith(">"):
            pid = line[1:].split(" # ")[0]
            mm = re.search(r"partial=(\d\d)", line)
            partial[pid] = mm.group(1) if mm else ""
meta["partial"] = meta.pid.map(partial).fillna("")
out = []
for lab, g in meta.groupby("label"):
    k = g[g.source == "kegg"].aa_len
    s = g[(g.source == "salterns") & (g.partial == "00")].aa_len
    if len(k) >= 5:
        use, src = k, "kegg"
    elif len(s) >= 5:
        use, src = s, "salterns_complete"
    else:   # UniProt-dominated classes contain PCR-amplicon fragments: median of the longest half
        a = g.aa_len.sort_values()
        use, src = a.iloc[len(a) // 2:], "all_longest_half"
    out.append(dict(label=lab, role=g.role.iloc[0], median_aa=float(use.median()), n_used=len(use), length_source=src,
                    n_kegg=len(k), n_salterns_complete=len(s)))
pd.DataFrame(out).to_csv("refs/label_lengths.tsv", sep="\t", index=False)
print(ri.groupby(["label", "genus_source"]).size().unstack(fill_value=0).to_string())
