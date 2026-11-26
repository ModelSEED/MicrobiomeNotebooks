#!/usr/bin/env python
"""Step 5b (v2): one row per DIAMOND-DB representative (v1 05b_rep_info.py, extended).
Taxonomy of members:
  * salterns members on rollout kept MAGs / co-assembly MAGs: GTDB R232 (meta/mags_three_way.tsv), as v1;
  * KEGG members: GTDB R232 lineage of the organism's genome (meta/kegg_org_gtdb.tsv, see 11b/11b2) -- v2 change, v1
    used the KEGG organism name; when the genome is not in R232: KEGG br08601 domain/phylum + organism genus (as v1);
  * UniProt members: UniProt (NCBI) lineage, as v1.
Representative genus (genus_source):
  'salterns_MAG' majority GTDB genus of binned salterns members; 'external' majority genus of KEGG/UniProt members;
  'transferred' salterns-only unbinned cluster of a counted class (target, marker, decoy_reported -- v1 did targets and
  reported classes only): best DIAMOND blastp hit (--sensitive) among (a) external proteins of counted labels,
  (b) binned salterns members of counted labels, (c) all proteins of the 1,324 first_pass MAGs (GTDB R232);
  genus kept at any identity as in v1 (transfer_pident recorded; the read tables bin by read identity).
Also: order (for the methanogen-lineage flag), genus_source_name (the reference's own name), methanogen_lineage /
anme flags, label lengths (refs/label_lengths.tsv, v1 rule). Writes refs/rep_info.tsv.gz.
Methanogen lineage = Archaea and (order, family or genus name starts with 'Methano', or genus Methermicoccus), minus
the methanotrophic/alkanotrophic ANME lineages (order Methanophagales; families Methanoperedenaceae,
Methanocomedenaceae, Methanogasteraceae; genera Methanoperedens, Methanocomedens, Methanogaster, Methanophaga,
Methanospirare), which are flagged 'anme'.
"""
import os, subprocess, gzip, re
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
DM = "/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond"
A = "/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin"
T = os.environ.get("T", "32")
os.chdir(W)
meta = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
m3 = pd.read_csv(f"{V1}/meta/mags_three_way.tsv", sep="\t", dtype=str).fillna("")
org = pd.read_csv("meta/kegg_org_gtdb.tsv", sep="\t", dtype=str).fillna("").set_index("org")


def ranks(lin):
    return dict(x.split("__", 1) for x in lin.split(";") if "__" in x)


m3r = {g: ranks(c) for g, c in zip(m3.genome, m3.classification)}
meta["order"] = ""
meta["genus_source_name"] = meta.genus
# salterns binned: order from GTDB
sb = (meta.source == "salterns") & (meta.mag != "")
meta.loc[sb, "order"] = [m3r.get(g, {}).get("o", "") for g in meta.mag[sb]]
# KEGG: GTDB R232 where the genome is in R232
kg = meta.source == "kegg"
codes = meta.pid[kg].str.split("|").str[2].str.split(":").str[0]
gl = codes.map(org.gtdb_lineage).fillna("")
has = gl != ""
idx = gl[has].index
rk = [ranks(x) for x in gl[has]]
for col, key in (("domain", "d"), ("phylum", "p"), ("order", "o"), ("family", "f"), ("genus", "g")):
    meta.loc[idx, col] = [r.get(key, "") for r in rk]
# UniProt: order from the lineage string
up = meta.source == "uniprot"
lin = pd.read_csv("refs/ext_meta.tsv", sep="\t", dtype=str).fillna("").set_index("pid").lineage
meta.loc[up, "order"] = [dict((m.group(2), m.group(1).strip()) for m in re.finditer(r"([^,]+?) \(([a-z ]+)\)", lin.get(p, ""))).get("order", "")
                         for p in meta.pid[up]]
meta.to_csv("refs/labeled_meta_tax.tsv.gz", sep="\t", index=False)

c99 = pd.read_csv("refs/clu/clusters99.tsv", sep="\t", header=None, names=["label", "rep", "member"], dtype=str)
m = c99.merge(meta.drop(columns=["label"]), left_on="member", right_on="pid", how="left")
assert m.pid.notna().all()


def majority(sub, col):
    """per rep: most frequent non-empty value of col among the rows of sub (ties: alphabetical)"""
    t = sub[sub[col] != ""].groupby(["rep", col]).size().reset_index(name="n")
    t = t.sort_values(["rep", "n", col], ascending=[True, False, True]).drop_duplicates("rep")
    return t.set_index("rep")[col]


def joined(sub, col):
    t = sub[sub[col] != ""][["rep", col]].drop_duplicates().sort_values(["rep", col])
    return t.groupby("rep")[col].agg(",".join)


# vectorized (v1 looped over representatives; too slow at the v2 size)
m["is_sal"] = m.source == "salterns"
m["is_bin"] = m.is_sal & (m.mag != "")
m["is_kegg"] = m.source == "kegg"
m["is_up"] = m.source == "uniprot"
sal, binned, ext = m[m.is_sal], m[m.is_bin], m[~m.is_sal]
grp = m.groupby("rep", sort=False)
ri = pd.DataFrame({"label": grp.label.first(), "role": grp.role.first(), "n_members": grp.size(),
                   "n_salterns": grp.is_sal.sum(), "n_salterns_binned": grp.is_bin.sum(),
                   "n_kegg": grp.is_kegg.sum(), "n_uniprot": grp.is_up.sum()})
ri["samples"] = joined(sal, "sample")
ri["mags"] = joined(binned, "mag")
ri["species_clusters"] = joined(binned, "species_cluster")
for col in ("genus", "family", "order", "phylum", "domain"):
    ri[f"{col}_mag"] = majority(binned, col)
    ri[f"{col}_ext"] = majority(ext, col)
ri["organism_ext"] = majority(ext, "organism")
e2 = ext.merge(ri.genus_ext.rename("gx"), left_on="rep", right_index=True)
ri["genus_source_name_ext"] = majority(e2[e2.genus == e2.gx], "genus_source_name")
ri = ri.fillna("").reset_index().rename(columns={"index": "rep"})
for c in ("n_members", "n_salterns", "n_salterns_binned", "n_kegg", "n_uniprot"):
    ri[c] = ri[c].astype(int)
print("reps", len(ri), flush=True)

# transfer
counted = ri.role.isin(["target", "marker", "decoy_reported"])
need = ri[counted & (ri.n_salterns_binned == 0) & (ri.n_kegg + ri.n_uniprot == 0)]
refs = meta[meta.role.isin(["target", "marker", "decoy_reported"]) & ((meta.source != "salterns") | (meta.mag != ""))]
os.makedirs("refs/tx", exist_ok=True)
open("refs/tx/q.ids", "w").write("\n".join(need.rep) + "\n")
open("refs/tx/r.ids", "w").write("\n".join(refs.pid) + "\n")
subprocess.run(f"{A}/seqkit grep -f refs/tx/q.ids db/survey.faa > refs/tx/q.faa && "
               f"{A}/seqkit grep -f refs/tx/r.ids refs/labeled.faa > refs/tx/r.faa && "
               f"sed 's/ # .*//' genes_firstpass/firstpass.faa >> refs/tx/r.faa && "
               f"nice -n 10 {DM} makedb --in refs/tx/r.faa -d refs/tx/r --threads {T} --quiet && "
               f"nice -n 10 {DM} blastp -q refs/tx/q.faa -d refs/tx/r -o refs/tx/hits.tsv -k 1 --sensitive --threads {T} --quiet "
               f"-t tmp --outfmt 6 qseqid sseqid pident length evalue bitscore", shell=True, check=True)
h = pd.read_csv("refs/tx/hits.tsv", sep="\t", header=None, names="q s pid len e bits".split(), dtype={"q": str, "s": str})
mi = meta.set_index("pid")
fpg = {}
for gname, r in m3r.items():
    if gname.startswith("first_pass|"):
        fpg[gname.split("|", 1)[1]] = r
tr = {}
for q, s, p in zip(h.q, h.s, h.pid):
    if s.startswith("fp|"):
        r = fpg.get(s.split("|")[1], {})
        tr[q] = (r.get("g", ""), r.get("f", ""), r.get("o", ""), r.get("p", ""), r.get("d", ""), p, "first_pass_MAG", r.get("g", ""))
    else:
        x = mi.loc[s]
        tr[q] = (x.genus, x.family, x.order, x.phylum, x.domain, p, x.label, x.genus_source_name)
E = ("",) * 8
for i, col in enumerate(["genus_transferred", "family_transferred", "order_transferred", "phylum_transferred",
                         "domain_transferred", "transfer_pident", "transfer_hit_label", "genus_source_name_transferred"]):
    ri[col] = [tr.get(r, E)[i] for r in ri.rep]


def pick(r):
    if r.n_salterns_binned > 0:
        return r.genus_mag, r.family_mag, r.order_mag, r.phylum_mag, r.domain_mag, "salterns_MAG", r.genus_mag
    if r.n_kegg + r.n_uniprot > 0:
        return r.genus_ext, r.family_ext, r.order_ext, r.phylum_ext, r.domain_ext, "external", r.genus_source_name_ext
    if r.genus_transferred or r.family_transferred or r.domain_transferred:
        return (r.genus_transferred, r.family_transferred, r.order_transferred, r.phylum_transferred, r.domain_transferred,
                "transferred", r.genus_source_name_transferred)
    return "", "", "", "", "", "none", ""


p = ri.apply(pick, axis=1, result_type="expand")
ri[["genus", "family", "order", "phylum", "domain", "genus_source", "genus_source_name"]] = p
ANME_O = {"Methanophagales"}
ANME_F = {"Methanoperedenaceae", "Methanocomedenaceae", "Methanogasteraceae"}
ANME_G = {"Methanoperedens", "Methanocomedens", "Methanogaster", "Methanophaga", "Methanospirare"}


def lineage_class(r):
    if r.domain != "Archaea":
        return ""
    if r.order in ANME_O or r.family in ANME_F or r.genus in ANME_G or "ANME" in f"{r.order}{r.family}{r.genus}":
        return "anme"
    if any(x.startswith("Methano") for x in (r.order, r.family, r.genus)) or r.genus == "Methermicoccus":
        return "methanogen"
    return "other_archaea"


ri["archaeal_lineage"] = ri.apply(lineage_class, axis=1)
ri.to_csv("refs/rep_info.tsv.gz", sep="\t", index=False)

# label lengths (v1 rule)
ln = {}
for line in open("refs/labeled.faa"):
    if line.startswith(">"):
        cur = line[1:].strip()
    else:
        ln[cur] = len(line.strip())
meta["aa_len"] = meta.pid.map(ln)
partial = {}
for base in sorted(set(meta[meta.source == "salterns"].set)):
    for line in open(f"{V1}/genes/{base}/{base}.faa"):
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
    else:
        a = g.aa_len.sort_values()
        use, src = a.iloc[len(a) // 2:], "all_longest_half"
    out.append(dict(label=lab, role=g.role.iloc[0], median_aa=float(use.median()), n_used=len(use), length_source=src,
                    n_kegg=len(k), n_salterns_complete=len(s)))
pd.DataFrame(out).to_csv("refs/label_lengths.tsv", sep="\t", index=False)
print(ri.groupby(["role", "genus_source"]).size().unstack(fill_value=0).to_string())
print(ri[ri.role.isin(["target", "marker"])].groupby("archaeal_lineage").size().to_string())
