#!/usr/bin/env python
"""KEGG organism code -> assembly accession (11b) -> GTDB R232 lineage (release232/taxonomy/gtdb_taxonomy.tsv, which
lists every R232 genome, representatives and members). Matched on the accession number without version; a GCA/GCF
pair shares the number. Writes meta/kegg_org_gtdb.tsv (org, tnum, accession, taxid, kegg_name, gtdb_accession,
gtdb_lineage, gtdb_genus, match).
Usage: python -I 11b2_map_org_gtdb.py <kegg_org_assembly.tsv> <gtdb_taxonomy.tsv> <out.tsv>"""
import sys, re
import pandas as pd
ka, gt, out = sys.argv[1:4]
k = pd.read_csv(ka, sep="\t", header=None, names=["org", "tnum", "accession", "taxid", "kegg_name"], dtype=str).fillna("")
k = k.drop_duplicates("org", keep="last")
g = {}
for line in open(gt):
    acc, lin = line.rstrip("\n").split("\t")
    num = acc.split("_", 2)[2].split(".")[0]          # RS_GCF_000016525.1 -> 000016525
    pref = acc.split("_")[1]                             # GCF / GCA
    g.setdefault(num, {})[pref] = (acc, lin)
rows = []
for _, r in k.iterrows():
    m = re.match(r"(GC[AF])_(\d+)\.\d+", r.accession)
    hit, how = None, "no_accession"
    if m:
        d = g.get(m.group(2), {})
        if m.group(1) in d:
            hit, how = d[m.group(1)], "same_prefix"
        elif d:
            hit, how = list(d.values())[0], "GCA_GCF_pair"
        else:
            how = "not_in_GTDB_R232"
    gen = re.search(r"g__([^;]*)", hit[1]).group(1) if hit else ""
    rows.append(dict(org=r.org, tnum=r.tnum, accession=r.accession, taxid=r.taxid, kegg_name=r.kegg_name,
                     gtdb_accession=hit[0] if hit else "", gtdb_lineage=hit[1] if hit else "", gtdb_genus=gen, match=how))
o = pd.DataFrame(rows)
o.to_csv(out, sep="\t", index=False)
print(o.match.value_counts().to_string())
kg = o.kegg_name.str.replace(r"^(Candidatus|uncultured) ", "", regex=True).str.split().str[0].str.strip("[]")
d = o[(o.gtdb_genus != "") & (kg != o.gtdb_genus)]
print("orgs whose GTDB genus differs from the KEGG name genus:", len(d), "of", (o.gtdb_genus != "").sum())
