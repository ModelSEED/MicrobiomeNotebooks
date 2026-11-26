"""NCBI taxonomy id for each of the 1,206 representatives, for TRANSYT (KBDLGenomeAnnotation TRANSYT_LOCAL needs a
tax_id; the KBDL genomes carry GTDB R232 taxonomy without one).
Rule: GTDB-Tk's closest reference (r232_closest_reference) -> its ncbi_taxid; otherwise the GTDB R232 species
representative whose gtdb_taxonomy shares the MAG's deepest named rank (first by accession) -> its ncbi_taxid, so the
NCBI lineage TRANSYT filters by is right down to that rank. Output: rep_taxids.tsv (mag, tax_id, method, source).
Usage: python3 make_taxids.py <mag_census.tsv> <reps.txt> <metadata dir>"""
import gzip, sys, csv
import pandas as pd
mc = pd.read_csv(sys.argv[1], sep="\t").set_index("mag")
reps = [l.strip() for l in open(sys.argv[2]) if l.strip()]
acc_tax, rep_by_prefix = {}, {}
for f in ("bac120_metadata_r232.tsv.gz", "ar53_metadata_r232.tsv.gz"):
    with gzip.open(f"{sys.argv[3]}/{f}", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a = r["accession"][3:]                      # RS_/GB_ prefix off
            acc_tax[a] = r["ncbi_taxid"]
            if r["gtdb_representative"] == "t":
                ranks = r["gtdb_taxonomy"].split(";")
                for k in range(1, len(ranks) + 1):
                    rep_by_prefix.setdefault(";".join(ranks[:k]), (a, r["ncbi_taxid"]))
out = []
for m in reps:
    row = mc.loc[m]; ref = str(row.r232_closest_reference)
    if ref not in ("nan", "") and ref[3:] in acc_tax:
        out.append((m, acc_tax[ref[3:]], "closest_reference", ref)); continue
    ranks = [x for x in str(row.r232_classification).split(";") if len(x) > 3]
    for k in range(len(ranks), 0, -1):
        hit = rep_by_prefix.get(";".join(ranks[:k]))
        if hit:
            out.append((m, hit[1], f"rank_{ranks[k-1][0]}", hit[0])); break
    else:
        out.append((m, "", "none", ""))
pd.DataFrame(out, columns=["mag", "tax_id", "method", "source"]).to_csv("rep_taxids.tsv", sep="\t", index=False)
d = pd.DataFrame(out, columns=["mag", "tax_id", "method", "source"])
print(len(d), d.method.value_counts().to_dict(), "empty:", (d.tax_id == "").sum())
