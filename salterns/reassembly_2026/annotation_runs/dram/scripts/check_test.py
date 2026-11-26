"""Sanity checks on the one-MAG DRAM test: python check_test.py ANNOT_DIR DISTILL_DIR. Exit 1 on failure."""
import sys, os, pandas as pd
a_dir, d_dir = sys.argv[1], sys.argv[2]
a = pd.read_csv(os.path.join(a_dir, "annotations.tsv"), sep="\t", index_col=0)
bad = []
print("genes", len(a), "columns", list(a.columns))
for c in ("ko_id", "kegg_hit", "pfam_hits", "cazy_ids", "peptidase_id", "bin_taxonomy", "bin_completeness"):
    if c not in a.columns:
        bad.append(f"missing column {c}")
if "bin_taxonomy" in a and not str(a.bin_taxonomy.iloc[0]).startswith("d__"):
    bad.append(f"bin_taxonomy not GTDB: {a.bin_taxonomy.iloc[0]}")
if "bin_completeness" in a:
    print("bin_completeness", a.bin_completeness.iloc[0], "bin_contamination", a.bin_contamination.iloc[0], "taxonomy", a.bin_taxonomy.iloc[0])
    if not a.bin_completeness.iloc[0] > 0: bad.append("bin_completeness is 0 (name mismatch with checkm table)")
if "ko_id" in a:
    f = a.ko_id.notna().mean(); print(f"fraction with KO {f:.3f}")
    if f < 0.2: bad.append(f"KO fraction {f:.3f} < 0.2")
    kos = set(",".join(a.ko_id.dropna()).split(","))
    for k, g in [("K00399", "mcrA"), ("K16178", "mtbB"), ("K16176", "mtmB"), ("K14083", "mttB"), ("K14082", "mtbA"), ("K16180", "pylB")]:
        print(f"  {g} {k}: {'yes' if k in kos else 'no'}")
    if "K00399" not in kos: bad.append("no mcrA in Methanolobus test MAG")
for c in ("pfam_hits", "cazy_hits", "peptidase_hit"):
    if c in a:
        s = a[c].dropna().astype(str)
        print(f"{c}: {len(s)} genes, e.g. {s.iloc[0][:90] if len(s) else '-'}")
        if s.str.contains("<html|<!DOCTYPE", case=False).any(): bad.append(f"HTML in {c}")
if "pfam_hits" in a and a.pfam_hits.notna().mean() < 0.2: bad.append("Pfam fraction < 0.2")
for f in ("trnas.tsv", "rrnas.tsv", "genes.faa", "genes.gff", "scaffolds.fna"):
    p = os.path.join(a_dir, f); print(f, os.path.exists(p) and os.path.getsize(p))
x = pd.read_excel(os.path.join(d_dir, "metabolism_summary.xlsx"), sheet_name=None)
print("distill sheets", list(x))
if "Salterns methane" not in x:
    bad.append("no 'Salterns methane' sheet")
else:
    s = x["Salterns methane"]; gcol = [c for c in s.columns if c not in ("gene_id", "gene_description", "module", "header", "subheader")]
    hit = s[s[gcol[0]] > 0]
    print("Salterns methane rows", len(s), "present:", ", ".join(f"{r.subheader}={r[gcol[0]]}" for _, r in hit.iterrows()))
    if len(s) != 65: bad.append(f"Salterns methane has {len(s)} rows, expected 65")
for f in ("genome_stats.tsv", "product.html", "product.tsv"):
    if not os.path.exists(os.path.join(d_dir, f)): bad.append(f"distill missing {f}")
print("TEST_CHECK", "FAIL: " + "; ".join(bad) if bad else "OK")
sys.exit(1 if bad else 0)
