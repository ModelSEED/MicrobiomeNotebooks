#!/usr/bin/env bash
# Reviewer check: DRAM 1.5 distill (incl. product.html with altair 6.2.2) on a tiny synthetic annotations.tsv.
source /scratch1/afreiburger/salterns_2026/dram/scripts/dram_env.sh
S=$W/review/synth; rm -rf $S; mkdir -p $S; cd $S
python -I - <<'PY'
import pandas as pd
rows = []
kos = {"gA": ["K00399", "K00401", "K00402", "K16178", "K14083", "K16176", "K00200"],
       "gB": ["K06163", "K06164", "K06165", "K06166", "K14083", "K21578", "K02586"]}
for g, ks in kos.items():
    for i, k in enumerate(ks, 1):
        rows.append({"": f"{g}_c1_{i}", "fasta": g, "scaffold": f"{g}_c1", "gene_position": i, "start_position": i*1000,
                     "end_position": i*1000+900, "strandedness": 1, "rank": "D", "ko_id": k, "kegg_hit": f"{k} test",
                     "bin_taxonomy": "d__Archaea;p__x" if g == "gA" else "d__Bacteria;p__y", "bin_completeness": 95.0, "bin_contamination": 1.0})
pd.DataFrame(rows).set_index("").to_csv("annotations.tsv", sep="\t")
PY
nice -n 10 DRAM.py distill -i annotations.tsv -o out --custom_distillate "$SHEET" > distill.log 2>&1; echo "distill rc=$?"
ls -la out
