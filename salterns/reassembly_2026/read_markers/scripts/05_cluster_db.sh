#!/usr/bin/env bash
# Step 5: resolve salterns K10944 proteins (pmoA vs amoA vs other CuMMO) by best hit to the labelled external CuMMO
# references; cluster every label separately (MMseqs2; targets/markers/reported classes at 99 % identity, decoys at 95 %,
# coverage >= 80 % of the shorter sequence); also 90 % clusters (all labels) for the hold-out validation;
# build the DIAMOND DB of the representatives.
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers
E=/scratch1/afreiburger/salterns_2026/envs
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
MM=$E/dram15/bin/mmseqs
DM=$E/bt_dastool/bin/diamond
PY=$A/python3
export TMPDIR=$W/tmp
T=${T:-16}
cd $W
mkdir -p refs/clu db
# 5a. K10944 resolution
awk -F'\t' 'NR>1 && $2=="K10944_unresolved"{print $1}' refs/labeled_meta.tsv > refs/clu/k10944_q.ids
awk -F'\t' 'NR>1 && $4!="salterns" && ($2=="pmoA"||$2=="amoA"||$2=="cuMMO_other"){print $1}' refs/labeled_meta.tsv > refs/clu/k10944_r.ids
$A/seqkit grep -f refs/clu/k10944_q.ids refs/labeled.faa > refs/clu/k10944_q.faa
$A/seqkit grep -f refs/clu/k10944_r.ids refs/labeled.faa > refs/clu/k10944_r.faa
nice -n 10 $DM makedb --in refs/clu/k10944_r.faa -d refs/clu/k10944_r --quiet
nice -n 10 $DM blastp -q refs/clu/k10944_q.faa -d refs/clu/k10944_r -o refs/clu/k10944_hits.tsv -k 1 --sensitive \
  --threads $T --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore || true
$PY -I - <<'EOF'
import pandas as pd
m = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
h = pd.read_csv("refs/clu/k10944_hits.tsv", sep="\t", header=None, names="q s pid len e bits".split())
lab = dict(zip(m.pid, m.label))
res = {q: lab[s] for q, s in zip(h.q, h.s)}
sel = m.label == "K10944_unresolved"
m.loc[sel, "how"] = m.loc[sel, "how"] + "+K10944_by_best_external_hit"
m.loc[sel, "label"] = [res.get(p, "cuMMO_other") for p in m.pid[sel]]
m.to_csv("refs/labeled_meta.tsv", sep="\t", index=False)
h.to_csv("refs/clu/k10944_resolution.tsv", sep="\t", index=False)
print("K10944 salterns proteins resolved:", sel.sum(), m.loc[sel, "label"].value_counts().to_dict())
EOF
# 5b. per-label clustering
$PY -I - <<'EOF'
import pandas as pd
m = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
p = pd.read_csv("meta/ko_panel.tsv", sep="\t", dtype=str).fillna("")
role = dict(zip(p.gene, p.role))
extra = {"amoA": "decoy_reported", "cuMMO_other": "decoy_reported", "mttB_family_other": "target", "euk_marker": "decoy",
         "SDIMO_other": "decoy", "GRE_other": "decoy", "Rieske_other": "decoy", "hepD": "decoy_reported"}
role.update(extra)
m["role"] = [role.get(l, "decoy" if l.endswith("_like") else "?") for l in m.label]
assert (m.role != "?").all(), m.label[m.role == "?"].unique()
m.to_csv("refs/labeled_meta.tsv", sep="\t", index=False)
with open("refs/clu/labels.tsv", "w") as fh:
    for lab, g in m.groupby("label"):
        open(f"refs/clu/{lab}.ids", "w").write("\n".join(g.pid) + "\n")
        fh.write(f"{lab}\t{g.role.iloc[0]}\t{len(g)}\n")
EOF
: > db/survey.faa; : > refs/clu/clusters99.tsv; : > refs/clu/clusters90.tsv
while IFS=$'\t' read lab role n; do
  id=0.99; [ "$role" = decoy ] && id=0.95
  $A/seqkit grep -f refs/clu/$lab.ids refs/labeled.faa > refs/clu/$lab.faa
  for pass in main c90; do
    pid=$id; [ $pass = c90 ] && pid=0.90
    rm -rf $TMPDIR/mm_$lab
    nice -n 10 $MM easy-cluster refs/clu/$lab.faa refs/clu/${lab}_$pass $TMPDIR/mm_$lab --min-seq-id $pid -c 0.8 --cov-mode 0 \
      --threads $T -v 1 > /dev/null
    rm -rf $TMPDIR/mm_$lab
  done
  awk -v l=$lab '{print l"\t"$0}' refs/clu/${lab}_main_cluster.tsv >> refs/clu/clusters99.tsv
  awk -v l=$lab '{print l"\t"$0}' refs/clu/${lab}_c90_cluster.tsv >> refs/clu/clusters90.tsv
  sed "s/^>\(.*\)$/>\1/" refs/clu/${lab}_main_rep_seq.fasta >> db/survey.faa
  rm -f refs/clu/${lab}_main_all_seqs.fasta refs/clu/${lab}_c90_all_seqs.fasta refs/clu/${lab}_c90_rep_seq.fasta
  echo "$lab $role $n -> $(grep -c '^>' refs/clu/${lab}_main_rep_seq.fasta) reps"
done < refs/clu/labels.tsv
nice -n 10 $DM makedb --in db/survey.faa -d db/survey --threads $T --quiet
md5sum db/survey.faa db/survey.dmnd > db/survey.md5
echo "$(date -u +%FT%TZ) survey DB: $(grep -c '^>' db/survey.faa) representatives"
