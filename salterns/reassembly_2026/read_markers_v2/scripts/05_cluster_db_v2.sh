#!/usr/bin/env bash
# Step 5 (v2): v1 05_cluster_db.sh plus
#  5a'. eukaryote resolution (v2): salterns proteins that carry a counted target label L for which eukaryote KEGG
#       references exist (class L_euk) and that are NOT on a kept / co-assembly MAG are compared (DIAMOND blastp
#       --sensitive) with the external references of {L, L_euk}; if the best hit is an L_euk reference the protein is
#       relabelled L_euk (decoy_reported). This mirrors v1's K10944 resolution and keeps algal / fungal genes of the
#       Dunaliella-rich waters out of the prokaryotic counts.
#  5b. per-label MMseqs2 clustering as v1 (99 % targets / markers / reported classes, 95 % decoys; c 0.8; 90 % for the
#      hold-out units); files are named by a filesystem-safe form of the label (refs/clu/labels.tsv: label, role, n, safe).
set -euo pipefail
W=/scratch1/afreiburger/salterns_2026/read_markers_v2
E=/scratch1/afreiburger/salterns_2026/envs
A=/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin
MM=$E/dram15/bin/mmseqs
DM=$E/bt_dastool/bin/diamond
PY=$A/python3
export TMPDIR=$W/tmp
T=${T:-32}
cd $W
mkdir -p refs/clu db
# 5a. K10944 resolution (v1)
awk -F'\t' 'NR>1 && $2=="K10944_unresolved"{print $1}' refs/labeled_meta.tsv > refs/clu/k10944_q.ids
awk -F'\t' 'NR>1 && $4!="salterns" && ($2=="pmoA"||$2=="amoA"||$2=="cuMMO_other"){print $1}' refs/labeled_meta.tsv > refs/clu/k10944_r.ids
$A/seqkit grep -f refs/clu/k10944_q.ids refs/labeled.faa > refs/clu/k10944_q.faa
$A/seqkit grep -f refs/clu/k10944_r.ids refs/labeled.faa > refs/clu/k10944_r.faa
nice -n 10 $DM makedb --in refs/clu/k10944_r.faa -d refs/clu/k10944_r --quiet
nice -n 10 $DM blastp -q refs/clu/k10944_q.faa -d refs/clu/k10944_r -o refs/clu/k10944_hits.tsv -k 1 --sensitive \
  --threads $T --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore || true
# 5a'. eukaryote resolution: queries / references
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
euk = {l[:-4] for l in m.label.unique() if l.endswith("_euk")}
q = m[(m.source == "salterns") & m.label.isin(euk) & (m.mag == "")]
r = m[(m.source != "salterns") & (m.label.isin(euk) | m.label.isin({e + "_euk" for e in euk}))]
open("refs/clu/euk_q.ids", "w").write("\n".join(q.pid) + "\n")
open("refs/clu/euk_r.ids", "w").write("\n".join(r.pid) + "\n")
print("eukaryote resolution: queries", len(q), "references", len(r), "labels", len(euk))
EOF
$A/seqkit grep -f refs/clu/euk_q.ids refs/labeled.faa > refs/clu/euk_q.faa
$A/seqkit grep -f refs/clu/euk_r.ids refs/labeled.faa > refs/clu/euk_r.faa
nice -n 10 $DM makedb --in refs/clu/euk_r.faa -d refs/clu/euk_r --quiet
nice -n 10 $DM blastp -q refs/clu/euk_q.faa -d refs/clu/euk_r -o refs/clu/euk_hits.tsv -k 25 -e 1e-5 --sensitive \
  --threads $T --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore
$PY -I - <<'EOF'
import pandas as pd
m = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
lab = dict(zip(m.pid, m.label))
h = pd.read_csv("refs/clu/euk_hits.tsv", sep="\t", header=None, names="q s pid len e bits".split())
h["ql"], h["sl"] = h.q.map(lab), h.s.map(lab)
h = h[(h.sl == h.ql) | (h.sl == h.ql + "_euk")].sort_values(["q", "bits"], ascending=[True, False]).drop_duplicates("q")
to_euk = set(h.q[h.sl == h.ql + "_euk"])
sel = m.pid.isin(to_euk)
m.loc[sel, "how"] = m.loc[sel, "how"] + "+eukaryote_by_best_external_hit"
m.loc[sel, "label"] = m.loc[sel, "label"] + "_euk"
m.to_csv("refs/labeled_meta.tsv", sep="\t", index=False)
h.to_csv("refs/clu/euk_resolution.tsv", sep="\t", index=False)
print("salterns proteins relabelled eukaryote:", sel.sum(), m.loc[sel, "label"].value_counts().head(30).to_dict())
EOF
# 5b. roles, per-label clustering
$PY -I - <<'EOF'
import re
import pandas as pd
m = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
p = pd.read_csv("meta/ko_panel_v2.tsv", sep="\t", dtype=str).fillna("")
role = dict(zip(p.gene, p.role))
role.update({"amoA": "decoy_reported", "cuMMO_other": "decoy_reported", "mttB_family_other": "target", "euk_marker": "decoy",
             "SDIMO_other": "decoy", "GRE_other": "decoy", "Rieske_other": "decoy"})
def r_of(l):
    if l in role:
        return role[l]
    if l.endswith("_like") or l.endswith("_other"):
        return "decoy"
    if l.endswith("_euk"):
        return "decoy_reported"
    return "?"
m["role"] = [r_of(l) for l in m.label]
assert (m.role != "?").all(), m.label[m.role == "?"].unique()
m.to_csv("refs/labeled_meta.tsv", sep="\t", index=False)
with open("refs/clu/labels.tsv", "w") as fh:
    for lab, g in m.groupby("label"):
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", lab).strip("_")
        open(f"refs/clu/{safe}.ids", "w").write("\n".join(g.pid) + "\n")
        fh.write(f"{lab}\t{g.role.iloc[0]}\t{len(g)}\t{safe}\n")
EOF
# per-label clustering, 8 labels at a time x 4 threads (v1 ran them one after another); results concatenated in
# labels.tsv order, so the output does not depend on the parallelism
export A MM TMPDIR
clu1() {
  lab="$1"; role="$2"; safe="$4"
  id=0.99; [ "$role" = decoy ] && id=0.95
  [ -f refs/clu/${safe}.done ] && return 0
  $A/seqkit grep -f refs/clu/$safe.ids refs/labeled.faa > refs/clu/$safe.faa 2> /dev/null
  for pass in main c90; do
    pid=$id; [ $pass = c90 ] && pid=0.90
    rm -rf $TMPDIR/mm_$safe
    nice -n 10 $MM easy-cluster refs/clu/$safe.faa refs/clu/${safe}_$pass $TMPDIR/mm_$safe --min-seq-id $pid -c 0.8 --cov-mode 0 \
      --threads 4 -v 1 > /dev/null
    rm -rf $TMPDIR/mm_$safe
  done
  rm -f refs/clu/${safe}_main_all_seqs.fasta refs/clu/${safe}_c90_all_seqs.fasta refs/clu/${safe}_c90_rep_seq.fasta refs/clu/$safe.faa
  touch refs/clu/${safe}.done
}
export -f clu1
tr '\t' '\037' < refs/clu/labels.tsv | xargs -P 8 -I{} bash -c 'IFS=$'"'"'\037'"'"' read -r a b c d <<< "{}"; clu1 "$a" "$b" "$c" "$d"'
: > db/survey.faa; : > refs/clu/clusters99.tsv; : > refs/clu/clusters90.tsv
while IFS=$'\t' read lab role n safe; do
  [ -f refs/clu/${safe}.done ] || { echo "clustering of $lab missing" >&2; exit 1; }
  awk -v l="$lab" 'BEGIN{OFS="\t"} {print l, $0}' refs/clu/${safe}_main_cluster.tsv >> refs/clu/clusters99.tsv
  awk -v l="$lab" 'BEGIN{OFS="\t"} {print l, $0}' refs/clu/${safe}_c90_cluster.tsv >> refs/clu/clusters90.tsv
  cat refs/clu/${safe}_main_rep_seq.fasta >> db/survey.faa
  echo "$lab $role $n -> $(grep -c '^>' refs/clu/${safe}_main_rep_seq.fasta) reps"
done < refs/clu/labels.tsv
nice -n 10 $DM makedb --in db/survey.faa -d db/survey --threads $T --quiet
md5sum db/survey.faa db/survey.dmnd > db/survey.md5
echo "$(date) survey DB: $(grep -c '^>' db/survey.faa) representatives"
