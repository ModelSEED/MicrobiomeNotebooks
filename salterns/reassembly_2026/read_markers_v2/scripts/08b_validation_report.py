#!/usr/bin/env python
"""Step 8b: per-label validation metrics at the chosen cut-offs (cutoffs.tsv), from the simulated hold-out reads.
  recall           = held-out reads of label L (>= 30 aa inside the CDS) assigned to L / all such reads
  recall by divergence of the held-out gene: its best blastp identity to the validation DB (val/holdout_vs_valdb.tsv)
                     nearest validation-DB reference OF THE SAME LABEL, binned >= 90, 80-90, 70-80, 60-70, < 60 % (0 = no hit)
  precision        = TP / (TP + reads of other held-out labels or of background genes assigned to L)
  fp_per_M_background = background reads assigned to L per million background reads
  uh_assigned      = reads simulated from KOfam-unassigned target homologs that were assigned to L (truth unknown)
Writes val/validation_per_label.tsv and val/validation_recall_by_divergence.tsv.
Usage: 08b_validation_report.py <cutoffs.tsv> <val dmnd tsv(.gz)> <tag>
"""
import sys, os
from collections import Counter, defaultdict
import pandas as pd

sys.argv, args = sys.argv[:1], sys.argv[1:]
import importlib.util
spec = importlib.util.spec_from_file_location("cl", os.path.join(os.path.dirname(os.path.abspath(__file__)), "08_classify.py"))
cl = importlib.util.module_from_spec(spec); spec.loader.exec_module(cl)
cut_p, dm_p, tag = args
cut = pd.read_csv(cut_p, sep="\t", dtype=str)
g = dict(zip(cut.param, cut.value))
emax, lmin, margin, pdef = float(g["EMAX"]), int(g["LMIN"]), float(g["MARGIN"]), float(g["PMIN_DEFAULT"])
pmin = {k[5:]: float(v) for k, v in g.items() if k.startswith("PMIN:")}
truth = pd.read_csv("val/sim_truth.tsv", sep="\t", dtype=str)
hg = pd.read_csv("val/holdout_genes.tsv", sep="\t", dtype=str)
hv = pd.read_csv("val/holdout_vs_valdb.tsv", sep="\t", header=None, names="q s pident len e bits".split(), dtype={"q": str})
hlab = dict(zip(hg.gene, hg.label))
hv = hv[[cl.LAB.get(x) == hlab.get(q) for q, x in zip(hv.q, hv.s)]]      # nearest reference of the same label
best_id = hv.groupby("q").pident.max()
hits = dict(cl.groups([dm_p]))
pred = {}
for r in truth.read:
    L, b = cl.decide(hits.get(r, []), emax, lmin, margin, pmin, pdef)
    pred[r] = L if L not in (None, "ambiguous", "low_identity") else ("unassigned" if L is None else L)
truth["pred"] = truth.read.map(pred)
truth["div"] = truth.gene.map(best_id).fillna(0)


def dbin(x):
    return ">=90" if x >= 90 else ("80-90" if x >= 80 else ("70-80" if x >= 70 else ("60-70" if x >= 60 else "<60")))


truth["div_bin"] = truth["div"].map(dbin)
nbg = (truth.source == "background").sum()
rows, rows2 = [], []
for lab in sorted(set(truth.truth) | set(truth.pred)):
    if cl.ROLE.get(lab) not in ("target", "marker", "decoy_reported"):
        continue
    t = truth[truth.truth == lab]
    tp = (t.pred == lab).sum()
    fph = ((truth.pred == lab) & (truth.truth != lab) & (truth.source == "holdout")).sum()
    fpb = ((truth.pred == lab) & (truth.source == "background")).sum()
    uh = ((truth.pred == lab) & (truth.source == "unassigned_homolog")).sum()
    n_genes = t.gene.nunique()
    rows.append(dict(label=lab, role=cl.ROLE.get(lab), n_heldout_genes=n_genes, n_reads=len(t), tp=tp,
                     recall=tp / len(t) if len(t) else float("nan"),
                     ambiguous=(t.pred == "ambiguous").sum() / len(t) if len(t) else float("nan"),
                     low_identity=(t.pred == "low_identity").sum() / len(t) if len(t) else float("nan"),
                     to_other_label=((t.pred != lab) & ~t.pred.isin(["unassigned", "ambiguous", "low_identity"])).sum(),
                     fp_holdout=fph, fp_background=fpb, precision=tp / (tp + fph + fpb) if tp + fph + fpb else float("nan"),
                     fp_per_M_background=1e6 * fpb / nbg, uh_assigned=uh,
                     median_divergence_identity=t.drop_duplicates("gene")["div"].median() if len(t) else float("nan")))
    for b, tt in t.groupby("div_bin"):
        rows2.append(dict(label=lab, div_bin=b, n_genes=tt.gene.nunique(), n_reads=len(tt), recall=(tt.pred == lab).mean()))
out = pd.DataFrame(rows)
out.to_csv(f"val/validation_per_label_{tag}.tsv", sep="\t", index=False, float_format="%.4g")
pd.DataFrame(rows2).to_csv(f"val/validation_recall_by_divergence_{tag}.tsv", sep="\t", index=False, float_format="%.4g")
# confusion of held-out reads
conf = truth[truth.source == "holdout"].groupby(["truth", "pred"]).size().reset_index(name="reads")
conf.to_csv(f"val/validation_confusion_{tag}.tsv", sep="\t", index=False)
uhc = truth[truth.source == "unassigned_homolog"].groupby(["gene_label", "pred"]).size().reset_index(name="reads")
uhc.to_csv(f"val/validation_unassigned_homologs_{tag}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(out.to_string())
