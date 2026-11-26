#!/usr/bin/env python
"""Step 8d (v2): for chosen labels, recall / precision at the v1 cut-offs with an extra per-label identity cut
PMIN in {0, 50, 60, 70, 75, 80} (as v1 chose its per-label cuts), from the simulated hold-out reads (k5 run).
Writes val/identity_cuts.tsv. Usage: 08d_identity_cuts.py <label> [<label> ...]"""
import sys, os
import pandas as pd
labs, sys.argv = sys.argv[1:], sys.argv[:1]
import importlib.util
spec = importlib.util.spec_from_file_location("cl", "scripts/08_classify.py"); cl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cl)
cut = pd.read_csv("cutoffs.tsv", sep="\t", dtype=str); g = dict(zip(cut.param, cut.value))
emax, lmin, margin, pdef = float(g["EMAX"]), int(g["LMIN"]), float(g["MARGIN"]), float(g["PMIN_DEFAULT"])
pmin0 = {k[5:]: float(v) for k, v in g.items() if k.startswith("PMIN:")}
truth = pd.read_csv("val/sim_truth.tsv", sep="\t", dtype=str)
T = dict(zip(truth.read, truth.truth)); SRC = dict(zip(truth.read, truth.source))
hits = dict(cl.groups(["val/sim_mid-sensitive_k5.tsv.gz"]))
rows = []
for lab in labs:
    for p in [0, 50, 60, 70, 75, 80]:
        pm = dict(pmin0); pm[lab] = p
        tp = fp = fn = uh = 0; fps = {}
        for r, t in T.items():
            L, b = cl.decide(hits.get(r, []), emax, lmin, margin, pm, pdef)
            if t == lab:
                tp += L == lab; fn += L != lab
            elif L == lab:
                if SRC[r] == "unassigned_homolog":
                    uh += 1
                else:
                    fp += 1; fps[t] = fps.get(t, 0) + 1
        rows.append(dict(label=lab, PMIN=p, n_true=tp + fn, recall=tp / (tp + fn) if tp + fn else float("nan"),
                         precision=tp / (tp + fp) if tp + fp else float("nan"), fp=fp, uh_assigned=uh,
                         fp_from=";".join(f"{k}:{v}" for k, v in sorted(fps.items(), key=lambda x: -x[1])[:3])))
    print(pd.DataFrame([r for r in rows if r["label"] == lab]).round(3).to_string(), flush=True)
pd.DataFrame(rows).to_csv("val/identity_cuts.tsv", sep="\t", index=False, float_format="%.4f")
