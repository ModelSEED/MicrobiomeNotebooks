#!/usr/bin/env python
"""Step 8 (v2; v1 script, v2 paths, grid reduced, see README): best-hit read classification with decoys.

Per read (each mate separately), from the DIAMOND blastx hits (top 5 references, E <= 1e-3):
  1. keep hits with E <= EMAX and aligned length >= LMIN aa, ignoring the over-long marker references listed in
     refs/excluded_reps.tsv (08a_exclude_reps.py);
  2. best = the highest-bitscore hit; its label L = the label of that reference's cluster;
  3. if any kept hit with a different label scores >= best - MARGIN bits and either label is a target / marker /
     reported class, the read is 'ambiguous' (not counted);
  4. else the read is assigned to L if pident(best) >= PMIN[L] (PMIN default PMIN_DEFAULT), otherwise 'low_identity'.
Cut-offs come from cutoffs.tsv (written after the validation, see README); the defaults below are only for the grid.

  sample   <cutoffs.tsv> <out_prefix> <dmnd.tsv.gz> [...]  -> <out_prefix>.counts.tsv, <out_prefix>.reads.tsv.gz
  validate <val/dmnd.tsv.gz> <val/sim_truth.tsv>           -> val/grid.tsv, val/per_label_*.tsv
"""
import sys, gzip, os, itertools
from collections import defaultdict, Counter
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
os.chdir(W)
c99 = pd.read_csv("refs/clu/clusters99.tsv", sep="\t", header=None, names=["label", "rep", "member"], dtype=str)
LAB = dict(zip(c99.rep, c99.label))
meta = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
ROLE = dict(zip(meta.label, meta.role))
COUNTED = {l for l, r in ROLE.items() if r in ("target", "marker", "decoy_reported")}
EXCL = set(pd.read_csv("refs/excluded_reps.tsv", sep="\t").rep) if os.path.exists("refs/excluded_reps.tsv") else set()


def groups(paths):
    for p in paths:
        op = gzip.open if p.endswith(".gz") else open
        cur, hs = None, []
        with op(p, "rt") as fh:
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if f[0] != cur:
                    if cur is not None:
                        yield cur, hs
                    cur, hs = f[0], []
                # bits, pident, alen, evalue, rep
                hs.append((float(f[11]), float(f[2]), int(f[3]), float(f[10]), f[1]))
            if cur is not None:
                yield cur, hs


def decide(hs, emax, lmin, margin, pmin, pmin_default):
    k = [h for h in hs if h[3] <= emax and h[2] >= lmin and h[4] not in EXCL]
    if not k:
        return None, None
    b = max(k)
    L = LAB[b[4]]
    comp = max((h[0] for h in k if LAB[h[4]] != L and (L in COUNTED or LAB[h[4]] in COUNTED)), default=None)
    if comp is not None and comp >= b[0] - margin:
        return "ambiguous", b
    if b[1] < pmin.get(L, pmin_default):
        return "low_identity", b
    return L, b


def cmd_sample(cut_p, prefix, paths):
    cut = pd.read_csv(cut_p, sep="\t", dtype=str)
    g = dict(zip(cut.param, cut.value))
    emax, lmin, margin, pdef = float(g["EMAX"]), int(g["LMIN"]), float(g["MARGIN"]), float(g["PMIN_DEFAULT"])
    pmin = {k[5:]: float(v) for k, v in g.items() if k.startswith("PMIN:")}
    cnt = Counter(); amb = Counter(); low = Counter()
    with gzip.open(prefix + ".reads.tsv.gz", "wt") as fo:
        fo.write("read\tmate\tlabel\trep\tpident\talen\tbits\n")
        for p in paths:
            mate = "1" if os.path.basename(p).replace(".tsv.gz", "").endswith("_1") else "2"
            for q, hs in groups([p]):
                L, b = decide(hs, emax, lmin, margin, pmin, pdef)
                if L is None:
                    continue
                if L == "ambiguous":
                    amb[LAB[b[4]]] += 1
                elif L == "low_identity":
                    low[LAB[b[4]]] += 1
                else:
                    cnt[L] += 1
                    if L in COUNTED:
                        fo.write(f"{q}\t{mate}\t{L}\t{b[4]}\t{b[1]}\t{b[2]}\t{b[0]}\n")
    labs = sorted(set(cnt) | set(amb) | set(low))
    pd.DataFrame({"label": labs, "role": [ROLE.get(l, "") for l in labs], "reads": [cnt[l] for l in labs],
                  "ambiguous_best": [amb[l] for l in labs], "low_identity_best": [low[l] for l in labs]}) \
        .to_csv(prefix + ".counts.tsv", sep="\t", index=False)


def cmd_validate(dm_p, truth_p, out_p="val/grid.tsv"):
    """Grid over EMAX x LMIN x MARGIN x PMIN (one global identity cut). Per counted label: n_true, tp, recall;
    false positives split by origin: fp_holdout (reads of another held-out label), fp_background (genes without any
    panel hit), and separately uh_assigned (reads of KOfam-unassigned target homologs, truth unknown);
    precision = tp / (tp + fp_holdout + fp_background)."""
    truth = pd.read_csv(truth_p, sep="\t", dtype=str)
    T = dict(zip(truth.read, truth.truth))
    SRC = dict(zip(truth.read, truth.source))
    hits = dict(groups([dm_p]))
    grid = []
    for emax, lmin, margin, pdef in itertools.product([1e-5], [30], [0, 2], [0, 50, 60, 70]):
        tp = Counter(); nt = Counter(); fph = Counter(); fpb = Counter(); uh = Counter()
        for r, t in T.items():
            L, b = decide(hits.get(r, []), emax, lmin, margin, {}, pdef)
            p = L if L not in (None, "ambiguous", "low_identity") else "unassigned"
            if t != "none":
                nt[t] += 1
                if p == t:
                    tp[t] += 1
                    continue
            if p == "unassigned":
                continue
            src = SRC[r]
            if src == "unassigned_homolog":
                uh[p] += 1
            elif src == "background":
                fpb[p] += 1
            else:
                fph[p] += 1
        for lab in sorted(set(nt) | set(fph) | set(fpb) | set(uh)):
            if ROLE.get(lab) not in ("target", "marker", "decoy_reported"):
                continue
            fp = fph[lab] + fpb[lab]
            grid.append(dict(EMAX=emax, LMIN=lmin, MARGIN=margin, PMIN=pdef, label=lab, role=ROLE.get(lab),
                             n_true=nt[lab], tp=tp[lab], fp_holdout=fph[lab], fp_background=fpb[lab], uh_assigned=uh[lab],
                             recall=tp[lab] / nt[lab] if nt[lab] else float("nan"),
                             precision=tp[lab] / (tp[lab] + fp) if tp[lab] + fp else float("nan")))
    pd.DataFrame(grid).to_csv(out_p, sep="\t", index=False)


if __name__ == "__main__":
    if sys.argv[1] == "sample":
        cmd_sample(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif sys.argv[1] == "validate":
        cmd_validate(*sys.argv[2:])
