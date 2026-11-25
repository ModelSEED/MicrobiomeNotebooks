#!/usr/bin/env python3
"""Threshold sensitivity, part 2a: can contigs shorter than the 2.5 kb binning floor be assigned to our MAG species as
well as >= 2.5 kb ones? (per Andrew 2026-11-25: "is 2.5kb a reasonable threshold, or would lowering it to 1kb better
capture the system?")
The KBase assemblies (narrative 261700, MEGAHIT app, min_contig_length = 2000 in every provenance checked) contain no
contig under 2 kb, so their real unbinned contigs fall in two bands, 2.0-2.5 kb and >= 2.5 kb. The 1.0-2.0 kb band is
tested here by ground-truth fragments only; real 1-2 kb contigs come from re-assemblies with a 1 kb floor
(reassembly/, a second pass of this script with --tag reasm).

One method for both bands, so the method cannot confound the threshold: minimap2 2.28, preset asm20 with base-level
alignment (-c), against all 13,164 MAGs. To keep the index small, it is inverted: the test contigs are the index
(a few Gbp) and the MAGs (~38 Gbp) are streamed as queries. Each alignment scores identity (PAF matches / block length,
gaps count as mismatches) and the aligned fraction of the TEST contig.

Assignment of a test contig: species cluster (census run 4) of the best hit with identity >= 0.95 and aligned
fraction >= 0.5 (also reported at >= 0.8). If a hit from a DIFFERENT cluster scores within 0.01 (identity x
fraction) of the best, the contig is "ambiguous" (conserved or mobile sequence). Otherwise it is "assigned", or
"none" when nothing passes.

A. Ground truth (upper bound, since real unbinned contigs are harder than binned sequence): fragments cut from MAGs in
   species clusters seen in >= 2 samples. One random MAG per cluster for GT_MAGS clusters; GT_PER_BAND fragments per
   band per MAG. Lengths: 1.0-2.0 kb uniform; 2.0-2.5 kb and >= 2.5 kb drawn from the real per-band length
   distribution of the selected assemblies (capped at 50 kb). Hits to the source MAG are dropped. Reported per band
   and per length bin: correct, wrong, ambiguous, none; precision = correct / (correct + wrong).
B. Real unbinned contigs (in no MAG of their own assembly) of N_ASM assemblies that have admitted MAGs.
   Selection: the two Bakta timing-test assemblies plus a seeded, habitat-stratified draw of downloaded assemblies of
   20-400 Mbp. Outcome per contig, under the same-sample abundance rule (unbinned/DESIGN.md):
     same_sample  = assigned to a cluster with a MAG from the SAME run (its abundance may merge there)
     cross_sample = assigned only to clusters without a MAG in this run (identifies; moves no abundance)
     ambiguous / none
   Reported per band by contigs, bp and abundance (contig length x MEGAHIT k-mer coverage).
Outputs: unbinned/threshold/assign_<tag>_{gt,real}_contigs.tsv.gz, _selection.tsv, _summary.json (tag "kbase" here).
Work files (large) in ~/gtdbtk_kbase/mm2/: all_mags.fna (prefixed `<mag>|<contig>`, reusable for the unbinned
pipeline), idx_<tag>.fna, hits_<tag>.paf.gz.
Usage: ~/Documents/py_venv/bin/python -I threshold_assign.py [--threads 16]   (resumable: finished steps are kept)
"""
import argparse, glob, gzip, json, os, random, re, subprocess, sys, time
from collections import defaultdict
import pandas as pd

MM2 = os.path.expanduser("~/Documents/EmilyKin/meta/mag_gene_ab/tools/minimap2-2.28_x64-linux/minimap2")
MAGS = os.path.expanduser("~/gtdbtk_kbase/genomes"); ASM = os.path.expanduser("~/gtdbtk_kbase/assemblies")
W = os.path.expanduser("~/gtdbtk_kbase/mm2"); os.makedirs(W, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "threshold"); os.makedirs(OUT, exist_ok=True)
CENSUS = os.path.join(HERE, "..", "census", "outputs")
TEST_ASM = ["SRR7963043", "SRR10127268"]
N_ASM, GT_MAGS, GT_PER_BAND, SEED = 20, 300, 10, 20261125
ID_MIN, AF_MIN, AF_STRICT, TIE = 0.95, 0.5, 0.8, 0.01
ap = argparse.ArgumentParser(); ap.add_argument("--threads", type=int, default=16); ap.add_argument("--tag", default="kbase")
a = ap.parse_args(); T = a.tag
log = lambda m: print(f"[{time.strftime('%F %T')}] {m}", flush=True)

cen = pd.read_csv(f"{CENSUS}/mag_census.tsv", sep="\t", usecols=["mag", "accession", "habitat", "species_cluster"])
clu = dict(zip(cen.mag, cen.species_cluster)); run_of = dict(zip(cen.mag, cen.accession))
runs_with_mag = defaultdict(set)                       # cluster -> runs that have a MAG of it
for m, c in clu.items(): runs_with_mag[c].add(run_of[m])
hab = cen.groupby("accession").habitat.first().to_dict()


BANDS = ("1.0-2.0", "2.0-2.5", "ge2.5")
band_of = lambda L: "ge2.5" if L >= 2500 else "2.0-2.5" if L >= 2000 else "1.0-2.0"
LBINS = [1000, 1500, 2000, 2500, 5000, 10000, 10**9]
LLAB = ["1.0-1.5", "1.5-2.0", "2.0-2.5", "2.5-5", "5-10", "ge10"]


def fasta(path):
    name = None; seq = []
    for l in open(path):
        if l[0] == ">":
            if name is not None: yield name, "".join(seq)
            name = l[1:].rstrip("\n"); seq = []
        else: seq.append(l.strip())
    if name is not None: yield name, "".join(seq)


# 1. all MAGs, prefixed (written once)
allm = f"{W}/all_mags.fna"
if not os.path.exists(allm):
    log("writing all_mags.fna")
    with open(allm + ".tmp", "w") as fo:
        for f in sorted(glob.glob(f"{MAGS}/*.fna")):
            m = os.path.basename(f)[:-4]
            for l in open(f):
                fo.write(f">{m}|{l[1:].split()[0]}\n" if l[0] == ">" else l)
    os.rename(allm + ".tmp", allm)

# 2. selection, real contigs and ground-truth fragments -> idx.fna
idx = f"{W}/idx_{T}.fna"; sel_tsv = f"{OUT}/assign_{T}_selection.tsv"
if not os.path.exists(idx):
    rnd = random.Random(SEED)
    binned = defaultdict(set)
    for f in glob.glob(f"{MAGS}/*.fna"):
        run = os.path.basename(f)[:-4].split("__")[-1]
        for l in open(f):
            if l[0] == ">": binned[run].add(l[1:].split()[0])
    pa = pd.read_csv(f"{OUT}/per_assembly.tsv", sep="\t")
    pool = pa[pa.has_mags & pa.bp.between(20e6, 400e6) & ~pa.run.isin(TEST_ASM)].copy()
    pool["habitat"] = pool.run.map(hab).fillna("unknown")
    by_h = {h: rnd.sample(list(g.run), len(g)) for h, g in pool.groupby("habitat")}
    pick = list(TEST_ASM)
    while len(pick) < N_ASM and any(by_h.values()):
        for h in sorted(by_h):
            if by_h[h] and len(pick) < N_ASM: pick.append(by_h[h].pop())
    lens = {}; meta = []
    with open(idx + ".tmp", "w") as fo:
        for run in pick:
            f = glob.glob(f"{ASM}/{run}_*.fna")[0]
            for name, s in fasta(f):
                cid = name.split()[0]; L = len(s)
                if cid in binned[run] or L < 1000: continue
                band = band_of(L); lens.setdefault(band, []).append(L)
                m = re.search(r"multi=([\d.]+)", name); cov = float(m.group(1)) if m else 0.0
                fo.write(f">R|{run}|{cid}\n{s}\n"); meta.append((f"R|{run}|{cid}", run, cid, band, L, cov))
        # ground truth: one MAG per cluster seen in >= 2 samples
        cl = pd.read_csv(f"{CENSUS}/clusters.tsv", sep="\t", usecols=["species_cluster", "n_samples"])
        multi = sorted(cl[cl.n_samples >= 2].species_cluster)
        members = defaultdict(list)
        for m, c in clu.items(): members[c].append(m)
        gt_clusters = rnd.sample(multi, min(GT_MAGS, len(multi)))
        for c in gt_clusters:
            m = rnd.choice(sorted(members[c]))
            ctg = [(n.split()[0], s) for n, s in fasta(f"{MAGS}/{m}.fna")]
            for band in BANDS:
                for i in range(GT_PER_BAND):
                    L = rnd.randrange(1000, 2000) if band == "1.0-2.0" else min(rnd.choice(lens[band]), 50000)
                    ok = [(n, s) for n, s in ctg if len(s) >= L]
                    if not ok: continue
                    n, s = rnd.choices(ok, weights=[len(s) for _, s in ok])[0]
                    st = rnd.randrange(0, len(s) - L + 1)
                    q = f"G|{m}|{band}|{i}"; fo.write(f">{q}\n{s[st:st + L]}\n")
                    meta.append((q, run_of[m], f"{n}:{st}", band, L, 0.0))
    pd.DataFrame(meta, columns=["id", "run", "contig", "band", "length", "multi"]).to_csv(f"{W}/idx_{T}_meta.tsv.gz", sep="\t", index=False)
    pd.DataFrame({"run": pick, "habitat": [hab.get(r, "unknown") for r in pick]}).to_csv(sel_tsv, sep="\t", index=False)
    os.rename(idx + ".tmp", idx)
    log(f"idx.fna: {len(pick)} assemblies, {sum(1 for x in meta if x[0][0] == 'R')} real contigs, "
        f"{sum(1 for x in meta if x[0][0] == 'G')} ground-truth fragments")

# 3. minimap2: test contigs indexed (one part), MAGs streamed
paf = f"{W}/hits_{T}.paf.gz"
if not os.path.exists(paf):
    log("minimap2 running")
    cmd = (f"set -o pipefail; {MM2} -x asm20 -c -N 200 -p 0.5 -I 64G -t {a.threads} {idx} {allm} 2> {W}/mm2_{T}.log "
           f"| gzip -1 > {paf}.tmp")
    rc = subprocess.run(["bash", "-c", cmd]).returncode
    if rc: sys.exit(f"minimap2 failed rc={rc}; see {W}/mm2_{T}.log")
    os.rename(paf + ".tmp", paf); log("minimap2 done")

# 4. assign
meta = pd.read_csv(f"{W}/idx_{T}_meta.tsv.gz", sep="\t").set_index("id")
best = defaultdict(dict)        # test id -> {cluster: (score, identity, af, mag)} for hits passing ID_MIN and AF_MIN
src = {i: i.split("|")[1] for i in meta.index if i[0] == "G"}
with gzip.open(paf, "rt") as fh:
    for l in fh:
        p = l.split("\t", 12)
        q, t = p[0], p[5]
        mag = q.split("|", 1)[0]
        if src.get(t) == mag: continue                 # ground truth: never its own MAG
        tl, ts, te, mt, bl = int(p[6]), int(p[7]), int(p[8]), int(p[9]), int(p[10])
        ident = mt / bl; af = (te - ts) / tl
        if ident < ID_MIN or af < AF_MIN: continue
        c = clu[mag]; sc = ident * af
        if sc > best[t].get(c, (0,))[0]: best[t][c] = (sc, ident, af, mag)


def call(t, af_min):
    h = {c: v for c, v in best.get(t, {}).items() if v[2] >= af_min}
    if not h: return "none", "", "", 0.0, 0.0
    r = sorted(h.items(), key=lambda kv: -kv[1][0]); c, (sc, ident, af, mag) = r[0]
    if len(r) > 1 and r[1][1][0] >= sc - TIE: return "ambiguous", c, mag, ident, af
    return "assigned", c, mag, ident, af


rows = []
for t, m in meta.iterrows():
    for lab, af_min in (("af50", AF_MIN), ("af80", AF_STRICT)):
        st, c, mag, ident, af = call(t, af_min)
        if t[0] == "G":
            truth = clu[t.split("|")[1]]
            st2 = "none" if st == "none" else "ambiguous" if st == "ambiguous" else ("correct" if c == truth else "wrong")
        else:
            st2 = st if st != "assigned" else ("same_sample" if m.run in runs_with_mag[c] else "cross_sample")
        rows.append((t, lab, m.run, m.band, m.length, m.multi, st2, c, mag, round(ident, 4), round(af, 3)))
d = pd.DataFrame(rows, columns=["id", "rule", "run", "band", "length", "multi", "outcome", "cluster", "best_mag", "identity", "aligned_frac"])
gt, real = d[d.id.str.startswith("G|")], d[d.id.str.startswith("R|")]
gt.to_csv(f"{OUT}/assign_{T}_gt_contigs.tsv.gz", sep="\t", index=False)
real.to_csv(f"{OUT}/assign_{T}_real_contigs.tsv.gz", sep="\t", index=False)
summ = {"params": {"minimap2": "2.28 -x asm20 -c -N 200 -p 0.5", "identity_min": ID_MIN, "af_min": [AF_MIN, AF_STRICT],
                   "tie": TIE, "seed": SEED, "assemblies": N_ASM, "gt_mags": GT_MAGS, "gt_per_band": GT_PER_BAND},
        "ground_truth": {}, "real": {}}
gt = gt.assign(lbin=pd.cut(gt.length, LBINS, right=False, labels=LLAB).astype(str))
for key, cols in (("ground_truth", ["rule", "band"]), ("ground_truth_by_length", ["rule", "lbin"])):
    summ.setdefault(key, {})
    for (rule, band), g in gt.groupby(cols):
        v = g.outcome.value_counts(); n = len(g)
        summ[key][f"{rule}|{band}"] = {"n": n, **{k: round(v.get(k, 0) / n * 100, 2) for k in ("correct", "wrong", "ambiguous", "none")},
                                       "precision": round(v.get("correct", 0) / max(1, v.get("correct", 0) + v.get("wrong", 0)) * 100, 2)}
real = real.assign(bp=real.length, abund=real.length * real.multi)
for (rule, band), g in real.groupby(["rule", "band"]):
    o = {"contigs": len(g)}
    for w in ("contigs", "bp", "abund"):
        tot = len(g) if w == "contigs" else g[w].sum()
        s = g.groupby("outcome").size() if w == "contigs" else g.groupby("outcome")[w].sum()
        o[w] = {k: round(float(s.get(k, 0)) / tot * 100, 2) for k in ("same_sample", "cross_sample", "ambiguous", "none")}
    summ["real"][f"{rule}|{band}"] = o
json.dump(summ, open(f"{OUT}/assign_{T}_summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
