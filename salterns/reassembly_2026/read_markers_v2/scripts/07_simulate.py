#!/usr/bin/env python
"""Step 7 (v2; v1 script with v2 paths): hold-out validation set for the read classifier.

* Hold-out unit = a 90 %-identity cluster of one label (refs/clu/clusters90.tsv). Per label, 20 % of its 90 % clusters
  (at least 1 when the label has >= 3 clusters, at most 60), random.Random(1).
* Validation DB = the production representatives (db/survey.faa) minus every representative whose 90 % cluster is
  held out -> val/val_survey.faa (+ .dmnd). So a held-out gene has no >= 90 %-identical relative in the validation DB.
* Simulated reads (148 bp, uniform start, random strand, 0.3 % substitutions) from one member per held-out cluster:
  - salterns member with a complete gene (partial=00): window = CDS +/- 150 bp of its contig (so reads straddle the
    gene ends as in real data);
  - else a KEGG member: CDS from KEGG REST ntseq (no flanks);
  - UniProt-only clusters are not simulated (no nucleotide sequence).
  Reads per gene = window length / 30 (about 5x), max 100. Truth = the gene's label if the read overlaps the CDS by
  >= 90 nt (30 aa), else 'none'.
* Background: 4,000 random salterns genes without any panel hit (E <= 1e-5), +/- 150 bp, truth 'none'.
* KOfam-unassigned homologs: up to 30 complete salterns genes per target KO that hit the target profile but pass no
  panel KO and are not in the DB (0.5-1 x threshold; refs/unassigned_homologs.tsv minus the '_like' decoys), truth
  'none' but evaluated separately (their true function is unknown).
Writes val/sim_reads.fa, val/sim_truth.tsv, val/holdout_genes.tsv (+ .faa for the identity-to-DB stratification).
"""
import os, sys, random, re, gzip, time, glob
from collections import defaultdict
import pandas as pd
import urllib.request

W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
R = "/scratch1/afreiburger/salterns_2026/rollout"
os.chdir(W)
os.makedirs("val", exist_ok=True)
rng = random.Random(1)
RL, FLANK, ERR = 148, 150, 0.003
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(p):
    op = gzip.open if p.endswith(".gz") else open
    name, seq = None, []
    with op(p, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = l[1:].rstrip("\n"), []
            else:
                seq.append(l.strip())
    if name:
        yield name, "".join(seq)


meta = pd.read_csv("refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
lab_of = dict(zip(meta.pid, meta.label))
src_of = dict(zip(meta.pid, meta.source))
role_of_label = dict(zip(meta.label, meta.role))
c90 = pd.read_csv("refs/clu/clusters90.tsv", sep="\t", header=None, names=["label", "rep", "member"], dtype=str)
c99 = pd.read_csv("refs/clu/clusters99.tsv", sep="\t", header=None, names=["label", "rep", "member"], dtype=str)
clu90_of = dict(zip(c90.member, c90.rep))

# ---- choose held-out 90 % clusters
held = set()
for lab, g in c90.groupby("label"):
    reps = sorted(g.rep.unique())
    if len(reps) < 3:
        continue
    k = min(60, max(1, round(0.2 * len(reps))))
    held.update(rng.sample(reps, k))
# ---- validation DB
keep = [r for r in c99.rep.unique() if clu90_of[r] not in held]
keep = set(keep)
n = 0
with open("val/val_survey.faa", "w") as fo:
    for name, s in read_fasta("db/survey.faa"):
        if name.split()[0] in keep:
            fo.write(f">{name}\n{s}\n"); n += 1
print("validation DB reps", n, "of", c99.rep.nunique(), "held-out 90% clusters", len(held), flush=True)

# ---- pick one simulated member per held-out cluster
gene_info = {}   # salterns pid -> (sample, contig, start, end, strand, partial)
def parse_prodigal_header(h):
    f = h.split(" # ")
    pid = f[0]
    st, en, strand = int(f[1]), int(f[2]), int(f[3])
    partial = re.search(r"partial=(\d\d)", h).group(1)
    return pid, st, en, strand, partial

members = defaultdict(list)
for rep, m in zip(c90.rep, c90.member):
    if rep in held:
        members[rep].append(m)
need_salterns = defaultdict(set)   # base -> pids we need coordinates for
for rep, ms in members.items():
    for m in ms:
        if src_of[m] == "salterns":
            base = "coasm520" if m.startswith("coasm|") else m.split("|")[0]
            need_salterns[base].add(m)
for base, pids in need_salterns.items():
    for h, _ in read_fasta(f"{V1}/genes/{base}/{base}.faa"):
        pid = h.split(" # ")[0]
        if pid in pids:
            p, st, en, strand, partial = parse_prodigal_header(h)
            gene_info[p] = (base, st, en, strand, partial)
choice = {}
kegg_needed = {}
for rep in sorted(members):
    ms = sorted(members[rep])
    sal = [m for m in ms if m in gene_info and gene_info[m][4] == "00"]
    kg = [m for m in ms if src_of[m] == "kegg"]
    if sal:
        choice[rep] = rng.choice(sal)
    elif kg:
        choice[rep] = rng.choice(kg); kegg_needed[choice[rep]] = choice[rep].split("|", 2)[2]
print("held-out clusters simulated", len(choice), "(salterns", sum(1 for v in choice.values() if v in gene_info), ")", flush=True)

# ---- KEGG ntseq
kegg_nt = {}
ids = sorted(kegg_needed.items())
for i in range(0, len(ids), 10):
    batch = ids[i:i + 10]
    url = "https://rest.kegg.jp/get/" + "+".join(g for _, g in batch) + "/ntseq"
    for t in range(5):
        try:
            txt = urllib.request.urlopen(url, timeout=120).read().decode()
            break
        except Exception:
            time.sleep(5 * (t + 1)); txt = ""
    cur, seq = None, {}
    for l in txt.splitlines():
        if l.startswith(">"):
            cur = l[1:].split()[0]; seq[cur] = []
        elif cur:
            seq[cur].append(l.strip())
    for pid, g in batch:
        if g in seq:
            kegg_nt[pid] = "".join(seq[g]).upper()
    time.sleep(0.35)
print("KEGG ntseq fetched", len(kegg_nt), "of", len(kegg_needed), flush=True)

# ---- background genes (no panel hit at all)
panel_hit = set()
for f in glob.glob(f"{V1}/hmm/*.tbl") + glob.glob("hmm_new/*.tbl"):   # v2: no hit to any of the 286 v2 profiles
    for l in open(f):
        if not l.startswith("#"):
            panel_hit.add(l.split()[0])
bg = {}
samples = sorted(os.listdir(f"{R}/reads"))
for s in rng.sample(samples, 10):
    hs = [h for h, _ in read_fasta(f"{V1}/genes/{s}/{s}.faa")]
    hs = [h for h in hs if h.split(" # ")[0] not in panel_hit and "partial=00" in h]
    for h in rng.sample(hs, 400):
        p, st, en, strand, partial = parse_prodigal_header(h)
        bg[p] = (s, st, en, strand, partial)

# ---- KOfam-unassigned homologs (hit a target profile, pass no panel KO): up to 30 complete genes per top target KO
uh_tab = pd.read_csv("refs/unassigned_homologs.tsv", sep="\t", dtype=str)
uh_tab = uh_tab[~uh_tab.set.str.startswith(("pfam_", "name_", "ext", "rev_")) & ~uh_tab.pid.isin(lab_of)]   # not the '_like' decoys
uh_pick = defaultdict(set)
for ko, g in uh_tab.groupby("top_target_ko"):
    for p in rng.sample(sorted(g.pid), min(60, len(g))):
        base = "coasm520" if p.startswith("coasm|") else p.split("|")[0]
        uh_pick[base].add(p)
uh = {}
for base, pids in uh_pick.items():
    for h, _ in read_fasta(f"{V1}/genes/{base}/{base}.faa"):
        pid = h.split(" # ")[0]
        if pid in pids and "partial=00" in h:
            q, st, en, strand, partial = parse_prodigal_header(h)
            uh[q] = (base, st, en, strand, partial)
uh_ko = dict(zip(uh_tab.pid, uh_tab.top_target_ko))
kept_uh = {}
for ko in sorted(set(uh_ko[p] for p in uh)):
    ps = sorted(p for p in uh if uh_ko[p] == ko)
    for p in ps[:30]:
        kept_uh[p] = uh[p]
uh = kept_uh
bg.update(uh)
print("unassigned homologs simulated", len(uh), flush=True)

# ---- contig windows for salterns genes
want = defaultdict(dict)    # fasta path -> contig -> list of pids
for p, (base, st, en, strand, partial) in list(gene_info.items()) + list(bg.items()):
    if p not in bg and p not in choice.values():
        continue
    if base == "coasm520":
        mag = p.split("|")[1]; contig = p.split("|")[2].rsplit("_", 1)[0]
        fp = f"/scratch1/afreiburger/salterns_2026/coasm_qc/mags_kept/{mag}.fa"
    else:
        contig = p.split("|", 1)[1].rsplit("_", 1)[0]
        fp = f"{R}/{base}/asm/contigs.fa"
    want[fp].setdefault(contig, []).append(p)
windows = {}
for fp, cs in want.items():
    for h, s in read_fasta(fp):
        c = h.split()[0]
        if c in cs:
            for p in cs[c]:
                base, st, en, strand, partial = gene_info.get(p) or bg.get(p)
                a, b = max(0, st - 1 - FLANK), min(len(s), en + FLANK)
                windows[p] = (s[a:b].upper(), st - 1 - a, en - a)   # window, cds start (0-based), cds end (excl)
print("windows", len(windows), flush=True)


def mutate(r):
    r = list(r)
    for i in range(len(r)):
        if rng.random() < ERR:
            r[i] = rng.choice([x for x in "ACGT" if x != r[i]])
    return "".join(r)


nread = 0
with open("val/sim_reads.fa", "w") as fr, open("val/sim_truth.tsv", "w") as ft, open("val/holdout_genes.tsv", "w") as fh:
    ft.write("read\tgene\tgene_label\ttruth\toverlap_nt\tsource\n")
    fh.write("gene\tlabel\trole\tsource\tcluster90\tlen_nt\n")
    todo = [(p, lab_of[p], "holdout") for p in choice.values()] + \
           [(p, ("homolog_of_" + uh_ko[p]) if p in uh else "none", "unassigned_homolog" if p in uh else "background") for p in bg]
    for p, lab, kind in todo:
        if p in windows:
            win, cs, ce = windows[p]
        elif p in kegg_nt:
            win, cs, ce = kegg_nt[p], 0, len(kegg_nt[p])
        else:
            continue
        if len(win) < RL:
            continue
        if kind == "holdout":
            fh.write(f"{p}\t{lab}\t{role_of_label.get(lab, '')}\t{src_of[p]}\t{clu90_of[p]}\t{ce - cs}\n")
        nr = min(100, max(3, len(win) // 30))
        for i in range(nr):
            a = rng.randint(0, len(win) - RL)
            r = win[a:a + RL]
            ov = max(0, min(a + RL, ce) - max(a, cs))
            if rng.random() < 0.5:
                r = rc(r)
            nread += 1
            rid = f"sim{nread}"
            fr.write(f">{rid}\n{mutate(r)}\n")
            truth = lab if (kind == "holdout" and ov >= 90) else "none"
            ft.write(f"{rid}\t{p}\t{lab}\t{truth}\t{ov}\t{kind}\n")
print("simulated reads", nread)
