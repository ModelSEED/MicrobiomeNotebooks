#!/usr/bin/env python
"""Step 11c (v2 deliverable 1): read-based genus profile from the 8 ribosomal-protein marker families, computed from the
v1 per-read tables (read_markers/results/per_sample/<s>.reads.tsv.gz; v1 files are only read).

Genus of a marker read = genus of its best-hit DB representative (v1 refs), resolved in tiers:
  1. salterns_MAG   : majority GTDB R232 genus of the cluster's binned salterns members (rollout kept MAGs and
                      co-assembly MAGs; meta/mags_three_way.tsv) -- the v1 10b rule. If the binned members have no
                      GTDB genus (novel genus) the read is 'unassigned'.
  2. kegg           : cluster without binned salterns members but with KEGG GENES members: majority GTDB R232 genus of
                      the members' genomes (KEGG organism -> assembly accession via KEGG REST gn:, 11b ->
                      release232 gtdb_taxonomy.tsv, accession number matched without version, GCA<->GCF). Members whose
                      genome is not in GTDB R232 keep the KEGG organism genus (NCBI name).
  3. transferred    : salterns-only cluster without any binned member (v1 genus_source 'none'): best DIAMOND blastp
                      hit (--sensitive) of the representative among marker proteins of the SAME family that carry a
                      genus (binned rollout / co-assembly MAG members, the 1,324 first_pass MAGs (11a), KEGG members);
                      the genus is taken only if the hit has >= 95 % identity, otherwise 'unassigned'.
  4. unassigned     : everything else.
genus_harmonized = GTDB R232 genus (GTDB suffixes such as _A kept); genus_source_name = the name as it came from the
reference (GTDB name for salterns / first_pass MAGs, KEGG organism genus for KEGG genes); Halomonas and Vreelandella
stay separate genera.

copies_per_GE (genomes per genome equivalent; depth-based as v1 10_summarize.py):
  depth(g, f) = reads(g, f) x r_s / L_eff(f)    for each of the 8 marker families f (L_eff = 3 x median aa + r_s - 180 + 1)
  copies_per_GE(g) = mean_f depth(g, f) / GE_depth_mean(s)        -> sums to 1 over all rows of a sample (incl. unassigned)
share_of_marker_reads = reads(g) / all marker reads of the sample (incl. the 'unassigned' row, so shares sum to 1);
share_of_assigned_marker_reads (extra column) = reads(g) / marker reads with a genus.
Outputs (results/): marker_genus_profile.tsv, marker_genus_profile_name_map.tsv, marker_rep_genus.tsv (per marker rep)
"""
import os, re, sys, gzip, subprocess
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
DM = "/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond"
A = "/scratch1/afreiburger/emilykin/processed/.snakemake_envs/544db00dd5c254ecfa7c6335967a046f_/bin"
KOLIST = "/scratch1/afreiburger/salterns_2026/kofam_db/ko_list"
os.chdir(W)
MARKERS = ["uL2", "uL11", "uS3", "uS11", "uS12", "uS9", "uL4", "uL16"]
MK_KO = {"K02886": "uL2", "K02867": "uL11", "K02982": "uS3", "K02948": "uS11", "K02950": "uS12", "K02996": "uS9",
         "K02926": "uL4", "K02930": "uL4", "K02878": "uL16", "K02866": "uL16"}
TRANSFER_MIN_PID = 95.0
samples = sorted(os.listdir("/scratch1/afreiburger/salterns_2026/rollout/reads"))


def gtdb_genus(lineage):
    m = re.search(r"g__([^;]*)", lineage or "")
    return m.group(1) if m else ""


def name_genus(org):
    w = org.replace("[", "").replace("]", "").split()
    if not w:
        return ""
    if w[0] in ("Candidatus", "uncultured") and len(w) > 1:
        return w[1]
    return w[0]


# ------------------------------------------------------------------ reference taxonomy
meta = pd.read_csv(f"{V1}/refs/labeled_meta.tsv", sep="\t", dtype=str).fillna("")
c99 = pd.read_csv(f"{V1}/refs/clu/clusters99.tsv", sep="\t", header=None, names=["label", "rep", "member"], dtype=str)
c99 = c99[c99.label.isin(MARKERS)]
m3 = pd.read_csv(f"{V1}/meta/mags_three_way.tsv", sep="\t", dtype=str).fillna("")
g_of_mag = dict(zip(m3.genome, m3.classification.map(gtdb_genus)))
org = pd.read_csv("meta/kegg_org_gtdb.tsv", sep="\t", dtype=str).fillna("").set_index("org")
mi = meta.set_index("pid")


def member_tax(pid):
    """-> (tier, genus_harmonized, genus_source_name) for one labelled member"""
    r = mi.loc[pid]
    if r.source == "salterns":
        if r.mag:
            g = g_of_mag.get(r.mag, r.genus)
            return "salterns_MAG", g, g
        return "salterns_unbinned", "", ""
    if r.source == "kegg":
        code = pid.split("|")[2].split(":")[0]
        kname = name_genus(r.organism)
        gg = org.at[code, "gtdb_genus"] if code in org.index else ""
        return "kegg", (gg if gg else kname), kname
    return "uniprot", r.genus, r.genus


rows = []
for rep, g in c99.groupby("rep"):
    tx = [member_tax(p) for p in g.member]
    lab = g.label.iloc[0]
    sal = [t for t in tx if t[0] == "salterns_MAG"]
    ext = [t for t in tx if t[0] in ("kegg", "uniprot")]
    if sal:
        named = [t for t in sal if t[1]]
        if named:
            gh = Counter(t[1] for t in named).most_common(1)[0][0]
            rows.append((rep, lab, "salterns_MAG", gh, gh, len(sal), np.nan, ""))
        else:
            rows.append((rep, lab, "salterns_MAG_genus_unnamed", "unassigned", "", len(sal), np.nan, ""))
    elif ext:
        gh = Counter(t[1] for t in ext if t[1]).most_common(1)[0][0] if any(t[1] for t in ext) else "unassigned"
        raw = Counter(t[2] for t in ext if t[1] == gh).most_common(1)[0][0] if gh != "unassigned" else ""
        rows.append((rep, lab, "kegg", gh, raw, len(ext), np.nan, ""))
    else:
        rows.append((rep, lab, "transfer_pending", "unassigned", "", 0, np.nan, ""))
rg = pd.DataFrame(rows, columns=["rep", "label", "tier", "genus_harmonized", "genus_source_name", "n_members_used",
                                 "transfer_pident", "transfer_hit"])
print(rg.tier.value_counts().to_string(), flush=True)

# ------------------------------------------------------------------ transfer for salterns-only unbinned clusters
os.makedirs("refs/mk_tx", exist_ok=True)
q = rg[rg.tier == "transfer_pending"]
open("refs/mk_tx/q.ids", "w").write("\n".join(q.rep) + "\n")
# reference proteins: labelled v1 marker members with a genus (binned salterns, KEGG) + first_pass marker proteins
ref_tax = {}
for pid, lab in zip(meta.pid, meta.label):
    if lab not in MARKERS:
        continue
    t = member_tax(pid)
    if t[0] in ("salterns_MAG", "kegg") and t[1]:
        ref_tax[pid] = (lab, t[1], t[2], t[0])
kl = pd.read_csv(KOLIST, sep="\t", dtype=str)
THR = {k: (float(t), st) for k, t, st in zip(kl.knum, kl.threshold, kl.score_type) if k in MK_KO}
best = {}
for l in open("hmm/firstpass_markers.tbl"):
    if l.startswith("#"):
        continue
    f = l.split()
    pid, ko, full, dom = f[0], f[2], float(f[5]), float(f[8])
    t, st = THR[ko]
    sc = full if st == "full" else dom
    if sc >= 0.5 * t and (pid not in best or full > best[pid][1]):
        best[pid] = (MK_KO[ko], full)
n_fp = 0
for pid, (lab, _) in best.items():
    genome = "first_pass|" + pid.split("|")[1]
    g = g_of_mag.get(genome, "")
    if g:
        ref_tax[pid] = (lab, g, g, "first_pass_MAG"); n_fp += 1
print("transfer refs:", len(ref_tax), "of which first_pass", n_fp, flush=True)
open("refs/mk_tx/r_v1.ids", "w").write("\n".join(p for p in ref_tax if not p.startswith("fp|")) + "\n")
open("refs/mk_tx/r_fp.ids", "w").write("\n".join(p for p in ref_tax if p.startswith("fp|")) + "\n")
subprocess.run(f"{A}/seqkit grep -f refs/mk_tx/q.ids {V1}/db/survey.faa > refs/mk_tx/q.faa && "
               f"{A}/seqkit grep -f refs/mk_tx/r_v1.ids {V1}/refs/labeled.faa > refs/mk_tx/r.faa && "
               f"{A}/seqkit grep -f refs/mk_tx/r_fp.ids genes_firstpass/firstpass.faa | sed 's/ # .*//' >> refs/mk_tx/r.faa && "
               f"nice -n 10 {DM} makedb --in refs/mk_tx/r.faa -d refs/mk_tx/r --quiet && "
               f"nice -n 10 {DM} blastp -q refs/mk_tx/q.faa -d refs/mk_tx/r -o refs/mk_tx/hits.tsv -k 50 -e 1e-5 --sensitive "
               f"--threads 8 --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore", shell=True, check=True)
h = pd.read_csv("refs/mk_tx/hits.tsv", sep="\t", header=None, names="q s pident len e bits".split(), dtype={"q": str, "s": str})
lab_of_rep = dict(zip(rg.rep, rg.label))
h = h[[ref_tax[s][0] == lab_of_rep[qq] for qq, s in zip(h.q, h.s)]]       # same marker family only
h = h.sort_values(["q", "bits", "pident"], ascending=[True, False, False]).drop_duplicates("q")
tr = {qq: (s, p) for qq, s, p in zip(h.q, h.s, h.pident)}
for i in rg.index[rg.tier == "transfer_pending"]:
    rep = rg.at[i, "rep"]
    if rep in tr:
        s, p = tr[rep]
        rg.at[i, "transfer_pident"] = p
        rg.at[i, "transfer_hit"] = s
        if p >= TRANSFER_MIN_PID:
            lab, gh, raw, src = ref_tax[s]
            rg.at[i, "tier"] = f"transferred_{src}"
            rg.at[i, "genus_harmonized"] = gh
            rg.at[i, "genus_source_name"] = raw
            continue
    rg.at[i, "tier"] = "unassigned_transfer_lt95"
rg.to_csv("results/marker_rep_genus.tsv", sep="\t", index=False)
print(rg.tier.value_counts().to_string(), flush=True)

# ------------------------------------------------------------------ per-sample profile
ge = pd.read_csv(f"{V1}/results/genome_equivalents.tsv", sep="\t", index_col=0)
LL = pd.read_csv(f"{V1}/refs/label_lengths.tsv", sep="\t").set_index("label").median_aa
cut = pd.read_csv(f"{V1}/cutoffs.tsv", sep="\t", dtype=str)
LMIN = int(dict(zip(cut.param, cut.value))["LMIN"])
R = rg.set_index("rep")
out, nm_rows = [], []
for s in samples:
    rd = pd.read_csv(f"{V1}/results/per_sample/{s}.reads.tsv.gz", sep="\t", dtype={"read": str, "mate": str})
    rd = rd[rd.label.isin(MARKERS)].join(R[["tier", "genus_harmonized", "genus_source_name"]], on="rep")
    assert rd.genus_harmonized.notna().all(), s
    cnt = pd.read_csv(f"{V1}/results/per_sample/{s}.counts.tsv", sep="\t").set_index("label").reads
    assert len(rd) == int(cnt.reindex(MARKERS).sum()), s           # every counted marker read is in the per-read table
    rs = ge.at[s, "read_len"]
    leff = {f: 3 * LL[f] + rs - 6 * LMIN + 1 for f in MARKERS}
    gdm = ge.at[s, "GE_depth_mean"]
    ntot = len(rd)
    nassigned = (rd.genus_harmonized != "unassigned").sum()
    for gh, g in rd.groupby("genus_harmonized"):
        per_f = g.groupby("label").size().reindex(MARKERS, fill_value=0)
        cpg = np.mean([per_f[f] * rs / leff[f] for f in MARKERS]) / gdm
        names = g[g.genus_source_name != ""].groupby("genus_source_name").size().sort_values(ascending=False)
        tiers = g.tier.str.replace(r"^transferred_.*", "transferred", regex=True).value_counts()
        out.append(dict(sample=s, genus_harmonized=gh, genus_source_name=";".join(names.index),
                        marker_reads=len(g), share_of_marker_reads=len(g) / ntot, copies_per_GE=cpg,
                        share_of_assigned_marker_reads=(len(g) / nassigned) if gh != "unassigned" else np.nan,
                        reads_best_hit_ge95=int((g.pident >= 95).sum()),
                        reads_tier_salterns_MAG=int(tiers.get("salterns_MAG", 0)), reads_tier_kegg=int(tiers.get("kegg", 0)),
                        reads_tier_transferred=int(tiers.get("transferred", 0)),
                        n_marker_families_with_reads=int((per_f > 0).sum())))
        for nmn, k in names.items():
            nm_rows.append((s, gh, nmn, k))
    print(s, ntot, round(nassigned / ntot, 3), flush=True)
p = pd.DataFrame(out).sort_values(["sample", "marker_reads"], ascending=[True, False])
chk = p.groupby("sample").copies_per_GE.sum()
assert (abs(chk - 1) < 1e-9).all(), chk
p.to_csv("results/marker_genus_profile.tsv", sep="\t", index=False, float_format="%.6g")
nm = pd.DataFrame(nm_rows, columns=["sample", "genus_harmonized", "genus_source_name", "reads"])
nmt = nm.groupby(["genus_harmonized", "genus_source_name"]).reads.sum().reset_index()
src = rg.groupby(["genus_harmonized", "genus_source_name"]).tier.agg(lambda x: ";".join(sorted(set(x)))).reset_index()
nmt = nmt.merge(src, on=["genus_harmonized", "genus_source_name"], how="left")
nmt["renamed"] = nmt.genus_harmonized != nmt.genus_source_name
nmt.sort_values("reads", ascending=False).to_csv("results/marker_genus_profile_name_map.tsv", sep="\t", index=False)
print("done", len(p), "rows")
