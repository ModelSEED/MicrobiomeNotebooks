#!/usr/bin/env python
"""Step 4 (v2): label proteins and assemble the read-survey protein reference. v1 04_build_db.py, extended:

  ext    <refs_raw> -> refs/ext.faa + refs/ext_meta.tsv        (KEGG GENES prokaryote + eukaryote sets, UniProt sets)
  label             -> refs/labeled.faa + refs/labeled_meta.tsv (after 03_panel_hmmsearch_v2.sh, salterns + ext)

Unchanged from v1: KOfam 2026-10-26 thresholds, DRAM best-KO rule (lowest full E among passing panel KOs), KEGG KO =
label of a KEGG gene, eukaryote genes of marker KOs -> euk_marker, K10944 by taxonomy, v1 UniProt family rules,
salterns rules (marker relaxed 0.5 x threshold; decoy-labelled salterns proteins kept only if they hit a target profile;
complete target homologs < 0.5 x threshold -> '<gene>_like'; 0.5-1 x 'twilight' homologs left out), mpnS/hepD rule.
v2 changes:
  * panel = meta/ko_panel_v2.tsv (286 KOs). Salterns hits = v1 hmm/<s>.tbl (62 v1 profiles) + hmm_new/<s>.<chunk>.tbl
    (224 new profiles), same -Z, so E-values are comparable.
  * eukaryote KEGG genes of counted target KOs (refs_raw/kegg_euk/) -> '<gene>_euk' (decoy_reported).
  * Rieske family rule (pfam_PF00848): 'Rieske_other' unless the top profile is any Rieske target (cntA, CMO, gbcA, msmA).
  * UniProt reviewed Pfam families of the generic targets (rev_PFxxxxx): member passes a panel KO -> that label; top
    profile is one of the family's target KOs but nothing passes -> left out (twilight); else '<family>_other' (decoy).
  * KOs without a KOfam threshold: mpnS / hepD keep the v1 rule (top panel profile, full E <= 1e-30, HMM cov >= 0.5);
    dmoA (K16967), phnY (K23996), phnZ (K05915) get a threshold (full score) from meta/synthetic_thresholds.tsv
    (04b: 1 + the highest score of any other panel KEGG gene on that profile) and are then treated like KOfam KOs.
"""
import os, sys, re, gzip, glob, csv
from collections import defaultdict
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers_v2"
V1 = "/scratch1/afreiburger/salterns_2026/read_markers"
R = "/scratch1/afreiburger/salterns_2026/rollout"
KOLIST = "/scratch1/afreiburger/salterns_2026/kofam_db/ko_list"
os.chdir(W)
panel = pd.read_csv("meta/ko_panel_v2.tsv", sep="\t", dtype=str).fillna("")
KO2GENE = dict(zip(panel.ko, panel.gene))
KO2ROLE = dict(zip(panel.ko, panel.role))
TARGET_KOS = set(panel.ko[panel.role == "target"])
MARKER_KOS = set(panel.ko[panel.role == "marker"])
kl = pd.read_csv(KOLIST, sep="\t", dtype=str)
THR = {k: (float(t) if t not in ("-", "") else None, st) for k, t, st in zip(kl.knum, kl.threshold, kl.score_type)}
NOTHR_V1 = {"K18049", "K12905"}                       # mpnS, hepD: v1 rule
SYN = {}
if os.path.exists("meta/synthetic_thresholds.tsv"):
    st_ = pd.read_csv("meta/synthetic_thresholds.tsv", sep="\t", dtype=str)
    SYN = {k: float(v) for k, v in zip(st_.ko, st_.min_full_score)}
# v2: the derived floors of dmoA / phnY / phnZ act as KOfam thresholds (full score), so these KOs compete in the
# best-KO rule (DmoA E9JFX9 also passes the LadA K20938 threshold; by E-value it is DmoA)
for _k, _v in SYN.items():
    THR[_k] = (_v, "full")
NOTHR_V2 = {k for k in TARGET_KOS if THR.get(k, (None,))[0] is None} - NOTHR_V1     # (empty once SYN is set)
PLEN = {}
for line in open("db/panel_v2.hmm"):
    if line.startswith("NAME"):
        cur = line.split()[1]
    elif line.startswith("LENG"):
        PLEN[cur] = int(line.split()[1])

AMO_AMOA = ("Nitrosomonas", "Nitrosospira", "Nitrosococcus", "Nitrosovibrio", "Nitrosolobus", "Nitrospira",
            "Nitrosopumilus", "Nitrososphaera", "Nitrosotalea", "Nitrosocaldus", "Nitrosarchaeum", "Nitrosotenuis",
            "Nitrosocosmicus", "Nitrosopelagicus", "Nitrosomarinus", "Cenarchaeum")
AMO_AMOA_HIGHER = ("Nitrosomonadaceae", "Nitrosococcaceae", "Nitrosopumilaceae", "Nitrososphaeraceae",
                   "Nitrosocaldaceae", "Nitrososphaerota", "Thaumarchaeota", "Nitrososphaeria")
AMO_PMOA = ("Methylococcus", "Methylomonas", "Methylobacter", "Methylomicrobium", "Methylosarcina", "Methylocaldum",
            "Methylothermus", "Methylohalobius", "Methylomarinum", "Methyloprofundus", "Methylovulum", "Methylosoma",
            "Methyloglobulus", "Methylogaea", "Methylomagnum", "Methylotuvimicrobium", "Methylosinus", "Methylocystis",
            "Methylocapsa", "Methylocella", "Methyloferula", "Methylacidiphilum", "Methylacidimicrobium",
            "Crenothrix", "Methylomirabilis", "Methylicorpusculum", "Methylospira", "Methyloterricola",
            "Methylococcaceae", "Methylomonadaceae", "Methylocystaceae", "Methylacidiphilaceae")
FAMILY_DEFAULT = {
    "pfam_PF02249": ("K00399", "mcrA"), "pfam_PF02745": ("K00399", "mcrA"),
    "pfam_PF02241": ("K00401", "mcrB"), "pfam_PF02783": ("K00401", "mcrB"), "pfam_PF02240": ("K00402", "mcrG"),
    "pfam_PF06007": ("K06163", "phnJ"), "pfam_PF12176": ("K04480", "mtaB"),
    "pfam_PF09505": ("K16178", "mtbB"), "pfam_PF05369": ("K16176", "mtmB"),
}
FAMILY_DECOY = {   # set -> (target KOs that must NOT be the top hit, decoy label for non-passing members)
    "pfam_PF02332": ({"K16157"}, "SDIMO_other"), "pfam_PF02901": ({"K20038"}, "GRE_other"),
    "pfam_PF00848": ({"K22443", "K00499", "K00479", "K16968"}, "Rieske_other"),
    "pfam_PF06253": ({"K14083"}, "mttB_family_other"),
    # v2: reviewed (Swiss-Prot) members of the generic targets' Pfam families
    "rev_PF00763": ({"K01491"}, "THF_DHG_CYH_other"), "rev_PF02882": ({"K01491"}, "THF_DHG_CYH_other"),
    "rev_PF00551": ({"K01433"}, "Formyl_trans_other"), "rev_PF01268": ({"K01938"}, "FTHFS_other"),
    "rev_PF00464": ({"K00600"}, "SHMT_other"), "rev_PF03313": ({"K01752"}, "SDH_other"),
    "rev_PF03315": ({"K01752"}, "SDH_other"), "rev_PF00389": ({"K00122"}, "2-Hacid_dh_other"),
    "rev_PF02826": ({"K00122"}, "2-Hacid_dh_other"), "rev_PF00909": ({"K03320"}, "Ammonium_transp_other"),
    "rev_PF00756": ({"K01070"}, "Esterase_other"), "rev_PF08240": ({"K00121", "K00148", "K17067"}, "MDR_other"),
    "rev_PF00107": ({"K00121", "K00148", "K17067"}, "MDR_other"), "rev_PF00296": ({"K04091", "K17228", "K16967", "K00320"}, "LLM_other"),
    "rev_PF00871": ({"K00925"}, "Acetate_kinase_other"), "rev_PF01515": ({"K00625", "K13788", "K15024", "K04020"}, "PTA_PTB_other"),
    "rev_PF00501": ({"K01895", "K20034"}, "AMP-binding_other"),
}


def genus_of_name(org):
    w = org.replace("[", "").replace("]", "").split()
    if not w:
        return ""
    if w[0] in ("Candidatus", "uncultured") and len(w) > 1:
        return w[1]
    return w[0]


def amo_label(taxstr):
    for g in AMO_PMOA:
        if g in taxstr:
            return "pmoA"
    for g in AMO_AMOA + AMO_AMOA_HIGHER:
        if g in taxstr:
            return "amoA"
    return "cuMMO_other"


def read_fasta(p, full_header=False):
    op = gzip.open if p.endswith(".gz") else open
    name, seq = None, []
    with op(p, "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = (l[1:].rstrip("\n") if full_header else l[1:].split()[0]), []
            else:
                seq.append(l.strip())
    if name:
        yield name, "".join(seq)


def br08601(raw):
    tax = {}
    A = B = C = D = ""
    for line in open(f"{raw}/kegg/br08601.txt"):
        if line.startswith("A"):
            A = re.sub(r"\s*\(\d+\)$", "", line[1:].strip())
        elif line.startswith("B"):
            B = re.sub(r"\s*\(\d+\)$", "", line[1:].strip()); C = D = ""
        elif line.startswith("C"):
            C = re.sub(r"\s*\(\d+\)$", "", line[1:].strip()); D = ""
        elif line.startswith("D"):
            D = re.sub(r"\s*\(\d+\)$", "", line[1:].strip())
        elif line.startswith("E"):
            f = line[1:].strip().split(None, 1)
            if len(f) == 2:
                tax[f[0]] = (A, B, C, D, f[1])
    return tax


# ----------------------------------------------------------------------------------------------------------- ext
def cmd_ext(raw):
    os.makedirs("refs", exist_ok=True)
    tax = br08601(raw)
    meta = []
    seen = set()
    with open("refs/ext.faa", "w") as fo:
        for sub in ("kegg", "kegg_euk"):
            for ko in panel.ko:
                fp = f"{raw}/{sub}/{ko}.faa"
                if not os.path.exists(fp):
                    continue
                for g, s in read_fasta(fp):
                    pid = f"kegg|{ko}|{g}"
                    if pid in seen or not s:
                        continue
                    seen.add(pid)
                    a, b, c, d, org = tax[g.split(":")[0]]
                    dom = b if a == "Prokaryotes" else "Eukaryota"
                    fo.write(f">{pid}\n{s.rstrip('*')}\n")
                    meta.append(dict(pid=pid, source="kegg", set=ko, domain=dom, phylum=c, family=d,
                                     genus=genus_of_name(org), organism=org, lineage=f"{a};{b};{c};{d}"))
        seen_acc = set()
        for fp in sorted(glob.glob(f"{raw}/uniprot/*.tsv")):
            st = os.path.basename(fp)[:-4]
            u = pd.read_csv(fp, sep="\t", dtype=str).fillna("")
            for _, r in u.iterrows():
                if st.startswith("rev_"):          # a Swiss-Prot entry can sit in two families of one target: keep once
                    if r.Entry in seen_acc:
                        continue
                    seen_acc.add(r.Entry)
                pid = f"up|{st}|{r.Entry}"
                lin = r["Taxonomic lineage"]
                ranks = dict((m.group(2), m.group(1).strip()) for m in re.finditer(r"([^,]+?) \(([a-z ]+)\)", lin))
                fo.write(f">{pid}\n{r.Sequence}\n")
                meta.append(dict(pid=pid, source="uniprot", set=st, domain=ranks.get("domain", ranks.get("superkingdom", "")),
                                 phylum=ranks.get("phylum", ""), family=ranks.get("family", ""),
                                 genus=ranks.get("genus", genus_of_name(r.Organism)), organism=r.Organism, lineage=lin))
    pd.DataFrame(meta).to_csv("refs/ext_meta.tsv", sep="\t", index=False)
    print("ext proteins", len(meta))


# --------------------------------------------------------------------------------------------------------- label
def parse_hmm(tbls, domtbls):
    """-> {pid: {ko: [fullE, fullScore, bestDom, cov]}}"""
    h = defaultdict(dict)
    for p in tbls:
        for l in open(p):
            if l.startswith("#"):
                continue
            f = l.split()
            h[f[0]][f[2]] = [float(f[4]), float(f[5]), 0.0, 0.0]
    spans = defaultdict(list)
    for p in domtbls:
        for l in open(p):
            if l.startswith("#"):
                continue
            f = l.split()
            pid, ko = f[0], f[3]
            if ko in h.get(pid, {}):
                h[pid][ko][2] = max(h[pid][ko][2], float(f[13]))
                spans[(pid, ko)].append((int(f[15]), int(f[16])))
    for (pid, ko), sp in spans.items():
        cov = set()
        for a, b in sp:
            cov.update(range(a, b + 1))
        h[pid][ko][3] = len(cov) / PLEN.get(ko, max(cov))
    return h


def salterns_hmm(base):
    t = [f"{V1}/hmm/{base}.tbl"] + sorted(glob.glob(f"hmm_new/{base}.[0-9][0-9].tbl"))
    d = [f"{V1}/hmm/{base}.domtbl"] + sorted(glob.glob(f"hmm_new/{base}.[0-9][0-9].domtbl"))
    assert len(t) == 5 and len(d) == 5, (base, t)
    return parse_hmm(t, d)


def kofam_call(hits):
    passing = []
    for ko, (e, sc, dom, cov) in hits.items():
        t, st = THR.get(ko, (None, ""))
        if t is None:
            continue
        if (sc if st == "full" else dom) >= t:
            passing.append((e, -sc, ko))
    top = max(hits, key=lambda k: hits[k][1]) if hits else None
    return (min(passing)[2] if passing else None), top


def special_label(hits, best, top):
    if best is None and top in NOTHR_V1 | NOTHR_V2:
        e, sc, dom, cov = hits[top]
        if e <= 1e-30 and cov >= 0.5 and (top not in SYN or sc >= SYN[top]):
            return KO2GENE[top]
    return None


def cmd_label():
    km = pd.read_csv("refs/ext_meta.tsv", sep="\t", dtype=str).fillna("").set_index("pid")
    status = pd.read_csv(f"{V1}/meta/contig_status.tsv.gz", sep="\t", dtype=str).fillna("")
    status["key"] = status["sample"] + "|" + status["contig"]
    st_of = dict(zip(status.key, zip(status.status, status.mag_id)))
    m3 = pd.read_csv(f"{V1}/meta/mags_three_way.tsv", sep="\t", dtype=str).fillna("").set_index("genome")
    out_meta, unassigned = [], []
    fo = open("refs/labeled.faa", "w")

    def emit(pid, seq, label, how, extra):
        fo.write(f">{pid}\n{seq.rstrip('*')}\n")
        d = dict(pid=pid, label=label, how=how)
        d.update(extra)
        out_meta.append(d)

    h = parse_hmm(["hmm_new/ext.tbl"], ["hmm_new/ext.domtbl"])
    for pid, seq in read_fasta("refs/ext.faa"):
        r = km.loc[pid]
        hits = h.get(pid, {})
        best, top = kofam_call(hits)
        taxs = f"{r.lineage};{r.organism}"
        extra = dict(source=r.source, set=r.set, domain=r.domain, phylum=r.phylum, family=r.family, genus=r.genus,
                     organism=r.organism, kofam_best=best or "", kofam_top=top or "", mag="", species_cluster="",
                     sample="", contig_status="")
        if r.source == "kegg":
            ko = r.set
            if KO2ROLE[ko] == "marker" and r.domain == "Eukaryota":
                emit(pid, seq, "euk_marker", "kegg_ko_eukaryote", extra)
            elif KO2ROLE[ko] == "target" and r.domain == "Eukaryota":
                emit(pid, seq, KO2GENE[ko] + "_euk", "kegg_ko_eukaryote", extra)
            elif r.domain == "Eukaryota":
                continue                     # (not expected: decoy KOs are fetched for prokaryotes only)
            elif ko == "K10944":
                emit(pid, seq, amo_label(taxs), "kegg_ko+taxonomy", extra)
            else:
                emit(pid, seq, KO2GENE[ko], "kegg_ko", extra)
            continue
        st = r.set
        if st in ("name_mpnS", "name_hepD"):
            emit(pid, seq, "mpnS" if st == "name_mpnS" else "hepD", "uniprot_name", extra)
            continue
        if st == "pfam_PF02461":
            emit(pid, seq, amo_label(taxs), "pfam_family+taxonomy", extra)
            continue
        sp = special_label(hits, best, top) if st.startswith("rev_") else None
        if best:
            emit(pid, seq, KO2GENE[best], "kofam_pass", extra)
        elif sp:                              # v2: reviewed members of a no-threshold KO (e.g. DmoA E9JFX9)
            emit(pid, seq, sp, "no_threshold_rule", extra)
        elif st in FAMILY_DEFAULT and top == FAMILY_DEFAULT[st][0] and hits[top][0] <= 1e-10:
            emit(pid, seq, FAMILY_DEFAULT[st][1], "family_default", extra)
        elif st in FAMILY_DECOY and top not in FAMILY_DECOY[st][0]:
            emit(pid, seq, FAMILY_DECOY[st][1], "family_decoy", extra)
        else:
            unassigned.append((pid, st, top or "", hits[top][1] if top else ""))

    inputs = sorted(os.listdir(f"{R}/reads")) + ["coasm520"]
    for base in inputs:
        h = salterns_hmm(base)
        need = {}
        for pid, hits in h.items():
            best, top = kofam_call(hits)
            has_target_hit = any(k in TARGET_KOS and v[0] <= 1e-5 for k, v in hits.items())
            lab, how = None, None
            if best:
                role = KO2ROLE[best]
                if role in ("target", "marker", "decoy_reported"):
                    lab, how = KO2GENE[best], "kofam_pass"
                elif role == "decoy" and has_target_hit:
                    lab, how = KO2GENE[best], "kofam_pass_decoy_homolog"
            else:
                sp = special_label(hits, best, top)
                if sp:
                    lab, how = sp, "no_threshold_rule"
                elif top in MARKER_KOS:
                    t, stype = THR[top]
                    e, sc, dom, cov = hits[top]
                    if (sc if stype == "full" else dom) >= 0.5 * t:
                        lab, how = KO2GENE[top], "marker_relaxed_0.5thr"
                if lab is None and has_target_hit:
                    tk = min((v[0], k) for k, v in hits.items() if k in TARGET_KOS)[1]
                    unassigned.append((pid, base, tk, hits[tk][1]))
                    t, stype = THR.get(tk, (None, ""))
                    if t is not None and (hits[tk][1] if stype == "full" else hits[tk][2]) < 0.5 * t:
                        lab, how = KO2GENE[tk] + "_like", "below_half_threshold_homolog"
            if lab:
                need[pid] = (lab, how, best or "", top or "")
        faa = f"{V1}/genes/{base}/{base}.faa"
        for hdr, seq in read_fasta(faa, full_header=True):
            pid = hdr.split()[0]
            if pid not in need:
                continue
            lab, how, best, top = need[pid]
            if how == "below_half_threshold_homolog" and "partial=00" not in hdr:
                continue
            if base == "coasm520":
                mag = pid.split("|")[1]
                g = "coassembly|" + mag
                sample, cst = "coasm", "kept_mag"
            else:
                sample, rest = pid.split("|", 1)
                contig = rest.rsplit("_", 1)[0]
                cst, mag = st_of[f"{sample}|{contig}"]
                g = "rollout|" + mag if mag else ""
            tx = m3.loc[g] if g and g in m3.index else None
            cls = tx.classification if tx is not None else ""
            rk = dict(x.split("__", 1) for x in cls.split(";") if "__" in x)
            extra = dict(source="salterns", set=base, domain=rk.get("d", ""), phylum=rk.get("p", ""),
                         family=rk.get("f", ""), genus=rk.get("g", ""), organism=rk.get("s", ""),
                         kofam_best=best, kofam_top=top, mag=g, species_cluster=tx.species_cluster if tx is not None else "",
                         sample=sample, contig_status=cst)
            if lab == "pmoA":
                lab = "K10944_unresolved"
            emit(pid, seq, lab, how, extra)
        print(base, len(need), flush=True)
    fo.close()
    md = pd.DataFrame(out_meta)
    md.to_csv("refs/labeled_meta.tsv", sep="\t", index=False)
    pd.DataFrame(unassigned, columns=["pid", "set", "top_target_ko", "score"]).to_csv("refs/unassigned_homologs.tsv", sep="\t", index=False)
    pd.set_option("display.max_rows", 1000)
    print(md.groupby(["label", "source"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "ext":
        cmd_ext(sys.argv[2])
    elif c == "label":
        cmd_label()
