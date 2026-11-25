#!/usr/bin/env python
"""Step 4: label proteins and assemble the read-survey protein reference.

  ext    <refs_raw> -> refs/ext.faa + refs/ext_meta.tsv            (external refs: KEGG GENES + UniProt families)
  status            -> meta/contig_status.tsv.gz                   (sample, contig, length, status, bin, mag_id)
  label             -> refs/labeled.faa + refs/labeled_meta.tsv    (after 03_panel_hmmsearch.sh on every input)

Labelling (KOfam 2026-10-26 thresholds from kofam_db/ko_list; DRAM best-KO rule):
  * a KO passes if its threshold is numeric and (score_type full: full-sequence score >= threshold;
    domain: best domain score >= threshold); the protein's KO = the passing KO with the lowest full E-value
    (ties: higher score).
  * KEGG GENES refs: the KEGG KO is the label (KEGG orthology, not KOfam). Eukaryote genes of marker KOs -> 'euk_marker'.
    K10944 -> pmoA / amoA / cuMMO_other by organism taxonomy (see AMO_* below).
  * UniProt family refs: KOfam passing KO; if none passes, a family default (FAMILY_DEFAULT) applies, else excluded.
  * salterns proteins (39 assemblies + 520 co-assembly MAGs): KOfam passing KO; markers: if nothing passes, the
    top-scoring panel KO is a marker KO and its score >= 0.5 x threshold -> that marker (archaeal ribosomal proteins
    fall below several KOfam thresholds). Decoy-labelled salterns proteins are kept only if they also hit a target
    profile at full E <= 1e-5 (i.e. they are homologs of a target). Proteins that hit a target profile but pass no
    panel KO are listed in refs/unassigned_homologs.tsv; those that are complete genes (partial=00) and score < 0.5 x
    the target threshold become a '<gene>_like' decoy class (e.g. ABC ATPases that are not PhnK/PhnL); the rest
    (0.5-1 x threshold, or partial genes) are excluded.
  * mpnS / hepD (KOfam threshold '-'): top-scoring panel KO is K18049 or K12905, full E <= 1e-30, HMM coverage >= 0.5,
    nothing else passes -> mpnS or hepD (whichever scores higher).
"""
import os, sys, re, gzip, glob, csv
from collections import defaultdict
import pandas as pd

W = "/scratch1/afreiburger/salterns_2026/read_markers"
R = "/scratch1/afreiburger/salterns_2026/rollout"
KOLIST = "/scratch1/afreiburger/salterns_2026/kofam_db/ko_list"
os.chdir(W)
panel = pd.read_csv("meta/ko_panel.tsv", sep="\t", dtype=str).fillna("")
KO2GENE = dict(zip(panel.ko, panel.gene))
KO2ROLE = dict(zip(panel.ko, panel.role))
TARGET_KOS = set(panel.ko[panel.role == "target"])
MARKER_KOS = set(panel.ko[panel.role == "marker"])
kl = pd.read_csv(KOLIST, sep="\t", dtype=str)
THR = {k: (float(t) if t not in ("-", "") else None, st) for k, t, st in zip(kl.knum, kl.threshold, kl.score_type)}
PLEN = {}
for line in open("db/panel.hmm") if os.path.exists("db/panel.hmm") else []:
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
FAMILY_DEFAULT = {  # uniprot set -> (KO whose top hit with full E <= 1e-10 gives the label, label if nothing passes)
    "pfam_PF02249": ("K00399", "mcrA"), "pfam_PF02745": ("K00399", "mcrA"),
    "pfam_PF02241": ("K00401", "mcrB"), "pfam_PF02783": ("K00401", "mcrB"), "pfam_PF02240": ("K00402", "mcrG"),
    "pfam_PF06007": ("K06163", "phnJ"), "pfam_PF12176": ("K04480", "mtaB"),
    "pfam_PF09505": ("K16178", "mtbB"), "pfam_PF05369": ("K16176", "mtmB"),
}
FAMILY_DECOY = {  # uniprot set -> (target KO that must NOT be the top hit, decoy label for non-passing members)
    "pfam_PF02332": ("K16157", "SDIMO_other"), "pfam_PF02901": ("K20038", "GRE_other"),
    "pfam_PF00848": ("K22443", "Rieske_other"), "pfam_PF06253": ("K14083", "mttB_family_other"),
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


# ----------------------------------------------------------------------------------------------------------- ext
def cmd_ext(raw):
    os.makedirs("refs", exist_ok=True)
    # KEGG organism taxonomy from br08601 (the gene lists are re-derived from the fetched .faa files, because two
    # fetch instances shared the marker KOs; see README)
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
    meta = []
    with open("refs/ext.faa", "w") as fo:
        seen = set()
        for ko in panel.ko:
            fp = f"{raw}/kegg/{ko}.faa"
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
        for fp in sorted(glob.glob(f"{raw}/uniprot/*.tsv")):
            st = os.path.basename(fp)[:-4]
            u = pd.read_csv(fp, sep="\t", dtype=str).fillna("")
            for _, r in u.iterrows():
                pid = f"up|{st}|{r.Entry}"
                lin = r["Taxonomic lineage"]
                ranks = dict((m.group(2), m.group(1).strip()) for m in re.finditer(r"([^,]+?) \(([a-z ]+)\)", lin))
                fo.write(f">{pid}\n{r.Sequence}\n")
                meta.append(dict(pid=pid, source="uniprot", set=st, domain=ranks.get("domain", ranks.get("superkingdom", "")),
                                 phylum=ranks.get("phylum", ""), family=ranks.get("family", ""),
                                 genus=ranks.get("genus", genus_of_name(r.Organism)), organism=r.Organism, lineage=lin))
    pd.DataFrame(meta).to_csv("refs/ext_meta.tsv", sep="\t", index=False)
    print("ext proteins", len(meta))


# -------------------------------------------------------------------------------------------------------- status
def cmd_status():
    rows = []
    mk = pd.read_csv(f"{R}/results/mags_kept.tsv", sep="\t", dtype=str)
    for s in sorted(os.listdir(f"{R}/reads")):
        lens = {}
        for l in open(f"{R}/{s}/asm/contigs.fa"):
            if l.startswith(">"):
                c = l[1:].split()[0]
                m = re.search(r"_length_(\d+)_", c)
                lens[c] = int(m.group(1)) if m else -1
        bins = pd.read_csv(f"{R}/{s}/keep/bins.tsv", sep="\t", dtype=str)
        c2 = {}
        for _, b in bins.iterrows():
            for l in open(b.fasta):
                if l.startswith(">"):
                    c = l[1:].split()[0]
                    assert c not in c2, (s, c)
                    c2[c] = (b.bin, b.mag_id if b.keep == "yes" else "", "kept_mag" if b.keep == "yes" else "other_bin")
        kept_ids = set(mk.mag_id[mk["sample"] == s])
        assert kept_ids == set(bins.mag_id[bins.keep == "yes"]), s
        for c, L in lens.items():
            b, m, st = c2.get(c, ("", "", "unbinned"))
            rows.append((s, c, L, st, b, m))
        assert all(c in lens for c in c2), s
    df = pd.DataFrame(rows, columns=["sample", "contig", "length", "status", "bin", "mag_id"])
    df.to_csv("meta/contig_status.tsv.gz", sep="\t", index=False)
    print(df.groupby(["sample", "status"]).length.sum().unstack().to_string())


# --------------------------------------------------------------------------------------------------------- label
def parse_hmm(base):
    """-> {pid: {ko: [fullE, fullScore, bestDom, cov]}}"""
    h = defaultdict(dict)
    for l in open(f"hmm/{base}.tbl"):
        if l.startswith("#"):
            continue
        f = l.split()
        h[f[0]][f[2]] = [float(f[4]), float(f[5]), 0.0, 0.0]
    spans = defaultdict(list)
    for l in open(f"hmm/{base}.domtbl"):
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


def kofam_call(hits):
    """hits {ko: [E, score, dom, cov]} -> (best passing KO or None, top-scoring KO)"""
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
    if best is None and top in ("K18049", "K12905"):
        e, sc, dom, cov = hits[top]
        if e <= 1e-30 and cov >= 0.5:
            return KO2GENE[top]
    return None


def cmd_label():
    km = pd.read_csv("refs/ext_meta.tsv", sep="\t", dtype=str).fillna("").set_index("pid")
    status = pd.read_csv("meta/contig_status.tsv.gz", sep="\t", dtype=str).fillna("")
    status["key"] = status["sample"] + "|" + status["contig"]
    st_of = dict(zip(status.key, zip(status.status, status.mag_id)))
    m3 = pd.read_csv("meta/mags_three_way.tsv", sep="\t", dtype=str).fillna("").set_index("genome")
    out_meta, unassigned = [], []
    fo = open("refs/labeled.faa", "w")

    def emit(pid, seq, label, how, extra):
        fo.write(f">{pid}\n{seq.rstrip('*')}\n")
        d = dict(pid=pid, label=label, how=how)
        d.update(extra)
        out_meta.append(d)

    # external
    h = parse_hmm("ext")
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
            elif ko == "K10944":
                emit(pid, seq, amo_label(taxs), "kegg_ko+taxonomy", extra)
            else:
                emit(pid, seq, KO2GENE[ko], "kegg_ko", extra)
            continue
        st = r.set
        if st in ("name_mpnS", "name_hepD"):
            lab = "mpnS" if st == "name_mpnS" else "hepD"
            emit(pid, seq, lab, "uniprot_name", extra)
            continue
        if st == "pfam_PF02461":
            emit(pid, seq, amo_label(taxs), "pfam_family+taxonomy", extra)
            continue
        if best:
            emit(pid, seq, KO2GENE[best] if KO2ROLE[best] != "marker" else KO2GENE[best], "kofam_pass", extra)
        elif st in FAMILY_DEFAULT and top == FAMILY_DEFAULT[st][0] and hits[top][0] <= 1e-10:
            emit(pid, seq, FAMILY_DEFAULT[st][1], "family_default", extra)
        elif st in FAMILY_DECOY and top != FAMILY_DECOY[st][0]:
            emit(pid, seq, FAMILY_DECOY[st][1], "family_decoy", extra)
        else:
            unassigned.append((pid, st, top or "", hits[top][1] if top else ""))

    # salterns: 39 assemblies + co-assembly MAGs
    inputs = sorted(os.listdir(f"{R}/reads")) + ["coasm520"]
    for base in inputs:
        h = parse_hmm(base)
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
                    lab, how = sp, "mpnS_hepD_rule"
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
                        # complete genes only (checked below): clearly below-threshold homolog -> '<gene>_like' decoy
                        lab, how = KO2GENE[tk] + "_like", "below_half_threshold_homolog"
            if lab:
                need[pid] = (lab, how, best or "", top or "")
        faa = f"genes/{base}/{base}.faa"
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
            if lab == "pmoA":   # K10944 salterns proteins: resolved later against the labelled external set
                lab = "K10944_unresolved"
            emit(pid, seq, lab, how, extra)
        print(base, len(need), flush=True)
    fo.close()
    md = pd.DataFrame(out_meta)
    md.to_csv("refs/labeled_meta.tsv", sep="\t", index=False)
    pd.DataFrame(unassigned, columns=["pid", "set", "top_target_ko", "score"]).to_csv("refs/unassigned_homologs.tsv", sep="\t", index=False)
    print(md.groupby(["label", "source"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "ext":
        cmd_ext(sys.argv[2])
    elif c == "status":
        cmd_status()
    elif c == "label":
        cmd_label()
