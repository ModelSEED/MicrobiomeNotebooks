#!/usr/bin/env python
"""Annotation-level look at the archaeal MePn consumers and the four TMA-route candidates.

Targets (taken from guilds_r232.csv, not hard-coded):
  * archaeal MePn MAGs   - archaea in either MePn set (phnJ carriers, or the hand-curated 2020 set)
  * TMA-route candidates - tma_producer == True (cntA, or genes "clustered with choline TMA-lyase")
  * context              - the only methanogen (mcrA), for its methylamine genes

Per MAG it collects R232 taxonomy + closest reference/ANI (GTDB-Tk 2.8.0), CheckM, contigs (MAG_contigs/),
DNA-weighted relative abundance per core, guild flags, 2020 modelled-set membership and flux, the C-P lyase reactions
of its base ModelSEED model (matched by compound ID: cpd25960 MePn, cpd26518 RPnTP, cpd26519 PRPn), then inventories
genes from the RAST genome JSON by explicit function-string rules (RULES below), with locations and
neighbourhoods. Sequence checks back the interpretation. All use the RAST protein translations, a local
Smith-Waterman alignment (BLOSUM62, gap open -11 / extend -1), and a null of 500 random proteins drawn from
the 310 MAGs (seeded):
  (1) positional candidates for genes RAST left unannotated inside a C-P lyase locus (the gene upstream of
      phnH, the gene downstream of phnM), scored against RAST-annotated PhnG/PhnN/PhnP/PhnM of other MAGs;
  (2) gene-by-gene comparison of the archaeal C-P lyase loci;
  (3) contig-edge fragments of one gene that may be a single gene split by the assembly (the contigs are stitched at
      their overlap and the ORF translated across it);
  (4) lone phnM genes outside a locus vs the operon PhnMs; archaeal PhnJ vs the closest bacterial PhnJ;
  (5) glycine-reductase protein B paralogs of the TMA candidates vs betaine-type (GrdH/GrdI) and
      glycine-type (GrdB/GrdE) B proteins of other MAGs.
Genes that RAST split into two ORFs at an in-frame Sec (UGA) or Pyl (UAG) codon are found by reading the first ORF on
through that codon (recoded_splits) and are counted as one gene.

Outputs (this folder):
  archaeal_mepn_and_tma_candidates.md         narrative + per-MAG tables
  archaeal_mepn_and_tma_candidates_genes.csv  one row per gene hit (coordinates 1-based, start <= stop)

Run (from this folder):  ~/Documents/py_venv/bin/python make_archaeal_mepn_tma_report.py   (about 1-2 min)
"""
import collections
import datetime
import json
import math
import random
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
COR = HERE.parent
S = COR.parent
RAST = S / "jplfaria/Cliff MAGs/RASTified_MAG_json"
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
ZERO = 1e-5                       # detection level used throughout the 2020 code
N_NULL = 500                      # random proteins in the alignment null
SEED = 7
OUT_MD = HERE / "archaeal_mepn_and_tma_candidates.md"
OUT_CSV = HERE / "archaeal_mepn_and_tma_candidates_genes.csv"


def key(name):
    m = PAT.search(str(name))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


# ---- gene rules: (panel, gene label, regex on the lower-cased RAST function, exclude regex) ---------------
# Panel "phn" = phosphonate genes, "pi" = Pi-starvation genes, "tma" = choline/carnitine/betaine/TMA genes and
# the non-TMA routes that compete for the same substrates, "methylamine" = methanogenic methylamine genes.
# A label says what the RAST string supports, nothing more.
RULES = [
    ("phn", "phnC", r"phosphonate abc transporter atp-binding protein|\bphnc\b", r"2-aminoethylphosphonate"),
    ("phn", "phnD", r"phosphonate abc transporter (phosphate-binding periplasmic component|substrate-binding)|\bphnd\b"
                    r"|abc-type phosphate/phosphonate transport system,? periplasmic component", r"2-aminoethylphosphonate|\bptx"),
    ("phn", "phnE", r"phosphonate abc transporter permease protein|\bphne[12]?\b", r"2-aminoethylphosphonate|\bptx"),
    ("phn", "ptxB/C (phosphite ABC)", r"\bptx[bc]\b", None),
    ("phn", "phnCDE-like ABC (cluster 12)", r"cluster 12, methionine/phosphonates", None),
    ("phn", "AEP ABC transporter", r"2-aminoethylphosphonate abc transporter", None),
    ("phn", "phnF", r"transcriptional regulator phnf", None),
    ("phn", "phnG", r"\bphng\b", None),
    ("phn", "phnH", r"\bphnh\b", None),
    ("phn", "phnI", r"\bphni\b", None),
    ("phn", "phnJ", r"c-p lyase|ec 4\.7\.1\.1\b", None),
    ("phn", "phnK", r"\bphnk\b", None),
    ("phn", "phnL", r"\bphnl\b", None),
    ("phn", "phnM", r"methylphosphonate 5-triphosphate diphosphatase|ec 3\.6\.1\.63", None),
    ("phn", "phnM-like", r"metal-dependent hydrolase involved in phosphonate metabolism", None),
    ("phn", "phnN", r"\bphnn\b|ec 2\.7\.4\.23", None),
    ("phn", "phnO-like", r"phosphonate utilization associated acetyltransferase|\bphno\b", None),
    ("phn", "phnP", r"\bphnp\b", None),
    ("phn", "phnA (YjdM-type)", r"^protein phna$|alkylphosphonate utilization operon protein phna", None),
    ("phn", "phnA (phosphonoacetate hydrolase)", r"phosphonoacetate hydrolase|ec 3\.11\.1\.2\b", None),
    ("phn", "phnW", r"2-aminoethylphosphonate:pyruvate aminotransferase|ec 2\.6\.1\.37|\bphnw\b", None),
    ("phn", "phnX", r"phosphonoacetaldehyde hydrolase|ec 3\.11\.1\.1\b|\bphnx\b", None),
    ("phn", "phnY", r"phosphonoacetaldehyde dehydrogenase|\bphny\b", None),
    ("phn", "phnZ", r"hydroxyethylphosphonate dioxygenase|\bphnz\b", None),
    ("phn", "mpnS", r"methylphosphonate synthase|\bmpns\b", None),
    ("phn", "pepM", r"phosphoenolpyruvate (phospho)?mutase|ec 5\.4\.2\.9\b", None),
    ("phn", "PEP-mutase family (CPEP mutase-like)", r"carboxyvinyl-carboxyphosphonate phosphorylmutase", None),
    ("phn", "ptxD (phosphite DH)", r"phosphonate dehydrogenase|phosphite dehydrogenase|ec 1\.20\.1\.1\b", None),
    ("pi", "pstS", r"\bpsts\b", None),
    ("pi", "pstC", r"\bpstc\b", None),
    ("pi", "pstA", r"\bpsta\b", None),
    ("pi", "pstB", r"\bpstb\b", None),
    ("pi", "phoU", r"\bphou\b", None),
    ("pi", "phoB", r"\bphob\b", None),
    ("pi", "phoR", r"\bphor\b", None),
    ("pi", "Pi-regulatory homolog (unspecified)", r"^phosphate regulatory protein homolog$", None),
    ("tma", "cntA", r"carnitine monooxygenase, oxygenase component|\bcnta\b|\byeaw\b", None),
    ("tma", "cntB", r"carnitine monooxygenase, reductase component|\bcntb\b|\byeax\b", None),
    ("tma", "cutC", r"choline trimethylamine-lyase|ec 4\.3\.99\.4", r"clustered with|activat"),   # not "cutc": copper CutC
    ("tma", "cutD", r"choline trimethylamine-lyase activat|\bcutd\b", None),
    ("tma", "cutC-cluster label (BMC gene)", r"clustered with choline trimethylamine-lyase", None),
    ("tma", "glycyl radical enzyme", r"pyruvate formate-lyase \(|formate acetyltransferase|formate c-acetyltransferase"
                                     r"|glycyl radical|glycerol dehydratase \(b12-independent\)"
                                     r"|hydroxyphenylacetate decarboxylase", r"activat"),
    ("tma", "GRE-activase-like (COG1180)", r"pyruvate-formate lyase-activating enzyme like|pyruvate formate-lyase activating", None),
    ("tma", "eutBC (ethanolamine ammonia-lyase)", r"ethanolamine ammonia-lyase", None),
    ("tma", "BMC / eut-pdu gene", r"ethanolamine utilization|propanediol utilization|polyhedral|microcompartment",
     r"clustered with choline"),
    ("tma", "BCCT transporter (betT/opuD)", r"\bbett\b|\bopud\b|\bbcct\b|\bcait\b|\bbetp\b|\bbetu\b|high-affinity choline uptake", None),
    ("tma", "choline permease LicB", r"choline permease licb", None),
    ("tma", "ProU/OpuA/Ota/Cho ABC", r"\bprov\b|\bprow\b|\bprox \(tc|\bopua[abc]\b|\bopub[abc]?\b"
                                     r"|\bopuc[abcd]?\b|\bota[abc]\b|\bcho[vwx]\b", None),
    ("tma", "betA/betB (choline -> betaine)", r"^choline dehydrogenase|betaine.aldehyde dehydrogenase|choline oxidase", None),
    ("tma", "betC (choline-sulfatase)", r"choline-sulfatase", None),
    ("tma", "grdA", r"reductase protein a @ selenocysteine|\bgrda\b", None),
    ("tma", "grdB/E (glycine reductase B)", r"glycine reductase component b", None),
    ("tma", "grdC/D (reductase C)", r"reductase component c", None),
    ("tma", "grdH/I (betaine reductase B)", r"betaine reductase component b", None),
    ("tma", "grdF/G (sarcosine reductase B)", r"sarcosine reductase component b", None),
    ("tma", "tmm", r"trimethylamine monooxygenase|\btmm\b", None),
    ("tma", "FMO, generic (Tmm family)", r"^flavin-containing monooxygenase$", None),
    ("tma", "torA/torZ (TMAO reductase)", r"trimethylamine-n-oxide reductase|trimethylamine n-oxide reductase|\btora\b|\btorz\b", None),
    ("tma", "TMA dehydrogenase", r"trimethylamine dehydrogenase|ec 1\.5\.8\.2\b", None),
    ("tma", "MttB superfamily, non-pyl", r"trimethylamine methyltransferase family protein|trimethylamine:corrinoid methyltransferase-like", None),
    ("tma", "corrinoid protein (unspecified)", r"^uncharacterized corrinoid protein$", None),
    ("tma", "betaine demethylation (gbcA/dgc)", r"\bgbc[ab]\b|glycine betaine demethylase|dimethylglycine (de)?(hydrogenase|demethylase|oxidase)|\bdgc[ab]\b", None),
    ("tma", "sarcosine oxidase/DH", r"sarcosine (oxidase|dehydrogenase)|monomeric sarcosine oxidase", None),
    ("tma", "BHMT (betaine-homocysteine MT)", r"betaine--homocysteine", None),
    ("tma", "carnitine 3-DH route (non-TMA)", r"carnitine 3-dehydrogenase|carnitinyl-coa", None),
    ("tma", "CaiB/BaiF-family CoA-transferase (generic)", r"caib/baif|l-carnitine dehydratase/bile acid", None),
    ("methylamine", "mttB (pyl, TMA MT)", r"\[trimethylamine--corrinoid protein\] co-methyltransferase|ec 2\.1\.1\.250"
                                          r"|^trimethylamine:corrinoid methyltransferase( / pyrrolysine-containing)?$", None),
    ("methylamine", "mtbB (pyl, DMA MT)", r"\[dimethylamine--corrinoid protein\] co-methyltransferase|ec 2\.1\.1\.249"
                                          r"|^dimethylamine:corrinoid methyltransferase( / pyrrolysine-containing)?$", None),
    ("methylamine", "mtmB (pyl, MMA MT)", r"\[methylamine--corrinoid protein\] co-methyltransferase|ec 2\.1\.1\.248", None),
    ("methylamine", "mttC/mtbC/mtmC (corrinoid)", r"(trimethylamine|dimethylamine|monomethylamine) methyltransferase corrinoid protein", None),
    ("methylamine", "mtbA (methylamine CoM MT)", r"methylamine-specific corrinoid protein\]:coenzyme m methyltransferase|ec 2\.1\.1\.247", None),
    ("methylamine", "pyl biosynthesis / pylS", r"pyrrolysine|pyrrolysyl", r"methyltransferase"),
    ("methylamine", "mcrA", r"methyl.coenzyme m reductase alpha", None),
]
RULES = [(p, g, re.compile(rx), re.compile(ex) if ex else None) for p, g, rx, ex in RULES]
CPL_CORE = ["phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"]      # C-P lyase complex + PhnM
POS_G, POS_X = "phnG (positional, unannotated)", "unannotated gene after phnM (positional)"
PHN_ORDER = (["phnC", "phnD", "phnE", "phnCDE-like ABC (cluster 12)", "AEP ABC transporter", "ptxB/C (phosphite ABC)",
              "phnF", "phnG", POS_G]
             + CPL_CORE[1:] + ["phnM-like", POS_X, "phnN", "phnO-like", "phnP", "phnA (YjdM-type)",
                               "phnA (phosphonoacetate hydrolase)", "phnW", "phnX", "phnY", "phnZ", "mpnS", "pepM",
                               "PEP-mutase family (CPEP mutase-like)", "ptxD (phosphite DH)"])
PI_ORDER = ["pstS", "pstC", "pstA", "pstB", "phoU", "phoB", "phoR", "Pi-regulatory homolog (unspecified)"]
TMA_ORDER = [g for p, g, _, _ in RULES if p == "tma"]
MA_ORDER = [g for p, g, _, _ in RULES if p == "methylamine"]
TMA_ANCHORS = {"cntA", "cntB", "cutC", "cutD", "cutC-cluster label (BMC gene)", "glycyl radical enzyme",
               "eutBC (ethanolamine ammonia-lyase)", "grdH/I (betaine reductase B)", "grdB/E (glycine reductase B)",
               "tmm", "torA/torZ (TMAO reductase)", "BCCT transporter (betT/opuD)"}
REF_FAMS = {"phnG": r"\bphng\b", "phnN": r"\bphnn\b|ec 2\.7\.4\.23", "phnP": r"\bphnp\b",
            "phnM": r"ec 3\.6\.1\.63", "phnJ": r"c-p lyase|ec 4\.7\.1\.1\b",
            "GrdH (betaine B, Sec)": r"betaine reductase component b beta",
            "GrdI (betaine B, alpha)": r"betaine reductase component b alpha",
            "GrdB (glycine B, Sec)": r"glycine reductase component b gamma",
            "GrdE (glycine B, alpha/beta)": r"glycine reductase component b beta"}
REF_FAMS = {k: re.compile(v) for k, v in REF_FAMS.items()}
# extra strings counted across the 310 MAGs (not gene-rule labels): partners for the betaine-demethylation route
X_FLAGS = {"x:gbcB": r"\bgbcb\b|glycine betaine demethylase subunit b",
           "x:dimethylglycine demethylation": r"dimethylglycine (de)?(hydrogenase|demethylase|oxidase)|\bdgc[ab]\b"}
X_FLAGS = {k: re.compile(v) for k, v in X_FLAGS.items()}
RECODE = {"TGA": ("Sec", "UGA", "selenocysteine"), "TAG": ("Pyl", "UAG", "pyrrolysine")}


def classify(fn):
    low = fn.lower()
    return [(p, g) for p, g, rx, ex in RULES if rx.search(low) and not (ex and ex.search(low))]


# ---- shared inputs ----------------------------------------------------------------------------------------
guilds = pd.read_csv(HERE / "guilds_r232.csv")
GI = guilds.set_index("mag")
rel = pd.read_csv(S / "data/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
dna = pd.read_csv(S / "data/SaltPondsDNA.csv").set_index("Sample")["DNA conc (µg/g soil)"]
PRE = sorted({"_".join(s.split("_")[:-1]) for s in dna.index})          # 12 cores
SITE = {p: ("restored" if p.startswith("SF2") else "reference" if p.startswith("R2A") else "unrestored") for p in PRE}
CORES = [p for s in ("unrestored", "restored", "reference") for p in PRE if SITE[p] == s]
AB = pd.DataFrame({p: (rel[f"{p}_D1"] * dna[f"{p}_D1"] + rel[f"{p}_D2"] * dna[f"{p}_D2"])
                   / (dna[f"{p}_D1"] + dna[f"{p}_D2"]) for p in PRE})   # DNA-weighted depth merge (raw basis)
checkm = pd.read_csv(S / "MAG_QC/checkm/summary/CheckM_summary_table.tsv", sep="\t")
checkm.index = checkm["Bin Name"].map(key)
gtdb = pd.concat([pd.read_csv(S / f"MAG_QC/gtdbtk_local_r232/gtdbtk.{d}.summary.tsv", sep="\t")
                  for d in ("ar53", "bac120")])
gtdb.index = gtdb["user_genome"].map(key)
fl_mepn = json.load(open(COR / "saltern_fluxes.json"))["0.07"]
fl_bet = json.load(open(COR / "saltern_fluxes_betaine.json"))
MEPN_MODEL = {key(k): v for k, v in fl_mepn["methane"].items()}       # 20 MAGs, 2020 MePn -> CH4 simulations
TMA_MODEL = {key(k): v for k, v in fl_bet["TMA"].items()}              # 15 MAGs, 2020 betaine -> TMA simulations
BET_UPTAKE = {key(k): v for k, v in fl_bet["betaine"].items()}
MEPN_DIR = S / "models/MePn_consumers"                                 # the curated 2020 MePn models
# C-P lyase chain in ModelSEED, matched by compound ID (the names are abbreviations such as "MePn [c0]")
CPL_CPDS = ("cpd25960", "cpd26518", "cpd26519")                        # MePn, RPnTP, PRPn
CPL_RXNS = {"rxn26456": "MePn → RPnTP", "rxn26457": "RPnTP → PRPn", "rxn26459": "PRPn → CH₄"}

arch_mepn = guilds[guilds.archaea & (guilds.mepn_consumer | guilds.mepn_2020)]
TARGETS = ([(r, "archaeal MePn") for _, r in arch_mepn.iterrows()]
           + [(r, "TMA candidate") for _, r in guilds[guilds.tma_producer].iterrows()]
           + [(r, "context: methanogen") for _, r in guilds[guilds.methanogen].iterrows()])


def load_genome(mag):
    d = json.load(open(RAST / f"{mag}.contigs__.RAST.json"))
    genes = []
    for f in d["features"]:
        c, st, sd, ln = f["location"][0]
        lo, hi = (st, st + ln - 1) if sd == "+" else (st - ln + 1, st)
        m = re.search(r"length_(\d+)", c)
        genes.append(dict(fid=f["id"].split(".")[-1], contig=c, contig_length=int(m[1]) if m else None,
                          start=lo, stop=hi, strand=sd,
                          function=" / ".join(f.get("functions") or [f.get("function", "")]).strip(),
                          protein=f.get("protein_translation", "").rstrip("*")))
    by = collections.defaultdict(list)
    for g in genes:
        by[g["contig"]].append(g)
    for c in by:
        by[c].sort(key=lambda g: g["start"])
        for i, g in enumerate(by[c]):
            g["idx"], g["n_on_contig"] = i, len(by[c])
            g["edge"] = g["start"] <= 30 or bool(g["contig_length"] and g["contig_length"] - g["stop"] <= 30)
    meta = {k: d.get(k) for k in ("domain", "genetic_code", "num_contigs", "dna_size", "gc_content")}
    return genes, by, meta


def contig_stats(mag):
    lens = sorted((len(r.seq) for r in SeqIO.parse(S / "MAG_contigs" / f"{mag}.contigs__.fna", "fasta")), reverse=True)
    run = 0
    for n50 in lens:
        run += n50
        if run >= sum(lens) / 2:
            break
    return dict(n_contigs=len(lens), size_bp=sum(lens), n50_bp=n50, longest_bp=lens[0])


def gap(a, b):
    """intergenic distance between two genes on one contig (negative = overlap)"""
    lo, hi = (a, b) if a["start"] <= b["start"] else (b, a)
    return hi["start"] - lo["stop"] - 1


# ---- alignment helpers --------------------------------------------------------------------------------------
AL = PairwiseAligner()
AL.mode = "local"
AL.substitution_matrix = substitution_matrices.load("BLOSUM62")
AL.open_gap_score, AL.extend_gap_score = -11, -1
AL_DNA = PairwiseAligner()
AL_DNA.mode, AL_DNA.match_score, AL_DNA.mismatch_score = "local", 1, -2
AL_DNA.open_gap_score, AL_DNA.extend_gap_score = -5, -2


def aln(a, b, al=AL):
    """local score, identity over aligned columns, aligned columns"""
    if not a or not b:
        return 0.0, 0.0, 0
    x = al.align(a, b)[0]
    ident = cols = 0
    for (a0, a1), (b0, b1) in zip(*x.aligned):
        cols += a1 - a0
        ident += sum(1 for i in range(a1 - a0) if a[a0 + i] == b[b0 + i])
    return x.score, ident / cols if cols else 0.0, cols


def null_scores(q, pool):
    s = sorted(AL.score(q, p) for p in pool)
    return s[int(0.99 * len(s))], s[-1]


def best_vs(q, refs, exclude_mag):
    best = (0.0, 0.0, 0, None)
    for m, p in refs:
        if m != exclude_mag and AL.score(q, p) > best[0]:
            sc, idn, cols = aln(q, p)
            best = (sc, idn, cols, m)
    return best


def contig_seqs(mag):
    return {r.id: str(r.seq) for r in SeqIO.parse(S / "MAG_contigs" / f"{mag}.contigs__.fna", "fasta")}


def readthrough(seq, g):
    """read an ORF on to the next in-frame stop: (its own stop codon, contig coordinate where the read-through ends)"""
    n = g["stop"] - g["start"] + 1
    s = seq[g["start"] - 1:] if g["strand"] == "+" else str(Seq(seq[:g["stop"]]).reverse_complement())
    for k in range(n, len(s) - 2, 3):
        if s[k:k + 3] in ("TAA", "TAG", "TGA"):
            ext = k + 3 - n
            return s[n - 3:n], (g["stop"] + ext if g["strand"] == "+" else g["start"] - ext)
    return s[n - 3:n], None


def recoded_splits(recs, by, hits):
    """genes RAST called as two ORFs because it read an in-frame Sec (UGA) or Pyl (UAG) codon as a stop: the labelled
    ORF says 'selenocysteine/pyrrolysine-containing', ends at that codon, and its in-frame read-through ends exactly at
    the stop of the next ORF on the same strand"""
    lab = {(h["contig"], h["idx"]): h["gene"] for h in hits}
    out = []
    for h in hits:
        cg = by[h["contig"]]
        j = h["idx"] + (1 if h["strand"] == "+" else -1)
        if not 0 <= j < len(cg) or cg[j]["strand"] != h["strand"]:
            continue
        stop, end = readthrough(recs[h["contig"]], h)
        if stop not in RECODE or RECODE[stop][2] not in h["function"].lower():
            continue
        c = cg[j]
        if end != (c["stop"] if h["strand"] == "+" else c["start"]):
            continue
        full = ((end - h["start"]) if h["strand"] == "+" else (h["stop"] - end)) + 1
        out.append(dict(n=h, c=c, c_label=lab.get((c["contig"], c["idx"]), ""), aa=RECODE[stop][0],
                        codon=RECODE[stop][1], full_aa=full // 3 - 1))
    return out


PYL_MT = ("mttB (pyl, TMA MT)", "mtbB (pyl, DMA MT)", "mtmB (pyl, MMA MT)")


def pyl_evidence(mag):
    """does a pyl-type methylamine methyltransferase label come with Pyl evidence? (a UAG-split ORF whose C-terminal
    ORF is unlabelled or carries the same label, and/or Pyl biosynthesis genes)"""
    genes, by, _ = load_genome(mag)
    hits = [dict(g, gene=lab) for g in genes for p, lab in classify(g["function"]) if p == "methylamine"]
    mt = [h for h in hits if h["gene"] in PYL_MT]
    sp = [s for s in recoded_splits(contig_seqs(mag), by, mt) if s["aa"] == "Pyl" and
          (s["c_label"] == s["n"]["gene"] or not s["c"]["function"])]
    return dict(n_orfs=len(mt), longest=max((len(h["protein"]) for h in mt), default=0), n_split=len(sp),
                full_aa={s["n"]["gene"]: s["full_aa"] for s in sp},
                n_machinery=sum(1 for h in hits if h["gene"] == "pyl biosynthesis / pylS"))


# ---- 1. one pass over all 310 RAST genomes: reference families, null pool, archaeal phn context ---------------
def scan_all():
    rng = random.Random(SEED)
    fam = collections.defaultdict(list)       # family -> [(MAG, protein)]
    pool, arch_ctx, glob = [], [], collections.defaultdict(set)   # glob: label -> MAGs (iterativeID) with it
    for mag in guilds["mag"]:
        genes, _, _ = load_genome(mag)
        row = GI.loc[mag]
        cnt = collections.Counter()
        for g in genes:
            low = g["function"].lower()
            for f, rx in REF_FAMS.items():
                if rx.search(low) and g["protein"]:
                    fam[f].append((mag, g["protein"]))
            for _, lab in classify(g["function"]):
                cnt[lab] += 1
                if lab == "torA/torZ (TMAO reductase)":     # DMSO-reductase genes with a dual TMAO-reductase label
                    cnt["x:TMAO reductase, " + ("DMSO-reductase dual label" if "dimethyl sulfoxide" in low else "TorA/TorZ-type")] += 1
            for lab, rx in X_FLAGS.items():
                if rx.search(low):
                    cnt[lab] += 1
            if len(g["protein"]) > 80 and rng.random() < 0.002:    # ~1e6 proteins -> ~2000, N_NULL kept
                pool.append((mag, g["protein"]))
        for lab in cnt:
            glob[lab].add(row.iterativeID_r232)
        if row.archaea:
            arch_ctx.append(dict(iterativeID=row.iterativeID_r232, mepn_consumer=row.mepn_consumer,
                                 mepn_2020=row.mepn_2020,
                                 **{g: cnt.get(g, 0) for g in ["phnC", "phnD", "phnE"] + CPL_CORE
                                    + ["phnM-like", "phnN", "phnP", "phnO-like"]}))
    rng.shuffle(pool)
    return fam, pool[:N_NULL], pd.DataFrame(arch_ctx), glob


# ---- 2. per-MAG inventory -----------------------------------------------------------------------------------
def cpl_contigs(hits):
    return sorted({h["contig"] for h in hits if h["gene"] in ("phnH", "phnI", "phnJ", "phnK", "phnL")})


def inventory(row, role, fam, pool):
    mag = row.mag
    genes, by, meta = load_genome(mag)
    panels = {"archaeal MePn": ("phn", "pi"), "TMA candidate": ("tma",),
              "context: methanogen": ("methylamine", "tma")}[role]
    hits = [dict(g, panel=p, gene=lab, evidence="RAST annotation")
            for g in genes for p, lab in classify(g["function"]) if p in panels]
    pos_checks = []
    if role == "archaeal MePn":           # positional candidates inside the C-P lyase locus only
        locus_contigs = cpl_contigs(hits)
        labelled = {(h["contig"], h["idx"]) for h in hits}
        pool_other = [p for m, p in pool if m != mag]
        for h in [h for h in hits if h["contig"] in locus_contigs and h["gene"] in ("phnH", "phnM", "phnM-like")]:
            cg = by[h["contig"]]
            up = h["gene"] == "phnH"         # phnG sits upstream of phnH; look downstream of phnM
            j = h["idx"] + (-1 if (h["strand"] == "+") == up else 1)
            if not 0 <= j < len(cg):
                continue
            c = cg[j]
            if c["strand"] != h["strand"] or (c["contig"], c["idx"]) in labelled:
                continue
            fams = ["phnG", "phnN", "phnP", "phnM"]
            n99, nmax = null_scores(c["protein"], pool_other)
            best = {f: best_vs(c["protein"], fam[f], mag) for f in fams}
            pc = dict(iterativeID=row.iterativeID_r232, position="upstream of phnH" if up else "downstream of phnM",
                      gene=c, gap_bp=gap(c, h), null99=n99, nullmax=nmax, best=best,
                      n_ref={f: sum(1 for m, _ in fam[f] if m != mag) for f in fams})
            pos_checks.append(pc)
            if up and best["phnG"][0] > nmax:
                hits.append(dict(c, panel="phn", gene=POS_G,
                                 evidence=f"inferred: gene upstream of phnH, best local score {best['phnG'][0]:.0f} vs "
                                          f"RAST PhnG of other MAGs (null max {nmax:.0f})"))
            elif not up:
                hits.append(dict(c, panel="phn", gene=POS_X,
                                 evidence=f"position only: best score vs RAST PhnN/PhnP/PhnG/PhnM "
                                          f"{max(b[0] for b in best.values()):.0f} (null max {nmax:.0f})"))
    return genes, by, meta, hits, pos_checks


def abundance(mag):
    a = AB.loc[mag]
    out = {p: a[p] for p in PRE}
    for s in ("unrestored", "restored", "reference"):
        out[f"mean_{s}"] = a[[p for p in PRE if SITE[p] == s]].mean()
    out["n_detected"] = int((a > ZERO).sum())
    out["max_core"], out["max"] = a.idxmax(), a.max()
    return out


def model_info(row):
    k = key(row.mag)
    out = {"mepn_set": k in MEPN_MODEL, "tma_set": k in TMA_MODEL, "mepn_ch4": MEPN_MODEL.get(k),
           "tma": TMA_MODEL.get(k), "bet": BET_UPTAKE.get(k), "tma_rxns": [], "base_cpl": None,
           "base_mepn_uptake": None, "curated_mepn_model": k in {key(p.name) for p in MEPN_DIR.glob("*.json")}}
    f = COR / "models" / f"{row.mag}.contigs__.RAST_gapfilled_betaine.json"
    if k in TMA_MODEL and f.exists():
        for r in json.load(open(f))["reactions"]:
            if r["metabolites"].get("cpd00441_c0", 0) > 0 and not r["id"].startswith(("EX_", "rxn_TMA_transport")):
                out["tma_rxns"].append((r["id"], r.get("name", ""), r.get("gene_reaction_rule", "") or "none"))
    base = COR / "models" / f"{row.mag}.contigs__.RAST.json"
    if base.exists():
        m = json.load(open(base))
        out["base_cpl"] = {r["id"].rsplit("_", 1)[0]: (r.get("gene_reaction_rule") or "") for r in m["reactions"]
                           if any(x.split("_")[0] in CPL_CPDS for x in r["metabolites"])}
        out["base_mepn_uptake"] = any(x["id"].startswith("cpd25960_e") for x in m["metabolites"])
    return out


def curated_mepn_summary():
    """what the curated 2020 MePn models (models/MePn_consumers) carry for MePn uptake and the C-P lyase step"""
    n, have = 0, collections.Counter()      # (step, reaction ID, name, has a gene rule) -> number of models
    for p in MEPN_DIR.glob("*.json"):
        n += 1
        for r in json.load(open(p))["reactions"]:
            mets, rid = r["metabolites"], r["id"].rsplit("_", 1)[0]
            step = ("uptake" if "cpd25960_e0" in mets and "cpd25960_c0" in mets else
                    "C-P lyase" if mets.get("cpd26519_c0", 0) < 0 else None)
            if step:
                have[step, rid, r.get("name", "").split(" [")[0], bool(r.get("gene_reaction_rule"))] += 1
    return n, have


# ---- 3. locus-level checks -----------------------------------------------------------------------------------
def cpl_locus(by, hits):
    """genes on the C-P lyase contigs, from the first to the last phn hit, +-1 gene"""
    locus = []
    for c in cpl_contigs(hits):
        idx = [h["idx"] for h in hits if h["panel"] == "phn" and h["contig"] == c]
        locus += by[c][max(0, min(idx) - 1): min(len(by[c]), max(idx) + 2)]
    return locus


def split_fragments(mag, hits):
    """same gene label on two contigs, both copies at a contig edge -> protein overlap + contig-end DNA alignment"""
    out, edge = [], collections.defaultdict(list)
    for h in hits:
        if h["panel"] == "phn" and h["edge"]:
            edge[h["gene"]].append(h)
    recs = None
    for lab, hs in edge.items():
        if len({h["contig"] for h in hs}) < 2:
            continue
        a, b = hs[0], hs[1]
        recs = recs or contig_seqs(mag)
        A, B, best = recs[a["contig"]], recs[b["contig"]], None
        # an end of one contig against the start of the other; B as given or reverse-complemented
        for rc in (False, True):
            Bo = str(Seq(B).reverse_complement()) if rc else B
            for left, right, a_end, b_end in ((A, Bo, "end", "start"), (Bo, A, "start", "end")):
                x = AL_DNA.align(left[-300:], right[:300])[0]
                if best is None or x.score > best["x"].score:
                    best = dict(x=x, left=left, right=right, a_left=left is A, rc=rc, a_end=a_end,
                                b_end=b_end if not rc else {"start": "end", "end": "start"}[b_end],
                                orient="reverse complement" if rc else "same strand")
        sc, idn, cols = aln(best["left"][-300:], best["right"][:300], AL_DNA)
        best.update(score=sc, ident=idn, cols=cols)
        out.append(dict(gene=lab, a=a, b=b, dna=best, **stitch_translate(a, b, best)))
    return out


def stitch_translate(a, b, d):
    """join the two contigs at the ends of their (ungapped) alignment and translate the ORF that spans both fragments.
    The overlap is taken once from each contig, to see which nucleotide mismatches change the protein."""
    x, left, right = d["x"], d["left"], d["right"]
    bl, br = x.aligned
    if len(bl) != 1:                                   # gapped overlap: no codon-by-codon comparison
        return dict(joined=None)
    off = len(left) - min(300, len(left))
    lx0, lx1, ry0, ry1 = off + bl[0][0], off + bl[0][1], br[0][0], br[0][1]
    nB = len(right) if d["a_left"] else len(left)

    def place(g, is_b):                                # fragment coordinates (1-based) in the joined sequence
        lo, hi, sd = g["start"], g["stop"], g["strand"]
        if is_b and d["rc"]:
            lo, hi, sd = nB - hi + 1, nB - lo + 1, {"+": "-", "-": "+"}[sd]
        if is_b == d["a_left"]:                        # this fragment's contig is the right-hand one
            lo, hi = lo + lx1 - ry1, hi + lx1 - ry1
        return lo, hi, sd

    (la, ha, sa), (lb, hb, sb) = place(a, False), place(b, True)
    lo, hi = min(la, lb), max(ha, hb)
    if sa != sb or (hi - lo + 1) % 3:
        return dict(joined=None)
    prot = {}
    for src, joined in (("left", left[:lx1] + right[ry1:]), ("right", left[:lx0] + right[ry0:])):
        nt = joined[lo - 1:hi]
        prot[src] = str(Seq(str(Seq(nt).reverse_complement()) if sa == "-" else nt).translate())
    if "*" in prot["left"][:-1] or "*" in prot["right"][:-1]:
        return dict(joined=None)
    joined_aa = len(prot["left"].rstrip("*"))
    codon = (lambda j: (j - (lo - 1)) // 3) if sa == "+" else (lambda j: (hi - 1 - j) // 3)
    mm = [codon(lx0 + i) for i in range(lx1 - lx0) if left[lx0 + i] != right[ry0 + i]]
    nonsyn = sorted({c for c in mm if prot["left"][c] != prot["right"][c]})
    names = (a["contig"], b["contig"]) if d["a_left"] else (b["contig"], a["contig"])
    return dict(joined=joined_aa, overlap_aa=len(a["protein"]) + len(b["protein"]) - joined_aa, n_mm=len(mm),
                n_mm_codons=len(set(mm)), n_mm_nonsyn=sum(1 for c in mm if c in nonsyn),
                aa_diffs=[(c + 1, prot["left"][c], prot["right"][c]) for c in nonsyn],
                left_contig=names[0], right_contig=names[1])


# ---- 4. main ------------------------------------------------------------------------------------------------
def main():
    fam, pool, arch_ctx, glob = scan_all()
    pool_p = [p for _, p in pool]
    R = {}
    for row, role in TARGETS:
        genes, by, meta, hits, pos = inventory(row, role, fam, pool)
        k = key(row.mag)
        g = gtdb.loc[k]
        R[row.iterativeID_r232] = dict(
            row=row, role=role, meta=meta, hits=hits, pos=pos, by=by, n_genes=len(genes),
            ab=abundance(row.mag), model=model_info(row), contigs=contig_stats(row.mag),
            gtdb=g.to_dict(), checkm=dict(completeness=checkm.loc[k, "Completeness"],
                                         contamination=checkm.loc[k, "Contamination"],
                                         lineage=checkm.loc[k, "Marker Lineage"]))
        v = R[row.iterativeID_r232]
        v["recoded"] = recoded_splits(contig_seqs(row.mag), by, hits)     # Sec/Pyl genes RAST split into two ORFs
        for s in v["recoded"]:
            how = f"one gene split by RAST at its {s['aa']} ({s['codon']}) codon"
            s["n"]["evidence"] += (f"; N-terminal ORF of {how}, continued by feature {s['c']['fid']}"
                                   + ("" if s["c_label"] else " (no rule label)"))
            for h in hits:
                if (h["contig"], h["idx"]) == (s["c"]["contig"], s["c"]["idx"]) and h["gene"] == s["n"]["gene"]:
                    h["evidence"] += f"; C-terminal ORF of {how}, after feature {s['n']['fid']}"
        if role == "archaeal MePn":
            v["cpl_locus"] = cpl_locus(by, hits)
            v["splits"] = split_fragments(row.mag, hits)
            v["phnJ_vs_bact"] = [(h, best_vs(h["protein"], [(m, p) for m, p in fam["phnJ"] if not GI.loc[m].archaea], row.mag))
                                 for h in hits if h["gene"] == "phnJ"]
        if role == "TMA candidate":
            v["grdB"] = []
            for h in hits:
                if h["gene"] in ("grdB/E (glycine reductase B)", "grdH/I (betaine reductase B)") and len(h["protein"]) >= 150:
                    v["grdB"].append((h, {f: best_vs(h["protein"], fam[f], row.mag)[0] for f in
                                          ("GrdH (betaine B, Sec)", "GrdI (betaine B, alpha)",
                                           "GrdB (glycine B, Sec)", "GrdE (glycine B, alpha/beta)")}))
    # gene-by-gene comparison of the archaeal C-P lyase loci
    loci = {i: v for i, v in R.items() if v["role"] == "archaeal MePn" and v["cpl_locus"]}
    synteny, ids = [], list(loci)
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            A, B = loci[ids[i]], loci[ids[j]]
            labA = {(h["contig"], h["idx"]): h["gene"] for h in A["hits"] if h["panel"] == "phn"}
            labB = {(h["contig"], h["idx"]): h["gene"] for h in B["hits"] if h["panel"] == "phn"}
            for gb in B["cpl_locus"]:
                if len(gb["protein"]) < 30:
                    continue
                (sc, idn, cols), ga = max(((aln(gb["protein"], ga["protein"]), ga) for ga in A["cpl_locus"]
                                           if len(ga["protein"]) >= 30), key=lambda t: t[0][0])
                n99, nmax = null_scores(gb["protein"], pool_p)
                synteny.append(dict(query=ids[j], subject=ids[i], q=gb, q_label=labB.get((gb["contig"], gb["idx"]), ""),
                                    s=ga, s_label=labA.get((ga["contig"], ga["idx"]), ""), score=sc, identity=idn,
                                    cols=cols, null99=n99, nullmax=nmax))
    # lone phnM-family genes in archaeal MePn MAGs without a C-P lyase locus, vs the operon PhnMs
    lone_m = []
    op_m = [(i, h) for i, v in loci.items() for h in v["hits"] if h["gene"] in ("phnM", "phnM-like")]
    for i, v in R.items():
        if v["role"] == "archaeal MePn" and not v["cpl_locus"]:
            for h in v["hits"]:
                if h["gene"] in ("phnM", "phnM-like"):
                    n99, nmax = null_scores(h["protein"], pool_p)
                    for oi, oh in op_m:
                        sc, idn, cols = aln(h["protein"], oh["protein"])
                        lone_m.append(dict(query=i, q=h, subject=oi, s=oh, score=sc, identity=idn, cols=cols,
                                           null99=n99, nullmax=nmax))
    df = write_csv(R)
    write_md(R, synteny, lone_m, arch_ctx, fam, df, glob)


def write_csv(R):
    rows = []
    for iid, v in R.items():
        r = v["row"]
        loc = {(g["contig"], g["idx"]) for g in v.get("cpl_locus", [])}
        for h in v["hits"]:
            rows.append(dict(mag=r.mag, iterativeID_r232=iid,
                             iterativeID_2020=r.iterativeID_2020 if isinstance(r.iterativeID_2020, str) else "",
                             role=v["role"], panel=h["panel"], gene=h["gene"], evidence=h["evidence"],
                             function=h["function"] or "(no function)", feature_id=f"{r.mag}.contigs__.RAST.{h['fid']}",
                             contig=h["contig"], contig_length=h["contig_length"], start=h["start"], stop=h["stop"],
                             strand=h["strand"], length_aa=len(h["protein"]), gene_on_contig=f"{h['idx'] + 1}/{h['n_on_contig']}",
                             at_contig_edge=h["edge"], in_CP_lyase_locus=(h["contig"], h["idx"]) in loc))
    df = pd.DataFrame(rows)
    order = {g: i for i, g in enumerate(PHN_ORDER + PI_ORDER + TMA_ORDER + MA_ORDER)}
    role_o = {"archaeal MePn": 0, "TMA candidate": 1, "context: methanogen": 2}
    df = (df.assign(_r=df.role.map(role_o), _o=df.gene.map(order).fillna(999))
            .sort_values(["_r", "iterativeID_r232", "_o", "contig", "start"]).drop(columns=["_r", "_o"]))
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV.name}: {len(df)} gene rows")
    return df


# ---- 5. markdown ---------------------------------------------------------------------------------------------
def node(c):
    m = re.match(r"(NODE_\d+)_length_(\d+)", c)
    return f"{m[1]} ({int(m[2]) / 1000:.1f} kb)" if m else c


def pct(x, d=3):
    """% of reads; values that would round to zero keep one significant figure"""
    if x is None or pd.isna(x):
        return "–"
    v = 100 * x
    if v == 0:
        return "0"
    return f"{v:.{d}f}" if v >= 10 ** -d else f"{v:.{int(-math.floor(math.log10(v)))}f}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def esc(s):
    return str(s).replace("|", "\\|")


def overview_rows(R, ids):
    rows = []
    for i in ids:
        v, r, g = R[i], R[i]["row"], R[i]["gtdb"]
        sp = r.species if isinstance(r.species, str) else f"{r.genus if isinstance(r.genus, str) else r.family} sp. (no species)"
        if isinstance(g.get("closest_genome_reference"), str):
            ref = f"{g['closest_genome_reference']}, ANI {g['closest_genome_ani']:.2f} %, AF {g['closest_genome_af']:.3f}"
        elif isinstance(g.get("warnings"), str) and "Closest reference" in g["warnings"]:
            ref = "no species; " + re.search(r"Closest reference [^;]+", g["warnings"])[0]
        elif isinstance(g.get("closest_placement_reference"), str):
            ref = (f"none within species radius; placement ref {g['closest_placement_reference']}, ANI "
                   f"{g['closest_placement_ani']:.2f} %, AF {g['closest_placement_af']:.3f}")
        else:
            ref = "none (placed by topology)"
        c, ck = v["contigs"], v["checkm"]
        rows.append([f"**{i}**", r.iterativeID_2020 if isinstance(r.iterativeID_2020, str) else "none", f"{r.phylum}; {r.family}",
                     f"*{sp}*", ref, f"{ck['completeness']:.2f} / {ck['contamination']:.2f}",
                     f"{c['n_contigs']} / {c['size_bp'] / 1e6:.2f} Mb / {c['n50_bp'] / 1000:.1f} kb"])
    return md_table(["MAG (R232)", "2020 ID", "Phylum; family", "R232 species", "Closest reference (GTDB-Tk)",
                     "CheckM compl. / contam. (%)", "Contigs / size / N50"], rows)


def abundance_table(R, ids):
    rows = []
    for i in ids:
        a = R[i]["ab"]
        rows.append([f"**{i}**"] + [pct(a[p]) for p in CORES] + [pct(a["mean_unrestored"]), pct(a["mean_restored"]),
                    pct(a["mean_reference"]), f"{a['n_detected']}/12"])
    return md_table(["MAG"] + CORES + ["mean unrestored", "mean restored", "mean reference", "detected"], rows)


def flags_table(R, ids):
    rows = []
    for i in ids:
        r, m = R[i]["row"], R[i]["model"]
        fl = lambda x: "yes" if x else "no"
        tma = "–"
        if m["tma_set"]:
            vals = sorted({round(x, 4) for x in m["tma"].values()})
            bet = sorted({round(x, 4) for x in m["bet"].values()})
            tma = f"TMA {', '.join(map(str, vals))} in all {len(m['tma'])} cores (betaine {', '.join(map(str, bet))})"
        rows.append([f"**{i}**", fl(r.mepn_consumer), fl(r.mepn_2020), fl(r.betaine_reducer), fl(r.betaine_2020),
                     fl(r.tma_producer), fl(m["mepn_set"]), fl(m["tma_set"]), tma])
    return md_table(["MAG", "MePn (phnJ)", "MePn 2020 set", "betaine (grdH)", "betaine 2020 set", "TMA candidate",
                     "in 2020 MePn model set (20)", "in 2020 betaine model set (15)", "2020 modelled flux"], rows)


def inventory_matrix(R, ids, order, df):
    sub = df[df.iterativeID_r232.isin(ids)]
    rows = []
    for g in order:
        cells = []
        for i in ids:
            s = sub[(sub.iterativeID_r232 == i) & (sub.gene == g)]
            nodes = ', '.join(sorted({node(c).split(' ')[0] for c in s.contig}))
            ct = s.evidence.str.extract(r"C-terminal ORF of one gene split by RAST at its (\w+) \((\w+)\) codon").dropna()
            if s.empty:
                cells.append("–")
            elif len(ct):         # count genes, not the two ORFs RAST made of a Sec/Pyl gene
                cells.append(f"{len(s) - len(ct)} ({len(s)} ORFs, split at {ct.iloc[0, 1]}/{ct.iloc[0, 0]}; {nodes})")
            else:
                cells.append(f"{len(s)} ({nodes})")
        if any(c != "–" for c in cells):
            rows.append([g] + cells)
    absent = [g for g in order if sub[sub.gene == g].empty]
    return md_table(["gene (rule label)"] + ids, rows), absent


def neighbourhood(v, contig, lo, hi, labels):
    cg = v["by"][contig]
    rows, prev = [], None
    for g in cg[max(0, lo): min(len(cg), hi + 1)]:
        lab = labels.get((g["contig"], g["idx"]), "")
        rows.append([g["idx"] + 1, f"{g['start']}–{g['stop']}", g["strand"], len(g["protein"]),
                     "–" if prev is None else gap(prev, g), f"**{lab}**" if lab else "",
                     esc(g["function"] or "(no function)")])
        prev = g
    head = f"`{contig}` ({cg[0]['n_on_contig']} genes, {cg[0]['contig_length']:,} bp)"
    return head + "\n\n" + md_table(["#", "start–stop", "strand", "aa", "gap (bp)", "rule label", "RAST function"], rows)


def write_md(R, synteny, lone_m, arch_ctx, fam, df, glob):
    A = [i for i, v in R.items() if v["role"] == "archaeal MePn"]
    T = [i for i, v in R.items() if v["role"] == "TMA candidate"]
    M = [i for i, v in R.items() if v["role"] == "context: methanogen"]
    P = lambda i, g: next((p for p in R[i]["pos"] if p["position"] == g), None)
    L = []
    w = L.append
    nat, hq, hr = (next(i for i in A if i.startswith(s)) for s in ("Natrinema", "Haloquadratum", "Halorubrum"))
    gN, gH = P(nat, "upstream of phnH"), P(hq, "upstream of phnH")
    xN, xH = P(nat, "downstream of phnM"), P(hq, "downstream of phnM")
    sG = next(s for s in synteny if s["q_label"] == POS_G)
    sJ = next(s for s in synteny if s["q_label"] == "phnJ")
    sM = next(s for s in synteny if s["q_label"] == "phnM")
    sX = next(s for s in synteny if s["q_label"] == POS_X)
    sO = next(s for s in synteny if s["q_label"] == "phnO-like")
    core_ids = [s for s in synteny if s["q_label"] in ("phnC", "phnD", "phnE", POS_G, "phnH", "phnI", "phnJ", "phnK",
                                                         "phnL", "phnM", POS_X, "phnO-like")]
    idr = (min(s["identity"] for s in core_ids), max(s["identity"] for s in core_ids))
    split = R[nat]["splits"][0] if R[nat]["splits"] else None
    hqI = next(h for h in R[hq]["hits"] if h["gene"] == "phnI")
    jb = {i: R[i]["phnJ_vs_bact"][0][1] for i in (nat, hq)}
    n_arch = len(arch_ctx)
    n_cpl = int(((arch_ctx[["phnH", "phnI", "phnJ", "phnK", "phnL"]].sum(axis=1)) > 0).sum())
    n_phnP = int((arch_ctx.phnP > 0).sum())
    n_phnCD = int(((arch_ctx.phnC + arch_ctx.phnD + arch_ctx.phnE) > 0).sum())
    n_m = int(((arch_ctx.phnM + arch_ctx["phnM-like"]) > 0).sum())
    hrM = [x for x in lone_m if x["query"] == hr]
    ab = {i: R[i]["ab"] for i in R}
    mod = {i: R[i]["model"] for i in R}
    caz = next(i for i in T if i.startswith("CAZFJG01"))
    j153 = next(i for i in T if i.startswith("J153"))
    uba = next(i for i in T if i.startswith("UBA11579"))
    t1 = next(i for i in T if i.startswith("T1Sed10"))
    hitsof = lambda i, g: [h for h in R[i]["hits"] if h["gene"] == g]
    cnt = lambda i, g: len(hitsof(i, g))
    grd = R[caz]["grdB"]
    grd_gly_wins = sum(1 for _, s in grd if max(s["GrdB (glycine B, Sec)"], s["GrdE (glycine B, alpha/beta)"])
                       > max(s["GrdH (betaine B, Sec)"], s["GrdI (betaine B, alpha)"]))
    meth = M[0] if M else None
    in310 = lambda lab: len(glob.get(lab, ()))
    nT = lambda lab: sum(1 for i in T if cnt(i, lab))
    grd_contigs = len({h["contig"] for h in hitsof(caz, "grdB/E (glycine reductase B)")})
    id2020 = R[caz]["row"].iterativeID_2020
    # genes, not ORFs: drop the C-terminal ORF of a Sec/Pyl gene that RAST split and that carries the same label
    n_split_c = lambda i, g: sum(1 for s in R[i]["recoded"] if s["n"]["gene"] == g and s["c_label"] == g)
    genes_of = lambda i, g: cnt(i, g) - n_split_c(i, g)
    edge_frag = lambda i, g: [h for h in hitsof(i, g) if h["edge"] and len(h["protein"]) < 150]
    # TMAO reductase labels across the 310 MAGs
    tor_t = sorted(glob.get("x:TMAO reductase, TorA/TorZ-type", ()))
    tor_d = sorted(set(glob.get("x:TMAO reductase, DMSO-reductase dual label", ())) - set(tor_t))
    gbcB = sorted(glob.get("x:gbcB", ()))
    # the J153.1 text below states that it has no GbcB and no dimethylglycine-demethylating enzyme
    assert j153 not in gbcB and j153 not in glob.get("x:dimethylglycine demethylation", ())
    # ... and that none of the four candidates carries a betaine or TMAO reductase
    assert not any(cnt(i, "torA/torZ (TMAO reductase)") or cnt(i, "grdH/I (betaine reductase B)") for i in T)
    # Pyl evidence behind each pyl-type methylamine methyltransferase label
    gi_ = guilds.set_index("iterativeID_r232")
    pylev = {o: pyl_evidence(gi_.loc[o, "mag"]) for o in sorted(set().union(*(glob.get(l, set()) for l in PYL_MT)))}
    pyl_ok = lambda o: pylev[o]["n_split"] > 0 or pylev[o]["n_machinery"] > 0
    lab_of = lambda i: {h["fid"]: h["gene"] for h in R[i]["hits"]}
    andjoin = lambda xs: ", ".join(xs[:-1]) + " and " + xs[-1] if len(xs) > 1 else "".join(xs)

    def gpr_text(i, rule):
        """a model gene rule with each CDS replaced by its rule label"""
        lab = lab_of(i)
        return re.sub(r"[^\s()]+?\.CDS\.(\d+)", lambda m: f"{lab.get(m[1], 'unlabelled')} CDS.{m[1]}", rule) if rule else "no gene rule"

    def cpl_text(i):
        rx = mod[i]["base_cpl"] or {}
        if not rx:
            return "no reaction with MePn, RPnTP or PRPn"
        return "; ".join(f"{r} ({CPL_RXNS.get(r, r)}; {gpr_text(i, rx[r])})" for r in sorted(rx))

    w("# Archaeal MePn consumers and TMA-route candidates: what the annotations support\n")
    w(f"Generated by `make_archaeal_mepn_tma_report.py` on {datetime.date.today().isoformat()}. "
      "Gene coordinates and RAST function strings for every hit are in `archaeal_mepn_and_tma_candidates_genes.csv`.\n")
    w("**Scope.** The archaea in either MePn set:")
    w(f"- the two *phnJ* carriers, {nat} and {hq};")
    w(f"- {hr}, which is in the hand-curated 2020 set only.\n")
    w(f"No other archaeal MAG is in either set. The report also covers the four TMA-route candidates ({', '.join(T)}). "
      f"{meth} is added as context, because it is the only methanogen that could turn TMA into CH₄.\n")
    w("Everything below comes from the RAST annotation (function strings and gene locations) and from sequence "
      "comparisons on the RAST protein translations. No profile-HMM search, phylogeny or experiment was run. "
      "Statements marked *inference* go beyond what the annotation states.\n")

    # ---------------- bottom line
    w("## Bottom line\n")
    w(f"1. **{nat} and {hq} both carry a C-P lyase locus that, at the annotation level, is complete enough to release CH₄ from MPn.**")
    w(f"   - Both loci hold the PhnCDE transporter and PhnH, I, J, K, L and M. The arrangement is the same in both: *phnDCE* is "
      f"divergent from *phnGHIJKLM*, followed by one unannotated gene.")
    w(f"   - RAST did not call PhnG in either MAG. The gene at the PhnG position is still a PhnG homologue (*inference*). "
      f"It scores {gN['best']['phnG'][0]:.0f} and {gH['best']['phnG'][0]:.0f} against RAST-annotated PhnG from other MAGs; "
      f"random proteins score at most {max(gN['nullmax'], gH['nullmax']):.0f}. The two copies are "
      f"{sG['identity']:.0%} identical to each other.")
    w(f"   - Across the two loci, matching genes are {idr[0]:.0%}–{idr[1]:.0%} identical; PhnJ is {sJ['identity']:.0%}.")
    w(f"   - In {hq} the whole locus is on one contig. In {nat} it is split at *phnI* across two contigs, whose ends overlap "
      f"by {split['dna']['cols']} nt at {split['dna']['ident']:.0%} identity.")
    w(f"   - What is missing: *phnN* in both (*phnP* is annotated elsewhere in both genomes). PhnN only recycles the ribose "
      f"product and is not needed to release CH₄.")
    w(f"2. **{hr} has no C-P lyase.**")
    w(f"   - Its 2020 membership rests on one lone *phnM* gene among unrelated genes. That gene is a real PhnM homologue, "
      f"{min(x['identity'] for x in hrM):.0%}–{max(x['identity'] for x in hrM):.0%} identical to the two operon PhnMs. "
      f"It also has generic PhnC-type and \"cluster 12\" ABC transporter genes and a PhnP-labelled metallo-β-lactamase.")
    w(f"   - PhnM alone cannot cleave the C–P bond. The *phnJ* rule is therefore right to exclude {hr}.")
    n_mod = sum(1 for i in A if mod[i]["mepn_set"])
    n_cur = sum(1 for i in A if mod[i]["curated_mepn_model"])
    w(f"3. **{'None' if not n_mod else n_mod} of the {len(A)} archaea {'was' if not n_mod else 'were'} modelled in 2020.**")
    w(f"   - {'None' if not n_mod else n_mod} is among the {len(MEPN_MODEL)} MAGs with MePn → CH₄ simulations "
      f"(`saltern_fluxes.json`), and {'none' if not n_cur else n_cur} has a model in `models/MePn_consumers/`.")
    w(f"   - Their base ModelSEED models (`correlational/models/*.RAST.json`) do contain C-P lyase reactions, found by the "
      f"compound IDs cpd25960 (MePn), cpd26518 (RPnTP) and cpd26519 (PRPn):")
    for i in (hq, nat, hr):
        rx = mod[i]["base_cpl"] or {}
        miss = [r for r in CPL_RXNS if r not in rx]
        txt = f"     - {i}: {cpl_text(i)}."
        if not miss and all(rx.values()):
            txt += " That is the whole intracellular MePn → CH₄ chain, each step with a gene rule."
        elif miss == ["rxn26457"] and hitsof(i, "phnM-like") and not hitsof(i, "phnM"):
            txt += (f" rxn26457 ({CPL_RXNS['rxn26457']}) is missing: RAST labels its PhnM \"{hitsof(i, 'phnM-like')[0]['function']}\" "
                    f"(CDS.{hitsof(i, 'phnM-like')[0]['fid']}), without EC 3.6.1.63.")
        elif len(rx) == 1:
            txt += f" Nothing else of the chain."
        w(txt)
    if not any(mod[i]["base_mepn_uptake"] for i in A):
        w(f"     - None of the {len(A)} has an MePn uptake reaction (no extracellular MePn, cpd25960_e0).")
    w(f"   - These base-model reactions may be how the hand-curated 2020 set came to include all {len(A)} archaea, {hr} "
      f"through {', '.join(sorted(mod[hr]['base_cpl'] or {}))} alone (*inference*, untested).")
    w(f"   - Both *phnJ* archaea live in the unrestored cores:")
    w(f"     - {nat}: mean {pct(ab[nat]['mean_unrestored'])} % of reads there, against {pct(ab[nat]['mean_restored'])} % "
      f"(restored) and {pct(ab[nat]['mean_reference'])} % (reference).")
    w(f"     - {hq}: {pct(ab[hq]['max'])} % in {ab[hq]['max_core']}, and ≤ "
      f"{pct(max(ab[hq][p] for p in PRE if p != ab[hq]['max_core']))} % in every other core.")
    w(f"4. **The annotations do not support TMA production by any of the four TMA-route candidates.**")
    zero310 = [n for n, l in (("CutC", "cutC"), ("CutD", "cutD"), ("CntB", "cntB")) if not in310(l)]
    zeroT = [n for n, l in (("a glycyl-radical enzyme", "glycyl radical enzyme"), ("Tmm", "tmm"),
                            ("TorA/TorZ", "torA/torZ (TMAO reductase)"), ("TMA dehydrogenase", "TMA dehydrogenase")) if not nT(l)]
    w(f"   - None of the 310 MAGs has {', '.join(zero310[:-1]) + ' or ' + zero310[-1] if len(zero310) > 1 else ''.join(zero310)} "
      f"annotated. None of the four candidates has {', '.join(zeroT[:-1]) + ' or ' + zeroT[-1]}.")
    bet = guilds[guilds.betaine_reducer]
    w(f"   - The TMA-forming enzymes the annotations do show are betaine reductase protein B (GrdH/GrdI), in the "
      f"{len(bet)} betaine reducers (a separate guild; table in section 2), and TMAO reductase: TorA/TorZ-type labels in "
      f"{len(tor_t)} MAGs, plus DMSO-reductase genes with a dual TMAO-reductase label in {len(tor_d)} more "
      f"({', '.join(tor_d)}). None of the four TMA candidates carries either enzyme.")
    w(f"   - The three \"clustered with choline TMA-lyase\" hits are shell or EutJ-like genes of bacterial microcompartment "
      f"(BMC) loci.")
    w(f"     - In {caz} the labelled gene shares its contig with ethanolamine ammonia-lyase (EutBC) and other Eut/Pdu "
      f"genes, which points to an ethanolamine, not choline, BMC (*inference*).")
    w(f"     - In {t1} EutBC is on a different contig from the labelled gene, in a bin with "
      f"{R[t1]['checkm']['contamination']:.1f} % contamination, so the BMC type of the labelled locus is unresolved.")
    w(f"     - In {uba} the locus runs off the contig edge.")
    w(f"   - {j153} has one CntA-labelled oxygenase, with no CntB-labelled gene in the genome and no reductase next to it. "
      f"Its other annotations describe routes that break choline, betaine and carnitine down without TMA, but by the same "
      f"standard these are incomplete too:")
    w(f"     - choline → betaine → dimethylglycine → sarcosine → glycine: GbcA is likewise a lone Rieske-oxygenase subunit "
      f"(no GbcB{', which is annotated only in ' + andjoin(gbcB) if gbcB else ''}), and no "
      f"dimethylglycine → sarcosine enzyme is annotated;")
    w("     - carnitine 3-dehydrogenase.")
    w(f"5. **{caz} (2020 ID {id2020}) is the only one of the four in the 2020 betaine models.**")
    w(f"   - Its modelled TMA flux of {', '.join(sorted({str(round(x, 4)) for x in mod[caz]['tma'].values()}))} comes from "
      f"`{mod[caz]['tma_rxns'][0][0]}` (\"{mod[caz]['tma_rxns'][0][1]}\"). That reaction has **no gene association**. "
      f"It equals the betaine uptake bound ({', '.join(sorted({str(round(-x, 4)) for x in mod[caz]['bet'].values()}))}, "
      f"set in the 2020 medium), which the solution reaches.")
    w(f"   - The genome encodes glycine-reductase protein B genes on {grd_contigs} contigs, but no betaine-type GrdH. Against the B proteins "
      f"of other MAGs, {grd_gly_wins} of its {len(grd)} B-protein paralogs (≥ 150 aa) score higher with glycine-type "
      f"GrdB/GrdE than with betaine-type GrdH/GrdI.")
    w(f"6. **Where the candidates live.**")
    w(f"   - The two candidates with the most TMA-like annotations, {j153} and {uba}, are reference-wetland organisms:")
    w(f"     - {j153}: mean {pct(ab[j153]['mean_reference'])} % in R2A, ≈ {pct(ab[j153]['mean_unrestored'], 2)} % elsewhere;")
    w(f"     - {uba}: mean {pct(ab[uba]['mean_reference'])} % in R2A, {pct(ab[uba]['mean_unrestored'], 4)} % in unrestored cores.")
    w(f"   - R2A is the site with net CH₄ uptake, so these two are unlikely TMA suppliers for methanogenesis in the "
      f"unrestored cores (*inference*).")
    if meth:
        w(f"7. **The TMA → CH₄ step itself is genomically present.**")
        w(f"   - {meth} encodes pyrrolysine-containing MttB, MtbB and MtmB (TMA, DMA and MMA methyltransferases), "
          f"their corrinoid proteins, MtbA and the pyrrolysine machinery.")
        w(f"   - The missing link is a demonstrated TMA source, not the methanogen.")
    w("")

    # ---------------- methods
    w("## How the inventory was made\n")
    w("**Gene rules.** The `RULES` regular expressions are matched against each RAST feature's function string, so a "
      "label is only as specific as the string. Some labels are broader:")
    w("- *phnCDE-like ABC (cluster 12)*, *phnP* and *CaiB/BaiF family* are generic;")
    w("- *FMO, generic*: Tmm is a flavin-containing monooxygenase, but so are many other enzymes;")
    w("- *MttB superfamily, non-pyl* is the bacterial non-pyrrolysine family. Its characterised members include a glycine "
      "betaine:corrinoid methyltransferase (Ticak et al. 2014, *PNAS* 111:E4668), which demethylates betaine rather than "
      "making TMA.\n")
    w("**Coordinates** are 1-based and inclusive, with start ≤ stop; strand is given separately. *gap* is the intergenic "
      "distance to the previous gene on the contig (negative = overlap). *gene_on_contig* (CSV) gives the gene's rank on "
      "its contig. *at_contig_edge* means the gene ends within 30 bp of a contig end, so it may be truncated.\n")
    w("**Recoded genes.** RAST reads an in-frame Sec (UGA) or Pyl (UAG) codon as a stop, so such a gene can appear as "
      "two ORFs. Where a \"selenocysteine-\" or \"pyrrolysine-containing\" ORF ends at that codon and its in-frame "
      "read-through ends at the stop of the next ORF, the two ORFs are counted as one gene; the CSV `evidence` column "
      "marks both.\n")
    w("**Sequence checks.** Local Smith–Waterman alignment of RAST protein translations (BLOSUM62, gap −11/−1).")
    w(f"- The null is {N_NULL} random proteins (> 80 aa) from the 310 MAGs: \"null 99th / max\" is that distribution's "
      f"99th percentile and maximum score.")
    w(f"- Reference sets are the RAST-annotated proteins of the 310 MAGs: PhnG n = {len(fam['phnG'])}, "
      f"PhnN n = {len(fam['phnN'])}, PhnP n = {len(fam['phnP'])}, PhnM n = {len(fam['phnM'])}, PhnJ n = {len(fam['phnJ'])}; "
      f"GrdH {len(fam['GrdH (betaine B, Sec)'])}, GrdI {len(fam['GrdI (betaine B, alpha)'])}, "
      f"GrdB {len(fam['GrdB (glycine B, Sec)'])}, GrdE {len(fam['GrdE (glycine B, alpha/beta)'])}. "
      "A query is never compared with genes from its own MAG.")
    w("- Raw scores are not E-values. A score well above the null maximum is good evidence of homology; a profile-HMM "
      "search would be the proper test.\n")
    w("**Abundance** is the relative abundance of reads in `data/Cliff_310MAG_relabund.tsv`, with the two depths merged "
      "by DNA weighting as in the R232 rerun. It is shown in % of reads. *Detected* means > 1e-5.\n")
    nB = sum(1 for v in R.values() if v["meta"]["domain"] == "B")
    nArB = sum(1 for v in R.values() if v["meta"]["domain"] == "B" and v["row"].archaea)
    w(f"**RAST domain.** {nB} of the {len(R)} genome objects examined here carry `domain: \"B\"` (genetic code 11), "
      f"including all {nArB} archaeal ones, so the archaea were annotated with bacterial settings. RAST k-mer calls "
      f"can also miss divergent archaeal homologues, as PhnG shows here.\n")

    # ---------------- archaeal
    w("## 1. Archaeal MePn consumers\n")
    w(overview_rows(R, A) + "\n")
    w("Guild flags and 2020 modelling:\n")
    w(flags_table(R, A) + "\n")
    w("Relative abundance per core (% of reads, DNA-weighted depth merge):\n")
    w(abundance_table(R, A) + "\n")
    tab, absent = inventory_matrix(R, A, PHN_ORDER + PI_ORDER, df)
    w("### Gene inventory (count, contigs)\n")
    w(tab + "\n")
    w(f"Not annotated in any of the three: {', '.join(absent)}.\n")

    w(f"### {nat}\n")
    r = R[nat]
    labN = {(h["contig"], h["idx"]): h["gene"] for h in r["hits"] if h["panel"] == "phn"}
    w(f"**Taxonomy and quality.**")
    w(f"- *{r['row'].species}*: ANI {r['gtdb']['closest_genome_ani']:.2f} % to {r['gtdb']['closest_genome_reference']}.")
    w(f"- CheckM {r['checkm']['completeness']:.1f} % complete; the bin is fragmented ({r['contigs']['n_contigs']} contigs, "
      f"N50 {r['contigs']['n50_bp'] / 1000:.1f} kb).\n")
    w(f"**Abundance.** Mean {pct(ab[nat]['mean_unrestored'])} % in the unrestored cores (maximum {pct(ab[nat]['max'])} % in "
      f"{ab[nat]['max_core']}), and at most {pct(max(ab[nat]['mean_restored'], ab[nat]['mean_reference']))} % on average "
      f"at the restored and reference sites.\n")
    w("**C-P lyase locus.** It spans two contigs:")
    w("- the first holds *phnM*-like – *phnL* – *phnK* – *phnJ* – the 3′ part of *phnI*, all on the − strand;")
    w("- the second holds the 5′ part of *phnI* – *phnH* – the unannotated PhnG candidate, also on −, with *phnD* – *phnC* "
      "– *phnE* divergent on +.\n")
    if split:
        w(f"The two *phnI* fragments are probably one gene split by the assembly:")
        w(f"- the fragment on `{node(split['a']['contig'])}` is {len(split['a']['protein'])} aa and sits at a contig edge;")
        w(f"- the fragment on `{node(split['b']['contig'])}` is {len(split['b']['protein'])} aa and also sits at a contig edge;")
        w(f"- stitching the contigs at the ends of their alignment and translating across the junction gives one "
          f"{split['joined']}-aa PhnI ORF, in which the two fragments overlap by {split['overlap_aa']} aa "
          f"({hq}'s intact PhnI is {len(hqI['protein'])} aa);")
        w(f"- the contig ends align over {split['dna']['cols']} nt at {split['dna']['ident']:.1%} identity "
          f"({split['a']['contig'].split('_')[1]} {split['dna']['a_end']} against "
          f"{split['b']['contig'].split('_')[1]} {split['dna']['b_end']}, {split['dna']['orient']}).\n")
        syn = split["n_mm"] - split["n_mm_nonsyn"]
        dif = "; ".join(f"residue {c}: {x} in the {split['left_contig'].split('_length')[0]} copy, {y} in the "
                        f"{split['right_contig'].split('_length')[0]} copy" for c, x, y in split["aa_diffs"])
        w(f"The nucleotide overlap is not exact. Of its {split['n_mm']} mismatches, {syn} are synonymous and "
          f"{split['n_mm_nonsyn']} {'changes' if split['n_mm_nonsyn'] == 1 else 'change'} the protein, so the "
          f"{split['overlap_aa']}-aa protein overlap differs at {len(split['aa_diffs'])} residue"
          f"{'' if len(split['aa_diffs']) == 1 else 's'} ({dif}). That fits contigs from two closely related haplotypes "
          f"in the bin (*inference*). Read pairs spanning the junction would settle it.\n")
    w(f"Two genes in the locus carry no RAST PhnG/PhnM label:")
    w(f"- The *phnM* position holds \"Metal-dependent hydrolase involved in phosphonate metabolism\", "
      f"{sM['identity']:.0%} identical to {hq}'s RAST-annotated PhnM (EC 3.6.1.63).")
    w(f"- The PhnG position holds a {len(gN['gene']['protein'])}-aa \"{gN['gene']['function'] or 'no function'}\". It "
      f"overlaps *phnH* by {-gN['gap_bp']} bp on the same strand and scores {gN['best']['phnG'][0]:.0f} against RAST PhnG "
      f"(null maximum {gN['nullmax']:.0f}; best against PhnN, PhnP or PhnM ≤ "
      f"{max(gN['best'][f][0] for f in ('phnN', 'phnP', 'phnM')):.0f}).\n")
    w(f"PhnJ is {jb[nat][1]:.0%} identical (over {jb[nat][2]} aa) to the closest bacterial PhnJ in the 310 MAGs "
      f"({GI.loc[jb[nat][3]].iterativeID_r232}).\n")
    w("**Other genes.**")
    w(f"- A *phnP*-labelled metallo-β-lactamase is on another contig.")
    w(f"- No *phnF* or *phnN*, and no phosphonatase (*phnW*/*phnX*), *phnY*/*phnZ*, *mpnS* or *pepM*.")
    w(f"- The \"carboxyvinyl-carboxyphosphonate phosphorylmutase\" (a PEP-mutase-family label) sits on an unrelated "
      f"contig with no phosphonate-biosynthesis neighbours. It is not evidence for MPn synthesis.")
    w(f"- Pi genes: *pstS* ×{cnt(nat, 'pstS')}, *pstC* ×{cnt(nat, 'pstC')}, *pstA* ×{cnt(nat, 'pstA')}, "
      f"*pstB* ×{cnt(nat, 'pstB')}, and {cnt(nat, 'Pi-regulatory homolog (unspecified)')} \"phosphate regulatory "
      f"protein homolog\". No *phoB*/*phoR* is annotated.\n")
    w("**Verdict.** All the subunits that make CH₄ from MPn are present (PhnG by inference): uptake (PhnCDE), "
      "RPnTP synthesis (PhnG/H/I/L), PhnM hydrolysis of RPnTP to PRPn, and the C–P lyase PhnJ; a role for PhnK was "
      "not established in vitro. The biochemistry is in Kamat, Williams & Raushel 2011, *Nature* 480:570.")
    w("Two caveats:")
    w("- the locus is split across contigs;")
    w("- whether haloarchaeal C-P lyase is expressed and active in situ is untested.\n")
    for c in sorted({g["contig"] for g in r["cpl_locus"]}, key=lambda c: -len([1 for k in labN if k[0] == c and labN[k] == "phnJ"])):
        idx = [g["idx"] for g in r["cpl_locus"] if g["contig"] == c]
        w(neighbourhood(r, c, min(idx), max(idx), labN) + "\n")

    w(f"### {hq}\n")
    r = R[hq]
    labH = {(h["contig"], h["idx"]): h["gene"] for h in r["hits"] if h["panel"] == "phn"}
    loc_c = sorted({g["contig"] for g in r["cpl_locus"]})
    gaps = [gap(a, b) for a, b in zip(r["cpl_locus"], r["cpl_locus"][1:]) if a["contig"] == b["contig"]
            and labH.get((a["contig"], a["idx"])) in CPL_CORE + [POS_G] and labH.get((b["contig"], b["idx"])) in CPL_CORE + [POS_G]]
    w(f"**Taxonomy and quality.**")
    w(f"- *{r['row'].species}*: ANI {r['gtdb']['closest_genome_ani']:.2f} % to {r['gtdb']['closest_genome_reference']}.")
    w(f"- CheckM {r['checkm']['completeness']:.1f} % complete; {r['contigs']['n_contigs']} contigs, "
      f"N50 {r['contigs']['n50_bp'] / 1000:.1f} kb.\n")
    w(f"**Abundance.** Almost all of it is in one core: {pct(ab[hq]['max'])} % in {ab[hq]['max_core']}, the core whose "
      f"D2 sample the bin came from. Every other core has ≤ {pct(max(ab[hq][p] for p in PRE if p != ab[hq]['max_core']))} %, "
      f"and it is detected in {ab[hq]['n_detected']}/12 cores.\n")
    w(f"**C-P lyase locus.** One contig, `{node(loc_c[0])}`, carries the whole locus.")
    w(f"- *phnD* – *phnC* – *phnE* – *phnO*-like acetyltransferase run on the − strand.")
    w(f"- Divergent on +: the PhnG candidate – *phnH* – *phnI* – *phnJ* – *phnK* – *phnL* – *phnM* – an unannotated gene.")
    w(f"- From the PhnG candidate to *phnM* the genes overlap or abut: gaps of {min(gaps)} to {max(gaps)} bp.\n")
    w(f"The PhnG candidate is a {len(gH['gene']['protein'])}-aa gene with no RAST function. It scores "
      f"{gH['best']['phnG'][0]:.0f} against RAST PhnG (null maximum {gH['nullmax']:.0f}) and is {sG['identity']:.0%} "
      f"identical to {nat}'s candidate. PhnJ is {jb[hq][1]:.0%} identical to the closest bacterial PhnJ "
      f"({GI.loc[jb[hq][3]].iterativeID_r232}).\n")
    w("**Other genes.**")
    w(f"- A second phosphonate-type ABC transporter (*phnC*/*phnD* plus a \"cluster 12\" permease) on another contig.")
    w(f"- *phnP* elsewhere.")
    w(f"- Pi genes: *pstS* ×{cnt(hq, 'pstS')}, *pstC* ×{cnt(hq, 'pstC')} (at a contig edge) and *phoU* ×{cnt(hq, 'phoU')}. "
      f"No *pstA*/*pstB* is annotated, which may reflect the 83.5 % completeness (*inference*).\n")
    w("**Verdict.** This is the most complete and best-assembled archaeal C-P lyase locus: CDE plus GHIJKLM on one "
      "contig, colinear with Natrinema.1. It is plausible as an MPn → CH₄ route.")
    w("- In this dataset the organism is abundant in only one unrestored core.\n")
    idx = [g["idx"] for g in r["cpl_locus"]]
    w(neighbourhood(r, loc_c[0], min(idx), max(idx), labH) + "\n")

    w(f"### {hr}\n")
    r = R[hr]
    labR = {(h["contig"], h["idx"]): h["gene"] for h in r["hits"] if h["panel"] == "phn"}
    hm = hitsof(hr, "phnM")[0]
    w(f"**Taxonomy and quality.**")
    w(f"- *Halorubrum* sp.; GTDB-Tk assigns no species. The placement neighbour is {r['gtdb']['closest_placement_reference']} "
      f"at {r['gtdb']['closest_placement_ani']:.2f} % ANI, which suggests an undescribed species (*inference*).")
    w(f"- CheckM {r['checkm']['completeness']:.1f} % / {r['checkm']['contamination']:.1f} %; N50 "
      f"{r['contigs']['n50_bp'] / 1000:.1f} kb.\n")
    w(f"**Abundance.** Mean {pct(ab[hr]['mean_unrestored'])} % in the unrestored cores (maximum {pct(ab[hr]['max'])} % in "
      f"{ab[hr]['max_core']}).\n")
    w("**phn genes.**")
    w(f"- The only C-P lyase-pathway gene is a lone *phnM* (EC 3.6.1.63; {len(hm['protein'])} aa) on "
      f"`{node(hm['contig'])}`, between unrelated genes.")
    w(f"- It is {', '.join(f'''{x['identity']:.0%} identical to {x['subject']}'s PhnM''' for x in hrM)} (null maximum "
      f"{hrM[0]['nullmax']:.0f}), so it is a genuine PhnM-family homologue.")
    w(f"- No PhnG, H, I, J, K or L is annotated, and no contig carries a C-P lyase cluster.")
    w(f"- The phosphonate-type ABC genes (one *phnC*-labelled ATPase plus a \"cluster 12\" permease) and the *phnP*-labelled "
      f"metallo-β-lactamase are generic.")
    w(f"- Pi genes: *pstS*, *pstC*, *pstA*, *pstB* ×2 each and *phoU* ×2, i.e. two Pst systems.\n")
    w(f"**Verdict.** Not an MPn → CH₄ organism by annotation.")
    w(f"- Without PhnJ and the PhnGHIL complex, PhnM has no substrate.")
    w(f"- At {r['checkm']['completeness']:.0f} % completeness, losing a whole seven-gene operon to incompleteness is "
      f"possible but unlikely (*inference*).")
    w(f"- The lone PhnM could be a relic of a lost operon or serve another function (*inference*, unverified).\n")
    w(neighbourhood(r, hm["contig"], hm["idx"] - 3, hm["idx"] + 3, labR) + "\n")

    w("### Gene-by-gene comparison of the two C-P lyase loci\n")
    w(f"Each gene of the {hq} locus was aligned with every gene of the {nat} locus; the best hit is shown.\n")
    rows = [[s["q"]["fid"], s["q_label"] or "–", esc(s["q"]["function"] or "(no function)")[:60], len(s["q"]["protein"]),
             s["s_label"] or "–", esc(s["s"]["function"] or "(no function)")[:50], f"{s['score']:.0f}",
             f"{s['identity']:.0%} / {s['cols']}", f"{s['null99']:.0f} / {s['nullmax']:.0f}"] for s in synteny]
    w(md_table([f"{hq} feature", "label", "RAST function", "aa", f"best {nat} hit: label", "its RAST function",
                "score", "identity / aligned aa", "null 99th / max"], rows) + "\n")
    w(f"- The *phnO*-like acetyltransferase of {hq} matches the gene in the same position in {nat} "
      f"({sO['identity']:.0%} identity), which RAST calls \"{sO['s']['function']}\".")
    w(f"- The genes after *phnM* are {sX['identity']:.0%} identical to each other. Neither scores above the null against "
      f"RAST PhnN, PhnP, PhnG or PhnM (best {max(b[0] for b in xH['best'].values()):.0f}; null maximum "
      f"{xH['nullmax']:.0f}). This is a conserved haloarchaeal gene of unknown function, not a recognisable PhnN.\n")

    w("### Context: phn genes in all 30 archaeal MAGs\n")
    w("Only the two *phnJ* carriers have any C-P lyase core subunit:")
    w(f"- **H–L subunits:** {n_cpl}/{n_arch} archaea have any of them.")
    w(f"- **PhnP labels:** {n_phnP}/{n_arch} archaea carry one.")
    w(f"- **PhnC/D/E labels:** {n_phnCD}/{n_arch} archaea carry one.")
    w(f"- **PhnM-family genes:** {n_m}/{n_arch} archaea carry one.\n")
    w("PhnP, PhnC/D/E and PhnM-family genes outside a locus are therefore not evidence of MePn use.\n")
    ctx = arch_ctx[(arch_ctx.drop(columns=["iterativeID", "mepn_consumer", "mepn_2020"]).sum(axis=1)) > 0]
    cols = ["phnC", "phnD", "phnE"] + CPL_CORE + ["phnM-like", "phnN", "phnP", "phnO-like"]
    w(md_table(["archaeal MAG", "phnJ set", "2020 set"] + cols,
               [[x["iterativeID"], "yes" if x["mepn_consumer"] else "", "yes" if x["mepn_2020"] else ""]
                + [x[c] or "" for c in cols] for _, x in ctx.iterrows()]) + "\n")
    w(f"The other {n_arch - len(ctx)} archaea have none of these labels. PhnG is annotated in "
      f"{int((arch_ctx.phnG > 0).sum())} of the {n_arch}; for the two *phnJ* carriers it was found by position.\n")

    # ---------------- TMA
    w("## 2. TMA-route candidates\n")
    w("These MAGs were flagged by `tma_producer` in `guilds_r232.csv`: either a CntA annotation, or a gene RAST calls "
      "\"... clustered with choline trimethylamine-lyase\".\n")
    w("Known TMA-forming routes and what each needs:")
    w("- **Choline → TMA:** the glycyl-radical enzyme CutC plus its activase CutD, usually inside a BMC (Craciun & "
      "Balskus 2012, *PNAS* 109:21307).")
    w("- **Carnitine → TMA:** the two-component Rieske oxygenase CntA/CntB (Zhu et al. 2014, *PNAS* 111:4268).")
    w("- **Betaine → TMA:** the selenoprotein betaine reductase, protein B GrdH/GrdI.")
    w("- **TMAO → TMA:** reduction by TorA/TorZ.\n")
    w(f"CutC belongs to the pyruvate formate-lyase (PFL) family of glycyl-radical enzymes. A CutC that RAST failed to "
      f"recognise would therefore most likely be labelled \"Pyruvate formate-lyase\": {in310('glycyl radical enzyme')} of "
      f"the 310 MAGs carry such a glycyl-radical enzyme label, but none of the four candidates does. Across the 310 MAGs, "
      f"the counts carrying a CutC, CutD, CntB, Tmm or TMA-dehydrogenase label are {in310('cutC')}, {in310('cutD')}, "
      f"{in310('cntB')}, {in310('tmm')} and {in310('TMA dehydrogenase')}.\n")
    bet = guilds[guilds.betaine_reducer]
    w(f"**For comparison: the {len(bet)} betaine reducers (*grdH*).** Betaine reductase turns glycine betaine into TMA "
      f"and acetyl phosphate. It is one of two TMA-forming enzymes annotated in this MAG set. The other is TMAO "
      f"reductase: TorA/TorZ-type labels in {len(tor_t)} MAGs ({', '.join(tor_t)}), and a dual "
      f"\"@ ... trimethylamine-N-oxide reductase\" label on DMSO-reductase genes in {len(tor_d)} more "
      f"({', '.join(tor_d)}). None of the four TMA candidates carries either enzyme. The betaine reducers were not "
      f"examined gene by gene here; for the pyl-MttB column, see the end of the {meth} section.\n")
    rows = []
    for _, b in bet.sort_values("iterativeID_r232").iterrows():
        a = abundance(b["mag"])
        rows.append([b["iterativeID_r232"], f"{b['phylum']}; {b['family']}", "yes" if b["betaine_2020"] else "",
                     "yes" if key(b["mag"]) in TMA_MODEL else "",
                     ("" if b["iterativeID_r232"] not in glob.get("mttB (pyl, TMA MT)", ()) else
                      "yes" if pyl_ok(b["iterativeID_r232"]) else "label only (no Pyl evidence)"),
                     pct(a["mean_unrestored"]), pct(a["mean_restored"]), pct(a["mean_reference"])])
    w(md_table(["MAG", "phylum; family", "2020 betaine set", "2020 betaine model", "pyl-MttB", "mean unrestored %",
                "mean restored %", "mean reference %"], rows) + "\n")
    w(overview_rows(R, T) + "\n")
    w(flags_table(R, T) + "\n")
    w("Relative abundance per core (% of reads):\n")
    w(abundance_table(R, T) + "\n")
    tab, absent = inventory_matrix(R, T, TMA_ORDER, df)
    w("### Gene inventory (count, contigs)\n")
    w(tab + "\n")
    w(f"Not annotated in any of the four: {', '.join(absent)}.\n")

    def tma_loci(i, anchors, win=3):
        r = R[i]
        lab = {(h["contig"], h["idx"]): h["gene"] for h in r["hits"]}
        out = []
        for c in sorted({h["contig"] for h in r["hits"] if h["gene"] in anchors}):
            idx = [h["idx"] for h in r["hits"] if h["contig"] == c and h["gene"] in anchors]
            out.append(neighbourhood(r, c, min(idx) - win, max(idx) + win, lab))
        return "\n\n".join(out)

    def tax_line(i):
        r = R[i]
        g = r["gtdb"]
        sp = r["row"].species if isinstance(r["row"].species, str) else "no species assigned"
        extra = f"; ANI {g['closest_genome_ani']:.2f} % to {g['closest_genome_reference']}" if isinstance(g.get("closest_genome_reference"), str) else ""
        return (f"{r['row']['phylum']}; {r['row']['order']}; {r['row']['family']}; *{sp}*{extra}. CheckM "
                f"{r['checkm']['completeness']:.1f} % / {r['checkm']['contamination']:.1f} %; {r['contigs']['n_contigs']} "
                f"contigs, N50 {r['contigs']['n50_bp'] / 1000:.1f} kb.")

    def ab_line(i):
        a = ab[i]
        return (f"Mean {pct(a['mean_unrestored'])} % (unrestored), {pct(a['mean_restored'])} % (restored), "
                f"{pct(a['mean_reference'])} % (reference); maximum {pct(a['max'])} % in {a['max_core']}; detected "
                f"in {a['n_detected']}/12 cores.")

    # UBA11579.1
    w(f"### {uba}\n")
    cc = hitsof(uba, "cutC-cluster label (BMC gene)")[0]
    gre = hitsof(uba, "GRE-activase-like (COG1180)")
    w(f"**Taxonomy and quality.** {tax_line(uba)}\n")
    w(f"**Abundance.** {ab_line(uba)} This is a reference-wetland organism.\n")
    w(f"**Evidence for the TMA flag.** \"{cc['function']}\" is gene {cc['idx'] + 1} of {cc['n_on_contig']} on "
      f"`{node(cc['contig'])}`, {cc['start'] - 1} bp from the start of the contig.")
    w(f"- Its neighbours are a PduA/PduJ-like shell protein and an acetaldehyde dehydrogenase (\"ethanolamine "
      f"utilization cluster\"), so the BMC locus runs off the contig.")
    w(f"- No glycyl-radical enzyme of any kind, and no EutBC, is annotated anywhere in the bin.")
    w(f"- The COG1180 activase-like gene ({len(gre)}) is on `{node(gre[0]['contig'])}`, with no BMC genes nearby.\n")
    w("**Other genes.**")
    w(f"- choline dehydrogenase ×{cnt(uba, 'betA/betB (choline -> betaine)')}, both at contig edges;")
    w(f"- choline-sulfatase ×{cnt(uba, 'betC (choline-sulfatase)')};")
    w(f"- one MttB-superfamily (non-pyl) gene and {cnt(uba, 'corrinoid protein (unspecified)')} unspecified corrinoid "
      f"proteins;")
    w("- no BCCT or ProU-type betaine transporter is annotated.\n")
    w("**Verdict.** Weak. The flag rests on the family label of a BMC accessory protein. The enzyme that defines the "
      "BMC type is not in the bin, which is 78 % complete. The annotation supports choline oxidation to betaine, not "
      "TMA release.\n")
    w(tma_loci(uba, TMA_ANCHORS) + "\n")

    # J153.1
    w(f"### {j153}\n")
    ca = hitsof(j153, "cntA")[0]
    w(f"**Taxonomy and quality.** {tax_line(j153)}\n")
    w(f"**Abundance.** {ab_line(j153)} Like {uba}, a reference-wetland organism.\n")
    w(f"**Evidence for the TMA flag.** One \"{ca['function']}\" ({len(ca['protein'])} aa) on `{node(ca['contig'])}`.")
    w(f"- Its neighbours are unrelated (glutamine synthetase, P-II, aspartate aminotransferase).")
    w(f"- CntB is annotated in {in310('cntB')} of the 310 MAGs; CntA in {in310('cntA')}.")
    fr = sorted((g for c in R[j153]["by"].values() for g in c
                 if g["function"].lower().startswith("flavodoxin reductases")), key=lambda g: g["contig"])
    vanb = [g for g in fr if "vanillate o-demethylase oxidoreductase" in g["function"].lower()]
    assert not [g for g in fr if g["contig"] == ca["contig"]]       # no flavin reductase on the CntA contig
    w(f"- CntA is a Rieske-type oxygenase. A k-mer label alone cannot separate it from related Rieske oxygenases that act "
      f"on other quaternary amines (*inference*). With no CntB-labelled gene and no reductase next to CntA, the "
      f"annotation does not support the TMA route.")
    w(f"- A partner elsewhere is not excluded. The genome has {len(fr)} \"Flavodoxin reductases (ferredoxin-NADPH "
      f"reductases) family 1\" genes on other contigs. One of them, {len(vanb[0]['protein'])} aa on "
      f"`{node(vanb[0]['contig'])}`, is also labelled \"Vanillate O-demethylase oxidoreductase\" (VanB type), the same "
      f"flavin-reductase plus [2Fe-2S] architecture as CntB (*inference*).\n")
    w("**Other genes.** The annotation describes aerobic catabolism of the same substrates by non-TMA routes, though "
      "not every step is annotated:")
    w(f"- choline uptake: ChoXWV, plus {cnt(j153, 'ProU/OpuA/Ota/Cho ABC')} ProU/OpuA/Ota/Cho ABC parts in all;")
    w(f"- choline → betaine: choline dehydrogenase;")
    w(f"- betaine demethylation: GbcA only. Its reductase partner GbcB is not annotated (in the 310 MAGs it is "
      f"annotated only in {andjoin(gbcB)}), so GbcA is as much a lone oxygenase subunit as CntA;")
    w(f"- no dimethylglycine dehydrogenase or demethylase is annotated, so the dimethylglycine → sarcosine step has no "
      f"annotated enzyme;")
    w(f"- sarcosine oxidase/dehydrogenase: ×{cnt(j153, 'sarcosine oxidase/DH')} genes, counting the α, β, γ and δ "
      f"subunits of heterotetrameric sarcosine oxidase separately;")
    w(f"- betaine–homocysteine methyltransferase;")
    w(f"- the carnitine 3-dehydrogenase route (\"carnitine 3-dehydrogenase / carnitine utilization associated "
      f"thioesterase\" and \"acyl-CoA synthetase clustered with carnitinyl-CoA dehydratase\");")
    w(f"- {cnt(j153, 'MttB superfamily, non-pyl')} MttB-superfamily (non-pyl) genes, each on a different contig, and "
      f"one generic FMO.\n")
    w("**Verdict.** Weak for TMA. By the same standard the non-TMA route is incomplete too: GbcA, like CntA, has no "
      "annotated reductase partner, and the dimethylglycine → sarcosine step is unannotated. The annotation leans "
      "toward choline, betaine and carnitine catabolism without TMA (choline dehydrogenase, sarcosine "
      "oxidase/dehydrogenase, the carnitine 3-dehydrogenase route), but supports neither oxygenase step on its own.\n")
    w(tma_loci(j153, {"cntA", "betaine demethylation (gbcA/dgc)", "carnitine 3-DH route (non-TMA)"}) + "\n")

    # CAZFJG01.2
    w(f"### {caz}\n")
    cc = hitsof(caz, "cutC-cluster label (BMC gene)")[0]
    w(f"**Taxonomy and quality.** {tax_line(caz)} In 2020 it was *Synergistales*.2.\n")
    w(f"**Abundance.** {ab_line(caz)} An unrestored-core organism.\n")
    w(f"**2020 model.** This MAG is in the 15-MAG betaine model set.")
    w(f"- Modelled TMA flux: {', '.join(sorted({str(round(x, 4)) for x in mod[caz]['tma'].values()}))} in each of the "
      f"{len(mod[caz]['tma'])} unrestored cores.")
    w(f"- That flux runs through `{mod[caz]['tma_rxns'][0][0]}` \"{mod[caz]['tma_rxns'][0][1]}\", which has gene "
      f"association: **{mod[caz]['tma_rxns'][0][2]}**. The reaction was added in curation, not derived from the "
      f"annotation.")
    w(f"- The flux equals the betaine uptake bound: the 2020 medium caps betaine uptake (`EX_cpd00540_e0`) at "
      f"{', '.join(sorted({str(round(-x, 4)) for x in mod[caz]['bet'].values()}))}, and the solution reaches it "
      f"(stored betaine flux {', '.join(sorted({str(round(x, 4)) for x in mod[caz]['bet'].values()}))}).\n")
    w(f"**Evidence for the TMA flag.** \"{cc['function']}\" sits on `{node(cc['contig'])}` inside a BMC locus that also "
      f"encodes ethanolamine ammonia-lyase heavy and light chains (EutBC), EutL, EutJ, PduL/PduU-like proteins and an "
      f"acetaldehyde/alcohol dehydrogenase.")
    w(f"- That is the gene content of an ethanolamine-utilisation BMC (*inference*); the TMA-lyase label belongs to the "
      f"shell-protein family.")
    w(f"- No glycyl-radical enzyme is annotated in the bin, which is {R[caz]['checkm']['completeness']:.0f} % complete.\n")
    w("**Betaine reductase.** No GrdH/GrdI (betaine-type protein B) is annotated. The bin encodes:")
    gB, gA = "grdB/E (glycine reductase B)", "grdA"
    ef = edge_frag(caz, gB)
    w(f"- {genes_of(caz, gB)} glycine-reductase protein B genes (GrdB/GrdE labels) on {grd_contigs} contigs. RAST "
      f"called them as {cnt(caz, gB)} ORFs, because it split {n_split_c(caz, gB)} Sec-containing GrdB into two labelled "
      f"ORFs at the UGA codon."
      + (f" The {genes_of(caz, gB)} include a fragment of {len(ef[0]['protein'])} aa at a contig edge "
         f"(`{node(ef[0]['contig']).split(' ')[0]}`);" if len(ef) == 1 else ";"))
    w(f"- GrdA ×{genes_of(caz, gA)}"
      + (f" (RAST split each at its Sec codon, giving {cnt(caz, gA)} ORFs)" if n_split_c(caz, gA) == genes_of(caz, gA)
         else f" ({cnt(caz, gA)} ORFs)" if n_split_c(caz, gA) else "") + " and one GrdC/D pair;")
    w(f"- betaine/choline transporters: OpuD (BCCT) and {cnt(caz, 'ProU/OpuA/Ota/Cho ABC')} ProU/Ota ABC parts.\n")
    w("Each B-protein paralog's best score against the B proteins of other MAGs:\n")
    w(md_table(["feature", "label", "aa", "vs GrdH (betaine, Sec)", "vs GrdI (betaine, α)", "vs GrdB (glycine, Sec)",
                "vs GrdE (glycine, α/β)"],
               [[h["fid"], h["gene"], len(h["protein"]), f"{s['GrdH (betaine B, Sec)']:.0f}", f"{s['GrdI (betaine B, alpha)']:.0f}",
                 f"{s['GrdB (glycine B, Sec)']:.0f}", f"{s['GrdE (glycine B, alpha/beta)']:.0f}"] for h, s in grd]) + "\n")
    w(f"{grd_gly_wins} of {len(grd)} paralogs are closer to the glycine-type than to the betaine-type B proteins.")
    w("- The \"glycine\" reference set is itself defined by RAST labels.")
    w("- Substrate specificity of protein B is only settled biochemically.")
    w("So this argues against, but does not rule out, a betaine-specific B protein (*inference*).\n")
    w("**Verdict.** No TMA-forming enzyme is annotated. The 2020 model's TMA output comes from a reaction without a "
      "gene, so that TMA flux is an assumption, not a genome prediction.\n")
    w(tma_loci(caz, {"cutC-cluster label (BMC gene)", "eutBC (ethanolamine ammonia-lyase)"}, 1) + "\n")
    w(tma_loci(caz, {"grdB/E (glycine reductase B)", "grdH/I (betaine reductase B)"}, 2) + "\n")

    # T1Sed10-7.1
    w(f"### {t1}\n")
    cc = hitsof(t1, "cutC-cluster label (BMC gene)")[0]
    eut = hitsof(t1, "eutBC (ethanolamine ammonia-lyase)")
    w(f"**Taxonomy and quality.** {tax_line(t1)}")
    w(f"- GTDB-Tk warns: \"{str(R[t1]['gtdb'].get('warnings')).strip()}\".")
    w(f"- With {R[t1]['checkm']['contamination']:.1f} % CheckM contamination, the bin is a mixture, and a single locus "
      f"cannot be tied to one organism.\n")
    w(f"**Abundance.** {ab_line(t1)}\n")
    w(f"**Evidence for the TMA flag.** \"{cc['function']}\" on `{node(cc['contig'])}` follows a PduA/PduJ-like shell "
      f"gene and a truncated acetaldehyde dehydrogenase at the contig start.")
    w(f"- A separate contig, `{node(eut[0]['contig'])}`, has EutA – EutB – EutC – EutL – acetaldehyde dehydrogenase, an "
      f"ethanolamine BMC core.")
    w(f"- No glycyl-radical enzyme, CutC/CutD, Grd, Tmm or Tor is annotated.")
    w(f"- Choline transport: BetT ×{cnt(t1, 'BCCT transporter (betT/opuD)')}.\n")
    assert all(e["contig"] != cc["contig"] for e in eut)
    w("**Verdict.** No TMA route is supported. The bin has ethanolamine BMC genes, but EutBC is on a different contig "
      "from the labelled gene, so the BMC type of the labelled locus is unresolved. The bin is also too contaminated to "
      "assign these genes to one genome.\n")
    w(tma_loci(t1, {"cutC-cluster label (BMC gene)", "eutBC (ethanolamine ammonia-lyase)"}, 2) + "\n")

    if meth:
        w(f"### Context: {meth}, the only methanogen\n")
        w(f"**Taxonomy and quality.** {tax_line(meth)}\n")
        w(f"**Abundance.** {ab_line(meth)}\n")
        tab, _ = inventory_matrix(R, [meth], MA_ORDER + TMA_ORDER, df)
        w(tab + "\n")
        w("**What this shows.** The genes for methanogenesis from TMA, DMA and MMA are present, including "
          "pyrrolysine-containing MttB, MtbB and MtmB.")
        w("- `guilds_r232.csv` column `ev_mtbB_mtmB` lists only the two corrinoid-protein strings for this MAG: its "
          "regular expression does not match the \"[...--corrinoid protein] Co-methyltransferase\" strings.")
        w("- That is a gap in the evidence column only. The `methanogen` flag uses *mcrA* and is unaffected.\n")
        others = sorted(set().union(*(glob.get(l, set()) for l in ("mttB (pyl, TMA MT)", "mtbB (pyl, DMA MT)",
                                                                   "mtmB (pyl, MMA MT)"))) - {meth})
        if others:
            gi = guilds.set_index("iterativeID_r232")
            ok, weak = [o for o in others if pyl_ok(o)], [o for o in others if not pyl_ok(o)]
            w(f"**Other MAGs with pyrrolysine-type methylamine methyltransferases.** {len(others)} bacterial MAGs carry "
              f"such labels, with no *mcrA*. In {len(ok)} of them ({', '.join(ok)}) the label is backed by Pyl evidence: "
              f"a labelled ORF ends at a UAG codon whose in-frame read-through runs into the next, unlabelled ORF (RAST "
              f"splitting the gene at the Pyl codon, as in {meth}), and/or Pyl-biosynthesis genes (PylS, PylB/C/D). "
              f"These are potential non-methanogenic TMA/DMA demethylators, i.e. competitors of {meth} for methylamines "
              f"(*inference*).")
            full = pylev[meth]["full_aa"].get("mttB (pyl, TMA MT)")
            for o in weak:
                w(f"- {o}'s labelled gene is a single {pylev[o]['longest']}-aa ORF"
                  + (f", close to the length of {meth}'s MttB read through its UAG ({full} aa)," if full else "")
                  + " so it is not a UAG-split fragment, and the bin has no Pyl-biosynthesis genes. It is probably a "
                    "non-pyl MttB-superfamily protein (*inference*); that family includes the MtgB glycine betaine "
                    "methyltransferase (Ticak et al. 2014).")
            w("")
            rows = []
            for o in others:
                m = gi.loc[o]
                labs = [l.split(" ")[0] for l in ("mttB (pyl, TMA MT)", "mtbB (pyl, DMA MT)", "mtmB (pyl, MMA MT)")
                        if o in glob.get(l, ())]
                a = abundance(m["mag"])
                e = pylev[o]
                ev = ", ".join(x for x in (f"{e['n_split']} UAG-split ORF{'s' if e['n_split'] > 1 else ''}" if e["n_split"] else "",
                                           f"{e['n_machinery']} Pyl-biosynthesis genes" if e["n_machinery"] else "") if x)
                rows.append([o, f"{m['phylum']}; {m['family']}", ", ".join(labs) + ("" if pyl_ok(o) else " (label only)"),
                             ev or "none: one continuous ORF, no Pyl-biosynthesis genes (probably non-pyl)",
                             "yes" if m["betaine_reducer"] else "", pct(a["mean_unrestored"]), pct(a["mean_restored"]),
                             pct(a["mean_reference"])])
            w(md_table(["MAG", "phylum; family", "pyl methyltransferases", "Pyl evidence", "betaine reducer (grdH)",
                        "mean unrestored %", "mean restored %", "mean reference %"], rows) + "\n")

    # ---------------- follow-up
    w("## What would test these readings\n")
    w("**MePn (Natrinema.1, Haloquadratum.4).**")
    w("1. **Confirm the gene calls.** Run profile-HMM searches (Pfam/TIGRFAM or KEGG KOfam) over the loci, especially "
      "the PhnG candidates and the gene after *phnM*. Then place PhnJ in a phylogeny with characterised bacterial PhnJ.")
    w("2. **Check that the locus is not a binning artefact.** Look for the same locus in the closest reference genomes "
      f"({R[nat]['gtdb']['closest_genome_reference']}, {R[hq]['gtdb']['closest_genome_reference']}). Map read pairs "
      f"across the {nat} *phnI* junction.")
    w("3. **Test expression.** Look for archaeal *phnJ* transcripts in the unrestored cores (metatranscriptome or "
      "RT-qPCR), and measure porewater Pi. In bacteria, C-P lyase is induced by Pi limitation; regulation in haloarchaea "
      "is not established (unverified).")
    w("4. **Test activity.** Amend unrestored cores with ¹³C-MPn, with and without added Pi, and measure ¹³CH₄. If an "
      "isolate of the reference species is available, test growth on MPn as the sole P source with CH₄ in the headspace.")
    ncur, have = curated_mepn_summary()

    def has_all(i):
        rx = mod[i]["base_cpl"] or {}
        miss = [r for r in CPL_RXNS if r not in rx]
        if not miss:
            return f"{i}'s base model already has {andjoin(list(CPL_RXNS))}" + (" with gene rules" if all(rx.values()) else "")
        m = hitsof(i, "phnM-like") or hitsof(i, "phnM")
        return (f"{i}'s lacks only {', '.join(miss)}" if len(miss) == 1 else f"{i}'s lacks {', '.join(miss)}") + (
            f", which its {m[0]['gene']} gene (CDS.{m[0]['fid']}) would carry" if miss == ["rxn26457"] and m else "")

    def gene_less(step):
        out = []
        for (st, rid, nm, g), n in sorted(have.items()):
            if st == step:
                k = f"all {ncur}" if n == ncur else f"{n} of {ncur}"
                name = f'"{nm}" ({rid}; {k})' if nm != rid else f"reaction, {rid} ({k})"
                out.append(f"{'' if g else 'a gene-less '}{name}")
        return "; ".join(out)

    w("5. **Modelling.** If MePn modelling is rerun:")
    w(f"   - {has_all(hq)}. {has_all(nat)}.")
    if not any(mod[i]["base_mepn_uptake"] for i in (hq, nat)):
        w(f"   - Neither has MePn uptake. The {ncur} curated models in `models/MePn_consumers/` add it as "
          f"{gene_less('uptake')} with an MePn exchange; here it could be tied to the PhnCDE genes. For the C-P lyase "
          f"step they use {gene_less('C-P lyase')}"
          + (", in place of rxn26459." if not any(k[1] == "rxn26459" for k in have) else "."))
    w(f"   - Drop {hr}.\n")
    w("**TMA (the four candidates).**")
    w("1. **Search the full assemblies, not just the bins.** Look for CutC/glycyl-radical domains and CntB-type "
      "reductases with HMMs, and extend the UBA11579.1 BMC contig in the assembly graph to find the locus's signature "
      "enzyme.")
    w("2. **Place the ambiguous genes phylogenetically.**")
    w("   - J153.1's CntA-like protein, alongside characterised CntA and other Rieske oxygenases;")
    w(f"   - {caz}'s GrdB paralogs, alongside characterised glycine-, sarcosine- and betaine-specific protein B.")
    w(f"   - Then compare with its 99.8 %-ANI reference genome, {R[caz]['gtdb']['closest_genome_reference']}.")
    w("3. **Incubations.** Add ¹³C-choline, ¹³C-carnitine and ¹³C-betaine to anoxic slurries of unrestored cores, with "
      "and without BES (a methanogenesis inhibitor, so TMA accumulates), and measure ¹³C-TMA and ¹³CH₄. Pair this with "
      "metatranscriptomics of *cutC*, *cntA* and *grd* genes.")
    w(f"4. **Modelling.** Until a TMA-forming gene is found, the betaine → TMA flux attributed to {caz} should be "
      "described as assumed.\n")

    w("## Caveats\n")
    w("- **RAST labels.** Function strings come from RAST k-mers (2020), and the archaea were annotated with domain "
      "\"B\". Missing labels for divergent archaeal genes are expected; PhnG is the clear case here.")
    w("- **Absence of a gene.** The bins are 78–100 % complete, so a missing gene is not proof the organism lacks it.")
    w("- **Alignment scores.** These are local-alignment scores against a random-protein null, not E-values, and "
      "positional homology is an inference.")
    w("- **Abundance.** It is relative read abundance (DNA-weighted depth merge), not activity, and it is bounded by "
      "binning: unbinned reads are not counted.")
    w("- **2020 sets.** The 2020 MePn and betaine sets were hand-curated; their gene criteria are not recorded beyond "
      "what `guilds_r232.csv` reconstructs.\n")
    OUT_MD.write_text("\n".join(L))
    print(f"wrote {OUT_MD.name}")


if __name__ == "__main__":
    main()
