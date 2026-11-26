"""Independent verifier: every taxonomy-bearing MAG file equals GTDB R232 by bin.

Reads the R232 summary from the prep results (not from the repo), the R220 lineage from git HEAD,
and for each file/column learns the HEAD convention (raw GTDB name vs suffix-stripped, s__ blank form)
then checks the working-tree cells against R232 with that same convention.
Read-only on the repo.
"""
import csv, io, json, re, subprocess, sys
from collections import Counter, defaultdict

REPO = "/home/freiburger/Documents/EmilyKin"
PREP = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PFX = dict(zip(RANKS, ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]))


def parse_lineage(s):
    out = {}
    for tok in s.split(";"):
        tok = tok.strip()
        for r, p in PFX.items():
            if tok.startswith(p):
                out[r] = tok[len(p):]
    return out


def load_summary(text):
    rows = list(csv.DictReader(io.StringIO(text.replace("\r\n", "\n")), delimiter="\t"))
    return {r["user_genome"]: r for r in rows}


r232 = load_summary(open(f"{PREP}/results/gtdbtk.r232.summary.tsv").read())
r220 = load_summary(subprocess.run(["git", "-C", REPO, "show", "HEAD:gtdbtk.bac120.summary.tsv"],
                                   capture_output=True, text=True, check=True).stdout)
assert len(r232) == 291 and set(r232) == set(r220)
L232 = {b: parse_lineage(r["classification"]) for b, r in r232.items()}
L220 = {b: parse_lineage(r["classification"]) for b, r in r220.items()}

strip = lambda x: re.sub(r"_[A-Z]+$", "", x) if x else x


def head_text(path):
    p = subprocess.run(["git", "-C", REPO, "show", f"HEAD:{path}"], capture_output=True, text=True)
    return p.stdout if p.returncode == 0 else None


def read_table(text, sep):
    return list(csv.DictReader(io.StringIO(text.replace("\r\n", "\n")), delimiter=sep))


SPECS = [
    # file, sep, bin column, {column: rank or 'lineage'}
    ("mag_quality_summary.tsv", "\t", "MAG", {r: r for r in RANKS}),
    ("mag_abundance_summary.tsv", "\t", "MAG", {r: r for r in RANKS}),
    ("mimag_hq_classification.tsv", "\t", "MAG", {"Phylum": "Phylum", "Genus": "Genus", "Species": "Species"}),
    ("Table_S7_MAG_functional_classification.csv", ",", "MAG", {"Genus": "Genus", "Family": "Family", "Phylum": "Phylum"}),
    ("mag_guild_classification.csv", ",", "MAG", {"Genus": "Genus"}),
    ("mag_percell_abundance.csv", ",", "MAG", {r: r for r in RANKS[1:]}),
    ("gene_ab_figure/data/taxonomy_labels.tsv", "\t", "MAG", {"Genus": "Genus", "Species": "Species"}),
    ("gene_ab_figure/offline_bundle/data/taxonomy_labels.tsv", "\t", "MAG", {"Genus": "Genus", "Species": "Species"}),
    ("gene_ab_figure/data/gene_copy_per_mag.tsv", "\t", "MAG", {"genus": "Genus"}),
    ("gene_ab_figure/data/nosz_clades.tsv", "\t", "mag", {"genus": "Genus"}),
    ("gene_ab_figure/data/nosz_cladeII_subclades.tsv", "\t", "mag", {"genus": "Genus"}),
    ("clade_classify/out/nosz_clades_updated.tsv", "\t", "mag", {"genus": "Genus"}),
    ("clade_classify/out/allbins_nosz_master.tsv", "\t", "MAG", {"genus": "Genus"}),
    ("gene_ab_figure/data/novel_bacteria_summary.tsv", "\t", "bin",
     {"phylum": "Phylum", "class": "Class", "order": "Order", "family": "Family", "genus": "Genus"}),
    ("gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase.tsv", "\t", "mag", {"genus": "Genus"}),
    ("gene_ab_figure/data/rpkm_long_MAG_gene_sample_phase_hq.tsv", "\t", "mag", {"genus": "Genus"}),
    ("gene_ab_figure/data/rpkm_long_MAG_gene_sample_291_preview.tsv", "\t", "mag", {"genus": "Genus"}),
    ("ko_gene_abundance/ko_rpkm_combined.tsv", "\t", "MAG", {"genus": "Genus"}),
    ("ko_gene_abundance/ko_rpkm_combined_ext.tsv", "\t", "MAG", {"genus": "Genus"}),
    ("ko_gene_abundance/ko_rpkm_combined_hq.tsv", "\t", "MAG", {"genus": "Genus"}),
    ("ko_gene_abundance/ko_rpkm_per_sample.tsv", "\t", "mag", {"genus": "Genus"}),
    ("ko_gene_abundance/ko_rpkm_per_sample_ext.tsv", "\t", "mag", {"genus": "Genus"}),
    ("ko_gene_abundance/ko_rpkm_per_sample_hq.tsv", "\t", "mag", {"genus": "Genus"}),
    ("phase_rpkm_chord/mag_category_assignments.csv", ",", "MAG", {"Genus": "Genus"}),
    ("phase_rpkm_chord/phase_mag_gene_rpkm_long.csv", ",", "MAG", {"Genus": "Genus"}),
    ("community_modeling/analysis/r13/annotations/taxonomy.tsv", "\t", "bin",
     {"gtdb_classification": "lineage", "phylum": "Phylum", "class": "Class", "order": "Order",
      "family": "Family", "genus": "Genus", "species": "Species"}),
    ("community_modeling/analysis/r13/annotations/panel/mag_panel.tsv", "\t", "bin",
     {"gtdb_phylum": "Phylum", "gtdb_class": "Class", "gtdb_family": "Family", "gtdb_genus": "Genus"}),
    ("community_modeling/analysis/r13/annotations/calls/summary.tsv", "\t", "bin",
     {"gtdb_phylum": "Phylum", "gtdb_order": "Order"}),
    ("fig1_tree/fig1_annotation.tsv", "\t", "tip", {r: r for r in RANKS}),
    ("Table_S6_skani_novelty.tsv", "\t", "MAG", {"mag_gtdb_taxonomy": "lineage"}),
]


def expected(lin_map, b, rank, conv):
    if rank == "lineage":
        return None
    v = lin_map[b].get(rank, "")
    if conv == "strip":
        v = strip(v)
    return v


def classify_conv(head_rows, bcol, col, rank):
    """Which convention reproduces HEAD against R220: raw or strip."""
    c = Counter()
    for row in head_rows:
        b = row.get(bcol)
        if b not in L220:
            continue
        cell = (row.get(col) or "").strip()
        raw = L220[b].get(rank, "")
        if cell == raw:
            c["raw"] += 1
        if cell == strip(raw):
            c["strip"] += 1
        c["n"] += 1
    if c["n"] == 0:
        return "raw", c
    return ("raw" if c["raw"] >= c["strip"] else "strip"), c


report = []
problems = 0
for path, sep, bcol, cols in SPECS:
    try:
        cur_text = open(f"{REPO}/{path}", encoding="utf-8").read()
    except FileNotFoundError:
        report.append(f"MISSING {path}")
        problems += 1
        continue
    cur = read_table(cur_text, sep)
    ht = head_text(path)
    head = read_table(ht, sep) if ht else []
    for col, rank in cols.items():
        if cur and col not in cur[0]:
            report.append(f"NOCOL {path} {col}")
            continue
        if rank == "lineage":
            bad = []
            n = 0
            for row in cur:
                b = row.get(bcol)
                if b not in r232:
                    continue
                n += 1
                if row[col].strip() != r232[b]["classification"]:
                    bad.append((b, row[col][:80]))
            report.append(f"{'OK ' if not bad else 'BAD'} {path} [{col}] lineage rows={n} mismatches={len(bad)} {bad[:3]}")
            problems += bool(bad)
            continue
        conv, hc = classify_conv(head, bcol, col, rank)
        n = 0
        bad = Counter()
        bad_ex = {}
        head_ok_bins = set()
        for row in head:
            b = row.get(bcol)
            if b in L220 and (row.get(col) or "").strip() == expected(L220, b, rank, conv):
                head_ok_bins.add(b)
        bad_headok = Counter()
        blank = 0
        for row in cur:
            b = row.get(bcol)
            if b not in L232:
                continue
            cell = (row.get(col) or "").strip()
            n += 1
            exp = expected(L232, b, rank, conv)
            if cell != exp:
                bad[b] += 1
                bad_ex[b] = (cell, exp)
                if b in head_ok_bins:
                    bad_headok[b] += 1
        tag = "OK " if not bad_headok else "BAD"
        problems += bool(bad_headok)
        report.append(
            f"{tag} {path} [{col}->{rank}] conv={conv} HEADvsR220 raw={hc['raw']}/{hc['n']} strip={hc['strip']}/{hc['n']} "
            f"| rows={n} mismatch_bins={len(bad)} (of which HEAD matched R220: {len(bad_headok)}) "
            f"ex={[ (b,)+bad_ex[b] for b in list(bad_headok)[:4]] or [(b,)+bad_ex[b] for b in list(bad)[:3]]}")

# JSON files
mas = json.load(open(f"{REPO}/mag_abundance_summary.json"))
bad = [(b, r) for b in mas for r in ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]
       if (mas[b]["taxonomy"].get(r) or "") != L232[b].get(r, "")]
report.append(f"{'OK ' if not bad else 'BAD'} mag_abundance_summary.json taxonomy mismatches={len(bad)} {bad[:5]}")
problems += bool(bad)
mt = json.load(open(f"{REPO}/leftover_reads_analysis/mag_taxonomy.json"))
bad = [(b, r) for b in mt for r in ["phylum", "class", "order", "family", "genus"]
       if (mt[b].get(r) or "") != L232[b].get(r.capitalize(), "")]
report.append(f"{'OK ' if not bad else 'BAD'} leftover_reads_analysis/mag_taxonomy.json mismatches={len(bad)} bins={len(mt)} {bad[:5]}")
problems += bool(bad)
mg = json.load(open(f"{REPO}/leftover_reads_analysis/mag_genera.json"))
exp_genera = sorted({L232[b].get("Genus") for b in mt if L232[b].get("Genus")})
got = sorted(mg) if isinstance(mg, list) else mg
report.append(f"mag_genera.json type={type(mg).__name__} len={len(mg)}; expected-from-R232 unique genera={len(exp_genera)}; "
              f"equal_as_set={set(mg) == set(exp_genera) if isinstance(mg, list) else 'n/a'}")
print("\n".join(report))
print("PROBLEMS", problems)
