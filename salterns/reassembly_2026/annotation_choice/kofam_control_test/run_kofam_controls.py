#!/usr/bin/env python
"""KOfam arm of the DRAM-vs-Bakta test, run locally on control MAGs (no DRAM, no Bakta installed here).

What it reproduces: DRAM 1.x's KOfam layer, i.e. prodigal -p meta gene calls on contigs >= min_contig,
`hmmsearch --domtblout` of KOfam profiles, and DRAM's `sig_scores` rule (score > KO threshold, using the
domain or full-sequence score as `score_type` says; KOs whose threshold is '-' are never called), then the
best KO per gene by full-sequence E-value. Only the ~67 marker KOs of this paper are searched
(../marker_ko_table.tsv), so a gene's "best KO" is best among the markers.

Truth: for the 2020 (original) control MAGs, the RAST genome JSONs that the models and guilds_r232.csv used.
Prodigal genes are matched to RAST CDSs by contig, strand and >= 50 % overlap of the shorter feature.

Inputs (read-only):
  salterns/MAG_contigs/<mag>.contigs__.fna                     original MAG contigs
  salterns/jplfaria/Cliff MAGs/RASTified_MAG_json/<mag>.contigs__.RAST.json
  reassembly_2026/work/checkm/<unit>/fasta/<bin>.fna            re-assembly MAGs
  KOFAM_DIR/markers.hmm (+ hmmpress), KOFAM_DIR/ko_list          KOfam release 2026-05/06 profiles + matching ko_list

Outputs (this folder): genes_vs_rast.tsv, marker_calls_by_mag.tsv, summary.txt
Run: ~/Documents/py_venv/bin/python run_kofam_controls.py <KOFAM_DIR> <WORK_DIR>
"""
import json
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
S = HERE.parents[2]                                   # salterns/
R = S / "reassembly_2026"
KOFAM, WORK = Path(sys.argv[1]), Path(sys.argv[2])
WORK.mkdir(parents=True, exist_ok=True)
BIN = Path.home() / "Documents/gtdbtk_venv/bin"        # HMMER 3.4, Prodigal 2.6.3
RAST = S / "jplfaria/Cliff MAGs/RASTified_MAG_json"

ORIG = {   # label -> 2020 MAG name (guilds_r232.csv `mag`)
    "Roseovarius.8": "Salt_Pond_MetaG_R2_B_D2_MG_DASTool_bins_metabat.38",
    "Albimonas.1": "Salt_Pond_MetaG_R2_restored_C_black_MG_DASTool_bins_metabat.8",
    "Marivita.2": "Salt_Pond_MetaG_R1_B_D2_MG_DASTool_bins_metabat.28",
    "Halomonas.1": "Salt_Pond_MetaG_R1_A_D1_MG_DASTool_bins.metabat.26",
    "Methanolobus.1": "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_concoct_out.21",
    "Halanaerobium.1": "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.45",
    "M55B127.1": "Salt_Pond_MetaG_R1_C_D2_MG_DASTool_bins_metabat.40",
    "DATFJG01.1": "Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_metabat.27",
    "Desulfobacter.2": None,          # filled from guilds_r232.csv below
    "Haloquadratum.4": "Salt_Pond_MetaG_R2_A_D2_MG_DASTool_bins_metabat.16",
    "Natrinema.1": "Salt_Pond_MetaG_R1_C_D1_MG_DASTool_bins_metabat.31",
    "J153.1": "Salt_Pond_MetaG_R2A_C_D1_MG_DASTool_bins_metabat.5",
}
NEW = {    # label -> re-assembly bin (MaxBin2 Roseovarius bins: contig-length test; best new Methanolobus)
    "new_Roseovarius_maxbin2_R1_A_D2": "maxbin2.Bin.048_sub.fasta__R1_A_D2",
    "new_Roseovarius_maxbin2_R1_B_D2": "maxbin2.Bin.032_sub.fasta__R1_B_D2",
    "new_Methanolobus_concoct_R1_sed": "concoct.bin.101.fasta__R1_sediment",
}
guilds = pd.read_csv(S / "correlational/r232/guilds_r232.csv")
ORIG["Desulfobacter.2"] = guilds.set_index("iterativeID_r232").loc["Desulfobacter.2", "mag"]
# every RAST-based guild member (22 MePn consumers, 9 betaine reducers, 1 methanogen, 4 TMA candidates)
for _, g in guilds[guilds.mepn_consumer | guilds.betaine_reducer | guilds.methanogen | guilds.tma_producer].iterrows():
    ORIG.setdefault(g.iterativeID_r232, g.mag)
# all re-assembly bins in Albimonas.1's species cluster (is its broken C-P lyase locus an assembly artefact?)
for b in ["concoct.bin.061.fasta__R1_C_D1", "concoct.bin.157_sub.fasta__R1_sediment", "concoct.bin.018.fasta__R2_B_D2",
          "metabat2.bin.062.fasta__R2_restored_C_black", "concoct.bin.058_sub.fasta__R2_sediment"]:
    NEW[f"new_Albimonas_{b.split('.fasta__')[0].replace('.', '')}_{b.split('__')[1]}"] = b
RAST_PATTERNS = re.compile(
    r"phosphon|C-P lyase|Phn[A-Z]\b|betaine|glycine reductase|sarcosine reductase|reductase (component|protein) [ABC]|"
    r"trimethylamine|dimethylamine|methylamine|methanol|methyl.coenzyme M|pyrrolysine|choline|carnitine|N-oxide|"
    r"Methylcobalamin|corrinoid", re.I)


def read_fasta(p):
    seqs, name = {}, None
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def unit_fasta(bin_name):
    hits = list((R / "work/checkm").glob(f"*/fasta/{bin_name}.fna"))
    assert len(hits) == 1, (bin_name, hits)
    return hits[0]


def call_genes(label, fasta, min_contig):
    """prodigal -p meta on contigs >= min_contig (DRAM default 2,500); protein ids <label>|<contig>_<n>."""
    out = WORK / f"{label}.min{min_contig}"
    faa, gff = Path(f"{out}.faa"), Path(f"{out}.gff")
    if faa.exists() and faa.stat().st_size:
        return faa, gff
    seqs = {k: v for k, v in read_fasta(fasta).items() if len(v) >= min_contig}
    flt = Path(f"{out}.fna")
    flt.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    subprocess.run(["nice", "-n", "10", str(BIN / "prodigal"), "-p", "meta", "-g", "11", "-q", "-i", str(flt),
                    "-a", f"{out}.raw.faa", "-f", "gff", "-o", str(gff)], check=True)
    with open(f"{out}.raw.faa") as fi, open(faa, "w") as fo:
        for line in fi:
            fo.write(f">{label}|{line[1:].split()[0]}\n" if line.startswith(">") else line)
    return faa, gff


def gene_coords(gff, label):
    co = {}
    for line in open(gff):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        n = re.search(r"ID=\d+_(\d+)", f[8])[1]
        co[f"{label}|{f[0]}_{n}"] = (f[0], int(f[3]), int(f[4]), f[6], "partial=00" not in f[8])
    return co


def rast_cds(mag):
    d = json.load(open(RAST / f"{mag}.contigs__.RAST.json"))
    out = []
    for x in d["features"]:
        loc = x["location"][0]
        contig, start, strand, length = loc[0], int(loc[1]), loc[2], int(loc[3])
        b, e = (start, start + length - 1) if strand == "+" else (start - length + 1, start)
        fn = "; ".join(x.get("functions") or [x.get("function", "")])
        out.append((x["id"].split(".")[-1], contig, b, e, strand, fn, len(x.get("protein_translation", ""))))
    return out


def main():
    ko_list = pd.read_csv(KOFAM / "ko_list", sep="\t", index_col=0)
    markers = pd.read_csv(HERE.parent / "marker_ko_table.tsv", sep="\t", index_col=0)
    hmm_len = {}
    for line in open(KOFAM / "markers.hmm"):
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.startswith("LENG"):
            hmm_len[name] = int(line.split()[1])
    jobs = [(lab, S / f"MAG_contigs/{m}.contigs__.fna", 2500) for lab, m in ORIG.items()]
    for lab, b in NEW.items():
        jobs += [(lab, unit_fasta(b), 2500), (lab, unit_fasta(b), 1000)]
    with ThreadPoolExecutor(8) as ex:
        called = list(ex.map(lambda j: (j, call_genes(*j)), jobs))
    allfaa = WORK / "all.faa"
    with open(allfaa, "w") as fo:
        for (lab, _, mc), (faa, _) in called:
            for line in open(faa):
                fo.write(line.replace(f">{lab}|", f">{lab}@{mc}|") if line.startswith(">") else line)
    dom = WORK / "markers.domtblout"
    subprocess.run(["nice", "-n", "10", str(BIN / "hmmsearch"), "--cpu", "8", "--domtblout", str(dom), "-o", "/dev/null",
                    str(KOFAM / "markers.hmm"), str(allfaa)], check=True)
    rows = []
    for line in open(dom):
        if line.startswith("#"):
            continue
        f = line.split()
        rows.append({"gene": f[0], "glen": int(f[2]), "ko": f[3], "full_evalue": float(f[6]), "full_score": float(f[7]),
                     "domain_score": float(f[13]), "hmm_from": int(f[15]), "hmm_to": int(f[16])})
    h = pd.DataFrame(rows)
    h["threshold"] = h.ko.map(ko_list.threshold)
    h["score_type"] = h.ko.map(ko_list.score_type)
    h["score"] = h.full_score.where(h.score_type == "full", h.domain_score)
    h["sig"] = (h.score_type.isin(["full", "domain"])) & (h.score > pd.to_numeric(h.threshold, errors="coerce"))
    # per gene: best significant KO (DRAM: lowest full E-value among significant), and best KO of any score
    best_any = h.sort_values(["gene", "full_evalue"]).drop_duplicates("gene").set_index("gene")
    sig = h[h.sig].sort_values(["gene", "full_evalue"]).drop_duplicates("gene").set_index("gene")
    hmm_cov = (h.groupby(["gene", "ko"]).apply(lambda x: (x.hmm_to.max() - x.hmm_from.min() + 1)) /
               h.groupby(["gene", "ko"]).ko.first().map(hmm_len))

    out = []
    for (lab, fasta, mc), (faa, gff) in called:
        coords = gene_coords(gff, lab)
        rast = rast_cds(ORIG[lab]) if lab in ORIG else []
        by_contig = defaultdict(list)
        for r in rast:
            by_contig[r[1]].append(r)
        seen_rast = set()
        for g, (contig, b, e, strand, partial) in coords.items():
            key = f"{lab}@{mc}|{g.split('|', 1)[1]}"
            match = None
            for r in by_contig.get(contig, []):
                if r[4] != strand:
                    continue
                ov = min(e, r[3]) - max(b, r[2]) + 1
                if ov > 0 and ov >= 0.5 * min(e - b + 1, r[3] - r[2] + 1):
                    match = r
                    break
            rast_fn = match[5] if match else ""
            if match:
                seen_rast.add(match[0])
            s = sig.loc[key] if key in sig.index else None
            a = best_any.loc[key] if key in best_any.index else None
            if s is None and a is None and not RAST_PATTERNS.search(rast_fn):
                continue
            out.append({
                "mag": lab, "min_contig": mc, "gene": key, "contig": contig, "begin": b, "end": e, "strand": strand,
                "partial_gene": partial, "aa_len": (e - b + 1) // 3 - 1,
                "kofam_ko": s.name and s.ko if s is not None else "",
                "kofam_gene": markers.gene.get(s.ko, "") if s is not None else "",
                "kofam_score": round(s.score, 1) if s is not None else None,
                "best_marker_ko_any_score": a.ko if a is not None else "",
                "best_marker_gene_any_score": markers.gene.get(a.ko, "") if a is not None else "",
                "best_score": round(a.score, 1) if a is not None else None,
                "best_threshold": a.threshold if a is not None else None,
                "best_hmm_cov": round(float(hmm_cov.get((key, a.ko), float("nan"))), 2) if a is not None else None,
                "rast_cds": match[0] if match else "", "rast_function": rast_fn,
                "rast_aa_len": match[6] if match else None,
            })
        for r in rast:   # RAST marker-like CDS with no prodigal gene at all
            if r[0] not in seen_rast and RAST_PATTERNS.search(r[5]):
                out.append({"mag": lab, "min_contig": mc, "gene": "", "contig": r[1], "begin": r[2], "end": r[3],
                            "strand": r[4], "rast_cds": r[0], "rast_function": r[5], "rast_aa_len": r[6]})
    t = pd.DataFrame(out)
    t.to_csv(HERE / "genes_vs_rast.tsv", sep="\t", index=False)

    # per MAG x marker: KOfam call count (DRAM rule) and RAST guild-rule evidence
    calls = t[t.kofam_ko.fillna("") != ""].groupby(["mag", "min_contig", "kofam_gene"]).size().unstack(fill_value=0)
    calls.to_csv(HERE / "marker_calls_by_mag.tsv", sep="\t")
    with open(HERE / "summary.txt", "w") as fh:
        fh.write(f"genes searched: {sum(1 for l in open(allfaa) if l.startswith('>'))}; hits {len(h)}; "
                 f"significant (DRAM rule) {int(h.sig.sum())}\n")
        for (lab, _, mc), (faa, _) in called:
            n = sum(1 for l in open(faa) if l.startswith(">"))
            fh.write(f"{lab} min_contig={mc}: {n} genes\n")
    print(open(HERE / "summary.txt").read())


if __name__ == "__main__":
    main()
