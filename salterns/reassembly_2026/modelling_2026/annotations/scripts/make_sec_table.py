#!/usr/bin/env python
"""Sec read-through rescue rows (KOfam tracks) for the 57 modelling representatives, with genomic ORF coordinates.

The per-MAG KOfam tables (annotations/kofam/<mag_id>.tsv.gz) are the pre-rescue Prodigal calls. The rescue
(annotation_runs/kofam*/scripts/sec_rescue.py) joins Prodigal pieces through an in-frame TGA (U) and re-scores the joined
protein; accepted junctions have status 'sec_readthrough' (new call) or 'confirmed' (a piece already had that KO).
Their KO (best_ko_all) replaces the calls of the Prodigal genes listed in `pieces`.

ORF coordinates (1-based, inclusive, stop codon included unless three_prime_partial):
  + strand: start = tga_lo - 3*ext_upstream_aa ; end = tga_hi + 3*ext_downstream_aa (+3 for the stop)
  - strand: end = tga_hi + 3*ext_upstream_aa ; start = tga_lo - 3*ext_downstream_aa (-3 for the stop)
Checked: when the upstream piece was called, its outer coordinate equals the ORF's; the ORF spans every piece.

usage: make_sec_table.py <reassembly_2026 dir> <annotations dir>
"""
import sys
from pathlib import Path

import pandas as pd

R, A = Path(sys.argv[1]), Path(sys.argv[2])
mags = pd.read_csv(A / "scripts/mags57.tsv", sep="\t")
src = {
    "rollout": R / "annotation_runs/kofam_rollout/results/sec_readthrough/sec_readthrough_junctions_scored.tsv",
    "coassembly": R / "annotation_runs/kofam/results/sec_readthrough/sec_readthrough_junctions_scored.tsv",
}
rows = []
for st, p in src.items():
    j = pd.read_csv(p, sep="\t")
    sub = mags[mags.set == st]
    m = j[j.mag.isin(sub.kofam_name)].merge(sub[["mag_id", "kofam_name"]], left_on="mag", right_on="kofam_name")
    m["source_table"] = str(p.relative_to(R.parent))
    rows.append(m)
d = pd.concat(rows, ignore_index=True)

stop = (~d.three_prime_partial.astype(bool)).astype(int) * 3
plus = d.strand == "+"
d["orf_start"] = (d.tga_lo - 3 * d.ext_upstream_aa).where(plus, d.tga_lo - 3 * d.ext_downstream_aa - stop)
d["orf_end"] = (d.tga_hi + 3 * d.ext_downstream_aa + stop).where(plus, d.tga_hi + 3 * d.ext_upstream_aa)
assert ((d.orf_end - d.orf_start + 1) == 3 * (d.joined_len + (stop // 3))).all()

# check against the Prodigal piece coordinates in the per-MAG KOfam tables
n_chk = 0
for mid, g in d.groupby("mag_id"):
    k = pd.read_csv(A / f"kofam/{mid}.tsv.gz", sep="\t", usecols=["gene", "contig", "start", "end", "strand"]).set_index("gene")
    for _, r in g.iterrows():
        for pc in str(r.pieces).split(";"):
            if pc in ("", "nan"):
                continue
            s, e = k.loc[pc, "start"], k.loc[pc, "end"]
            assert r.orf_start <= s and e <= r.orf_end, (mid, r.junction_id, pc, s, e, r.orf_start, r.orf_end)
            n_chk += 1
        if r.up_called and isinstance(r.up_gene, str) and r.up_gene:
            s, e = k.loc[r.up_gene, "start"], k.loc[r.up_gene, "end"]
            assert (r.orf_start == s) if r.strand == "+" else (r.orf_end == e), (mid, r.junction_id)

d["accepted"] = d.status.isin(["sec_readthrough", "confirmed"])
d["ko"] = d.best_ko_all.where(d.accepted, "")
cols = ["junction_id", "mag_id", "kofam_name", "contig", "contig_len", "strand", "orf_start", "orf_end", "tga_lo", "tga_hi",
        "u_pos", "joined_len", "three_prime_partial", "up_gene", "down_gene", "pieces", "pieces_best_ko", "status",
        "accepted", "ko", "called_marker", "joined_score_called", "piece_max_score_called", "best_ko_all",
        "best_ko_all_score", "best_ko_all_threshold", "secis_model", "secis_score", "secis_evalue", "source_table"]
d = d[cols].sort_values(["mag_id", "contig", "orf_start"])
d.to_csv(A / "kofam/sec_readthrough_57.tsv", sep="\t", index=False)
print(f"{len(d)} junctions in {d.mag_id.nunique()} MAGs; {n_chk} piece coordinates checked")
print(d.groupby("status").size().to_string())
print(d[d.accepted].groupby(["called_marker"]).size().to_string())
