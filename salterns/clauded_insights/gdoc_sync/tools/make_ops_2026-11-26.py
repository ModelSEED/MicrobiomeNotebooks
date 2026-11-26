#!/usr/bin/env python
"""Text suggestions for the Google Doc, 2026-11-26: the metagenome sections on the 2026 data (read survey v2, the
138-KO panel, Figures 1, S4 and S6, Tables S1 and S4) and the small v2 shifts in text suggested on 2026-11-25.
Writes ops for text_delta.py. Long new texts are taken from the review docx (its tracked insertions, with their
italics and subscripts as markup) so that the Doc and the review carry the same words; short edits are listed here.

  make_ops_2026-11-26.py <review.docx> <read_doc.json> <out ops.json>
Already in the Doc and therefore not here: the Figure S6 caption and image (2026-11-25), P33 in Andrew's wording,
and P96 '80%' -> '79%' (the 2026-11-26 test edit).
"""
import json
import re
import sys
from pathlib import Path

from docx import Document

sys.path.insert(0, str(Path(__file__).resolve().parent))
import review_to_suggestions as R  # noqa: E402
import text_delta as T  # noqa: E402


def markup(chars):
    """(char, state, fmt) -> markup string: *italic*, ^sup^, ~sub~."""
    out, cur, buf = [], None, ""
    def flush():
        nonlocal buf
        if buf:
            w = {"i": "*", "superscript": "^", "subscript": "~"}.get(cur, "")
            out.append(f"{w}{buf}{w}")
        buf = ""
    for ch, _, f in chars:
        k = "i" if f.get("i") else (f.get("va") or None)
        if k != cur:
            flush(); cur = k
        buf += ch
    flush()
    return "".join(out)


def review_inserted(review, prefix, whole_paragraph=False):
    """Markup of the Claude-inserted text of the review paragraph whose inserted text starts with prefix."""
    hits = []
    for p in Document(review).element.body.iter(R.W + "p"):
        chars, _ = R.tokens(p)
        ins = [c for c in chars if c[1] == "cins"]
        txt = "".join(c[0] for c in ins).lstrip("\t")
        if txt.startswith(prefix):
            hits.append(ins)
    assert len(hits) == 1, (prefix, len(hits))
    return markup(hits[0])


def main(review, read_p, out):
    doc = T.load(read_p)
    paras = T.paragraphs(doc)
    ops = []

    def op(key, anchor, kind, new, old=None, expect=None):
        o = {"key": key, "anchor": anchor, "kind": kind, "new": new}
        if old is not None:
            o["old"] = old
        if expect:
            o["expect"] = expect
        ops.append(o)

    # Methods
    op("p28_amplicons", "DNA was extracted and 16S rRNA gene amplicons", "replace", "amplicons were processed",
       old="amplicons and metagenomes were processed")
    op("p28_metagenomes", "DNA was extracted and 16S rRNA gene amplicons", "insert_before",
       "The 39 metagenomes, 26 of sediment and 13 of water (Table S4), were re-assembled from the raw reads as "
       "described below. ", old="Sediments were sent to the UC Davis")
    op("p30_panel", "The 24 metagenomes were quality filtered and trimmed with BBtools", "replace_all",
       review_inserted(review, "To build on the findings of Zhou et al. (2022), who reported differences in the"))
    m6_doc = T.find_paragraph(paras, "To count the guild genes independently of assembly and binning")
    m6_old = T.proposed(m6_doc)[0].rstrip("\n")
    op("m6_read_survey", "To count the guild genes independently of assembly and binning", "replace",
       review_inserted(review, "To count the methylated-compound and guild genes"), old=m6_old)
    op("m7_assembly", "To count the guild genes independently of assembly and binning", "append_paragraph",
       review_inserted(review, "The panel was also counted in the assembled genes."))
    edits = [e for e in json.load(open(Path(__file__).resolve().parents[2] / "r232_review_edits_2026-11-26.json"))["edits"]
             if e.get("origin") == "metagenome_2026_11_26"]
    byold = {(e["paragraph"], e["old"]): e for e in edits}
    p50 = "Metabolomic data were analyzed using a combination of R"
    e = byold[(50, "Differences in gene abundances among the three habitat types were tested with Kruskal-Wallis tests.")]
    op("p50_kw", p50, "replace", e["new"], old=e["old"])
    e = byold[(50, "Only metabolites present in at least three samples")]
    op("p50_s4", p50, "insert_before", e["new"], old=e["old"])
    e = byold[(50, "Correlations between Bray-Curtis dissimilarities were calculated with Mantel tests.")]
    op("p50_mantel", p50, "replace", e["new"], old=e["old"])
    # Results
    op("p57_figure1", "were present in at least one metagenome (Table S1)", "replace_all",
       review_inserted(review, "Of the 108 methylated-compound KOs in the panel"))
    rs = "Counting the guild genes in all reads"
    op("rs_methanogenium", rs, "replace", "(6.9%)", old="(6.5%)")
    op("rs_methanococcoides", rs, "replace", "(6.1%)", old="(5.9%)")
    op("rs_roseovarius", rs, "replace", "(39%)", old="(37%)")
    op("rs_markers", rs, "replace", "as did 42%", old="as did 43%")
    p62 = "were positively correlated with CH4 flux while glucose"
    for key, oldtxt in (("p62_figS4", "TMA, trehalose, acetate, and MMA were significantly positively correlated with "
                                      "mcrABG while sucrose, 2-oxoisocaproate, fructose, glucose, and proline were the "
                                      "most negatively correlated with mcrA (Figure S4)."),
                        ("p62_mantel", "(Mantel test, r = 0.78, p < 0.001)")):
        e = byold[(62, oldtxt)]
        op(key, p62, "replace", e["new"], old=e["old"])
    e = [x for x in edits if x["paragraph"] == 67][0]
    op("p67_figS6", "Several aerobic and anaerobic taxa from the hypersaline unrestored", "replace", e["new"], old=e["old"])
    syn = "Figure 10 sets the hypothesized methane sources side by side"
    op("p78_grdI_genus", syn, "replace", "DATFJG01", old="Halanaerobium")
    op("p78_tma_range", syn, "replace", "2–10)", old="2–9)")
    # Discussion
    e = [x for x in edits if x["paragraph"] == 81][0]
    op("p81_bc", "Despite taxonomic and metabolomic dissimilarity being positively correlated", "replace", e["new"],
       old=e["old"])
    e = [x for x in edits if x["paragraph"] == 93][0]
    op("p93_genera", "Several aerobic bacteria were isolated from the hypersaline salterns", "replace", e["new"],
       old=e["old"])
    op("p99_mpn_share", "Taken together, these results attribute", "replace", "1–7%", old="1–8%")
    op("p102_folders", "MAG quality assessments, GTDB R232 classifications", "replace",
       "salterns/correlational/final585); the read survey and the gene panel are in "
       "salterns/reassembly_2026/read_markers_v2 and salterns/reassembly_2026/community_ko, the re-assemblies and bins "
       "in KBase narrative 275645, and salterns/FIGURES.md lists the script behind each figure and table.",
       old="salterns/correlational/final585).")
    # captions
    op("cap_fig1", "Relative abundances of genes (z-score transformed DESeq2", "replace_all",
       review_inserted(review, "Figure 1. Relative abundance of genes involved in the cycling"))
    op("cap_figS4", "Heatmap of Spearman correlation coefficients between metabolite", "replace_all",
       review_inserted(review, "Figure S4. Spearman correlations between metabolite concentrations"))
    op("cap_tableS1", "List of KOs involved in methylated compound cycling", "replace_all",
       review_inserted(review, "Table S1. The 138 KEGG orthologs"))
    op("cap_tableS4", "Metadata, including IMG/M and NCBI accessions, for the 24 metagenomes", "replace_all",
       review_inserted(review, "Table S4. The 39 metagenomes of this study"))
    # every op must match the review's words
    for o in ops:
        t = T.parse(o["new"])[0]
        assert "*" not in t and "~" not in t, o["key"]
    json.dump(ops, open(out, "w"), indent=1, ensure_ascii=False)
    print(f"{len(ops)} ops -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
