#!/usr/bin/env python
"""Table S1 (2026): the methylated-compound gene panel, with read- and assembly-based abundance statistics.

Sheets
  Table S1         one row per panel KO (138): identity and panel status; the read survey (24 sediment metagenomes:
                   reads, copies per 100 genome equivalents per restoration type, Kruskal-Wallis with BH FDR over all
                   detected panel KOs, Dunn letters, and the survey's validation on held-out genes); the assembly-based
                   KOfam profile (same statistics, assembled fraction, Spearman rho with the read survey); and the
                   species of the final MAG set that carry the KO.
  Original list    the 39 KOs of the original Table S1 and what became of each.
  Notes            column definitions.
Inputs: ko_panel_methylated.tsv, table_s1_original_39.tsv, results/panel_stats_{read,assembly}.tsv (panel_stats.py),
the read survey named in ../read_survey_current.json (validation/, meta/merged_v2.json), assembly_based/results/.
Writes ../../Table S1 - methylated-compound gene panel (2026).xlsx (salterns/).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
R26 = os.path.dirname(HERE)
SAL = os.path.dirname(R26)
sys.path.insert(0, HERE)
import profiles as PR  # noqa: E402

OUT = os.path.join(SAL, "Table S1 - methylated-compound gene panel (2026).xlsx")
TYPES = [("unrestored", "U"), ("restored", "R"), ("reference", "Ref")]
STATUS = {"S1": "original Table S1", "S1-relabel": "original Table S1 (compound or process relabeled)",
          "S1-replaced": "replaces a mis-assigned original KO", "new": "added (2026)",
          "FigS4": "Figure S4 set (methanogenesis pathways)"}


def letters(r):
    if pd.isna(r.get("q")) or r.get("q") >= 0.05:
        return ""
    return "/".join(str(r[f"letter_{t}"]) for t, _ in TYPES)


def read_validation():
    cur = json.load(open(os.path.join(R26, "read_survey_current.json")))
    root = os.path.join(R26, os.path.dirname(cur["results"].rstrip("/")))
    vf = os.path.join(root, "validation", "validation_per_label_final.tsv")
    if not os.path.exists(vf):
        return None, {}, cur
    v = pd.read_csv(vf, sep="\t").set_index("label")
    merged = {}
    mf = os.path.join(root, "validation", "validation_merged_classes_v2.tsv")
    if os.path.exists(mf):
        for _, r in pd.read_csv(mf, sep="\t").iterrows():
            for m in str(r["members"]).split(";"):
                merged.setdefault(m, []).append(r["label"])
    return v, merged, cur


def main():
    P = PR.panel().set_index("ko")
    T = pd.DataFrame(index=P.index)
    T["KO"] = P.index
    T["Gene"] = P.gene
    T["KOfam definition"] = P.kofam_definition
    T["Compound"] = P.compound
    T["Compound class (Figure 1)"] = P.compound_class
    T["Process"] = P.process
    T["Type"] = P.ptype
    T["Figure"] = np.where(P.status == "FigS4", "S4", "1")
    T.loc[[k for k in P.index if P.loc[k, "gene"] in {"mttB", "mttC", "mtbB", "mtbC", "mtmB", "mtmC", "mtbA", "mtaB",
                                                         "mtaC", "mtaA", "mtgB", "mtsA", "mcrA", "mcrB", "mcrG", "phnG",
                                                         "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM"}], "Figure"] = "1, S4"
    T["Panel status"] = P.status.map(STATUS)
    T["Specificity"] = P.specificity
    T["Note"] = P.note.fillna("")

    # ---- read survey
    rs = os.path.join(HERE, "results", "panel_stats_read.tsv")
    val, merged, cur = read_validation()
    if os.path.exists(rs):
        r = pd.read_csv(rs, sep="\t").set_index("ko")
        T["Reads: reads (24 sediment)"] = r.n_24
        T["Reads: samples detected (of 24)"] = r.n_samples_detected
        for t, ab in TYPES:
            T[f"Reads: copies per 100 GE, {t}"] = r[f"mean_per100GE_{t}"]
        T["Reads: KW H"], T["Reads: KW p"], T["Reads: KW FDR q"] = r.H, r.p, r.q
        T["Reads: Dunn letters (U/R/Ref)"] = r.apply(letters, axis=1)
        if val is not None:
            lab = P.label
            T["Reads: validation recall"] = [val.recall.get(l, np.nan) for l in lab]
            T["Reads: validation precision"] = [val.precision.get(l, np.nan) for l in lab]
            notes = []
            for k in P.index:
                l = P.loc[k, "label"]
                n = []
                if l not in val.index or pd.isna(val.recall.get(l)) or val.n_reads.get(l, 0) == 0:
                    n.append("not testable (no held-out reference)")
                else:
                    if val.n_reads[l] >= 50 and (val.recall[l] < 0.80 or val.precision[l] < 0.90):
                        n.append("flagged (recall < 0.80 or precision < 0.90)")
                if l in merged:
                    n.append("reported with " + ", ".join(merged[l]))
                notes.append("; ".join(n))
            T["Reads: validation note"] = notes
        T.attrs["read_version"] = cur.get("version")

    # ---- assembly-based
    a = pd.read_csv(os.path.join(HERE, "results", "panel_stats_assembly.tsv"), sep="\t").set_index("ko")
    T["Assembly: genes (24 sediment)"] = a.n_24
    T["Assembly: samples detected (of 24)"] = a.n_samples_detected
    for t, ab in TYPES:
        T[f"Assembly: copies per 100 GE, {t}"] = a[f"mean_per100GE_{t}"]
    T["Assembly: KW FDR q"] = a.q
    T["Assembly: Dunn letters (U/R/Ref)"] = a.apply(letters, axis=1)
    if os.path.exists(rs):                         # concordance with the read survey (24 sediment samples)
        rd = PR.wide(PR.read_profile(), "copies_per_GE", PR.SED, kos=list(P.index))
        asm = PR.assembly_profile()
        aw = asm.pivot(index="ko", columns="sample", values="copies_per_GE_readGE").reindex(index=P.index,
                                                                                           columns=PR.SED)
        rho, frac = [], []
        for k in P.index:
            x, y = rd.loc[k], aw.loc[k]
            ok = x.notna() & y.notna()
            rho.append(spearmanr(x[ok], y[ok])[0] if ok.sum() > 2 and x[ok].nunique() > 1 and y[ok].nunique() > 1
                       else np.nan)
            frac.append(float(y[ok].sum() / x[ok].sum()) if ok.any() and x[ok].sum() > 0 else np.nan)
        T["Assembly vs reads: assembled fraction (24 sediment)"] = frac
        T["Assembly vs reads: Spearman rho (24 sediment)"] = rho

    # ---- MAG carriers (final MAG set, species representatives)
    m = pd.read_csv(os.path.join(HERE, "assembly_based", "results", "panel_mag_attribution_summary.tsv"),
                    sep="\t").set_index("ko")
    T["MAGs: species carrying the KO"] = m.n_species_rep_carries.reindex(T.index).fillna(0).astype(int)
    T["MAGs: archaeal species"] = m.n_archaeal_species_rep.reindex(T.index).fillna(0).astype(int)
    T["MAGs: most abundant carriers (genus, mean fraction of reads in the 12 unrestored samples)"] = \
        m.top5_rep_carriers.reindex(T.index).fillna("")

    # ---- original list
    O = pd.read_csv(os.path.join(HERE, "table_s1_original_39.tsv"), sep="\t")
    def fate(r):
        if r.in_new_panel:
            st = P.loc[r.ko, "status"] if r.ko in P.index else ""
            return "kept" + (" (relabeled: " + str(P.loc[r.ko, "note"]) + ")" if st == "S1-relabel" else "")
        rep = P[P.note.fillna("").str.contains(r.ko)]
        return ("replaced by " + ", ".join(f"{k} {g}" for k, g in zip(rep.index, rep.gene))) if len(rep) else "dropped"
    O["In the 2026 panel"] = O.apply(fate, axis=1)
    O = O.rename(columns={"ko": "KO", "label_in_figure1": "Label in the original Figure 1",
                          "category_in_figure1": "Category in the original Figure 1", "kofam_definition": "KOfam definition"})
    O = O.drop(columns=["in_new_panel"])

    notes = pd.DataFrame({"Column": [
        "Panel status", "GE", "Reads: …", "Reads: KW FDR q", "Reads: Dunn letters", "Reads: validation …",
        "Assembly: …", "Assembly vs reads", "MAGs: …"], "Definition": [
        "Where each KO comes from: the original Table S1 list (39 KOs; six were mis-assigned KOs and are replaced, "
        "see 'Original list'), KOs added in 2026 to cover methylated-compound pathways the original list missed, and "
        "the Figure S4 methanogenesis-pathway set.",
        "Genome equivalent: the mean depth of eight universal single-copy ribosomal-protein genes; copies per GE ≈ "
        "copies of the gene per genome in the community.",
        "Read survey: every quality-filtered read compared with DIAMOND blastx against the panel's reference proteins "
        "and their paralogs (decoys); a read counts when its best hit is a panel reference above the identity cut-off. "
        "Means are over the samples of each restoration type (unrestored R1 + R2, n = 12; restored SF2, n = 6; "
        "reference R2A, n = 6).",
        "Kruskal–Wallis across the three restoration types, Benjamini–Hochberg FDR over all detected panel KOs "
        "(the same family as Figure 1).",
        "Dunn's test (Holm-adjusted) when FDR q < 0.05: restoration types sharing a letter do not differ; 'a' is the "
        "highest mean rank.",
        "Simulated reads from held-out reference genes: recall (reads of the gene recovered) and precision (reads "
        "called as the gene that are truly from it). Genes that cannot be told apart in reads are also reported as a "
        "merged class.",
        "KOfam (hmmsearch with each KO's threshold; best-scoring passing KO per gene) on all 12.06 million genes "
        "assembled from the 39 metagenomes; copies per GE from the contig depths of the genes over the depth of the "
        "assembled ribosomal-protein genes.",
        "Assembled fraction: the gene's assembly-based copies over its read-based copies, both normalized by the "
        "read-based GE; Spearman rho between the two across the 24 sediment samples.",
        "Species of the final MAG set (585 species) whose representative carries the KO (KOfam); the most abundant "
        "carriers are listed by genus with their mean fraction of reads in the 12 unrestored samples."]})

    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        T.to_excel(xw, sheet_name="Table S1", index=False)
        O.to_excel(xw, sheet_name="Original list", index=False)
        notes.to_excel(xw, sheet_name="Notes", index=False)
        from openpyxl.styles import Alignment, Font
        for name, df in (("Table S1", T), ("Original list", O), ("Notes", notes)):
            ws = xw.sheets[name]
            ws.freeze_panes = "C2" if name == "Table S1" else "A2"
            for c in ws[1]:
                c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical="top")
            for i, col in enumerate(df.columns, 1):
                w = min(60, max(10, int(df[col].astype(str).str.len().quantile(0.9)) + 2)) if len(df) else 12
                ws.column_dimensions[ws.cell(1, i).column_letter].width = w if name != "Notes" or i == 1 else 110
                if df[col].dtype.kind == "f":
                    fmt = "0.000" if ("rho" in col or "fraction" in col or "recall" in col or "precision" in col) \
                        else ("0.0E+00" if col.endswith(("p", "q")) else "0.000")
                    for row in ws.iter_rows(min_row=2, min_col=i, max_col=i):
                        row[0].number_format = fmt
            if name == "Notes":
                for row in ws.iter_rows(min_row=2):
                    for c in row:
                        c.alignment = Alignment(wrap_text=True, vertical="top")
    print("wrote", OUT, T.shape, "read survey:", T.attrs.get("read_version", "not included"))


if __name__ == "__main__":
    main()
