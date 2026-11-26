"""Map the x-axis labels of main-DOCX Figure 4 (Emily's 5-class heatmap, image4.png; labels transcribed from the
image) to MAG bins through the map eras, then to the R232 ID.  Era priority: the submission workbook S7 at HEAD
(6a04588 map, the workbook Emily plots from), then 408cb39/HEAD, ffe668c, 924fa18."""
import csv, json, subprocess, io, os
import openpyxl
R = os.path.expanduser("~/Documents/EmilyKin")
P = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan"
LAB = {
 "Complete": "Ca_Accumulibacter.6 Thauera.34 Ancalomicrobiaceae.1_m Ca_Accumulibacter.8 Ca_Competibacter.49 Accumulibacter.2_m Ca_Competibacter.26 alphaI_cluster.1 Giesbergeria.1_m Azospira.3 Gammaproteobacteria.1_m Zoogloea.11 Ca_Competibacter.32",
 "N2O sink": "Ca_Accumulibacter.23 Azonexus.12 Terrimonas.22 Ca_Leptovillus.11 Thermoanaerobaculia.1_m Ignavibacteria.21 Desulfobacillus.1_m Terrimonas.3 Azonexus.1_m Ca_Sarcinithrix.2 Accumulibacter.1_m Anaerolineae.2 Flavobacteriales.7 Rhodothermaceae.1 Thermoanaerobaculia.2_m",
 "N2O producer": "Ca_Competibacter.57 Lysobacter.2 Pseudoxanthomonas.23 Bdellovibrio.13",
 "Upstream (NO3/NO2)": "Competibacter.1_m Pseudoxanthomonas_A.1_m Rhodospirillales.3 Phototrophicaceae.4_m Ignavibacteria.5 Competibacter_A.1_m Brevundimonas.4",
 "Non-denitrifier": "Ca_Competibacter.48 Gracilibacteria.30 Gracilibacteria.27 Micavibrionaceae.1_m Ca_Competibacter.109 Altimarinus.1_m Dokdonella.14 Ca_Falkowbacteria.7",
}
ws = openpyxl.load_workbook(f"{R}/MAG_target_gene_filled_UPDATED (1).pre_r232.bak.xlsx", read_only=True)["Table S7_functional"]
rows = list(ws.iter_rows(values_only=True))
hdr = rows[0]
s7 = {}
cls = {}
for r in rows[1:]:
    if r[0]:
        s7.setdefault(r[1], []).append(r[0]); cls[r[0]] = r[hdr.index("5 class denitrification ")] if "5 class denitrification " in hdr else None
eras = {}
for e in ["408cb39", "ffe668c", "924fa18"]:
    m = {}
    for r in csv.DictReader(open(f"{P}/docs/work/mag_summary_{e}.tsv"), delimiter="\t"):
        m.setdefault(r["mag_iterativeID"], []).append(r["MAG"])
    eras[e] = m
r232 = json.load(open(f"{R}/mag_iterativeID_old_to_new.json"))
gt = {r["user_genome"]: r["classification"] for r in csv.DictReader(open(f"{R}/gtdbtk.bac120.summary.tsv", newline=""), delimiter="\t")}
out = []
for c, labs in LAB.items():
    for lab in labs.split():
        bins = s7.get(lab); era = "workbook(6a04588)"
        if not bins:
            for e in ["408cb39", "ffe668c", "924fa18"]:
                if eras[e].get(lab):
                    bins = eras[e][lab]; era = e; break
        b = bins[0] if bins and len(bins) == 1 else ("AMBIGUOUS:" + ",".join(bins) if bins else "NOT FOUND")
        new = r232.get(b, "")
        genus = gt[b].split(";")[5][3:] if b in gt else ""
        out.append([c, lab, b, era, cls.get(b, ""), new, "changes" if new and new != lab else "same", genus])
w = csv.writer(open("fig4_label_map.tsv", "w", newline=""), delimiter="\t", lineterminator="\n")
w.writerow(["class_panel", "label_in_figure", "bin", "resolved_via", "workbook_5class", "R232_ID", "status", "R232_genus"])
w.writerows(out)
for o in out: print("\t".join(map(str, o)))
print(sum(o[6] == "changes" for o in out), "of", len(out), "labels change")
