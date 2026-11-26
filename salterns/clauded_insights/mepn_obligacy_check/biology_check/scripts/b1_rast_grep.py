"""B1. RAST function strings (2020 RASTified MAG JSON) for the roles that matter, across all 24 MePn-consumer MAGs.
Synonym-tolerant regexes, so a '0 hits' is not an artefact of wording."""
import sys
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import re, json
from common_bio import *

PAT = {
 "HSK_any (homoserine kinase / 2.7.1.39 / ThrB / ThrH)": r"homoserine kinase|2\.7\.1\.39|\bthrB\b|\bthrH\b|homoserine phosphotransferase",
 "SerB_PSP (phosphoserine phosphatase 3.1.3.3)": r"phosphoserine phosphatase|3\.1\.3\.3\b",
 "Hom (homoserine dehydrogenase)": r"homoserine dehydrogenase|1\.1\.1\.3\b",
 "ThrC (threonine synthase)": r"threonine synthase|4\.2\.3\.1\b",
 "LtaE (threonine aldolase)": r"threonine aldolase|4\.1\.2\.48|4\.1\.2\.5\b",
 "Tdh/Kbl (Thr dehydrogenase / 2-amino-3-ketobutyrate CoA ligase)": r"threonine 3-dehydrogenase|threonine dehydrogenase|1\.1\.1\.103|2-amino-3-ketobutyrate|glycine C-acetyltransferase|2\.3\.1\.29",
 "GlyA (SHMT)": r"serine hydroxymethyltransferase|2\.1\.2\.1\b",
 "GCS (glycine cleavage T/H/P/L)": r"glycine cleavage|glycine dehydrogenase \(decarboxylating|aminomethyltransferase|1\.4\.4\.2|2\.1\.2\.10",
 "SGAT/AGT (serine/alanine-glyoxylate aminotransferase)": r"glyoxylate aminotransferase|2\.6\.1\.45|2\.6\.1\.44|glyoxylate transaminase",
 "ICL/MS (isocitrate lyase / malate synthase)": r"isocitrate lyase|4\.1\.3\.1\b|malate synthase|2\.3\.3\.9",
 "EMC (ethylmalonyl-CoA pathway)": r"crotonyl-CoA (carboxylase|reductase)|ethylmalonyl-CoA|methylsuccinyl-CoA|mesaconyl|methylfumaryl|malyl-CoA lyase|beta-methylmalyl|malyl-CoA thioesterase|4\.1\.3\.24|1\.3\.1\.85|5\.4\.99\.63",
 "PurU (10-formyl-THF deformylase)": r"formyltetrahydrofolate deformylase|3\.5\.1\.10\b|\bpurU\b",
 "Fhs (formate--THF ligase)": r"formate--tetrahydrofolate ligase|formyltetrahydrofolate synthetase|6\.3\.4\.3",
 "FolD / MTHFD": r"methylenetetrahydrofolate dehydrogenase|methenyltetrahydrofolate cyclohydrolase|1\.5\.1\.5|3\.5\.4\.9",
 "MetF (CH2-THF reductase)": r"methylenetetrahydrofolate reductase|1\.5\.1\.20|1\.5\.1\.54",
 "FDH (formate dehydrogenase)": r"formate dehydrogenase",
 "Formaldehyde (FrmA/FrmB/Gfa/FDH-GSH-indep)": r"S-\(hydroxymethyl\)glutathione|S-formylglutathione|glutathione-dependent formaldehyde|formaldehyde dehydrogenase|formaldehyde-activating",
 "Hut formyl outlets (HutG/HutF/FTCD)": r"N-formylglutamate|formiminoglutamase|formimidoylglutamase|formimidoylglutamate|formiminotransferase|formiminoglutamate",
 "Prs (ribose-phosphate pyrophosphokinase)": r"ribose-phosphate pyrophosphokinase|phosphoribosyl ?pyrophosphate synth|2\.7\.6\.1\b|PRPP synth",
 "PhnN (ribose 1,5-bisphosphokinase)": r"phnN|ribose 1,5-bisphosph|2\.7\.4\.23",
 "PhnP (PRcP phosphodiesterase)": r"phnP|1,2-cyclic phosph|3\.1\.4\.55",
 "C-P lyase core PhnGHIJKLM": r"\bPhn[GHIJKLM]\b|C-P lyase|methylphosphonate 5-(tri)?phosphate",
 "Amn (AMP nucleosidase)": r"AMP nucleosidase|3\.2\.2\.4\b",
 "PNP / DeoD / MTAP": r"purine nucleoside phosphorylase|purine-nucleoside phosphorylase|2\.4\.2\.1\b|methylthioadenosine phosphorylase|2\.4\.2\.28",
 "DeoB (phosphopentomutase)": r"phosphopentomutase|5\.4\.2\.7",
 "Adenosine kinase": r"adenosine kinase|2\.7\.1\.20\b",
 "APRT (adenine phosphoribosyltransferase)": r"adenine phosphoribosyltransferase|2\.4\.2\.7\b",
 "Pfs/MtnN (MTA/SAH nucleosidase, 5'-dA)": r"methylthioadenosine nucleosidase|adenosylhomocysteine nucleosidase|3\.2\.2\.9\b|3\.2\.2\.16|5'-deoxyadenosine|MTA/SAH",
}

g = guilds()
rows = []
out = {}
for mid in model_ids():
    gj = rast(mid)
    fns = [(short(f["id"]), fn) for f in gj["features"] for fn in (f.get("functions") or [])]
    rec = {}
    for lab, rx in PAT.items():
        r = re.compile(rx, re.I)
        rec[lab] = sorted({(a, b) for a, b in fns if r.search(b)})
    out[iid(mid)] = {"model": mid, "completeness": g[mag_of(mid)]["completeness"],
                     "contamination": g[mag_of(mid)]["contamination"], "n_features": len(gj["features"]),
                     "hits": {k: [f"{a}: {b}" for a, b in v] for k, v in rec.items()}}

json.dump(out, open(VB / "b1_rast_grep.json", "w"), indent=1)
# compact count table
labs = list(PAT)
abbrev = ["HSK", "SerB", "Hom", "ThrC", "LtaE", "Tdh/Kbl", "GlyA", "GCS", "SGAT", "ICL/MS", "EMC", "PurU", "Fhs",
          "FolD", "MetF", "FDH", "HCHO", "Hut", "Prs", "PhnN", "PhnP", "CPlyase", "Amn", "PNP", "DeoB", "AdoK", "APRT", "Pfs"]
print("iterativeID".ljust(22), "compl  cont ", " ".join(a.rjust(5) for a in abbrev))
for k in sorted(out, key=lambda x: (x not in OBLIGATE, x)):
    v = out[k]
    print(("*" if k in OBLIGATE else " ") + k.ljust(21), f"{float(v['completeness']):5.1f} {float(v['contamination']):5.1f}",
          " ".join(str(len(v["hits"][l])).rjust(5) for l in labs))
print()
for k in ["Roseovarius.9", "Roseovarius.11", "Albimonas.1"]:
    print("=====", k, out[k]["model"], "completeness", out[k]["completeness"], "contamination", out[k]["contamination"])
    for l in labs:
        print("  ", l)
        for h in out[k]["hits"][l]:
            print("       ", h)
