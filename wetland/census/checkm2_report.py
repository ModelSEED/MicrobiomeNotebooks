#!/usr/bin/env python3
"""Save the CheckM2 results for all wetland MAGs as a KBase report in narrative 272312 (per Andrew, 2026-11-15:
"run CheckM2 and save those reports").

CheckM2 1.1.0 (DB uniref100.KO.1 v3) was run on poplar (census/run_census_compute.sh) and re-run with --allmodels
--resume (same gene calls and DIAMOND hits) so both completeness models are reported. This script writes, under
<out dir>:
  checkm2_quality_report_allmodels.tsv   CheckM2's own table (Name = MAG name) + iterativeID column
  checkm2_vs_checkm1.tsv                 per MAG: iterativeID, CheckM2 (general, specific, chosen model, contamination,
                                         tier), KBase CheckM v1 (kb_Msuite), DAS Tool SCG, GTDB R232 lineage
and, with --apply, uploads them plus an HTML summary and saves a visible KBaseReport.Report
'CheckM2_all_MAGs_report_<UTC date>' in narrative 272312. MAG objects are not modified.
Usage: KB_AUTH_TOKEN=... python -I -B checkm2_report.py <checkm2 quality_report.tsv> <census outputs dir> <out dir> [--apply]
"""
import html, os, sys, tempfile, time, zipfile
from collections import Counter
import pandas as pd

sys.path.insert(0, "/home/freiburger/gtdbtk_kbase")
sys.path.insert(0, "/home/freiburger/gtdbtk_kbase/checkm_push")
import kbase_mags as k
from writeback_checkm import blob_upload, BLOB_CANON

QR, CD, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
APPLY = "--apply" in sys.argv
os.makedirs(OUT, exist_ok=True)
q = pd.read_csv(QR, sep="\t")
ids = pd.read_csv(os.path.join(CD, "mag_iterativeIDs.tsv"), sep="\t").set_index("mag")
cen = pd.read_csv(os.path.join(CD, "mag_census.tsv"), sep="\t").set_index("mag")
assert len(q) == 12859 and set(q.Name) == set(ids.index), "CheckM2 table must cover exactly the 12,859 MAGs"
need = {"Completeness_General", "Completeness_Specific", "Completeness_Model_Used", "Contamination"}
assert need <= set(q.columns), f"not an --allmodels report: missing {need - set(q.columns)}"
q.insert(1, "iterativeID", q.Name.map(ids.iterativeID))
f1 = os.path.join(OUT, "checkm2_quality_report_allmodels.tsv")
q.to_csv(f1, sep="\t", index=False)


def chosen(r):
    return r.Completeness_Specific if str(r.Completeness_Model_Used).startswith("Neural") else r.Completeness_General


def tier(c, x):
    return "high" if c > 90 and x < 5 else "medium" if c >= 50 and x < 10 else "low"


rows = []
for r in q.itertuples(index=False):
    c = cen.loc[r.Name]
    comp = chosen(r)
    rows.append({"mag": r.Name, "iterativeID": r.iterativeID, "kbase_ref": c.kbase_ref,
                 "checkm2_completeness": comp, "checkm2_completeness_general": r.Completeness_General,
                 "checkm2_completeness_specific": r.Completeness_Specific, "checkm2_model": r.Completeness_Model_Used,
                 "checkm2_contamination": r.Contamination, "checkm2_tier": tier(comp, r.Contamination),
                 "checkm1_completeness": c.checkm1_completeness, "checkm1_contamination": c.checkm1_contamination,
                 "checkm1_tier": tier(c.checkm1_completeness, c.checkm1_contamination),
                 "das_completeness": c.das_completeness, "das_contamination": c.das_contamination,
                 "genome_size": r.Genome_Size, "contig_n50": r.Contig_N50, "r232_classification": c.r232_classification})
cmp_ = pd.DataFrame(rows)
f2 = os.path.join(OUT, "checkm2_vs_checkm1.tsv")
cmp_.to_csv(f2, sep="\t", index=False)
t2, t1 = Counter(cmp_.checkm2_tier), Counter(cmp_.checkm1_tier)
phy = cmp_.r232_classification.str.split(";").str[1].str[3:]
cpr = phy.eq("Patescibacteriota")
dpann = phy.isin({"Nanobdellota", "Aenigmatarchaeota", "Micrarchaeota", "Nanohalarchaeota", "SpSt-1190", "Altiarchaeota"})
pas = lambda c, x: float(((c > 70) & (x < 10)).mean() * 100)
summary = {
    "MAGs": len(cmp_), "CheckM2 tiers": dict(t2), "CheckM v1 tiers": dict(t1),
    "model used": dict(Counter(cmp_.checkm2_model)),
    "CPR (n)": int(cpr.sum()),
    "CPR >70/<10 CheckM v1 -> CheckM2 (%)": (round(pas(cmp_.checkm1_completeness[cpr], cmp_.checkm1_contamination[cpr]), 1),
                                            round(pas(cmp_.checkm2_completeness[cpr], cmp_.checkm2_contamination[cpr]), 1)),
    "DPANN (n)": int(dpann.sum()),
    "DPANN >70/<10 CheckM v1 -> CheckM2 (%)": (round(pas(cmp_.checkm1_completeness[dpann], cmp_.checkm1_contamination[dpann]), 1),
                                              round(pas(cmp_.checkm2_completeness[dpann], cmp_.checkm2_contamination[dpann]), 1)),
}
print(summary)
if not APPLY:
    print(f"wrote {f1}, {f2}. Dry run; re-run with --apply to save the KBase report.")
    sys.exit(0)

stage = tempfile.mkdtemp(prefix="ck2_report_", dir=OUT)
body = "".join(f"<tr><td>{html.escape(str(kk))}</td><td>{html.escape(str(v))}</td></tr>" for kk, v in summary.items())
page = f"""<!doctype html><html><head><meta charset="utf-8"><title>CheckM2 (all wetland MAGs)</title>
<style>body{{font-family:sans-serif;margin:16px}}td{{border:1px solid #ccc;padding:3px 8px}}</style></head><body>
<h2>CheckM2 1.1.0 on all {len(cmp_)} MAGs of narrative {k.WSID}</h2>
<p>CheckM2 1.1.0 (DIAMOND DB uniref100.KO.1), run on poplar with --allmodels; both completeness models are reported and
the model CheckM2 chose is used for the tier (high: &gt;90% complete and &lt;5% contaminated; medium: &ge;50 / &lt;10;
low: otherwise). Compared with KBase's CheckM v1 (kb_Msuite, reduced tree). MAG objects were not modified.</p>
<table>{body}</table><p>Per-MAG tables are attached (checkm2_vs_checkm1.tsv, checkm2_quality_report_allmodels.tsv).</p>
</body></html>"""
hdir = os.path.join(stage, "html")
os.makedirs(hdir)
open(os.path.join(hdir, "index.html"), "w").write(page)
hz = os.path.join(stage, "output.zip")
with zipfile.ZipFile(hz, "w", zipfile.ZIP_DEFLATED) as z:
    z.write(os.path.join(hdir, "index.html"), "index.html")
links = []
for p, d in ((f2, "Per MAG: iterativeID, CheckM2 (both models, chosen, contamination, tier), CheckM v1, DAS Tool SCG, "
                  "GTDB R232 lineage, KBase object"),
             (f1, "CheckM2 1.1.0 --allmodels quality_report.tsv (Name = MAG name) with the iterativeID column")):
    hid, node = blob_upload(p, os.path.basename(p))
    links.append({"handle": hid, "description": d, "name": os.path.basename(p), "label": "", "URL": f"{BLOB_CANON}/node/{node}"})
hid, node = blob_upload(hz, "output.zip")
rep = {"text_message": (f"CheckM2 1.1.0 for {len(cmp_)} MAGs: {t2['high']} high, {t2['medium']} medium, {t2['low']} low "
                        f"quality (CheckM v1: {t1['high']}/{t1['medium']}/{t1['low']})."),
       "warnings": [], "objects_created": [], "file_links": links,
       "html_links": [{"handle": hid, "description": "CheckM2 summary", "name": "index.html", "label": "index.html",
                       "URL": f"{BLOB_CANON}/node/{node}"}], "direct_html_link_index": 0}
name = f"CheckM2_all_MAGs_report_{time.strftime('%Y%m%d', time.gmtime())}"
prov = [{"time": time.strftime("%Y-%m-%dT%H:%M:%S+0000", time.gmtime()), "script": "census/checkm2_report.py",
         "description": "CheckM2 1.1.0 (--allmodels) on all wetland MAGs, run on poplar; summary vs CheckM v1"}]
info = k.rpc(k.WS, "Workspace.save_objects", [{"id": k.WSID, "objects": [
    {"type": "KBaseReport.Report", "name": name, "data": rep, "provenance": prov}]}], tries=1)[0][0]
back = k.rpc(k.WS, "Workspace.get_objects2", [{"objects": [{"ref": k.upa(info)}]}])[0]["data"][0]["data"]
print(f"Saved {name} as {k.upa(info)}; read back {len(back['file_links'])} files + {len(back['html_links'])} html")
