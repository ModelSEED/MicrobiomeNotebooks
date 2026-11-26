"""Independent re-derivation of Table S6 buckets from the R232 skani output.

Written from the Table S6 rules (best hit = highest ANI among hits with AF_query >= 15 %; if none,
highest ANI overall -> low_AF_unresolved; no hit -> <85; buckets >=95 / 85-95 / <85), not from the
builder's code. Reference taxonomy from the GTDB R232 bac120 file in gtdb_downloads (the critic's
source), cross-checked against the gtdbtk_data copy. Read-only on the repo.
"""
import csv, gzip, io, os, re, subprocess
from collections import Counter

REPO = "/home/freiburger/Documents/EmilyKin"
PREP = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep"
SK = f"{REPO}/reprocess_291/skani_novelty_r232/mag_vs_gtdb_r232_skani.tsv"

tax = {}
with gzip.open(f"{PREP}/gtdb_downloads/r232/bac120_taxonomy_r232.tsv.gz", "rt") as fh:
    for line in fh:
        a, t = line.rstrip("\n").split("\t")[:2]
        tax[a[3:]] = t
tax2 = {}
for line in open(os.path.expanduser("~/Documents/gtdbtk_data/release232/taxonomy/gtdb_taxonomy.tsv")):
    a, t = line.rstrip("\n").split("\t")[:2]
    tax2[re.sub(r"^(RS_|GB_)", "", a)] = t

r232 = {r["user_genome"]: r for r in csv.DictReader(open(f"{PREP}/results/gtdbtk.r232.summary.tsv"), delimiter="\t")}

hits = {}
qfiles = set()
for r in csv.DictReader(open(SK), delimiter="\t"):
    q = os.path.basename(r["Query_file"])
    qfiles.add(os.path.dirname(r["Query_file"]))
    mag = q[:-3] if q.endswith(".fa") else q
    acc = re.search(r"(GC[AF]_\d+\.\d+)", os.path.basename(r["Ref_file"])).group(1)
    hits.setdefault(mag, []).append((float(r["ANI"]), float(r["Align_fraction_query"]), float(r["Align_fraction_ref"]), acc))

mags = sorted(r232)
cls = {}
best = {}
for m in mags:
    hs = hits.get(m, [])
    el = [h for h in hs if h[1] >= 15.0]
    if el:
        b = max(el, key=lambda h: h[0])
        c = ">=95" if b[0] >= 95 else ("85-95" if b[0] >= 85 else "<85")
    elif hs:
        b = max(hs, key=lambda h: h[0])
        c = "low_AF_unresolved"
    else:
        b = None
        c = "<85"
    cls[m] = c
    best[m] = b

print("query dirs:", qfiles)
print("MAGs with any hit:", len(hits), "of", len(mags), "; hits:", sum(len(v) for v in hits.values()),
      "; max hits per MAG:", max(len(v) for v in hits.values()))
print("re-derived buckets:", Counter(cls.values()))

cur = {r["MAG"]: r for r in csv.DictReader(open(f"{REPO}/Table_S6_skani_novelty.tsv", newline=""), delimiter="\t")}
assert set(cur) == set(mags), (len(cur), len(mags))
print("committed-table buckets:", Counter(r["ANI_class"] for r in cur.values()))
bad = []
for m in mags:
    r = cur[m]
    b = best[m]
    exp = dict(ANI_class=cls[m])
    if b:
        exp.update(skani_best_ref=b[3], skani_ANI=f"{b[0]:.2f}", skani_AF_query=f"{b[1]:.2f}", skani_AF_ref=f"{b[2]:.2f}",
                   skani_best_taxonomy=tax.get(b[3], ""))
    else:
        exp.update(skani_best_ref="", skani_ANI="", skani_AF_query="", skani_AF_ref="", skani_best_taxonomy="")
    exp["mag_gtdb_taxonomy"] = r232[m]["classification"]
    for k, v in exp.items():
        if r[k] != v:
            bad.append((m, k, r[k][:70], v[:70]))
print("cell mismatches vs re-derivation:", len(bad), bad[:10])
miss_tax = [m for m in mags if best[m] and best[m][3] not in tax]
print("best refs absent from bac120_taxonomy_r232:", miss_tax)
dis = [m for m in mags if best[m] and tax.get(best[m][3]) != tax2.get(best[m][3])]
print("bac120 vs gtdbtk_data taxonomy disagreements for best refs:", dis)

# HEAD comparison and change table
head = {r["MAG"]: r for r in csv.DictReader(io.StringIO(subprocess.run(
    ["git", "-C", REPO, "show", "HEAD:Table_S6_skani_novelty.tsv"], capture_output=True, text=True).stdout.replace("\r\n", "\n")), delimiter="\t")}
print("HEAD buckets:", Counter(r["ANI_class"] for r in head.values()))
trans = Counter((head[m]["ANI_class"], cls[m]) for m in mags)
print("bucket transitions HEAD->R232:", sorted(trans.items(), key=lambda x: -x[1]))
for m in ["coasm_bin.476", "coasm_bin.535", "CAN_5_bin.70", "coasm_bin.260", "coasm_bin.185"]:
    print(m, "HEAD:", head[m]["ANI_class"], head[m]["skani_ANI"], head[m]["skani_AF_query"], head[m]["skani_best_taxonomy"].split(";")[-1],
          "| R232:", cur[m]["ANI_class"], cur[m]["skani_ANI"], cur[m]["skani_AF_query"], cur[m]["skani_best_taxonomy"].split(";")[-1],
          "| all hits:", sorted(hits.get(m, []), reverse=True)[:3])
# sentence check
sent_bad = 0
for m in mags:
    s = cur[m]["Novelty_interpretation"]
    c = cls[m]
    key = {">=95": "Same species", "85-95": "Same genus, novel species", "<85": ("Novel at genus", "No skani hit"),
           "low_AF_unresolved": "Unresolved by ANI"}[c]
    if not s.startswith(key if isinstance(key, str) else key):
        sent_bad += 1
print("sentence/class inconsistencies:", sent_bad)
# relation to GTDB-Tk species calls
gt_sp = {m: r232[m]["classification"].split("s__")[-1] for m in mags}
agree = Counter()
for m in mags:
    has_sp = bool(gt_sp[m])
    agree[(cls[m], has_sp)] += 1
print("bucket x GTDB-Tk species assigned:", sorted(agree.items()))
inc = [(m, cls[m], gt_sp[m], r232[m]["closest_genome_ani"]) for m in mags if (cls[m] == ">=95") != bool(gt_sp[m])]
print("bucket >=95 vs GTDB-Tk species disagreements:", len(inc))
for x in inc:
    print("   ", x, "| skani best:", best[x[0]])
