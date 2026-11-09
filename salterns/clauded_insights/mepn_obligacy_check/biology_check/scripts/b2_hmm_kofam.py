"""B2. Independent of RAST wording: KOfam HMMs (local profiles.tar.gz, extracted subset) searched against the RAST
protein translations of all 24 MePn-consumer MAGs with pyhmmer (cpus=8). KOfam adaptive thresholds from ko_list,
honouring score_type (full = sequence score, domain = best-domain score). Also records sub-threshold hits >= 40% of
the threshold for the KOs that decide the argument, so 'absent' is graded."""
import sys, json, csv
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
from kos import KOS, ALL
import pyhmmer

KO_LIST = Path("/home/freiburger/Documents/EmilyKin/kofam_new/ko_list")
thr = {}
with open(KO_LIST) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        if r["knum"] in ALL:
            thr[r["knum"]] = (float(r["threshold"]) if r["threshold"] not in ("-", "") else None, r["score_type"], r["definition"])

alpha = pyhmmer.easel.Alphabet.amino()
seqs = []
owner = {}
for mid in model_ids():
    g = rast(mid)
    for f in g["features"]:
        p = f.get("protein_translation")
        if not p:
            continue
        name = f"{iid(mid)}|{short(f['id'])}"
        fn = "; ".join(f.get("functions") or [])
        s = pyhmmer.easel.TextSequence(name=name.encode(), description=fn.encode()[:500], sequence=p.replace("*", ""))
        seqs.append(s.digitize(alpha))
        owner[name] = (iid(mid), fn)
print("proteins", len(seqs))
hmms = []
for k in ALL:
    with pyhmmer.plan7.HMMFile(str(VB / "kofam_subset" / f"{k}.hmm")) as hf:
        hmms.append(hf.read())
Z = len(seqs)
res = []
for top in pyhmmer.hmmsearch(hmms, seqs, cpus=8, E=1e-3, Z=Z):
    ko = top.query.name.decode() if isinstance(top.query.name, bytes) else top.query.name
    t, stype, defin = thr[ko]
    for h in top:
        name = h.name.decode() if isinstance(h.name, bytes) else h.name
        full = h.score
        dom = max(d.score for d in h.domains) if len(h.domains) else float("-inf")
        sc = full if stype == "full" else dom
        res.append({"ko": ko, "mag": owner[name][0], "gene": name.split("|")[1], "rast_function": owner[name][1],
                    "score_used": round(sc, 1), "score_type": stype, "threshold": t,
                    "pass": (t is not None and sc >= t), "frac_of_threshold": (round(sc / t, 2) if t else None),
                    "evalue": h.evalue})
json.dump(res, open(VB / "b2_hmm_hits.json", "w"), indent=0)
# Matrix of passing hits
mags = sorted({iid(m) for m in model_ids()}, key=lambda x: (x not in OBLIGATE, x))
passing = {}
for r in res:
    if r["pass"]:
        passing.setdefault((r["mag"], r["ko"]), []).append(r["gene"])
for grp, kos in KOS.items():
    print(f"\n### {grp}")
    print("MAG".ljust(22) + " ".join(k.rjust(7) for k in kos))
    for m in mags:
        print(("*" if m in OBLIGATE else " ") + m.ljust(21) + " ".join(str(len(passing.get((m, k), []))).rjust(7) for k in kos))
    for k in kos:
        print(f"   {k}: thr {thr[k][0]} ({thr[k][1]}) {thr[k][2]}")
print("\n### Decisive KOs: every hit >= 0.4 x threshold in the 3 obligate MAGs (and best hit in others)")
DEC = ["K00872", "K02204", "K02203", "K00948", "K01433", "K01243", "K01244", "K00772", "K00759", "K01938", "K01241",
       "K03783", "K03784", "K01839", "K00856", "K01637", "K01638", "K08691", "K14446", "K14447", "K01458"]
for k in DEC:
    print(f"-- {k} {thr[k][2]} (thr {thr[k][0]}, {thr[k][1]})")
    for m in mags:
        hs = sorted([r for r in res if r["ko"] == k and r["mag"] == m], key=lambda r: -r["score_used"])
        if not hs:
            print(f"   {m:22s} none (E<=1e-3)")
            continue
        best = hs[0]
        tag = "PASS" if best["pass"] else "sub"
        shown = [h for h in hs if h["frac_of_threshold"] is not None and h["frac_of_threshold"] >= 0.4][:3] if m in OBLIGATE else [best]
        if not shown:
            shown = [best]
        for h in shown:
            print(f"   {m:22s} {('PASS' if h['pass'] else 'sub '):4s} {h['gene']:10s} score {h['score_used']:7.1f} ({h['frac_of_threshold']}x thr)  RAST: {h['rast_function'][:90]}")
