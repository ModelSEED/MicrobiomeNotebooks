#!/usr/bin/env python
"""Analyse the 9-MAG bakta_test run (fetched from poplar).

Compares, per MAG and per marker:
  * Prodigal (-p meta) proteins + KOfam            ("prodigal_kofam")
  * Bakta CDS proteins + KOfam                     ("bakta_kofam")
  * Bakta's own KEGG db_xrefs (PSC/IPS/UPS only)   ("bakta_xref")
  * eggNOG-mapper KEGG_ko on Bakta proteins        ("eggnog")
KOfam rule = DRAM's / run_kofam_controls.py: score (full or domain, per ko_list score_type) > threshold; KOs with
threshold '-' never called; per gene the significant KO with the lowest full-sequence E-value.

Usage: analyze_bakta_test.py <RES_DIR fetched from poplar> <OUT_DIR> <marker_ko_table.tsv>
RES_DIR layout = /scratch1/afreiburger/salterns_2026/bakta_test (mags.tsv, labels.txt, bakta/, prodigal/, kofam/all.domtblout.gz,
eggnog/, logs/) plus ko_list.
"""
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

RES, OUT, MARK = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_colwidth", 70)

labels = [l.strip() for l in open(RES / "labels.txt") if l.strip()]
mags = pd.read_csv(RES / "mags.tsv", sep="\t", index_col=0)
markers = pd.read_csv(MARK, sep="\t", index_col=0)
ko_list = pd.read_csv(RES / "ko_list", sep="\t", index_col=0)
thr = pd.to_numeric(ko_list.threshold, errors="coerce")
stype = ko_list.score_type

# ---------------------------------------------------------------- KOfam hits
cols = ["gene", "glen", "ko", "qlen", "full_e", "full", "dom", "hf", "ht"]
rows = []
# keep only rows that can matter: any domain above the KO's threshold (a call), or any hit to a marker KO (for the
# "best marker KO of any score" columns). Everything else is noise under hmmsearch's default E<=10 reporting.
THR, ST, MK = thr.to_dict(), stype.to_dict(), set(markers.index)
n_all = 0
with gzip.open(RES / "kofam/all.domtblout.gz", "rt") as fh:
    for line in fh:
        n_all += 1
        f = line.split(None, 22)
        ko, full, dom = f[3], float(f[7]), float(f[13])
        t = THR.get(ko)
        sc = full if ST.get(ko) == "full" else dom
        if ko in MK or (t == t and t is not None and sc > t):
            rows.append((f[0], int(f[2]), ko, int(f[5]), float(f[6]), full, dom, int(f[15]), int(f[16])))
h = pd.DataFrame(rows, columns=cols)
del rows
print(f"KOfam domtblout: {n_all:,} rows read, {len(h):,} kept (significant or marker KO)")
h["thr"] = h.ko.map(thr)
h["stype"] = h.ko.map(stype)
h["score"] = h.full.where(h.stype == "full", h.dom)
h["sig"] = h.stype.isin(["full", "domain"]) & (h.score > h.thr)
# per gene x KO: max score over domains, min full E, hmm span
gk = h.groupby(["gene", "ko"]).agg(full_e=("full_e", "min"), score=("score", "max"), sig=("sig", "any"), thr=("thr", "first"),
                                   hf=("hf", "min"), ht=("ht", "max"), qlen=("qlen", "first"), glen=("glen", "first")).reset_index()
gk["hmm_cov"] = (gk.ht - gk.hf + 1) / gk.qlen
best_sig = gk[gk.sig].sort_values(["gene", "full_e"]).drop_duplicates("gene").set_index("gene")
best_any = gk.sort_values(["gene", "full_e"]).drop_duplicates("gene").set_index("gene")
gk_idx = gk.set_index(["gene", "ko"])
print(f"KOfam: {len(h):,} domain rows, {h.gene.nunique():,} genes with any hit, {len(best_sig):,} genes with a significant KO")


def kofam_call(gid):
    return best_sig.ko.get(gid, "")


def ko_score(gid, ko):
    try:
        r = gk_idx.loc[(gid, ko)]
        return round(float(r.score), 1)
    except KeyError:
        return None


# ---------------------------------------------------------------- eggNOG
egg = {}
ep = RES / "eggnog/bakta_test.emapper.annotations"
if ep.exists():
    lines = [x.rstrip("\n").split("\t") for x in open(ep) if not x.startswith("##") and x.strip()]
    e = pd.DataFrame(lines[1:], columns=lines[0]).rename(columns={"#query": "query"})
    for _, r in e.iterrows():
        kos = [] if r.KEGG_ko in ("-", None) or pd.isna(r.KEGG_ko) else [k.replace("ko:", "") for k in r.KEGG_ko.split(",")]
        egg[r["query"]] = dict(kos=kos, name=r.Preferred_name, desc=r.Description, ec=r.EC, og=r.max_annot_lvl)
    print(f"eggNOG: {len(e):,} Bakta proteins annotated, {sum(1 for v in egg.values() if v['kos']):,} with KEGG_ko")
else:
    print("eggNOG annotations not present yet")

# ---------------------------------------------------------------- Prodigal genes
def prodigal_genes(l):
    out = {}
    for line in open(RES / f"prodigal/{l}.gff"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        n = re.search(r"ID=\d+_(\d+)", f[8])[1]
        gid = f"prodigal|{l}|{f[0]}_{n}"
        b, e_, s = int(f[3]), int(f[4]), f[6]
        out[gid] = dict(contig=f[0], begin=b, end=e_, strand=s, stop3=e_ if s == "+" else b, aa_len=(e_ - b + 1) // 3 - 1,
                        partial=re.search(r"partial=(\d\d)", f[8])[1])
    return out


# ---------------------------------------------------------------- Bakta features
def bakta_features(l):
    d = json.load(open(RES / f"bakta/{l}/{l}.json"))
    seqmap = {s["id"]: s.get("orig_id", s["id"]) for s in d["sequences"]}
    out = {}
    for f in d["features"]:
        if f["type"] not in ("cds", "sorf"):
            continue
        gid = f"bakta|{l}|{f['locus']}"
        xr = f.get("db_xrefs", [])
        kos = sorted({x.split(":", 1)[1] for x in xr if x.startswith("KEGG:")})
        ecs = sorted({x.split(":", 1)[1] for x in xr if x.startswith("EC:")})
        # Bakta 1.12 json: 'psc' is present for any PSC alignment, with valid=True only at >=90 % id / >=80 % cov
        if f.get("ups") or f.get("ips"):
            src = "IPS/UPS"
        elif (f.get("psc") or {}).get("valid"):
            src = "PSC"
        elif f.get("pscc"):
            src = "PSCC"
        elif f.get("psc"):
            src = "PSC<50%"
        else:
            src = "none"
        aa = f.get("aa", "")
        exc = f.get("exception") or f.get("exceptions") or ""
        b, e_ = int(f["start"]), int(f["stop"])
        s = f["strand"]
        out[gid] = dict(type=f["type"], contig=seqmap.get(f["sequence"], f["sequence"]), begin=b, end=e_, strand=s,
                        stop3=e_ if s == "+" else b, aa_len=len(aa), product=f.get("product", ""), gene=f.get("gene") or "",
                        hypothetical=bool(f.get("hypothetical")), kegg=",".join(kos), ec=",".join(ecs), source=src,
                        has_U="U" in aa, exception=json.dumps(exc) if exc else "",
                        pseudo=bool(f.get("pseudo") or f.get("pseudogene")), edge=bool(f.get("edge")),
                        psc_id=(f.get("psc") or {}).get("uniref90_id", ""), pscc_id=(f.get("pscc") or {}).get("uniref50_id", ""),
                        db_xrefs=";".join(xr))
    return d, out


def match_genes(bk, pr):
    """Bakta -> Prodigal: same contig+strand+3' end; else >= 50 % overlap of the shorter (same strand)."""
    by_stop = {(v["contig"], v["strand"], v["stop3"]): k for k, v in pr.items()}
    by_contig = defaultdict(list)
    for k, v in pr.items():
        by_contig[(v["contig"], v["strand"])].append((v["begin"], v["end"], k))
    m = {}
    for k, v in bk.items():
        hit = by_stop.get((v["contig"], v["strand"], v["stop3"]))
        how = "stop"
        if hit is None:
            how = "overlap"
            cands = []
            for b, e_, pk in by_contig.get((v["contig"], v["strand"]), []):
                ov = min(e_, v["end"]) - max(b, v["begin"]) + 1
                if ov > 0:
                    cands.append((ov, pk, b, e_))
            ok = [c for c in cands if c[0] >= 0.5 * min(c[3] - c[2] + 1, v["end"] - v["begin"] + 1)]
            if ok:
                hit = ",".join(c[1] for c in sorted(ok, key=lambda c: c[2]))
        m[k] = (hit or "", how if hit else "none")
    return m


def parse_time(p):
    if not p.exists():
        return {}
    t = p.read_text()
    g = lambda pat: re.search(pat, t)
    wall = g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)")[1]
    parts = [float(x) for x in wall.split(":")]
    wall_s = parts[0] * 3600 + parts[1] * 60 + parts[2] if len(parts) == 3 else parts[0] * 60 + parts[1]
    return dict(wall_min=round(wall_s / 60, 1), cpu_min=round((float(g(r"User time \(seconds\): (\S+)")[1]) +
                float(g(r"System time \(seconds\): (\S+)")[1])) / 60, 1),
                max_rss_gb=round(int(g(r"Maximum resident set size \(kbytes\): (\d+)")[1]) / 1e6, 2),
                exit=int(g(r"Exit status: (\d+)")[1]))


MARKER_KOS = list(markers.index)
per_mag, marker_rows, gene_rows, all_bakta = [], [], [], []
for l in labels:
    d, bk = bakta_features(l)
    pr = prodigal_genes(l)
    mm = match_genes(bk, pr)
    bcds = {k: v for k, v in bk.items() if v["type"] == "cds"}
    n_b, n_p = len(bcds), len(pr)
    kb = {k: kofam_call(k) for k in bk}
    kp = {k: kofam_call(k) for k in pr}
    xr = {k: v["kegg"].split(",") if v["kegg"] else [] for k, v in bk.items()}
    eg = {k: egg.get(f"{l}|{k.split('|')[2]}", {}).get("kos", []) for k in bk}
    src = pd.Series([v["source"] for v in bcds.values()]).value_counts()
    # agreement on Bakta CDS that have a KOfam KO
    agree = defaultdict(int)
    for k in bcds:
        kf = kb[k]
        if not kf:
            continue
        agree["kofam"] += 1
        if xr[k]:
            agree["xref_any"] += 1
            agree["xref_same"] += kf in xr[k]
        if eg[k]:
            agree["egg_any"] += 1
            agree["egg_same"] += kf in eg[k]
    r = dict(label=l, taxon=mags.taxon.get(l, ""), bakta_cds=n_b, bakta_sorf=sum(v["type"] == "sorf" for v in bk.values()),
             prodigal_genes=n_p,
             bakta_stop_matched=sum(1 for k in bcds if mm[k][1] == "stop"),
             bakta_cds_with_U=sum(v["has_U"] for v in bcds.values()),
             bakta_cds_with_exception=sum(bool(v["exception"]) for v in bcds.values()),
             bakta_pseudo=sum(v["pseudo"] for v in bk.values()),
             frac_bakta_kofam=round(sum(bool(kb[k]) for k in bcds) / n_b, 3),
             frac_prodigal_kofam=round(sum(bool(v) for v in kp.values()) / n_p, 3),
             frac_bakta_xref_KO=round(sum(bool(xr[k]) for k in bcds) / n_b, 3),
             frac_eggnog_KO=round(sum(bool(eg[k]) for k in bcds) / n_b, 3) if egg else None,
             frac_bakta_EC=round(sum(bool(v["ec"]) for v in bcds.values()) / n_b, 3),
             frac_bakta_nonhypothetical=round(sum(not v["hypothetical"] for v in bcds.values()) / n_b, 3),
             frac_IPS_UPS=round(src.get("IPS/UPS", 0) / n_b, 3), frac_PSC=round(src.get("PSC", 0) / n_b, 3),
             frac_PSCC_only=round(src.get("PSCC", 0) / n_b, 3), frac_weak_PSC_only=round(src.get("PSC<50%", 0) / n_b, 3),
             frac_no_hit=round(src.get("none", 0) / n_b, 3),
             kofam_KO_also_in_bakta_xref=f"{agree['xref_same']}/{agree['xref_any']}",
             kofam_KO_also_in_eggnog=f"{agree['egg_same']}/{agree['egg_any']}" if egg else "",
             distinct_KO_prodigal_kofam=len({v for v in kp.values() if v}),
             distinct_KO_bakta_kofam=len({kb[k] for k in bcds if kb[k]}),
             distinct_KO_bakta_xref=len({x for k in bcds for x in xr[k]}),
             distinct_KO_eggnog=len({x for k in bcds for x in eg[k]}) if egg else None)
    r.update({f"bakta_{k}": v for k, v in parse_time(RES / f"logs/bakta_{l}.time").items()})
    r.update({f"prodigal_{k}": v for k, v in parse_time(RES / f"logs/prodigal_{l}.time").items()})
    per_mag.append(r)
    for ko in MARKER_KOS:
        marker_rows.append(dict(label=l, ko=ko, gene=markers.gene[ko],
                                prodigal_kofam=sum(v == ko for v in kp.values()),
                                bakta_kofam=sum(kb[k] == ko for k in bk),
                                bakta_xref=sum(ko in xr[k] for k in bk),
                                eggnog=sum(ko in eg[k] for k in bk) if egg else None))
    # gene-level table for marker-relevant Bakta CDS (any marker KO >= 100 score, or a marker KO call by any method)
    for k, v in bk.items():
        pids = mm[k][0].split(",") if mm[k][0] else []
        mk_any = best_any.ko.get(k, "")
        rel = (kb[k] in MARKER_KOS or any(x in MARKER_KOS for x in xr[k] + eg[k])
               or (mk_any in MARKER_KOS and best_any.score.get(k, 0) >= 100)
               or any(kp.get(p, "") in MARKER_KOS for p in pids))
        all_bakta.append(dict(label=l, locus=k.split("|")[2], kofam=kb[k], xref=",".join(xr[k]), eggnog=",".join(eg[k]),
                              product=v["product"], source=v["source"]))
        if not rel:
            continue
        gene_rows.append(dict(
            label=l, bakta_locus=k.split("|")[2], contig=v["contig"], begin=v["begin"], end=v["end"], strand=v["strand"],
            bakta_aa=v["aa_len"], has_U=v["has_U"], exception=v["exception"], bakta_gene=v["gene"], bakta_product=v["product"],
            bakta_source=v["source"], bakta_xref_KO=",".join(xr[k]), bakta_EC=v["ec"], psc=v["psc_id"], pscc=v["pscc_id"],
            bakta_kofam_KO=kb[k], bakta_kofam_gene=markers.gene.get(kb[k], "") if kb[k] else "",
            bakta_kofam_score=round(float(best_sig.score[k]), 1) if kb[k] else None,
            bakta_best_any_KO=mk_any, bakta_best_any_score=round(float(best_any.score[k]), 1) if mk_any else None,
            bakta_best_any_thr=best_any.thr.get(k) if mk_any else None,
            bakta_best_any_hmm_cov=round(float(best_any.hmm_cov[k]), 2) if mk_any else None,
            eggnog_KO=",".join(eg[k]), eggnog_name=egg.get(f"{l}|{k.split('|')[2]}", {}).get("name", ""),
            eggnog_desc=egg.get(f"{l}|{k.split('|')[2]}", {}).get("desc", ""),
            prodigal_match=mm[k][1], prodigal_genes=",".join(p.split("|")[2] for p in pids),
            prodigal_aa=",".join(str(pr[p]["aa_len"]) for p in pids),
            prodigal_kofam_KO=",".join(kp.get(p, "") or "-" for p in pids),
            prodigal_best_any=",".join(f"{best_any.ko.get(p, '-')}:{round(float(best_any.score[p]), 1) if p in best_any.index else '-'}"
                                       for p in pids),
            K21579_score_bakta=ko_score(k, "K21579"),
            K21579_score_prodigal=",".join(str(ko_score(p, "K21579")) for p in pids)))
    # Prodigal marker calls with no Bakta CDS counterpart at all
    matched_p = {p for k in bk for p in (mm[k][0].split(",") if mm[k][0] else [])}
    for p, v in pr.items():
        if p not in matched_p and kp[p] in MARKER_KOS:
            gene_rows.append(dict(label=l, bakta_locus="(none)", contig=v["contig"], begin=v["begin"], end=v["end"],
                                  strand=v["strand"], prodigal_match="prodigal-only", prodigal_genes=p.split("|")[2],
                                  prodigal_aa=str(v["aa_len"]), prodigal_kofam_KO=kp[p],
                                  K21579_score_prodigal=str(ko_score(p, "K21579"))))

pm = pd.DataFrame(per_mag).set_index("label")
mk = pd.DataFrame(marker_rows)
gr = pd.DataFrame(gene_rows)
pm.to_csv(OUT / "per_mag_summary.tsv", sep="\t")
gr.to_csv(OUT / "marker_genes_detail.tsv", sep="\t", index=False)
wide = mk.pivot_table(index=["ko", "gene"], columns="label", values=["prodigal_kofam", "bakta_kofam", "bakta_xref", "eggnog"],
                      aggfunc="sum", fill_value=0)
mk.to_csv(OUT / "marker_calls_long.tsv", sep="\t", index=False)
# compact per marker: "P/B/X/E" = prodigal_kofam / bakta_kofam / bakta_xref / eggnog counts
mk["PBXE"] = mk.apply(lambda r: f"{r.prodigal_kofam}/{r.bakta_kofam}/{r.bakta_xref}/{'' if pd.isna(r.eggnog) else int(r.eggnog)}", axis=1)
compact = mk.pivot(index=["ko", "gene"], columns="label", values="PBXE")[labels]
nz = mk.groupby("ko")[["prodigal_kofam", "bakta_kofam", "bakta_xref", "eggnog"]].sum()
compact = compact[compact.index.get_level_values(0).isin(nz[nz.sum(axis=1) > 0].index)]
compact.to_csv(OUT / "marker_calls_PBXE.tsv", sep="\t")
tot = mk.groupby("label")[["prodigal_kofam", "bakta_kofam", "bakta_xref", "eggnog"]].sum().loc[labels]
tot.to_csv(OUT / "marker_calls_totals.tsv", sep="\t")
pd.DataFrame(all_bakta).to_csv(OUT / "bakta_cds_ko_calls.tsv.gz", sep="\t", index=False)

print("\n=== per MAG ===")
print(pm.T.to_string())
print("\n=== marker calls per MAG (prodigal_kofam/bakta_kofam/bakta_xref/eggnog) ===")
print(compact.to_string())
print("\n=== marker totals ===")
print(tot.to_string())
