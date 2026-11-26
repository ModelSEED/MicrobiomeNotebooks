#!/usr/bin/env python
"""Read-only check of how the salterns scripts read the 310 RAST genome JSONs.

For every MAG in correlational/r232/guilds_r232.csv it reads
  jplfaria/Cliff MAGs/RASTified_MAG_json/<mag>.contigs__.RAST.json
and compares four ways of getting the MAG's function strings:
  trunc  - rerun_correlational_r232.build_guilds (lines 107, 112): raw-text regex, stops each list at the first ']'
  aware  - make_guild_membership_comparison.rast_functions (lines 76-77): raw-text, string-aware
  jfeat  - json.load, features[].functions (what make_archaeal_mepn_tma_report / rerun_betaine_grdH_r232 use)
  jall   - json.load, features[].functions + cdss[].functions
For each guild pattern (rerun_correlational_r232.GENES, copied verbatim) it reports whether the MAG-level hit set
differs between methods, and it rebuilds the ev_* columns from 'trunc' to compare with the committed guilds_r232.csv.
Writes only to the directory given by --out.
"""
import argparse
import collections
import json
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

S = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns")
RAST = S / "jplfaria/Cliff MAGs/RASTified_MAG_json"
GENES = {   # verbatim from correlational/r232/rerun_correlational_r232.py lines 52-60
    "phnJ": r"c-p lyase|ec 4\.7\.1\.1\b",
    "phn_accessory": r"methylphosphonate 5-triphosphate (synthase|diphosphatase)|phnn|ec 2\.7\.4\.23",
    "grdH": r"betaine reductase component b",
    "cutC_cluster": r"choline trimethylamine.lyase",
    "cntA": r"carnitine monooxygenase",
    "mcrA": r"methyl.coenzyme m reductase",
    "mtbB_mtmB": r"(dimethylamine|monomethylamine)\s*(:corrinoid )?methyltransferase",
}
EXTRA = {   # verbatim from correlational/r232/make_guild_membership_comparison.py lines 48-58
    "grd_shared": r"glycine/sarcosine/betaine reductase|glycine reductase component b|^grdx protein",
    "methylamine_corrinoid": r"methylamine--corrinoid protein\] co-methyltransferase|ec 2\.1\.1\.2(47|48|49|50)\b"
                             r"|(trimethylamine|dimethylamine|monomethylamine) methyltransferase corrinoid protein"
                             r"|methylcobalamin:coenzyme m methyltransferase",
    "phn_other": r"\bphn[a-z]\d?\b|phosphonate.*(abc|transport)|(abc|transport).*phosphonate",
}
TRUNC_BLOCK = re.compile(r'"functions"\s*:\s*\[(.*?)\]', re.S)                                   # build_guilds l.107
AWARE_BLOCK = re.compile(r'"functions"\s*:\s*\[((?:\s*"(?:[^"\\]|\\.)*"\s*,?)*)\s*\]', re.S)     # rast_functions l.76


def one(mag):
    text = (RAST / f"{mag}.contigs__.RAST.json").read_text(errors="ignore")
    trunc = {s for blk in TRUNC_BLOCK.findall(text) for s in re.findall(r'"([^"]*)"', blk)}
    aware = {s for blk in AWARE_BLOCK.findall(text) for s in re.findall(r'"((?:[^"\\]|\\.)*)"', blk)}
    d = json.loads(text)
    feats, cdss = d.get("features", []), d.get("cdss", [])
    jfeat = {s for f in feats for s in f.get("functions", [])}
    jcds = {s for c in cdss for s in c.get("functions", [])}
    jall = jfeat | jcds
    cds_by_id = {c["id"]: c for c in cdss}
    n_feat_cds_mismatch = sum(1 for f in feats if f.get("cdss") and cds_by_id.get(f["cdss"][0], {}).get("functions", [""])
                              != (f.get("functions") or [""]))
    nfun = collections.Counter(len(f.get("functions", [])) for f in feats)
    strings = [s for f in feats for s in f.get("functions", [])]
    st = dict(
        mag=mag, n_features=len(feats), n_cdss=len(cdss), n_noncoding=len(d.get("non_coding_features", [])),
        n_features_no_functions_key=sum(1 for f in feats if "functions" not in f),
        n_features_multi_entry=sum(v for k, v in nfun.items() if k >= 2), max_entries=max(nfun) if nfun else 0,
        n_feat_cds_function_mismatch=n_feat_cds_mismatch,
        n_unique_functions_jall=len(jall - {""}), n_unique_trunc=len(trunc - {""}), n_unique_aware=len(aware - {""}),
        n_dropped_by_trunc=len(jall - trunc), n_aware_minus_jall=len(aware - jall),
        n_trunc_not_in_jall=len(trunc - jall),
        n_str_slash=sum(" / " in s for s in strings), n_str_at=sum(" @ " in s for s in strings),
        n_str_semicolon=sum("; " in s for s in strings), n_str_hash=sum("#" in s for s in strings),
        n_str_bracket=sum(("[" in s or "]" in s) for s in strings), n_str_quote=sum('"' in s for s in strings),
        n_str_backslash=sum("\\" in s for s in strings), n_str_nonascii=sum(any(ord(c) > 127 for c in s) for s in strings),
        n_raw_unicode_escapes=text.count("\\u"),
        domain=d.get("domain"), genetic_code=d.get("genetic_code"), num_contigs=d.get("num_contigs"),
        ontologies=";".join(f"{k}:{len(v)}" for k, v in (d.get("ontologies_present") or {}).items()),
        ontology_events=";".join(f"{e.get('id')}|{e.get('method')}|{e.get('method_version')}|{e.get('timestamp')}"
                                 for e in d.get("ontology_events", [])),
        genome_id=d.get("id"),
    )
    hits = {}
    for name, pat in {**GENES, **{"x_" + k: v for k, v in EXTRA.items()}}.items():
        rx = re.compile(pat)
        for meth, fs in (("trunc", trunc), ("aware", aware), ("jall", jall), ("jfeat", jfeat)):
            hits[(name, meth)] = sorted(s for s in fs if rx.search(s.lower()))
    dropped_examples = sorted(jall - trunc)[:5]
    return st, hits, dropped_examples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    g = pd.read_csv(S / "correlational/r232/guilds_r232.csv").fillna("")
    mags = list(g.mag)
    with Pool(min(a.procs, 16)) as pool:
        res = pool.map(one, mags, chunksize=4)
    stats = pd.DataFrame([r[0] for r in res])
    stats.to_csv(a.out / "rast_parsing_per_mag.csv", index=False)
    # method agreement per pattern at the MAG level (presence) and string level
    rows, dropped_any = [], collections.Counter()
    for (st, hits, drop) in res:
        for s in drop:
            dropped_any[s] += 1
    for name in list(GENES) + ["x_" + k for k in EXTRA]:
        pres = {m: {st["mag"] for st, h, _ in res if h[(name, m)]} for m in ("trunc", "aware", "jall", "jfeat")}
        str_diff = sum(1 for st, h, _ in res if h[(name, "trunc")] != h[(name, "jall")])
        rows.append(dict(pattern=name, n_mags_trunc=len(pres["trunc"]), n_mags_aware=len(pres["aware"]),
                         n_mags_jall=len(pres["jall"]), n_mags_jfeat=len(pres["jfeat"]),
                         mags_jall_not_trunc=";".join(sorted(g.set_index("mag").loc[list(pres["jall"] - pres["trunc"]), "iterativeID_r232"])),
                         n_mags_string_set_differs_trunc_vs_jall=str_diff))
    agree = pd.DataFrame(rows)
    agree.to_csv(a.out / "rast_parsing_pattern_agreement.csv", index=False)
    # distinct RAST strings each pattern matched (json-parsed, all features + CDSs), with MAG counts
    vocab = []
    iid = g.set_index("mag")["iterativeID_r232"].to_dict()
    for name in list(GENES) + ["x_" + k for k in EXTRA]:
        c = collections.defaultdict(set)      # string -> MAGs whose json-parsed functions contain it
        lost = collections.defaultdict(set)   # string -> MAGs where the truncating parser misses it
        for st, h, _ in res:
            tr = set(h[(name, "trunc")])
            for s in h[(name, "jall")]:
                c[s].add(st["mag"])
                if s not in tr:
                    lost[s].add(st["mag"])
        for s, ms in sorted(c.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            vocab.append(dict(pattern=name, regex=(GENES.get(name) or EXTRA[name[2:]]), rast_string=s, n_mags=len(ms),
                              n_mags_missed_by_trunc_parser=len(lost.get(s, ())),
                              mags=";".join(sorted(iid[m] for m in ms)) if len(ms) <= 30 else f"{len(ms)} MAGs"))
    pd.DataFrame(vocab).to_csv(a.out / "guild_pattern_matched_strings.csv", index=False)
    # rebuild ev_* columns from 'trunc' and compare with the committed guilds_r232.csv
    by = {st["mag"]: h for st, h, _ in res}
    mism = {}
    for name in GENES:
        col = f"ev_{name}"
        rebuilt = g.mag.map(lambda m: " | ".join(by[m][(name, "trunc")]))
        mism[col] = int((rebuilt != g[col].astype(str)).sum())
    flags = {"mepn_consumer": "phnJ", "betaine_reducer": "grdH", "methanogen": "mcrA"}
    flag_mism = {c: int((g.mag.map(lambda m: bool(by[m][(p, 'trunc')])) != g[c].astype(bool)).sum()) for c, p in flags.items()}
    tma = g.mag.map(lambda m: bool(by[m][("cutC_cluster", "trunc")] or by[m][("cntA", "trunc")]))
    flag_mism["tma_producer"] = int((tma != g["tma_producer"].astype(bool)).sum())
    summary = dict(
        n_mags=len(mags),
        ev_column_mismatches_vs_committed=mism, flag_mismatches_vs_committed=flag_mism,
        totals={k: int(stats[k].sum()) for k in stats.columns if k.startswith("n_")},
        mags_with_dropped_strings=int((stats.n_dropped_by_trunc > 0).sum()),
        max_dropped_in_one_mag=int(stats.n_dropped_by_trunc.max()),
        domain_counts=stats.domain.value_counts().to_dict(),
        genetic_code_counts=stats.genetic_code.value_counts().to_dict(),
        ontology_events=stats.ontology_events.value_counts().to_dict(),
        most_common_dropped=dropped_any.most_common(15),
        any_trunc_string_not_in_json=int(stats.n_trunc_not_in_jall.sum()),
        any_aware_string_not_in_json=int(stats.n_aware_minus_jall.sum()),
    )
    json.dump(summary, open(a.out / "rast_parsing_summary.json", "w"), indent=1, default=str)
    print(json.dumps(summary, indent=1, default=str)[:6000])
    print(agree.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
