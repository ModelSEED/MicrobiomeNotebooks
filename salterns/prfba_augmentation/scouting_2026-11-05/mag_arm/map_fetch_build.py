#!/usr/bin/env python
"""prFBA MAG-arm smoke test, steps 3-5 (after mag_skani_hits.py search + select_references.py).

  3. map each selected GTDB accession -> BV-BRC genome_id with mag_skani_hits.load_acc_map
     (not wired into the shipped pipeline; this is the manual opt-in);
  4. fetch the selected genome's functions with the SHIPPED fast builder
     (build_synthetic_genomes.fast_build_genome / build_genome), with one change: a
     browser-like User-Agent installed on urllib (the shipped _api_json sends none -> 403);
     PGFams with prefetch_gene_families.fetch_pgfams (shipped, has its own UA);
  5. merge with the SHIPPED build_for_asv -> create_synthetic_genome, read the per-function
     probability it writes, and compare the reference function set with the MAG's own RAST
     annotation (exact strings, ModelSEED-style roles, lenient roles) + marker genes.
  Also: failure-mode probes of the unchanged code, and a prototype MAG(p=1)+reference merge.

Writes only under OUT (scratch). Never writes to prFBA's caches: genome cache dir is scratch,
no gene-provider is constructed, run with -B so no __pycache__ lands in the repos.
"""
import argparse
import collections
import copy
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request

PRFBA = "/home/freiburger/Documents/prFBA"
sys.path.insert(0, PRFBA)
import mag_skani_hits as M            # noqa: E402  load_acc_map
import select_references as S         # noqa: E402
import prefetch_gene_families as PF   # noqa: E402
import build_synthetic_genomes as BSG  # noqa: E402  (inserts codiffusion_bioreactor on sys.path)

ACC_MAP = f"{PRFBA}/reference_maps/gtdb_r226_acc_to_bvbrc.tsv.gz"
RAST_DIR = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/jplfaria/Cliff MAGs/RASTified_MAG_json"
UA = "Mozilla/5.0 prFBA"

# ---- marker genes (RAST / BV-BRC product-name patterns, lower-case). The first five reuse
# correlational/r232/rerun_correlational_r232.py:52-60; the rest are added here. ----
MARKERS = {
    "PhnJ (C-P lyase, EC 4.7.1.1)": r"c-p lyase|ec 4\.7\.1\.1\b|\bphnj\b",
    "PhnG": r"\bphng\b",
    "PhnH": r"\bphnh\b",
    "PhnI": r"\bphni\b",
    "PhnK": r"\bphnk\b",
    "PhnL": r"\bphnl\b",
    "PhnM": r"\bphnm\b|methylphosphonate 5-triphosphate diphosphatase|ec 3\.6\.1\.63\b",
    "GrdH (betaine reductase B, EC 1.21.4.4)": r"betaine reductase component b|ec 1\.21\.4\.4\b|\bgrdh\b",
    "McrA (EC 2.8.4.1)": r"methyl.coenzyme m reductase[^/@;]*alpha|\bmcra\b",
    "Mcr alpha/beta/gamma (EC 2.8.4.1)": r"methyl.coenzyme m reductase (alpha|beta|gamma) subunit|ec 2\.8\.4\.1\b",
    "CutC (choline TMA-lyase)": r"choline trimethylamine.lyase|ec 4\.3\.99\.4\b",
    "CntA (carnitine monooxygenase)": r"carnitine monooxygenase|\bcnta\b",
    "TorA/TorZ (TMAO reductase)": r"trimethylamine.n.oxide reductase|tmao reductase|\btor[az]\b|ec 1\.7\.2\.3\b",
    # pyrrolysine MttB/MtbB/MtmB: RAST names them '[X--corrinoid protein] Co-methyltransferase (EC 2.1.1.250/249/248)'
    "MttB (EC 2.1.1.250)": r"\[?trimethylamine-+corrinoid protein\]? co-methyltransferase|ec 2\.1\.1\.250\b|trimethylamine:corrinoid methyltransferase|\bmttb\b",
    "MtbB (EC 2.1.1.249)": r"\[?dimethylamine-+corrinoid protein\]? co-methyltransferase|ec 2\.1\.1\.249\b|dimethylamine:corrinoid methyltransferase|\bmtbb\b",
    "MtmB (EC 2.1.1.248)": r"\[methylamine-+corrinoid protein\]? co-methyltransferase|^methylamine-+corrinoid protein co-methyltransferase|ec 2\.1\.1\.248\b|monomethylamine:corrinoid methyltransferase|\bmtmb\b",
    "MttB superfamily, non-pyl homolog (RAST 'Trimethylamine methyltransferase family')": r"trimethylamine methyltransferase family",
}
# strings a pattern catches that are NOT the enzyme named by the marker
EXCLUDE = {
    "McrA (EC 2.8.4.1)": r"thioamidation|restriction|5-methylcytosine",
    "CutC (choline TMA-lyase)": r"activating|clustered with",
    "CntA (carnitine monooxygenase)": r"reductase",
    "TorA/TorZ (TMAO reductase)": r"\btor[cdesrt]\b|cytochrome c|chaperone|sensor|regulator",
}


def marker_hits(funcs):
    out = {}
    for name, pat in MARKERS.items():
        hits = sorted(f for f in funcs if re.search(pat, f.lower()))
        ex = EXCLUDE.get(name)
        if ex:
            hits = [f for f in hits if not re.search(ex, f.lower())]
        out[name] = hits
    return out


# ---- normalisation levels for the set comparison ----
_SPLIT = re.compile(r"\s*;\s+|\s+[@/]\s+")   # ModelSEED role split: '; ', ' @ ', ' / '


def roles(fn):
    fn = re.sub(r"\s*#.*$", "", fn)           # drop RAST '# comment' / '## comment'
    return {re.sub(r"\s+", " ", r).strip().lower() for r in _SPLIT.split(fn) if r.strip()}


def lenient(role):
    role = re.sub(r"\(\s*(ec|tc)[^)]*\)", "", role)   # drop EC / TC numbers
    return re.sub(r"[^a-z0-9]", "", role)


def is_hypo(fn):
    return fn.strip().lower() in ("hypothetical protein", "")


def mag_rast_functions(mag):
    g = json.load(open(os.path.join(RAST_DIR, f"{mag}.RAST.json")))
    funcs = collections.Counter()
    for f in g["features"]:
        for fn in f.get("functions") or []:
            if fn and fn.strip():
                funcs[fn] += 1
    return g, funcs


def compare_sets(ref_funcs, mag_funcs):
    """ref_funcs, mag_funcs: sets of function strings."""
    ref_nh = {f for f in ref_funcs if not is_hypo(f)}
    mag_nh = {f for f in mag_funcs if not is_hypo(f)}
    absent_exact = ref_nh - mag_nh
    mag_roles = set().union(*(roles(f) for f in mag_nh)) if mag_nh else set()
    mag_len = {lenient(r) for r in mag_roles}
    absent_role = {f for f in ref_nh if not roles(f) <= mag_roles}
    absent_len = {f for f in ref_nh if not {lenient(r) for r in roles(f)} <= mag_len}
    ref_roles = set().union(*(roles(f) for f in ref_nh)) if ref_nh else set()
    return dict(
        n_ref_distinct_functions=len(ref_funcs), n_ref_hypothetical=len(ref_funcs - ref_nh),
        n_ref_nonhypothetical=len(ref_nh), n_mag_distinct_functions=len(mag_funcs),
        n_mag_nonhypothetical=len(mag_nh),
        n_absent_exact=len(absent_exact), n_absent_role=len(absent_role),
        n_absent_role_lenient=len(absent_len),
        n_ref_roles=len(ref_roles), n_ref_roles_absent=len(ref_roles - mag_roles),
        n_mag_nonhypo_absent_from_ref_exact=len(mag_nh - ref_nh),
    ), absent_exact, absent_role, absent_len


# ---- prototype: MAG's own RAST annotation (p=1, observed) + reference functions (inferred) ----
def merge_mag_with_reference(mag_genome, ref_synth, ref_sources, mag_id, *, inferred_weight=1.0, match="role"):
    """MAG(p=1) + reference merge, KBase-Genome shaped (what load_genomes_to_modelseedpy reads).

    * every MAG RAST feature is kept as-is (its own id, sequences, functions) with
      probability 1.0 and provenance 'observed_in_mag' (or 'observed_in_mag+reference' when the
      selected references also carry the function); p_reference records the reference support;
    * every non-hypothetical reference function NOT observed in the MAG becomes one new feature
      with probability = reference prevalence k/n (exactly create_synthetic_genome's value)
      * inferred_weight, provenance 'inferred_from_reference', and the carrying reference ids.
    match = 'exact' (string identity, the merge's own key) or 'role' (ModelSEED role split)
    decides when a reference function counts as already observed in the MAG."""
    ref = {}
    for f in ref_synth.get("features", []):
        for fn in f.get("functions") or []:
            ref[fn] = {"p": float(f.get("probability", 0.0)), "aliases": f.get("aliases", [])}
    carriers = collections.defaultdict(list)
    for g in ref_sources:
        for fn in {fn for ft in g.get("features", []) for fn in (ft.get("functions") or [])}:
            carriers[fn].append(str(g.get("id")))
    role_p = {}
    for fn, v in ref.items():
        for r in roles(fn):
            role_p[r] = max(role_p.get(r, 0.0), v["p"])

    def ref_support(fn):
        if fn in ref:
            return ref[fn]["p"]
        rs = roles(fn)
        return min(role_p[r] for r in rs) if rs and all(r in role_p for r in rs) else 0.0

    mag_fns = {fn for f in mag_genome["features"] for fn in (f.get("functions") or []) if fn and fn.strip()}
    mag_roles = set().union(*(roles(fn) for fn in mag_fns)) if mag_fns else set()

    def observed(fn):
        if match == "exact":
            return fn in mag_fns
        rs = roles(fn)
        return bool(rs) and rs <= mag_roles

    out = copy.deepcopy({k: v for k, v in mag_genome.items() if k not in ("features", "cdss", "mrnas")})
    out.update(id=f"{mag_id}.prfba_merged", features=[], cdss=[], mrnas=[], source="prFBA-MAG-merge",
               warnings=list(mag_genome.get("warnings", [])) +
               ["MAG RAST features at probability 1 + reference-inferred functions (prFBA MAG arm)"])
    cds_by_id = {c["id"]: c for c in mag_genome.get("cdss", [])}
    counts = collections.Counter()
    for f in mag_genome["features"]:
        fns = [fn for fn in (f.get("functions") or []) if fn and fn.strip()]
        nf = copy.deepcopy(f)
        # KBase {ontology: {term: [event]}} -> {ontology: [terms]}, the list form the loader expects
        nf["ontology_terms"] = {k: (sorted(v) if isinstance(v, dict) else v)
                                for k, v in (f.get("ontology_terms") or {}).items()}
        sup = max((ref_support(fn) for fn in fns), default=0.0)
        nf["probability"] = 1.0
        nf["p_reference"] = sup
        nf["provenance"] = ("observed_in_mag_unannotated" if not fns else
                            "observed_in_mag+reference" if sup > 0 else "observed_in_mag")
        out["features"].append(nf)
        for cid in f.get("cdss", []):
            if cid in cds_by_id:
                c = copy.deepcopy(cds_by_id[cid])
                c.update(probability=1.0, provenance=nf["provenance"], ontology_terms=nf["ontology_terms"])
                out["cdss"].append(c)
        counts[nf["provenance"]] += 1
    n = 0
    for fn in sorted(ref):
        if is_hypo(fn) or observed(fn):
            continue
        n += 1
        fid = f"{mag_id}.ref_{n}"
        p = ref[fn]["p"]
        nf = {"id": fid, "type": "gene", "functions": [fn], "aliases": ref[fn]["aliases"],
              "location": [], "dna_sequence": "", "protein_translation": "",
              "ontology_terms": {"RAST": [fn]}, "cdss": [f"{fid}.CDS"],
              "probability": round(p * inferred_weight, 6), "p_reference": p,
              "provenance": "inferred_from_reference", "reference_genomes": carriers.get(fn, [])}
        out["features"].append(nf)
        c = dict(nf, id=f"{fid}.CDS", type="CDS", parent_gene=fid)
        c.pop("cdss", None)
        out["cdss"].append(c)
        counts["inferred_from_reference"] += 1
    out["feature_counts"] = {"gene": len(out["features"]), "CDS": len(out["cdss"])}
    return out, dict(counts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    cache_dir = os.path.join(out, "genome_cache")
    build_dir = os.path.join(out, "build")
    merged_dir = os.path.join(out, "merged")
    for d in (cache_dir, build_dir, merged_dir, os.path.join(out, "compare")):
        os.makedirs(d, exist_ok=True)
    res = {"timings_s": {}, "probes": {}, "mags": {}}
    T = res["timings_s"]

    sel = json.load(open(a.selection))
    sel.pop("_meta", None)
    hits = json.load(open(a.hits))

    # ---- probe A: the 16S default knobs (what build_synthetic_genomes*.py --hits apply) ----
    aud16, sum16 = S.build_audit(hits)                    # knobs=None -> DEFAULT_KNOBS (16S)
    audani, sumani = S.build_audit(hits, knobs=json.load(open(f"{PRFBA}/ani_knobs.json")))
    res["probes"]["selection_16S_default_knobs"] = {"summary": sum16, "per_mag": {
        m: {"tier": r["tier"], "n_selected": r["n_selected"], "flags": r["flags"]} for m, r in aud16.items()}}
    res["probes"]["selection_ani_knobs_summary"] = sumani

    # ---- step 3: accession -> BV-BRC genome_id ----
    t0 = time.perf_counter()
    acc = M.load_acc_map(ACC_MAP)
    T["load_acc_map"] = round(time.perf_counter() - t0, 2)
    res["acc_map_rows"] = len(acc)

    # ---- probe B: unchanged _api_json (no User-Agent) ----
    t0 = time.perf_counter()
    try:
        BSG._api_json("genome", "eq(genome_id,1293054.3)&select(genome_id)&limit(1)", retries=1)
        res["probes"]["api_json_no_UA"] = "OK (no 403)"
    except urllib.error.HTTPError as e:
        res["probes"]["api_json_no_UA"] = f"HTTPError {e.code}: {e.read()[:40]!r}"
    except Exception as e:  # noqa: BLE001
        res["probes"]["api_json_no_UA"] = f"{type(e).__name__}: {e}"
    T["probe_no_UA"] = round(time.perf_counter() - t0, 2)

    # install a browser-like UA for every later urllib call (_api_json uses bare urlopen)
    opener = urllib.request.build_opener()
    opener.addheaders = [("User-Agent", UA)]
    urllib.request.install_opener(opener)

    # ---- probe C: the fast builder fed an UNMAPPED GTDB accession (what the shipped build does) ----
    t0 = time.perf_counter()
    g_acc = BSG.fast_build_genome("GCA_020055655.1")      # direct call: nothing is cached
    res["probes"]["fast_build_genome_on_GTDB_accession"] = {
        "genome_id": "GCA_020055655.1", "n_features": len(g_acc["features"]),
        "scientific_name": g_acc["scientific_name"], "raised": False}
    T["probe_accession_build"] = round(time.perf_counter() - t0, 2)

    BSG._patch_taxonomy_viz()
    for mag, rec in sel.items():
        R = res["mags"][mag] = {"tier": rec["tier"], "n_selected": rec["n_selected"], "selected": []}
        if not rec["n_selected"]:
            R["status"] = "abstain"
            continue
        rec_m = copy.deepcopy(rec)
        for s in rec_m["selected"]:
            row = acc.get(s["genome_id"]) or {}
            info = {k: row.get(k, "") for k in (
                "gtdb_accession", "bvbrc_genome_id", "match_type", "gtdb_ncbi_taxid", "ncbi_genome_category",
                "genome_status", "length_ratio", "bvbrc_label_taxid_consistent", "taxid_consistent",
                "in_pgfam_cache", "cds", "candidates")}
            info.update(selected_accession=s["genome_id"], identity=s["identity"], species=s["species"])
            R["selected"].append(info)
            s["gtdb_accession"] = s["genome_id"]
            s["genome_id"] = row.get("bvbrc_genome_id") or s["genome_id"]
        gids = [s["genome_id"] for s in rec_m["selected"]]

        # step 4a: fast builder (shipped code, UA patched) -> cached genome dict in scratch
        for gid in gids:
            t0 = time.perf_counter()
            g = BSG.build_genome(gid, cache_dir, "fast")
            dt = time.perf_counter() - t0
            feats = g["features"]
            pat = [f for f in feats if dict((x[0], x[1]) for x in f.get("aliases", [])).get("PATRIC_id")]
            R.setdefault("fetch", {})[gid] = {
                "seconds": round(dt, 2), "scientific_name": g["scientific_name"], "taxonomy": g["taxonomy"],
                "n_features_with_product": len(feats), "n_PATRIC_features": len(pat),
                "n_RefSeq_features_no_patric_id": len(feats) - len(pat),
                "n_distinct_functions_all": len({fn for f in feats for fn in f["functions"]}),
                "n_distinct_functions_PATRIC": len({fn for f in pat for fn in f["functions"]}),
                "n_pgfams_from_features": len({f["pgfam_id"] for f in feats if f.get("pgfam_id")}),
            }
            t0 = time.perf_counter()
            pg = PF.fetch_pgfams(gid)
            R["fetch"][gid]["n_pgfams_fetch_pgfams"] = len(pg)
            R["fetch"][gid]["fetch_pgfams_seconds"] = round(time.perf_counter() - t0, 2)

        # step 5: shipped merge (cwd = build dir: create_synthetic_genome writes genome_objects/ there)
        here = os.getcwd()
        os.chdir(build_dir)
        for d in ("genome_objects", "taxonomies", "provenance", "ASVset_taxonomies"):
            os.makedirs(d, exist_ok=True)
        t0 = time.perf_counter()
        entry = BSG.build_for_asv(mag, rec_m, genome_cache_dir=cache_dir, builder="fast")
        T[f"build_for_asv:{mag}"] = round(time.perf_counter() - t0, 2)
        synth = json.load(open(f"{mag}.json"))
        os.chdir(here)
        R["build_entry"] = {k: entry.get(k) for k in ("status", "n_genomes", "genomes", "n_features",
                                                         "scientific_name", "gene_union")}
        R["synthetic_probability_counts"] = dict(collections.Counter(
            f.get("probability") for f in synth["features"]))

        # comparison with the MAG's own RAST annotation
        mag_g, mag_counts = mag_rast_functions(mag)
        mag_funcs = set(mag_counts)
        src = json.load(open(os.path.join(cache_dir, f"{gids[0]}.json")))
        ref_all = {fn for f in synth["features"] for fn in f["functions"]}
        ref_pat = {fn for f in src["features"] if dict((x[0], x[1]) for x in f.get("aliases", [])).get("PATRIC_id")
                   for fn in f["functions"]}
        R["mag_rast"] = {"n_features": len(mag_g["features"]),
                         "n_features_with_function": sum(1 for f in mag_g["features"] if f.get("functions")),
                         "ontology_events": mag_g.get("ontology_events")}
        for label, rf in (("reference_as_built (PATRIC+RefSeq products)", ref_all),
                          ("reference_PATRIC_only (RASTtk products)", ref_pat)):
            stats, ab_ex, ab_role, ab_len = compare_sets(rf, mag_funcs)
            mk_ref, mk_mag = marker_hits(rf), marker_hits(mag_funcs)
            mk = {}
            for name in MARKERS:
                mk[name] = {"in_reference": mk_ref[name], "in_mag_rast": mk_mag[name],
                            "absent_from_mag_exact": sorted(set(mk_ref[name]) & ab_ex),
                            "absent_from_mag_role": sorted(set(mk_ref[name]) & ab_role),
                            "absent_from_mag_role_lenient": sorted(set(mk_ref[name]) & ab_len)}
            # absent functions that exist only as RefSeq (non-PATRIC) product strings
            stats["n_absent_exact_refseq_only_strings"] = len(ab_ex - ref_pat)
            stats["n_absent_role_refseq_only_strings"] = len(ab_role - ref_pat)
            stats["n_absent_role_lenient_refseq_only_strings"] = len(ab_len - ref_pat)
            R[label] = {"stats": stats, "markers": mk}
            tag = "built" if label.startswith("reference_as") else "patric"
            with open(os.path.join(out, "compare", f"{mag}.absent_{tag}.json"), "w") as fh:
                json.dump({"absent_exact": sorted(ab_ex), "absent_role": sorted(ab_role),
                           "absent_role_lenient": sorted(ab_len)}, fh, indent=1)

        # prototype merge (role-level match; inferred weight left at 1.0 = shipped k/n)
        t0 = time.perf_counter()
        merged, counts = merge_mag_with_reference(mag_g, synth, [json.load(open(os.path.join(cache_dir, f"{g}.json"))) for g in gids], mag, inferred_weight=1.0, match="role")
        T[f"merge_prototype:{mag}"] = round(time.perf_counter() - t0, 2)
        with open(os.path.join(merged_dir, f"{mag}.prfba_merged.json"), "w") as fh:
            json.dump(merged, fh)
        R["merge_prototype_counts"] = counts
        R["merge_prototype_probability_counts"] = dict(collections.Counter(
            f.get("probability") for f in merged["features"]))

    with open(os.path.join(out, "compare", "smoke_results.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res["timings_s"], indent=1))
    print(json.dumps(res["probes"], indent=1)[:3000])


if __name__ == "__main__":
    main()
