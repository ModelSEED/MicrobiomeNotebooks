#!/usr/bin/env python3
"""Determine whether MAGs derivable from the wetland samples in
Sample_metadata_wetland_only.csv are present in SPIRE and/or MGnify.

The CSV is keyed by SRA run (SRR...) and BioProject (PRJNA...), not by MAG.
- SPIRE indexes MAGs by INSDC *biosample* (derived_from_sample), so we first
  crosswalk run -> biosample via the ENA Portal API, then join against SPIRE's
  bulk genome-metadata table.
- MGnify is queried per BioProject via its REST API (the accession= filter
  matches the BioProject form).

All network results are cached under ./db_membership/cache so the script is
cheap to re-run / resume after interruption.
"""
import csv
import gzip
import io
import json
import os
import sys
import time
import urllib.request
import urllib.error
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(HERE, "Sample_metadata_wetland_only.csv")
OUT_DIR = os.path.join(HERE, "db_membership")
CACHE_DIR = os.path.join(OUT_DIR, "cache")
SPIRE_GZ = os.path.join(CACHE_DIR, "spire_v1_genome_metadata.tsv.gz")
ENA_CACHE = os.path.join(CACHE_DIR, "ena_bioproject.json")
MGNIFY_CACHE = os.path.join(CACHE_DIR, "mgnify_studies.json")
MGNIFY_RUN_CACHE = os.path.join(CACHE_DIR, "mgnify_runs.json")

SPIRE_URL = "https://black.embl.de/~fullam/spire/metadata/spire_v1_genome_metadata.tsv.gz"
ENA_URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
           "?accession={acc}&result=read_run"
           "&fields=run_accession,sample_accession,secondary_sample_accession,study_accession")
MGNIFY_URL = "https://www.ebi.ac.uk/metagenomics/api/v1/studies?accession={acc}"
MGNIFY_RUN_URL = "https://www.ebi.ac.uk/metagenomics/api/v1/runs/{acc}"
MGNIFY_RUN_ANALYSES_URL = "https://www.ebi.ac.uk/metagenomics/api/v1/runs/{acc}/analyses"

os.makedirs(CACHE_DIR, exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def http_get(url, accept=None, retries=4, timeout=60):
    """GET with simple exponential backoff. Returns bytes or raises."""
    last = None
    url = urllib.parse.quote(url, safe=":/?&=,~%")  # tolerate stray chars
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "wetland-mag-check/1.0")
            if accept:
                req.add_header("Accept", accept)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None  # legitimate "not found"
            last = e
        except Exception as e:  # noqa: BLE001
            last = e
        time.sleep(1.5 * (attempt + 1))
    raise last


# ---------------------------------------------------------------------------
# 1. Parse the CSV
# ---------------------------------------------------------------------------
def parse_csv():
    runs = []           # (run, bioproject)
    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            run = (row.get("SRA run") or "").strip()
            proj = (row.get("Bioproject") or "").strip()
            if run:
                runs.append((run, proj))
    return runs


# ---------------------------------------------------------------------------
# 2. Download + index SPIRE genome metadata -> sample -> [genome rows]
# ---------------------------------------------------------------------------
def load_spire_index():
    if not os.path.exists(SPIRE_GZ):
        log(f"Downloading SPIRE genome metadata (~92 MB) -> {SPIRE_GZ}")
        data = http_get(SPIRE_URL, timeout=600)
        with open(SPIRE_GZ, "wb") as fh:
            fh.write(data)
    log("Indexing SPIRE genome metadata by derived_from_sample ...")
    sample_to_mags = {}
    with gzip.open(SPIRE_GZ, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        n = 0
        for row in reader:
            n += 1
            samp = row.get("derived_from_sample", "")
            if not samp:
                continue
            sample_to_mags.setdefault(samp, []).append(
                (row.get("genome_id", ""), row.get("classification", ""),
                 row.get("completeness", ""), row.get("contamination", ""))
            )
    log(f"  SPIRE: indexed {n} MAG rows across {len(sample_to_mags)} source samples")
    return sample_to_mags


# ---------------------------------------------------------------------------
# 3. ENA crosswalk: bioproject -> {run: {biosample, secondary, study}}
# ---------------------------------------------------------------------------
def load_cache(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def save_cache(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def ena_crosswalk(bioprojects):
    cache = load_cache(ENA_CACHE)
    todo = [b for b in bioprojects if b and b not in cache]
    log(f"ENA crosswalk: {len(cache)} cached, {len(todo)} bioprojects to fetch")
    for i, proj in enumerate(todo, 1):
        raw = http_get(ENA_URL.format(acc=proj))
        entry = {}
        if raw:
            text = raw.decode("utf-8", "replace").splitlines()
            if len(text) > 1:
                hdr = text[0].split("\t")
                idx = {h: j for j, h in enumerate(hdr)}
                for line in text[1:]:
                    parts = line.split("\t")
                    run = parts[idx["run_accession"]]
                    entry[run] = {
                        "biosample": parts[idx["sample_accession"]],
                        "secondary": parts[idx["secondary_sample_accession"]],
                        "study": parts[idx["study_accession"]],
                    }
        cache[proj] = entry
        if i % 25 == 0:
            log(f"  ENA {i}/{len(todo)} ...")
            save_cache(ENA_CACHE, cache)
        time.sleep(0.15)
    save_cache(ENA_CACHE, cache)
    return cache


# ---------------------------------------------------------------------------
# 4. MGnify: bioproject -> {present, mgys, samples_count}
# ---------------------------------------------------------------------------
def mgnify_lookup(bioprojects):
    cache = load_cache(MGNIFY_CACHE)
    todo = [b for b in bioprojects if b and b not in cache]
    log(f"MGnify: {len(cache)} cached, {len(todo)} bioprojects to fetch")
    for i, proj in enumerate(todo, 1):
        raw = http_get(MGNIFY_URL.format(acc=proj), accept="application/json")
        rec = {"present": False, "mgys": "", "samples_count": 0}
        if raw:
            try:
                d = json.loads(raw)
                data = d.get("data", [])
                if data:
                    rec["present"] = True
                    rec["mgys"] = ";".join(x["id"] for x in data)
                    rec["samples_count"] = sum(
                        int(x["attributes"].get("samples-count") or 0) for x in data)
            except Exception:  # noqa: BLE001
                pass
        cache[proj] = rec
        if i % 25 == 0:
            log(f"  MGnify {i}/{len(todo)} ...")
            save_cache(MGNIFY_CACHE, cache)
        time.sleep(0.15)
    save_cache(MGNIFY_CACHE, cache)
    return cache


def mgnify_run_lookup(run_accessions):
    """Direct per-run check: is the run in MGnify, and what analysis types?
    404 => not in MGnify. If present, fetch analyses to report experiment types
    (assembly analyses are what yield MAGs)."""
    cache = load_cache(MGNIFY_RUN_CACHE)
    todo = [r for r in run_accessions if r and r not in cache]
    log(f"MGnify per-run: {len(cache)} cached, {len(todo)} runs to fetch")
    for i, run in enumerate(todo, 1):
        raw = http_get(MGNIFY_RUN_URL.format(acc=run), accept="application/json")
        rec = {"present": False, "analysis_types": ""}
        if raw:
            rec["present"] = True
            araw = http_get(MGNIFY_RUN_ANALYSES_URL.format(acc=run),
                            accept="application/json")
            if araw:
                try:
                    ad = json.loads(araw)
                    types = sorted({a["attributes"].get("experiment-type", "")
                                    for a in ad.get("data", [])})
                    rec["analysis_types"] = ";".join(t for t in types if t)
                except Exception:  # noqa: BLE001
                    pass
        cache[run] = rec
        if i % 50 == 0:
            log(f"  MGnify per-run {i}/{len(todo)} ...")
            save_cache(MGNIFY_RUN_CACHE, cache)
        time.sleep(0.12)
    save_cache(MGNIFY_RUN_CACHE, cache)
    return cache


# ---------------------------------------------------------------------------
# 5. Merge + write results
# ---------------------------------------------------------------------------
def main():
    runs = parse_csv()
    # Only real INSDC BioProject accessions are queryable; rows tagged
    # "This study" carry non-INSDC run IDs (newly sequenced) and cannot be
    # matched by accession.
    bioprojects = sorted({p for _, p in runs if p.startswith("PRJ")})
    skipped = sorted({p for _, p in runs if p and not p.startswith("PRJ")})
    log(f"CSV: {len(runs)} run rows, {len(bioprojects)} INSDC bioprojects "
        f"({len(skipped)} non-INSDC project labels skipped: {skipped})")

    spire = load_spire_index()
    ena = ena_crosswalk(bioprojects)
    mgnify = mgnify_lookup(bioprojects)
    all_runs = sorted({r for r, _ in runs if r.startswith(("SRR", "ERR", "DRR"))})
    mgnify_runs = mgnify_run_lookup(all_runs)

    # run -> biosample lookup from ENA cache
    run_to_bs = {}
    for proj, runs_d in ena.items():
        for run, info in runs_d.items():
            run_to_bs[run] = info

    out_path = os.path.join(OUT_DIR, "mag_database_membership.csv")
    spire_hit_runs = mgnify_study_hit_runs = mgnify_run_hit = 0
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sra_run", "bioproject", "biosample",
                    "spire_mag_count", "spire_mag_ids", "spire_match_key",
                    "mgnify_study_present", "mgnify_study", "mgnify_samples_count",
                    "mgnify_run_present", "mgnify_run_analysis_types"])
        for run, proj in runs:
            info = run_to_bs.get(run, {})
            bs = info.get("biosample", "")
            # SPIRE: try biosample, secondary, then raw run id
            match_key, mags = "", []
            for key in (bs, info.get("secondary", ""), run):
                if key and key in spire:
                    match_key, mags = key, spire[key]
                    break
            mg = mgnify.get(proj, {})
            mr = mgnify_runs.get(run, {})
            if mags:
                spire_hit_runs += 1
            if mg.get("present"):
                mgnify_study_hit_runs += 1
            if mr.get("present"):
                mgnify_run_hit += 1
            w.writerow([
                run, proj, bs,
                len(mags), ";".join(m[0] for m in mags), match_key,
                mg.get("present", False), mg.get("mgys", ""),
                mg.get("samples_count", 0),
                mr.get("present", False), mr.get("analysis_types", ""),
            ])

    log("=" * 60)
    log(f"RESULTS written to {out_path}")
    log(f"  runs total:                   {len(runs)}")
    log(f"  runs with SPIRE MAG(s):       {spire_hit_runs}")
    log(f"  runs whose BioProject is in MGnify: {mgnify_study_hit_runs}")
    log(f"  runs directly present in MGnify:    {mgnify_run_hit}")
    bp_in_mgnify = sum(1 for b in bioprojects if mgnify.get(b, {}).get("present"))
    log(f"  bioprojects in MGnify:        {bp_in_mgnify}/{len(bioprojects)}")
    log("=" * 60)


if __name__ == "__main__":
    sys.exit(main())
