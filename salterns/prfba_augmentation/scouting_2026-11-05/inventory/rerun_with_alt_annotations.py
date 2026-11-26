#!/usr/bin/env python
"""Re-run the annotation-dependent salterns R232 analyses against an alternative annotation source,
writing ONLY to --out (an audit hook refuses any write under the protected repositories).

Annotation source (--annotations):
  * a directory holding one KBase-genome JSON per MAG named <mag>.contigs__.RAST.json (the layout of
    jplfaria/Cliff MAGs/RASTified_MAG_json; a list-wrapped genome as in its converted/ subfolder also works), or
  * a JSON file {name: [function strings]}; names may be the 310 'mag' strings, iterativeID_r232 labels, or any
    form key() understands (.contigs, .contigs__, bins./bins_, zero-padded maxbin numbers).
  --add-functions FILE unions {name: [function strings]} on top of --annotations (e.g. prFBA-added functions).
Parser (--parser):
  trunc  raw-text regex of rerun_correlational_r232.build_guilds (lines 107, 112): a function string containing ']'
         and every later string of the same list are dropped. Only valid for genome-JSON directories.
  json   json.load of features[].functions + cdss[].functions (identical, on the 310 RAST genomes, to the
         string-aware parser of make_guild_membership_comparison.rast_functions).
Parts (--parts, comma list):
  guilds   -> <out>/guilds_r232.csv                       (same columns and rules as build_guilds)
  corr     -> guild_methane_correlations_r232.csv, methanogen_archaea_by_core_r232.csv,
              mag_methane_correlations_r232.csv, guild_methane_{raw,dnanorm}_r232.png
  s11      -> figureS11_r232.png / .order.json            (figureS11_r232.py, run patched; needs 'guilds')
  membership -> guild_membership_2020_vs_r232.csv / .md   (make_guild_membership_comparison.py, run patched;
              its count asserts at lines 265-268 and 371 are kept, so it stops if membership changes)
  archaeal -> archaeal_mepn_and_tma_candidates.md / _genes.csv (make_archaeal_mepn_tma_report.py, patched; targets
              from <out>/guilds_r232.csv but genes/coordinates/proteins from the ORIGINAL RAST JSONs + MAG_contigs)
  modelling -> modelling_consistent_basis_r232.json, figure7_r232.png/.svg/_stats.json, modelling_numbers_r232.json
              (clauded_insights scripts, patched; stored 2020 fluxes subset by <out>/guilds_r232.csv; asserts 17/9 kept)
  prefilter_abund -> prefilter_sensitivity/abund_* (needs 'corr' outputs already in --out)
Run:  ~/Documents/py_venv/bin/python -B rerun_with_alt_annotations.py --out DIR [--annotations SRC]
          [--add-functions FILE] [--parser json|trunc] [--parts guilds,corr,s11,membership]
"""
import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True                         # never write __pycache__ into the repo

SALTERNS = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns")
R232 = SALTERNS / "correlational/r232"
ORIG_RAST = SALTERNS / "jplfaria/Cliff MAGs/RASTified_MAG_json"
HOME = Path.home() / "Documents"
PROTECTED = [HOME / "MicrobiomeNotebooks", HOME / "prFBA", HOME / "codiffusion_bioreactor"] + \
            [p for p in HOME.glob("ModelSEEDpy*")]
# outputs of rerun_correlational_r232.py that do not depend on annotations; figureS11_r232.py reads two of them
STATIC_R232_INPUTS = ["C_within_r232.csv", "glasso_direct_edges_r232.csv"]


def _install_write_guard(out):
    out = out.resolve()

    def protected(p):
        try:
            rp = Path(os.fsdecode(p)).resolve()
        except Exception:
            return False
        if rp == out or out in rp.parents:
            return False
        return any(rp == q or q in rp.parents for q in PROTECTED)

    def hook(event, args):
        if event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
            mode, flags = (args[1] if len(args) > 1 else None), (args[2] if len(args) > 2 else 0)
            writing = (isinstance(mode, str) and any(c in mode for c in "wax+")) or \
                      (isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND))
            if writing and protected(args[0]):
                raise PermissionError(f"write guard: refusing to write {args[0]}")
        elif event in ("os.mkdir", "os.remove", "os.rename", "os.rmdir", "shutil.rmtree", "os.unlink", "os.symlink"):
            targets = {"os.symlink": args[1:2], "os.rename": args[:2]}.get(event, args[:1])   # symlink: dst only
            for a in targets:
                if isinstance(a, (str, bytes, os.PathLike)) and protected(a):
                    if event == "os.mkdir" and Path(os.fsdecode(a)).exists():
                        continue          # mkdir(exist_ok=True) on an existing folder changes nothing
                    raise PermissionError(f"write guard: refusing {event} on {a}")
    sys.addaudithook(hook)


def _load_rr(out):
    sys.path.insert(0, str(R232))
    import rerun_correlational_r232 as rr   # module level only reads inputs (lines 72-96)
    sys.path.remove(str(R232))
    rr.HERE = out                           # every writer in the module uses HERE / "<file>"
    return rr


def _functions_from_dir(d, mags, parser):
    trunc_block = re.compile(r'"functions"\s*:\s*\[(.*?)\]', re.S)           # build_guilds line 107
    out = {}
    for mag in mags:
        path = Path(d) / f"{mag}.contigs__.RAST.json"
        text = path.read_text(errors="ignore")
        if parser == "trunc":
            out[mag] = {s for blk in trunc_block.findall(text) for s in re.findall(r'"([^"]*)"', blk)}  # line 112
        else:
            g = json.loads(text)
            g = g[0] if isinstance(g, list) else g
            out[mag] = {s for sec in ("features", "cdss") for f in g.get(sec, []) for s in (f.get("functions") or [])}
    return out


def _functions_from_json(path, rr):
    raw = json.load(open(path))
    by_key = {rr.key(m): m for m in rr.MAGS}
    by_id = {rr.ids[m]: m for m in rr.MAGS}
    out, unknown = {}, []
    for name, funcs in raw.items():
        mag = name if name in rr.ids else by_id.get(name) or by_key.get(rr.key(name))
        if mag is None:
            unknown.append(name)
            continue
        out.setdefault(mag, set()).update(funcs)
    if unknown:
        raise SystemExit(f"{len(unknown)} names in {path} match none of the 310 MAGs, e.g. {unknown[:3]}")
    return out


def build_guilds_from(rr, funcs_by_mag):
    """rerun_correlational_r232.build_guilds (lines 106-146) with the per-MAG function set passed in."""
    import pandas as pd
    evidence = {g: {} for g in rr.GENES}
    for mag in rr.MAGS:
        funcs = funcs_by_mag.get(mag, set())
        for g, pat in rr.GENES.items():
            hits = sorted(s for s in funcs if re.search(pat, s.lower()))
            if hits:
                evidence[g][mag] = hits
    # ---- from here verbatim (lines 117-146) ----
    MAGS, S, key, taxa_of, ids, ids_2020, AB = rr.MAGS, rr.S, rr.key, rr.taxa_of, rr.ids, rr.ids_2020, rr.AB
    GENES, HERE = rr.GENES, rr.HERE
    mepn_2020, bet_2020 = rr.marimo_set("mepe"), rr.marimo_set("betaine")
    checkm = pd.read_csv(S / "MAG_QC/checkm/summary/CheckM_summary_table.tsv", sep="\t")
    checkm.index = checkm["Bin Name"].map(key)
    rows = []
    for mag in MAGS:
        t, k = taxa_of[mag], key(mag)
        rows.append({
            "mag": mag, "iterativeID_r232": ids[mag], "iterativeID_2020": ids_2020.get(mag, ""),
            "domain": t["Kingdom"], "phylum": t["Phylum"], "order": t["Order"], "family": t["Family"],
            "genus": t["Genus"], "species": t["Species"],
            "completeness": checkm.loc[k, "Completeness"], "contamination": checkm.loc[k, "Contamination"],
            "mepn_consumer": mag in evidence["phnJ"], "mepn_2020": k in mepn_2020,
            "betaine_reducer": mag in evidence["grdH"], "betaine_2020": k in bet_2020,
            "tma_producer": mag in evidence["cutC_cluster"] or mag in evidence["cntA"],
            "methanogen": mag in evidence["mcrA"], "archaea": t["Kingdom"] == "Archaea",
            "methanosarcinales": t["Order"] == "Methanosarcinales",
            **{f"ev_{g}": " | ".join(evidence[g].get(mag, [])) for g in GENES},
            "mean_rel_abund_12cores": AB["raw"].loc[mag].mean(),
        })
    g = pd.DataFrame(rows)
    g.to_csv(HERE / "guilds_r232.csv", index=False)
    print("guild sizes:", {c: int(g[c].sum()) for c in ("mepn_consumer", "mepn_2020", "betaine_reducer",
          "betaine_2020", "tma_producer", "methanogen", "archaea", "methanosarcinales")})
    return g


def run_patched(script, subs, workdir):
    """exec a module-level script with exact text substitutions (each must occur exactly once), cwd = workdir."""
    src = Path(script).read_text()
    for old, new in subs:
        n = src.count(old)
        if n != 1:
            raise SystemExit(f"{script}: expected 1 occurrence of {old!r}, found {n}")
        src = src.replace(old, new)
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        exec(compile(src, str(script), "exec"), {"__name__": "__main__", "__file__": str(script)})
    finally:
        os.chdir(cwd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--annotations", default=str(ORIG_RAST))
    ap.add_argument("--add-functions", type=Path)
    ap.add_argument("--parser", choices=("json", "trunc"), default="json")
    ap.add_argument("--parts", default="guilds,corr")
    a = ap.parse_args()
    out = a.out.resolve()
    if any(out == q or q in out.parents for q in PROTECTED):
        raise SystemExit(f"--out {out} is inside a protected repository")
    out.mkdir(parents=True, exist_ok=True)
    _install_write_guard(out)
    parts = set(a.parts.split(","))
    t0 = time.time()
    rr = _load_rr(out)
    timing = {"import_rerun_correlational_r232_s": round(time.time() - t0, 1)}
    if parts & {"guilds", "corr", "s11", "membership", "archaeal", "modelling", "prefilter_abund"}:
        t = time.time()
        src = Path(a.annotations)
        if src.is_dir():
            funcs = _functions_from_dir(src, rr.MAGS, a.parser)
        else:
            if a.parser == "trunc":
                raise SystemExit("--parser trunc needs a genome-JSON directory")
            funcs = _functions_from_json(src, rr)
        if a.add_functions:
            for mag, extra in _functions_from_json(a.add_functions, rr).items():
                funcs.setdefault(mag, set()).update(extra)
        missing = [m for m in rr.MAGS if m not in funcs]
        if missing:
            print(f"WARNING: {len(missing)} MAGs have no function set (treated as empty), e.g. {missing[:2]}")
        timing["read_annotations_s"] = round(time.time() - t, 1)
        t = time.time()
        g = build_guilds_from(rr, funcs)
        timing["guilds_s"] = round(time.time() - t, 1)
    if "corr" in parts:
        t = time.time()
        corr = rr.guild_correlations(g)
        rr.per_mag_correlations(g)
        for basis in rr.AB:
            rr.plot_guilds(g, corr, basis)
        timing["corr_s"] = round(time.time() - t, 1)
    if "s11" in parts:
        t = time.time()
        for f in STATIC_R232_INPUTS:                     # annotation-independent inputs, linked from the repo
            if not (out / f).exists():
                (out / f).symlink_to(R232 / f)
        run_patched(R232 / "figureS11_r232.py",
                    [('HERE = Path(__file__).resolve().parent', f'HERE = Path({str(out)!r})')], out)
        timing["s11_s"] = round(time.time() - t, 1)
    if "membership" in parts:
        t = time.time()
        run_patched(R232 / "make_guild_membership_comparison.py",
                    [('HERE = Path(__file__).resolve().parent', f'HERE = Path({str(out)!r})'),
                     ('COR = HERE.parent', f'COR = Path({str(R232.parent)!r})')], out)
        timing["membership_s"] = round(time.time() - t, 1)
    if "archaeal" in parts:      # gene-level report; needs coordinates + proteins, so it always reads ORIG_RAST
        t = time.time()
        run_patched(R232 / "make_archaeal_mepn_tma_report.py",
                    [('HERE = Path(__file__).resolve().parent', f'HERE = Path({str(out)!r})'),
                     ('COR = HERE.parent', f'COR = Path({str(R232.parent)!r})')], out)
        timing["archaeal_s"] = round(time.time() - t, 1)
    if "modelling" in parts:     # stored-flux subsets by guild + Figure 7 (clauded_insights)
        t = time.time()
        CI = SALTERNS / "clauded_insights"
        g_csv = repr(str(out / "guilds_r232.csv"))
        run_patched(CI / "modelling_consistent_basis_r232.py",
                    [('HERE = Path(__file__).resolve().parent', f'HERE = Path({str(out)!r})'),
                     ('S = HERE.parent', f'S = Path({str(SALTERNS)!r})'),
                     ('GUILDS = pd.read_csv(CORR / "r232" / "guilds_r232.csv")', f'GUILDS = pd.read_csv({g_csv})')], out)
        run_patched(CI / "figure7_r232.py",
                    [('HERE = Path(__file__).resolve().parent', f'HERE = Path({str(out)!r})')], out)
        run_patched(CI / "modelling_numbers_r232.py",
                    [('HERE = Path(__file__).resolve().parent', f'HERE = Path({str(out)!r})'),
                     ('S = HERE.parent', f'S = Path({str(SALTERNS)!r})'),
                     ('G = pd.read_csv(CORR / "r232" / "guilds_r232.csv")', f'G = pd.read_csv({g_csv})'),
                     ('json.load(open(HERE / "figure7_stats.json"))', f'json.load(open({str(CI / "figure7_stats.json")!r}))')],
                    out)
        timing["modelling_s"] = round(time.time() - t, 1)
    if "prefilter_abund" in parts:   # needs the corr outputs in --out; imports the already-patched rr module
        t = time.time()
        run_patched(R232 / "prefilter_sensitivity_abundance.py",
                    [('HERE = Path(__file__).resolve().parent', f'HERE = Path({str(out)!r})')], out)
        timing["prefilter_abund_s"] = round(time.time() - t, 1)
    json.dump({"annotations": a.annotations, "add_functions": str(a.add_functions or ""), "parser": a.parser,
               "parts": sorted(parts), "timing": timing}, open(out / "run_manifest.json", "w"), indent=1)
    print("timing:", timing)


if __name__ == "__main__":
    main()
