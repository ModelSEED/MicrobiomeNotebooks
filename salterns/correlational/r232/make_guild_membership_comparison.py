#!/usr/bin/env python
"""Before/after membership of the methane-relevant guilds: the 2020 modelling sets vs the GTDB R232 gene-based guilds.

One row per (guild, MAG) for every MAG in
    MePn consumers     2020: the 27 hand-curated MAGs (`mepe` in ../analysis_2_marimo.py)
                       R232: phnJ carriers (C-P lyase, EC 4.7.1.1)
    Betaine reducers   2020: the 15 MAGs with betaine models (`betaine` in ../analysis_2_marimo.py)
                       R232: grdH carriers (betaine reductase component B, EC 1.21.4.4)
    Methanogen         R232 only: mcrA carriers                               (no 2020 set)
    TMA-route cand.    R232 only: cntA, or genes "clustered with choline TMA-lyase" (no 2020 set)

Membership, taxonomy, CheckM and the RAST evidence strings come from guilds_r232.csv (rerun_correlational_r232.py), so
the two files cannot disagree. Extra RAST strings (shared Grd components, methylamine methyl-transfer genes, other phn
and phosphonate-transporter genes) are read from the same RAST JSONs. The parser here is a string-aware version of the
one in rerun_correlational_r232.py: that one stops each "functions" list at the first ']', so it drops strings such as
"[Trimethylamine--corrinoid protein] Co-methyltransferase (EC 2.1.1.250)". Checked over all 310 MAGs: the truncation
changes none of the guilds_r232.csv marker hits (phnJ, phn accessory, grdH, cutC cluster, cntA, mcrA, mtbB_mtmB) and
none of the shared-Grd hits; in this report it matters only for Methanolobus.1's methylamine Co-methyltransferases.
The betaine TMA FVA regimes (fixed by the optimum or not) are read from
../betaine_grdH_r232/fva_tma_by_core_r232.csv (rerun_betaine_grdH_r232.py).

"Modelled in 2020" means the MAG has a key in ../saltern_fluxes.json ['0.07']['methane'] (MePn) or in
../saltern_fluxes_betaine.json ['TMA'] (betaine; also checked for the TMA candidates). The mean modelled flux is the
mean over the six unrestored cores of the FBA exchange flux exactly as stored in those files (EX_cpd01024_e0 for CH4,
EX_cpd00441_e0 for TMA, biomass fixed at the DNA-normalised abundance); the marimo converts to areal units only later.
Why a 2020 MePn member has no flux is read from the recorded output of ../analysis(2).ipynb (the model-loading cell
prints "skip <MAG>"; the growth check prints "<MAG> doesn't grow").

Relative abundance per site type = plain mean of the 24 samples grouped by site (unrestored R1+R2: 12 samples;
restored SF2: 6; reference R2A: 6), in % of reads (same definition as Table S5).

Writes guild_membership_2020_vs_r232.csv and guild_membership_2020_vs_r232.md next to this script.
Run:  ~/Documents/py_venv/bin/python make_guild_membership_comparison.py
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
COR = HERE.parent
S = COR.parent
RAST = S / "jplfaria/Cliff MAGs/RASTified_MAG_json"
PAT = re.compile(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool_bins[._](concoct_out|metabat|maxbin)\.(\d+)")
SITES = ("unrestored", "restored", "reference")
EXTRA = {   # supporting RAST strings beyond the marker genes already in guilds_r232.csv
    "grd_shared": r"glycine/sarcosine/betaine reductase|glycine reductase component b|^grdx protein",
    # Methylamine methyl-transfer: the Co-methyltransferases MttB/MtbB/MtmB (EC 2.1.1.250/249/248), their corrinoid
    # proteins, and the methyl-corrinoid:CoM methyltransferases (EC 2.1.1.247, MA0779). guilds_r232.csv's ev_mtbB_mtmB
    # column holds only the di-/monomethylamine corrinoid proteins, not the B-subunit methyltransferases it is named after.
    "methylamine_corrinoid": r"methylamine--corrinoid protein\] co-methyltransferase|ec 2\.1\.1\.2(47|48|49|50)\b"
                             r"|(trimethylamine|dimethylamine|monomethylamine) methyltransferase corrinoid protein"
                             r"|methylcobalamin:coenzyme m methyltransferase",
    # phn genes outside the phn_accessory pattern (PhnP, PhnF, ...) and phosphonate ABC transporter components
    "phn_other": r"\bphn[a-z]\d?\b|phosphonate.*(abc|transport)|(abc|transport).*phosphonate",
}
PHN_ACC = r"methylphosphonate 5-triphosphate (synthase|diphosphatase)|phnn|ec 2\.7\.4\.23"   # = phn_accessory, guilds_r232
PHN_LIST = "other C-P lyase pathway genes searched (PhnG/H/I/L/M/N)"


def key(name):
    """Bin identity that ignores '.contigs', 'bins.' vs 'bins_' and zero padding (maxbin.047)."""
    m = PAT.search(str(name))
    return f"{m[1]}|{m[2]}|{int(m[3])}" if m else None


def site_type(sample):
    return "restored" if sample.startswith("SF2") else "reference" if sample.startswith("R2A") else "unrestored"


def rast_functions(mag):
    text = (RAST / f"{mag}.contigs__.RAST.json").read_text(errors="ignore")
    # string-aware: a ']' inside a function string (e.g. "[Methylamine--corrinoid protein] ...") does not end the list
    blocks = re.findall(r'"functions"\s*:\s*\[((?:\s*"(?:[^"\\]|\\.)*"\s*,?)*)\s*\]', text, re.S)
    return {s for blk in blocks for s in re.findall(r'"((?:[^"\\]|\\.)*)"', blk)}


def split(ev):
    return [x for x in str(ev).split(" | ") if x and x != "nan"]


def phn_short(s):
    m = re.search(r"\b(Phn[A-Z])\b", s)
    return m[1] if m else "PhnM [EC 3.6.1.63]" if "EC 3.6.1.63" in s else s


def phn_other_short(strings):
    """Phn symbols outside the PhnG/H/I/L/M/N search, plus phosphonate ABC transporter genes without a symbol."""
    other = [s for s in strings if not re.search(PHN_ACC, s.lower())]
    syms = sorted({"Phn" + m[1].upper() for s in other for m in [re.search(r"\bphn([a-z]\d?)\b", s, re.I)] if m})
    tr = any(not re.search(r"\bphn[a-z]\d?\b", s, re.I) for s in other)
    return syms + (["(methionine/)phosphonate ABC transporter genes"] if tr else [])


def methyl_short(strings):
    """Methylamine methyl-transfer strings, with the '@ pyrrolysine-containing' variants folded into one entry."""
    base = sorted({s.replace(" @ pyrrolysine-containing", "") for s in strings})
    pyl = any("pyrrolysine-containing" in s for s in strings)
    return "; ".join(base) + (" (pyrrolysine-containing copies included)" if pyl else "")


def grd_components(strings):
    low = " ".join(strings).lower()
    named = [("protein A", "reductase protein a"), ("component C", "component c chain"), ("GrdX", "grdx protein"),
             ("glycine-reductase component B (EC 1.21.4.2)", "ec 1.21.4.2")]
    return [lab for lab, pat in named if pat in low]


def sortkey(i):
    base, _, n = str(i).rpartition(".")
    return (base, int(n)) if n.isdigit() else (str(i), 0)


# ---- inputs ------------------------------------------------------------------------------------------------------
g = pd.read_csv(HERE / "guilds_r232.csv")
g["iterativeID_2020"] = g["iterativeID_2020"].fillna("")
for c in [c for c in g.columns if c.startswith("ev_")] + ["species", "genus", "family"]:
    g[c] = g[c].fillna("")
g["key"] = g["mag"].map(key)
qc = pd.read_csv(S / "MAG_QC/comparison_vs_2020.csv").set_index("mag_key")
ids_cmp = pd.read_csv(COR / "iterativeIDs_r232_vs_2020.csv").set_index("mag")

rel = pd.read_csv(S / "data/Cliff_310MAG_relabund.tsv", sep="\t", index_col=0)
rel = rel.loc[[i for i in rel.index if key(i)]]                 # drops the unbinned_reads row
assert len(rel) == 310
rel.index = rel.index.map(key)
site_mean = rel.T.groupby(rel.columns.map(site_type)).mean().T * 100          # % of reads

flux_mepn = json.load(open(COR / "saltern_fluxes.json"))
flux_bet = json.load(open(COR / "saltern_fluxes_betaine.json"))
FLUX = {   # name -> {MAG key: {core: flux}}
    "CH4_0.07": {key(k): v for k, v in flux_mepn["0.07"]["methane"].items()},
    "CH4_0.33": {key(k): v for k, v in flux_mepn["0.33"]["methane"].items()},
    "TMA": {key(k): v for k, v in flux_bet["TMA"].items()},
}
assert all(len(v) == 6 for d in FLUX.values() for v in d.values()), "expected 6 unrestored cores per modelled MAG"

mepn_model_files = {key(p.name) for p in (S / "models/MePn_consumers").glob("*.json")}
bet_model_files = {key(p.name): p for p in (COR / "models").glob("*_gapfilled_betaine.json")}
nb_skip, nb_nogrow = set(), set()                         # recorded output of the 2020 notebook
for cell in json.load(open(COR / "analysis(2).ipynb"))["cells"]:
    for out in cell.get("outputs", []):
        text = "".join(out.get("text", []))
        nb_skip |= {key(m) for m in re.findall(r"^skip (\S+)", text, re.M) if key(m)}
        nb_nogrow |= {key(m) for m in re.findall(r"^(\S+) doesn't grow", text, re.M)}

# betaine reductase reaction (rxn17220_c0) gene association in each betaine model
bet_gpr, bet_in_base = {}, {}
for k, p in bet_model_files.items():
    rx = {r["id"]: r for r in json.load(open(p))["reactions"]}
    bet_gpr[k] = rx["rxn17220_c0"].get("gene_reaction_rule", "") if "rxn17220_c0" in rx else None
    base = Path(str(p).replace("_gapfilled_betaine", "_betaine"))
    bet_in_base[k] = base.exists() and any(r["id"] == "rxn17220_c0" for r in json.load(open(base))["reactions"])

# TMA FVA at the lactate optimum, per MAG and core, from the betaine re-run (../betaine_grdH_r232/)
fva = pd.read_csv(COR / "betaine_grdH_r232/fva_tma_by_core_r232.csv")
fva["key"] = fva["mag"].map(key)
assert set(fva.key) == set(FLUX["TMA"]) and (fva.groupby("key").size() == 6).all()
assert all(abs(r.TMA_existing - FLUX["TMA"][r.key][r.core]) < 1e-9 for r in fva.itertuples()), "FVA file vs stored TMA"
fva_by = fva.groupby("key").agg(n_alt=("regime", lambda s: int((s == "degenerate_alternative_optima").sum())),
                                forced0=("regime", lambda s: bool((s == "zero_forced").all())),
                                fva_min_max=("TMA_fva_min", "max"))

# ---- rows --------------------------------------------------------------------------------------------------------
GUILDS = [   # guild label, 2020 column, R232 column
    ("MePn consumer", "mepn_2020", "mepn_consumer"),
    ("Betaine reducer", "betaine_2020", "betaine_reducer"),
    ("Methanogen", None, "methanogen"),
    ("TMA-route candidate", None, "tma_producer"),
]
FLAGS = {"mepn_2020": "MePn 2020 set", "mepn_consumer": "MePn (phnJ)", "betaine_2020": "betaine 2020 set",
         "betaine_reducer": "betaine (grdH)", "methanogen": "methanogen (mcrA)", "tma_producer": "TMA-route candidate"}


def reason(guild, r, extra):
    if guild == "MePn consumer":
        acc = sorted({phn_short(s) for s in split(r.ev_phn_accessory)})
        acc_txt = ", ".join(acc) if acc else "none"
        if r.ev_phnJ:
            return f"phnJ present (C-P lyase, EC 4.7.1.1); {PHN_LIST}: {acc_txt}"
        oth = phn_other_short(extra["phn_other"])
        return (f"no phnJ (C-P lyase) annotated; {PHN_LIST}: {acc_txt}"
                + (f"; also annotated: {', '.join(oth)}" if oth else ""))
    if guild == "Betaine reducer":
        comps = grd_components(extra["grd_shared"])
        shared = f"other Grd proteins annotated: {', '.join(comps)}" if comps else "no other Grd proteins annotated"
        if r.ev_grdH:
            subs = sorted({"alpha" if "alpha" in s else "beta" for s in split(r.ev_grdH)})
            return f"grdH present (betaine reductase component B, EC 1.21.4.4; {'+'.join(subs)} subunit); {shared}"
        return f"no betaine-specific component B (grdH, EC 1.21.4.4) annotated; {shared}"
    if guild == "Methanogen":
        n_mcr = len(split(r.ev_mcrA))
        return (f"mcrA present ({n_mcr} methyl-CoM reductase strings incl. alpha subunit EC 2.8.4.1); "
                f"methylamine methyl-transfer genes: {methyl_short(extra['methylamine_corrinoid']) or 'none'}")
    if r.ev_cntA:
        return f"cntA present ({r.ev_cntA}); CutC itself not annotated"
    return f"only a gene 'clustered with choline TMA-lyase' ({r.ev_cutC_cluster}); CutC itself not annotated"


def model_status(guild, k, in_2020):
    if guild == "MePn consumer":
        if k in FLUX["CH4_0.07"]:
            return "simulated (models/MePn_consumers/)"
        if k not in mepn_model_files:
            note = "notebook printed 'skip'" if k in nb_skip else "not reported in notebook output"
            return f"not simulated: no model file in models/MePn_consumers/ ({note})"
        if k in nb_nogrow:
            return "not simulated: model has zero biomass on the 2020 lactate/MePn medium (notebook: \"doesn't grow\")"
        return "not simulated: reason not recorded"
    if k in FLUX["TMA"]:
        return "simulated (correlational/models/*_gapfilled_betaine.json)"
    return "no model" + ("" if in_2020 else " (not in any 2020 modelling set)")


rows = []
for guild, c20, c232 in GUILDS:
    sub = g[(g[c20] if c20 else False) | g[c232]]
    for r in sub.itertuples():
        k = r.key
        in20 = bool(getattr(r, c20)) if c20 else False
        in232 = bool(getattr(r, c232))
        change = ("kept" if in20 and in232 else "dropped" if in20 else "added") if c20 else "new guild (no 2020 set)"
        funcs = (rast_functions(r.mag) if guild in ("Betaine reducer", "Methanogen")
                 or (guild == "MePn consumer" and not r.ev_phnJ) else set())
        extra = {n: sorted(s for s in funcs if re.search(p, s.lower())) for n, p in EXTRA.items()}
        marker = {"MePn consumer": r.ev_phnJ, "Betaine reducer": r.ev_grdH, "Methanogen": r.ev_mcrA,
                  "TMA-route candidate": r.ev_cntA or r.ev_cutC_cluster}[guild]
        phn_oth = [s for s in extra["phn_other"] if not re.search(PHN_ACC, s.lower())]
        support = {"MePn consumer": " | ".join(split(r.ev_phn_accessory) + phn_oth),
                   "Betaine reducer": " | ".join(extra["grd_shared"]),
                   "Methanogen": " | ".join(sorted(set(extra["methylamine_corrinoid"] + split(r.ev_mtbB_mtmB)))),
                   "TMA-route candidate": ""}[guild]
        if guild == "MePn consumer":
            fl, fvar = FLUX["CH4_0.07"].get(k), "CH4 exchange EX_cpd01024_e0 (MePn:Pi >= 0.07)"
        else:
            fl, fvar = FLUX["TMA"].get(k), "TMA exchange EX_cpd00441_e0 (betaine model)"
        fl33 = FLUX["CH4_0.33"].get(k) if guild == "MePn consumer" else None
        gpr = bet_gpr.get(k) if k in FLUX["TMA"] else None
        rows.append({
            "guild": guild, "mag": r.mag, "iterativeID_2020": r.iterativeID_2020, "iterativeID_r232": r.iterativeID_r232,
            "domain": r.domain, "phylum": r.phylum, "family": r.family, "genus": r.genus, "species": r.species,
            "gtdb_2020": qc.loc[k, "gtdb_2020"], "taxon_changed_vs_2020": bool(ids_cmp.loc[r.mag, "taxon_changed"]),
            "in_2020": in20, "in_r232": in232, "change": change, "reason": reason(guild, r, extra),
            "marker_evidence": marker, "supporting_evidence": support,
            "completeness": r.completeness, "contamination": r.contamination,
            "mimag_medium": bool(r.completeness >= 50 and r.contamination < 10),
            **{f"mean_relabund_pct_{s}": site_mean.loc[k, s] for s in SITES},
            "also_in": "; ".join(FLAGS[c] for c in FLAGS if c not in (c20, c232) and getattr(r, c)),
            "modelled_2020": fl is not None, "model_status_2020": model_status(guild, k, in20),
            "flux_variable": fvar if fl is not None else "",
            "mean_modelled_flux": np.mean(list(fl.values())) if fl else np.nan,
            "mean_modelled_flux_CH4_ratio0.33": np.mean(list(fl33.values())) if fl33 else np.nan,
            "betaine_reductase_rxn17220_gpr": ("(empty)" if gpr == "" else gpr) if gpr is not None else "",
        })
ORDER = {gl: i for i, (gl, _, _) in enumerate(GUILDS)}
CHG = {"kept": 0, "dropped": 1, "added": 2}
rows.sort(key=lambda d: (ORDER[d["guild"]], CHG.get(d["change"], 3), sortkey(d["iterativeID_r232"])))
out = pd.DataFrame(rows)
out.to_csv(HERE / "guild_membership_2020_vs_r232.csv", index=False, float_format="%.6g")

# ---- checks ------------------------------------------------------------------------------------------------------
n = out.groupby("guild", observed=True).size().to_dict()
assert n == {"MePn consumer": 27, "Betaine reducer": 15, "Methanogen": 1, "TMA-route candidate": 4}, n
assert out[out.guild == "MePn consumer"].modelled_2020.sum() == len(FLUX["CH4_0.07"]) == 20
assert out[out.guild == "Betaine reducer"].modelled_2020.sum() == len(FLUX["TMA"]) == 15
assert (out.change == "added").sum() == 0
assert not out[["phylum", "iterativeID_r232", "completeness"]].isna().any().any()


# ---- markdown ----------------------------------------------------------------------------------------------------
def md_table(df):
    cells = df.astype(str).values.tolist()
    return "\n".join(["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(c.replace("|", "/") for c in row) + " |" for row in cells])


def ids(df):
    return ", ".join(sorted(df.iterativeID_r232, key=sortkey))


def pct(x):
    return f"{x:.3f}" if x >= 0.001 else f"{x:.1e}"


m = out[out.guild == "MePn consumer"]
b = out[out.guild == "Betaine reducer"]
meth = out[out.guild == "Methanogen"]
tma = out[out.guild == "TMA-route candidate"]
tot = site_mean.sum()            # all 310 MAGs, % of reads (the rest is unbinned)

sets = [("MePn 2020 set", m[m.in_2020]), ("MePn 2020 set, simulated subset", m[m.in_2020 & m.modelled_2020]),
        ("MePn R232 (phnJ)", m[m.in_r232]), ("MePn dropped (no phnJ)", m[m.change == "dropped"]),
        ("Betaine 2020 set (all simulated)", b[b.in_2020]), ("Betaine R232 (grdH)", b[b.in_r232]),
        ("Betaine dropped (no grdH)", b[b.change == "dropped"]),
        ("Methanogen (mcrA), R232 only", meth), ("TMA-route candidates, R232 only", tma),
        ("All 310 MAGs (reference)", None)]
sum_rows = []
for lab, df in sets:
    vals = tot if df is None else df[[f"mean_relabund_pct_{s}" for s in SITES]].sum().set_axis(list(SITES))
    sum_rows.append([lab, "310" if df is None else str(len(df))] + [pct(vals[s]) for s in SITES])
sum_tab = pd.DataFrame(sum_rows, columns=["Set", "n MAGs", "Unrestored (%)", "Restored (%)", "Reference (%)"])

xt_rows = []
for kept in (True, False):
    for sim in (True, False):
        cell = m[(m.in_r232 == kept) & (m.modelled_2020 == sim)]
        xt_rows.append(["phnJ carrier (kept)" if kept else "no phnJ (dropped)",
                        "simulated in 2020" if sim else "not simulated", str(len(cell)), ids(cell) or "–"])
xt_tab = pd.DataFrame(xt_rows, columns=["R232 status", "2020 simulation", "n", "MAGs"])

unsim = m[m.in_2020 & ~m.modelled_2020]
unsim_tab = pd.DataFrame([[r.iterativeID_r232, r.domain, "kept" if r.in_r232 else "dropped",
                           pct(r.mean_relabund_pct_unrestored), r.model_status_2020] for r in unsim.itertuples()],
                         columns=["MAG (R232 ID)", "Domain", "phnJ guild", "Unrestored rel. abund. (%)", "Why no flux"])
unsim_share = unsim.mean_relabund_pct_unrestored.sum() / m[m.in_2020].mean_relabund_pct_unrestored.sum()


def change_table(df, extra_cols=()):
    return md_table(pd.DataFrame([[r.iterativeID_r232, r.iterativeID_2020 or "(none)", r.change, r.reason,
                                   f"{r.completeness:.1f} / {r.contamination:.1f}", "yes" if r.mimag_medium else "no",
                                   pct(r.mean_relabund_pct_unrestored), *[getattr(r, c) for c in extra_cols]]
                                  for r in df.itertuples()],
                                 columns=["R232 ID", "2020 ID", "Change", "Reason (RAST)", "CheckM compl. / contam.",
                                          "MIMAG ≥ medium", "Unrestored (%)", *extra_cols]))


def nophnj_share(ratio):
    col = "mean_modelled_flux" if ratio == 0.07 else "mean_modelled_flux_CH4_ratio0.33"
    sim = m[m.modelled_2020]
    return sim.loc[~sim.in_r232, col].sum() / sim[col].sum()


bd = b[b.change == "dropped"]
bet_sim = b.assign(flux=b.mean_modelled_flux.round(4))
zero_tma = bet_sim[bet_sim.mean_modelled_flux.abs() < 1e-12]
cap_tma = bet_sim[(bet_sim.mean_modelled_flux - 0.1).abs() < 1e-9]
gpr_empty = (b.betaine_reductase_rxn17220_gpr == "(empty)").sum()
n_in_base = sum(bet_in_base.get(k, False) for k in b.mag.map(key))


def grd_list(df):
    return "; ".join(f"{'with' if flag else 'without'} *grdH*: " + ids(df[df.in_r232 == flag])
                     for flag in (True, False) if (df.in_r232 == flag).any())


fam20 = b.gtdb_2020.str.extract(r"o__([^;]*);f__([^;]*)").apply(lambda x: x[1] or f"no family; order {x[0]}", axis=1)
cap_fams = ", ".join(f"{i} ({f})" for i, f in zip(cap_tma.iterativeID_r232, fam20[cap_tma.index]))
halan = b[b.genus == "Halanaerobium"]
bk = b.assign(key=b.mag.map(key)).set_index("key")
bk = bk.join(fva_by)
fva_fixed = bk[(bk.n_alt == 0) & ~bk.forced0]
fva_forced = bk[bk.forced0]
fva_alt = bk[bk.n_alt > 0]
fva_alt_min0 = fva_alt[(fva_alt.n_alt == 6) & (fva_alt.fva_min_max == 0)]


def fva_note(mag):
    r = bk.loc[key(mag)]
    if r.forced0:
        return "forced to 0 (FVA maximum 0) in all 6 cores"
    if r.n_alt == 0:
        return "fixed by the optimum in all 6 cores"
    if r.n_alt == 6 and r.fva_min_max == 0:
        return "not unique: FVA minimum 0 in all 6 cores"
    return f"not unique in {r.n_alt} of 6 cores"


halan_alt0 = halan[halan.mag.map(key).isin(fva_alt_min0.index)]   # Halanaerobium whose TMA could be 0 in every core
assert len(halan_alt0) == 1, "text below says 'makes' (singular)"

md = f"""# Guild membership before (2020) and after (GTDB R232 gene rules)

Generated by `make_guild_membership_comparison.py` from `guilds_r232.csv`, `../saltern_fluxes*.json`,
`../analysis(2).ipynb` (recorded output), `../../data/Cliff_310MAG_relabund.tsv`,
`../betaine_grdH_r232/fva_tma_by_core_r232.csv` and the RAST JSONs.
Per-MAG detail is in `guild_membership_2020_vs_r232.csv` (one row per guild and MAG, {len(out)} rows).

## Counts

| Guild | 2020 set | R232 rule | R232 set | Kept | Dropped | Added |
|---|---|---|---|---|---|---|
| MePn consumer | 27 hand-curated MAGs | *phnJ* (C–P lyase, EC 4.7.1.1) | {m.in_r232.sum()} | {(m.change == 'kept').sum()} | {(m.change == 'dropped').sum()} | {(m.change == 'added').sum()} |
| Betaine reducer | 15 MAGs with betaine models | *grdH* (betaine reductase component B, EC 1.21.4.4) | {b.in_r232.sum()} | {(b.change == 'kept').sum()} | {(b.change == 'dropped').sum()} | {(b.change == 'added').sum()} |
| Methanogen | none | *mcrA* | {len(meth)} | – | – | – |
| TMA-route candidate | none | *cntA*, or genes "clustered with choline TMA-lyase" | {len(tma)} | – | – | – |

- "2020" names the pre-R232 sets used for modelling, not a build date: the betaine models and
  `saltern_fluxes_betaine.json` were made in 2026 (file dates), from the 15-MAG set listed in `analysis_2_marimo.py`.
- **Both gene-based sets are strict subsets of the 2020 sets.** Every *phnJ* carrier among the 310 MAGs was already
  in the 27, and every *grdH* carrier was already in the 15, so nothing was added; the R232 rules only remove members.
- The three archaeal members of the 2020 MePn set ({ids(m[m.domain == 'Archaea'])}) have no 2020 iterativeID because
  the 2020 pipeline assigned none to archaea; 2020 code fell back to the bin name for them. Their 2020 GTDB
  classification exists (`gtdb_2020` column, from `MAG_QC/comparison_vs_2020.csv`).
- {ids(tma[tma.in_2020 | tma.also_in.str.contains('betaine 2020')])} is both a dropped 2020 betaine member and a
  TMA-route candidate, so it appears in two rows.

## Changes and reasons

### MePn consumers: {(m.change == 'dropped').sum()} dropped
{change_table(m[m.change == 'dropped'], ('model_status_2020',)).replace('| model_status_2020 |', '| 2020 simulation |')}

"PhnM" is the conventional symbol for the RAST function "Alpha-D-ribose 1-methylphosphonate 5-triphosphate
diphosphatase (EC 3.6.1.63)"; RAST itself does not name the gene. The searched genes are the `phn_accessory` pattern
of `guilds_r232.csv` (PhnG/H/I/L synthase subunits, PhnM, PhnN); for the dropped members, "also annotated" lists the other Phn-named
genes and phosphonate ABC transporter genes in the same RAST genome. Membership depends only on *phnJ*. The full RAST strings
are in the CSV (`supporting_evidence`).

The {(m.change == 'kept').sum()} kept members all carry *phnJ*: {ids(m[m.change == 'kept'])}.

### Betaine reducers: {(b.change == 'dropped').sum()} dropped
{change_table(b[b.change == 'dropped'])}

The {(b.change == 'kept').sum()} kept members all carry *grdH*: {ids(b[b.change == 'kept'])}.
Of the {len(bd)} dropped members, {bd.supporting_evidence.str.contains('reductase protein A').sum()} carry Grd protein A,
{bd.supporting_evidence.str.contains('component C chain').sum()} carry component C (both shared by the glycine, sarcosine and
betaine reductases) and {bd.supporting_evidence.str.contains('EC 1.21.4.2').sum()} carry a glycine-reductase component B
(EC 1.21.4.2); none has the betaine-specific component B. Read as annotations, they look like glycine (Stickland)
reducers rather than betaine reducers. This is an inference: a missing *grdH* could also be a gap in the bin, though
all {len(bd)} are at least {bd.completeness.min():.1f} % complete.

### Methanogen and TMA-route candidates (new in R232; no 2020 set)
{md_table(pd.DataFrame([[r.guild, r.iterativeID_r232, r.reason, f"{r.completeness:.1f} / {r.contamination:.1f}",
                         pct(r.mean_relabund_pct_unrestored), pct(r.mean_relabund_pct_restored),
                         pct(r.mean_relabund_pct_reference), r.model_status_2020]
                        for r in pd.concat([meth, tma]).itertuples()],
                       columns=["Guild", "R232 ID", "Evidence (RAST)", "CheckM compl. / contam.", "Unrestored (%)",
                                "Restored (%)", "Reference (%)", "2020 model"]))}

The methylamine Co-methyltransferases (EC 2.1.1.248–250; the MtmB, MtbB and MttB subunits) are read here directly
from the RAST JSON. The `ev_mtbB_mtmB` column of `guilds_r232.csv` does not hold them, despite its name: it holds only
the di- and monomethylamine corrinoid proteins.

## Summed relative abundance per site type, before vs after

Sum over set members of each MAG's mean relative abundance (% of reads; plain mean of the 24 samples grouped by
site: unrestored R1 + R2 = 12 samples, restored SF2 = 6, reference R2A = 6). This is the Table S5 definition; it
differs from the DNA-weighted 12-core merge used for the correlations. Binned MAGs make up only part of the reads
(last row); the rest is unbinned.

{md_table(sum_tab)}

## Which 2020 MePn members were actually simulated (27 vs 20)

`saltern_fluxes.json` has 20 MAG keys (identical for both MePn:Pi ratios and all four variables). The other
{len(unsim)} members of the 27 were never simulated, for two recorded reasons:

{md_table(unsim_tab)}

- **No model file (3):** these are exactly the three archaeal members. `analysis_2_marimo.py` loads MePn models only
  from `models/MePn_consumers/`, which holds 24 files; the model-loading cell of `analysis(2).ipynb` printed
  `skip` for these three. Base RAST models for all three do exist in `correlational/models/`, so they were simply
  never curated into `models/MePn_consumers/`. Why is not recorded (that they are archaea and the curated models
  were built from bacterial templates is a plausible inference, not verified).
- **Zero growth (4):** the models load but their biomass optimum is 0 on the 2020 lactate/MePn medium, so the
  growth-check cell put them in `dontGrow` (notebook output: "doesn't grow"), and the simulation loop skips
  `dontGrow`. A one-off re-run of that growth check on the 24 files (2026-10-26; cobra 0.29.1, Gurobi, same
  medium; not part of this script) gave an optimum of 0 for exactly these four and 4.0–5.2 for the other 20.
  The notebook's recorded output of that cell differs for three growers: it shows 0.002, which matches what
  Roseovarius.9, Roseovarius.11 and Albimonas.1 give (0.0018–0.0021) when MePn uptake is bounded at 0.001, the value
  the substrate-addition cell writes into `media`. So the recorded check probably ran after that cell (its execution
  count is not recorded). The four zero-growth MAGs are the same under either bound.
- The unsimulated members held {unsim_share:.0%} of the 2020 set's summed unrestored abundance.

Dropped from the gene-based guild is not the same as unsimulated in 2020:

{md_table(xt_tab)}

So the 2020 model-predicted methane included {((~m.in_r232) & m.modelled_2020).sum()} MAGs without *phnJ*, and left out
{(m.in_r232 & ~m.modelled_2020).sum()} *phnJ* carriers, including both archaeal ones. The {((~m.in_r232) & m.modelled_2020).sum()} non-*phnJ* MAGs supply
{nophnj_share(0.07):.1%} of the summed modelled CH4 flux over the six unrestored cores at MePn:Pi 0.07
({nophnj_share(0.33):.1%} at 0.33), a small share of the modelled total. These shares use the exchange fluxes as
stored; converting to areal units weights the cores differently, so they are not directly comparable to areal-basis
figures (that conversion basis is superseded; see `../../clauded_insights/modelling_consistent_basis_r232.json`).

## Betaine models: TMA output does not follow *grdH*

- All 15 betaine-set MAGs were simulated (`saltern_fluxes_betaine.json` has 15 keys). Its `methane` entries are
  empty dicts: methane is not simulated in the betaine models, only TMA.
- In {gpr_empty} of 15 gap-filled betaine models the betaine reductase reaction (`rxn17220_c0`) has no gene
  association, and in {n_in_base} of 15 it is already present in the un-gapfilled `*_betaine.json`. Modelled TMA
  therefore does not depend on whether the genome carries *grdH*.
- Mean TMA flux is 0 for {len(zero_tma)} MAGs ({grd_list(zero_tma)}) and equals the 0.1 betaine-uptake bound for
  {len(cap_tma)} ({grd_list(cap_tma)}).
- These stored values are one solver's solution, which is not always unique. By the TMA FVA at the lactate optimum in the re-run
  (`../betaine_grdH_r232/fva_tma_by_core_r232.csv`), TMA is fixed in all 6 cores only for {ids(fva_fixed)}, and
  forced to 0 (FVA maximum 0) for {ids(fva_forced)}. The other {len(fva_alt)} ({ids(fva_alt)}) have alternative
  optima in at least one core; for {ids(fva_alt_min0)} the FVA minimum is 0 in all 6 cores.
- Against the text: draft P76 says four MAGs, "notably including *Halanaerobium*", make no TMA, which matches the
  {len(zero_tma)} zero-TMA MAGs above ({len(set(zero_tma.index) & set(halan.index))} of the {len(halan)} *Halanaerobium*). The R232 review copy committed
  in 574e52a (2026-10-21) had a comment on P76 saying all three *Halanaerobium* make none; by this file only
  {ids(halan[halan.mean_modelled_flux.abs() < 1e-12])} makes none
  ({"; ".join(f"{r.iterativeID_r232} averages {r.mean_modelled_flux:.3f}, {fva_note(r.mag)}" for r in halan[halan.mean_modelled_flux.abs() >= 1e-12].itertuples())}).
  So {ids(halan_alt0)} makes TMA only as the solver's pick and could equally make none. That comment is not in the
  review copy revised on 2026-10-26: its tracked rewrite of P76 (`../../clauded_insights/r232_review_edits_2026-10-26.json`)
  says Halanaerobium.2 is forced to make no TMA and that the TMA of other models, including Halanaerobium.1 and
  BM714.4, is not unique, which agrees with the FVA above.
  Draft P76 also calls the six MAGs at the maximum "Dethiosulfatibacteraceae and Acidaminobacteraceae"; by 2020 GTDB
  families they are {cap_fams}.
- The re-run of the betaine models on the 9 *grdH* carriers that the review comment asked for is in
  `../betaine_grdH_r232/` (see its README). With the original settings it reproduces the stored fluxes of the 9,
  finds no annotation-supported reaction or transporter missing from the two zero-TMA *grdH* models, and supplies
  the FVA used above. Its model-vs-measured fit statistics are not repeated here: they use the conversion bases that
  are superseded (see `../../clauded_insights/modelling_consistent_basis_r232.json`).

## Other files that touch this comparison, and what each lacks

| File | What it has | What it lacks |
|---|---|---|
| `guilds_r232.csv` | Per-MAG 2020 and R232 flags for both guilds, RAST marker strings, CheckM, R232 taxonomy | No change/reason column, no per-site-type abundance (only a 12-core DNA-weighted mean), nothing on modelling |
| `Table S5 - MAG taxonomy, quality and guilds (GTDB R232).xlsx` | Same flags for all 310 MAGs, MIMAG flag, per-site-type means (identical definition to this file), correlations | No change/reason, no accessory *phn* strings, nothing on modelling; no generating script found in the repo |
| `README.md` (this folder) and Table 4 of the R232 review docx (`clauded_insights/make_r232_review_docx.py`) | Set sizes and the IDs of the dropped members | No per-member quality, abundance or modelling; no set-level abundance before/after |
| `../iterativeIDs_r232_vs_2020.csv`, `../../MAG_QC/comparison_vs_2020.csv` | 2020 vs R232 labels, taxonomy and CheckM per MAG | No guild membership |
| `guild_methane_correlations_r232.csv` | Methane correlations for the 2020 and R232 sets side by side | Set level only, no membership |
| `../betaine_grdH_r232/` | Betaine models re-simulated for the 9 *grdH* carriers, TMA FVA per MAG and core, per-model diagnostics (origin and GPR of `rxn17220_c0`, annotated *grd* genes, O2 dependence) | Betaine guild only; no per-site-type abundance or membership comparison; its fit statistics use the conversion bases being replaced |
"""
(HERE / "guild_membership_2020_vs_r232.md").write_text(md)
print(out.groupby(["guild", "change"], observed=True).size().to_string())
print(sum_tab.to_string(index=False))
print("wrote", HERE / "guild_membership_2020_vs_r232.csv", "and .md")
