"""RAST-annotation support for template reactions (local ModelSEED GramNeg V6 template, the family fliu used:
gapfill.ipynb cell 62/72 loaded 'GramNegModelTemplateV6_test' from KBase ws 12218, which is not available locally).
A template reaction is 'RAST-supported' in a MAG when a role of one of its template complexes equals (after
modelseedpy normalize_role) one of the MAG's RAST functions (RASTified_MAG_json 'functions', split on '; ', ' / ',
' @ ' as fliu's builder did, build_models.ipynb cell 4)."""
import json
import re
import sys
from functools import lru_cache

sys.dont_write_bytecode = True
sys.path.insert(0, "/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/analysis")
import vcommon as V  # noqa: E402

TPL_V6 = V.TPL_DIR / "v6.0" / "GramNegModelTemplateV6.json"


def split_functions(fn):
    return [x for x in re.split("; | / | @", fn) if x]


@lru_cache(maxsize=None)
def template_json(path=str(TPL_V6)):
    return json.load(open(path))


@lru_cache(maxsize=None)
def template_reaction_roles(path=str(TPL_V6)):
    """{template rxn id (e.g. rxn00692_c): {'type':..., 'direction':..., 'complexes': {cpx: [role names]}}}"""
    t = template_json(path)
    roles = {r["id"]: r["name"] for r in t["roles"]}
    cpx = {c["id"]: [roles.get(cr["templaterole_ref"].split("/")[-1], cr["templaterole_ref"]) for cr in c["complexroles"]]
           for c in t["complexes"]}
    out = {}
    for r in t["reactions"]:
        out[r["id"]] = {"type": r["type"], "direction": r["direction"], "GapfillDirection": r.get("GapfillDirection"),
                        "name": r["name"],
                        "complexes": {ref.split("/")[-1]: cpx.get(ref.split("/")[-1], []) for ref in r["templatecomplex_refs"]}}
    return out


def normalize(s):
    from modelseedpy.core.msgenome import normalize_role
    return normalize_role(s)


@lru_cache(maxsize=None)
def mag_functions(model_id):
    """{normalized role: [(feature id, original function string)]} from the RASTified MAG genome JSON."""
    g = json.load(open(V.rast_path(model_id)))
    out = {}
    for f in g["features"]:
        for fn in f.get("functions", []) or []:
            for p in split_functions(fn):
                out.setdefault(normalize(p), []).append((f["id"].split("RAST.")[-1], fn))
    return out


def support(model_id, template_rxn_id):
    """RAST evidence for one template reaction in one MAG: {role: [genes]} for the matching roles, plus all roles."""
    tr = template_reaction_roles().get(template_rxn_id)
    if tr is None:
        return {"in_template": False}
    fns = mag_functions(model_id)
    hits, allroles = {}, []
    for cid, rl in tr["complexes"].items():
        for role in rl:
            allroles.append(role)
            n = normalize(role)
            if n in fns:
                hits[role] = sorted({g for g, _ in fns[n]})
    return {"in_template": True, "template_type": tr["type"], "template_direction": tr["direction"],
            "template_roles": sorted(set(allroles)), "rast_supported_roles": hits, "rast_supported": bool(hits)}


def grep_functions(model_id, pattern):
    """MAG features whose RAST function matches a regex (case-insensitive)."""
    g = json.load(open(V.rast_path(model_id)))
    rx = re.compile(pattern, re.I)
    return sorted({(f["id"].split("RAST.")[-1], fn) for f in g["features"] for fn in (f.get("functions") or []) if rx.search(fn)})
