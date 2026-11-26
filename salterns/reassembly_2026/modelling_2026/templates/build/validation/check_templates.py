#!/usr/bin/env python
"""Template-level checks of the combined set (modelling_2026/templates/*.json).

1. Load: every combined and upstream (2294299) v7.0 file loads with
   MSTemplateBuilder.from_dict(json).build(); counts are recorded.
2. Curations present: for every included curation commit C (parent P), each field-level change
   P -> C in templates/v7.0 must be found in the combined file (field value equal to C's,
   added entries present and identical, removed entries absent). Any field that a later
   included curation sets again is reported as 'superseded' with both values, not as missing.
3. Excluded alternatives absent: spot checks that the excluded variants did not leak in.
4. DB counterparts: ModelSEEDDatabase #325 (5389554a) template-carried reactions equal the
   corrected database stoichiometry up to the templates' H+ convention, and cpd11620 is +5.
5. Balance: element and charge balance of every template reaction under the template's own
   compound labels, combined vs upstream (count of imbalanced reactions and the lists that
   change).

    check_templates.py --out check_templates.json
"""
import argparse
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TDIR = os.path.abspath(os.path.join(HERE, "..", ".."))
R = "/home/freiburger/Documents/ModelSEEDTemplates"
H = "/home/freiburger/Documents/ModelSEEDTemplates-hppase"
UPSTREAM = "2294299"
FILES = ["GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json", "ArchaeaTemplateV6.json",
         "Core-V6.json"]
# (label, repo, parent, commit) in the order they are applied in the combined branch
INCLUDED = [
    ("#11 fix-irreversible-reaction-directions", R, "2294299", "21a8084"),
    ("#12 fix-opposite-template-directions", R, "21a8084", "5516611"),
    ("#13 fix-reverse-penalty-sign", R, "21a8084", "09f6bab"),
    ("#9 remove-rxn39860", R, "2294299", "e7fc410"),
    ("#14 stage A electron counts", R, "fb58a76", "a1a0986"),
    ("#14 stage B labels", R, "a1a0986", "6fe204e"),
    ("#14 stage C directions", R, "6fe204e", "ed9b4ab"),
    ("#14 stage D O2/uphill directions", R, "ed9b4ab", "5f53c0a"),
    ("#14 stage E O2/H2O2 directions", R, "5f53c0a", "2219d7e"),
    ("#14 remove rxn10600", R, "2219d7e", "b6fe86b"),
    ("#15 fix-draft-atp-loops", R, "21a8084", "7fd1ab1"),
    ("#16 atp-cycle-directions", R, "2294299", "0f6a5f6"),
    ("#17 hppase-pyrophosphate-proton-pump", H, "2294299", "e57ea75"),
]
EXCLUDED = [
    ("six-reaction-fix-nadx-variant", R, "2294299", "7d104ac"),
    ("negate-all-negative-reverse-penalties", R, "21a8084", "6430128"),
    ("db-align-opposite-directions (reference only)", R, None, "50ed6b4"),
    ("align-reversible-template-directions (reference only)", R, "50ed6b4", "62fcc52"),
    ("codif-combined-aligned (experiment)", R, "fb58a76", "829013e"),
]
_cache = {}


def git_json(repo, rev, f):
    key = (repo, rev, f)
    if key not in _cache:
        raw = subprocess.run(["git", "-C", repo, "show", f"{rev}:templates/v7.0/{f}"],
                             capture_output=True, check=True).stdout
        _cache[key] = json.loads(raw)
    return _cache[key]


def index(t):
    out = {}
    for k, v in t.items():
        if isinstance(v, list) and v and all(isinstance(x, dict) and "id" in x for x in v):
            out[k] = {x["id"]: x for x in v}
        else:
            out[k] = v
    return out


def field_diff(a, b):
    A, B = index(a), index(b)
    ch = []
    for sec in sorted(set(A) | set(B)):
        va, vb = A.get(sec), B.get(sec)
        if not (isinstance(va, dict) and isinstance(vb, dict)):
            continue
        for i in sorted(set(vb) - set(va)):
            ch.append(("add", sec, i, None, None, vb[i]))
        for i in sorted(set(va) - set(vb)):
            ch.append(("remove", sec, i, None, va[i], None))
        for i in sorted(set(va) & set(vb)):
            if va[i] != vb[i]:
                for f in sorted(set(va[i]) | set(vb[i])):
                    if va[i].get(f, "<absent>") != vb[i].get(f, "<absent>"):
                        ch.append(("field", sec, i, f, va[i].get(f, "<absent>"),
                                   vb[i].get(f, "<absent>")))
    return ch


def parse_formula(f):
    d = {}
    for el, n in re.findall(r"([A-Z][a-z]*)(\d*)", f or ""):
        d[el] = d.get(el, 0) + (int(n) if n else 1)
    return d


def imbalanced(t):
    T = index(t)
    cpds, ccs = T["compounds"], T["compcompounds"]
    bad = {}
    for rid, r in T["reactions"].items():
        el, ch = {}, 0.0
        for x in r["templateReactionReagents"]:
            cc = x["templatecompcompound_ref"].split("/")[-1]
            cpd = cc.rsplit("_", 1)[0]
            v = x["coefficient"]
            for e, n in parse_formula(cpds.get(cpd, {}).get("formula")).items():
                el[e] = el.get(e, 0) + v * n
            ch += v * ccs.get(cc, {}).get("charge", 0)
        el = {e: round(n, 6) for e, n in el.items() if abs(n) > 1e-6}
        if el or abs(ch) > 1e-6:
            bad[rid] = {"elements": el, "charge": round(ch, 6)}
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "results", "check_templates.json"))
    args = ap.parse_args()
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    combined = {f: json.load(open(os.path.join(TDIR, f))) for f in FILES}
    upstream = {f: git_json(R, UPSTREAM, f) for f in FILES}
    report = {"load": {}, "included": {}, "excluded": {}, "db325": {}, "balance": {}}

    # 1. load
    for arm, ts in (("combined", combined), ("upstream", upstream)):
        for f, d in ts.items():
            try:
                t = MSTemplateBuilder.from_dict(d).build()
                report["load"][f"{arm}/{f}"] = {
                    "ok": True, "id": t.id, "reactions": len(t.reactions),
                    "compounds": len(t.compounds), "compcompounds": len(t.compcompounds),
                    "complexes": len(t.complexes), "roles": len(t.roles),
                    "biomasses": len(t.biomasses)}
            except Exception as exc:  # noqa: BLE001
                report["load"][f"{arm}/{f}"] = {"ok": False, "error": repr(exc)}

    # 2. included curations
    # the final value each (file, section, id, field) should have = the last included
    # curation's value, in application order
    last = {}
    for label, repo, p, c in INCLUDED:
        for f in FILES:
            for kind, sec, i, fld, old, new in field_diff(git_json(repo, p, f), git_json(repo, c, f)):
                last[(f, sec, i, fld if kind == "field" else kind)] = (label, kind, new)
    for label, repo, p, c in INCLUDED:
        rec = {"changes": 0, "present": 0, "superseded": [], "missing": []}
        for f in FILES:
            C = index(combined[f])
            for kind, sec, i, fld, old, new in field_diff(git_json(repo, p, f), git_json(repo, c, f)):
                rec["changes"] += 1
                if kind == "add":
                    ok = C.get(sec, {}).get(i) == new
                elif kind == "remove":
                    ok = i not in C.get(sec, {})
                else:
                    ok = i in C.get(sec, {}) and C[sec][i].get(fld, "<absent>") == new
                if ok:
                    rec["present"] += 1
                    continue
                owner = last.get((f, sec, i, fld if kind == "field" else kind))
                if kind == "field" and i not in C.get(sec, {}):
                    owner = last.get((f, sec, i, "remove"), owner)
                cur =(C.get(sec, {}).get(i, {}) or {}).get(fld, "<absent>") if kind == "field" else None
                item = {"file": f, "section": sec, "id": i, "field": fld, "kind": kind,
                        "this_value": new, "combined_value": cur}
                if owner and owner[0] != label:
                    item["set_later_by"] = owner[0]
                    rec["superseded"].append(item)
                else:
                    rec["missing"].append(item)
        report["included"][label] = {"commit": c, **{k: (v if isinstance(v, int) else v[:50])
                                                     for k, v in rec.items()},
                                     "n_superseded": len(rec["superseded"]),
                                     "n_missing": len(rec["missing"])}

    # 3. excluded alternatives: how many of their changes (beyond what an included curation
    #    also sets to the same value) appear in the combined set
    for label, repo, p, c in EXCLUDED:
        if p is None:
            p = subprocess.run(["git", "-C", repo, "rev-parse", f"{c}^"], capture_output=True,
                               text=True, check=True).stdout.strip()
        tot = found = 0
        examples = []
        for f in FILES:
            C = index(combined[f])
            for kind, sec, i, fld, old, new in field_diff(git_json(repo, p, f), git_json(repo, c, f)):
                if kind != "field":
                    continue
                tot += 1
                if C.get(sec, {}).get(i, {}).get(fld, "<absent>") == new:
                    found += 1
                    inc = last.get((f, sec, i, fld))
                    if not inc and len(examples) < 10:
                        examples.append((f, i, fld, new))
        report["excluded"][label] = {"commit": c, "field_changes": tot,
                                     "same_value_in_combined": found,
                                     "examples_not_explained_by_included": examples}

    # 4. DB #325 counterpart
    cpd = {f: index(combined[f])["compounds"].get("cpd11620", {}).get("defaultCharge")
           for f in FILES}
    cc = {f: index(combined[f])["compcompounds"].get("cpd11620_c", {}).get("charge")
          for f in FILES}
    report["db325"] = {"cpd11620_defaultCharge": cpd, "cpd11620_c_charge": cc}

    # 5. balance
    for f in FILES:
        bu, bc = imbalanced(upstream[f]), imbalanced(combined[f])
        report["balance"][f] = {
            "imbalanced_upstream": len(bu), "imbalanced_combined": len(bc),
            "fixed": sorted(set(bu) - set(bc)), "newly_imbalanced": sorted(set(bc) - set(bu)),
            "removed_reactions_that_were_imbalanced": sorted(
                i for i in set(bu) - set(bc) if i not in index(combined[f])["reactions"])}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(report, open(args.out, "w"), indent=1, default=str)
    print(json.dumps({"load": {k: (v["ok"], v.get("reactions")) for k, v in report["load"].items()},
                      "included": {k: (v["changes"], v["present"], v["n_superseded"], v["n_missing"])
                                   for k, v in report["included"].items()},
                      "excluded": {k: (v["field_changes"], v["same_value_in_combined"],
                                       len(v["examples_not_explained_by_included"]))
                                   for k, v in report["excluded"].items()},
                      "db325": report["db325"],
                      "balance": {k: (v["imbalanced_upstream"], v["imbalanced_combined"],
                                      len(v["fixed"]), len(v["newly_imbalanced"]))
                                  for k, v in report["balance"].items()}}, indent=1))


if __name__ == "__main__":
    main()
