#!/usr/bin/env python
"""Classify the reaction and compound edits of ModelSEEDDatabase #324 (e8446e3e) and #325
(5389554a) against the combined templates (branch salterns-2026-combined, read from
modelling_2026/templates/*.json), under the templates' OWN compound labels.

For each edited reaction carried by GramNeg V7 (GramPos V7 and Archaea V6 give the same
result): is the template's current form balanced; does the database-corrected form, mapped
into template compartments and oriented like the template, balance under template labels; and
are all its compounds in the template. Output: db_prs_vs_templates.tsv next to this script.
"""
import csv, io, json, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tdiff import index                                       # noqa: E402
from dbprop import tsv, changed_files, parse_stoich, tmpl_stoich, to_tmpl_keys, norm  # noqa: E402
TDIR = os.path.abspath(os.path.join(HERE, "..", ".."))

def parse_formula(f):
    d = {}
    for el, n in re.findall(r"([A-Z][a-z]*)(\d*)", f or ""):
        d[el] = d.get(el, 0) + (int(n) if n else 1)
    return d

def balance(st, cpds, ccs):
    el, ch, missing = {}, 0.0, []
    for (cpd, comp), v in st.items():
        if cpd not in cpds or f"{cpd}_{comp}" not in ccs:
            missing.append(f"{cpd}_{comp}"); continue
        for e, n in parse_formula(cpds[cpd].get("formula")).items():
            el[e] = el.get(e, 0) + v * n
        ch += v * ccs[f"{cpd}_{comp}"]["charge"]
    el = {e: round(n, 6) for e, n in el.items() if abs(n) > 1e-6}
    return el, round(ch, 6), missing

def fmt(b):
    el, ch, miss = b
    s = "balanced" if not el and abs(ch) < 1e-6 else f"imbalanced {el} charge {ch:+g}"
    return s + (f"; missing {miss}" if miss else "")

t = index(json.load(open(os.path.join(TDIR, "GramNegModelTemplateV7.json"))))
rx, cpds, ccs = t["reactions"], t["compounds"], t["compcompounds"]
rows = []
for pr, rev in (("#324", "e8446e3e"), ("#325", "5389554a")):
    for fn in changed_files(rev):
        if "reaction_" not in fn:
            continue
        a, b = tsv(rev + "^", fn), tsv(rev, fn)
        for i in sorted(b):
            if i not in a or a[i]["stoichiometry"] == b[i]["stoichiometry"]:
                continue
            row = {"pr": pr, "reaction": i, "name": b[i]["name"][:60]}
            if i + "_c" not in rx:
                row.update(in_template="no"); rows.append(row); continue
            ts = tmpl_stoich(rx[i + "_c"])
            bef = to_tmpl_keys(parse_stoich(a[i]["stoichiometry"]))
            aft = to_tmpl_keys(parse_stoich(b[i]["stoichiometry"]))
            common = [k for k in ts if k in bef and k[0] != "cpd00067"]
            s = -1 if common and ts[common[0]] * bef[common[0]] < 0 else 1
            aft = {k: s * v for k, v in aft.items()}; bef = {k: s * v for k, v in bef.items()}
            same = ("db_after" if norm(ts) == norm(aft) else "db_after_except_H+" if norm(ts, True) == norm(aft, True)
                    else "db_before" if norm(ts) == norm(bef) else "db_before_except_H+" if norm(ts, True) == norm(bef, True)
                    else "neither")
            r = rx[i + "_c"]
            row.update(in_template="yes", template_dir=f"{r['direction']}/{r['GapfillDirection']}",
                       complexes=len(r["templatecomplex_refs"]), template_matches=same,
                       template_balance=fmt(balance(ts, cpds, ccs)),
                       db_corrected_balance_under_template_labels=fmt(balance(aft, cpds, ccs)))
            rows.append(row)
cols = ["pr", "reaction", "name", "in_template", "template_dir", "complexes", "template_matches",
        "template_balance", "db_corrected_balance_under_template_labels"]
with open(os.path.join(HERE, "db_prs_vs_templates.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
for r in rows:
    if r["pr"] == "#324":
        print(r["reaction"], r.get("in_template"), r.get("template_matches"), "|", r.get("template_balance"), "|", r.get("db_corrected_balance_under_template_labels"))
