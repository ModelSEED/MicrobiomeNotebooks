"""Compare ModelSEEDDatabase PR changes (#324 e8446e3e, #325 5389554a) with template contents."""
import csv, io, json, subprocess, sys, re
from tdiff import load, index, V7
DB = "/home/freiburger/Documents/ModelSEEDDatabase"
R = "/home/freiburger/Documents/ModelSEEDTemplates"
csv.field_size_limit(10**9)

def tsv(rev, path):
    raw = subprocess.run(["git", "-C", DB, "show", f"{rev}:{path}"], capture_output=True, check=True).stdout.decode()
    return {r["id"]: r for r in csv.DictReader(io.StringIO(raw), delimiter="\t", quoting=csv.QUOTE_NONE)}

def changed_files(rev):
    out = subprocess.run(["git", "-C", DB, "show", "--name-only", "--format=", rev], capture_output=True, check=True).stdout.decode().split()
    return [f for f in out if f.endswith(".tsv") and ("reaction_" in f or "compound_" in f)]

def parse_stoich(s):
    d = {}
    if not s: return d
    for part in s.split(";"):
        f = part.split(":")
        coef, cpd, comp = float(f[0]), f[1], int(f[2])
        d[(cpd, comp)] = d.get((cpd, comp), 0) + coef
    return d

def tmpl_stoich(r):
    d = {}
    for x in r["templateReactionReagents"]:
        cc = x["templatecompcompound_ref"].split("/")[-1]
        cpd, comp = cc.rsplit("_", 1)
        d[(cpd, comp)] = d.get((cpd, comp), 0) + x["coefficient"]
    return d

def norm(d, nohplus=False):
    return {k: v for k, v in d.items() if abs(v) > 1e-9 and not (nohplus and k[0] == "cpd00067")}

def to_tmpl_keys(dbst):
    # DB compartment index 0 -> the reaction's own compartment (c), 1 -> e
    return {(c, "c" if comp == 0 else "e"): v for (c, comp), v in dbst.items()}

def run(rev, label, tmpl_rev):
    files = changed_files(rev)
    rch, cch = {}, {}
    for f in files:
        a, b = tsv(rev + "^", f), tsv(rev, f)
        for i in b:
            if i in a:
                if "reaction_" in f and a[i]["stoichiometry"] != b[i]["stoichiometry"]:
                    rch[i] = (a[i], b[i])
                if "compound_" in f and (a[i]["charge"] != b[i]["charge"] or a[i]["formula"] != b[i]["formula"]):
                    cch[i] = (a[i], b[i])
    print(f"===== {label}: {len(rch)} reaction stoichiometries, {len(cch)} compounds changed")
    out = {"reactions": {}, "compounds": {}}
    for f in V7:
        t = index(load(R, tmpl_rev, f"templates/v7.0/{f}"))
        rx, cpds = t["reactions"], t["compounds"]
        cats = {"tmpl=db_before": [], "tmpl=db_after": [], "tmpl=db_after_exceptH": [], "tmpl=db_before_exceptH": [], "other": [], "absent": []}
        for i, (a, b) in sorted(rch.items()):
            tid = i + "_c"
            if tid not in rx: cats["absent"].append(i); continue
            ts = tmpl_stoich(rx[tid]); before = to_tmpl_keys(parse_stoich(a["stoichiometry"])); after = to_tmpl_keys(parse_stoich(b["stoichiometry"]))
            sign = None
            for s in (1, -1):
                bb = {k: s * v for k, v in before.items()}; aa = {k: s * v for k, v in after.items()}
                if norm(ts) == norm(bb): cats["tmpl=db_before"].append(i); sign = s; break
                if norm(ts) == norm(aa): cats["tmpl=db_after"].append(i); sign = s; break
                if norm(ts, True) == norm(aa, True): cats["tmpl=db_after_exceptH"].append(i); sign = s; break
                if norm(ts, True) == norm(bb, True): cats["tmpl=db_before_exceptH"].append(i); sign = s; break
            else:
                cats["other"].append(i)
            out["reactions"].setdefault(i, {})[f] = {"template": {f"{k[0]}_{k[1]}": v for k, v in norm(ts).items()},
                                                     "db_before": a["stoichiometry"], "db_after": b["stoichiometry"]}
        print(f" {f}: " + ", ".join(f"{k}={len(v)}" for k, v in cats.items()))
        for k in ("tmpl=db_before", "tmpl=db_before_exceptH", "other", "tmpl=db_after_exceptH"):
            if cats[k]: print(f"    {k}: {cats[k]}")
        cc = {"same_as_db_after": [], "same_as_db_before": [], "other": [], "absent": []}
        for i, (a, b) in sorted(cch.items()):
            if i not in cpds: cc["absent"].append(i); continue
            tc = cpds[i]; tval = (str(tc.get("defaultCharge")), tc.get("formula"))
            if tval == (b["charge"], b["formula"]): cc["same_as_db_after"].append(i)
            elif tval == (a["charge"], a["formula"]): cc["same_as_db_before"].append(i)
            else: cc["other"].append((i, tval, (a["charge"], a["formula"]), (b["charge"], b["formula"])))
        print(f"  compounds: " + ", ".join(f"{k}={len(v)}" for k, v in cc.items()))
        for k in ("same_as_db_before",):
            if cc[k]: print(f"    {k}: {cc[k][:70]}")
        for x in cc["other"][:6]: print("    other:", x)
    return rch, cch, out

if __name__ == "__main__":
    tmpl_rev = sys.argv[1] if len(sys.argv) > 1 else "b6fe86b"
    r1 = run("e8446e3e", "#324 charge-imbalance-corrections", tmpl_rev)
    r2 = run("5389554a", "#325 ferredoxin-one-electron", tmpl_rev)
    json.dump({"324": r1[2], "325": r2[2]}, open(f"dbprop_{tmpl_rev}.json", "w"), indent=1)
