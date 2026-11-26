#!/usr/bin/env python3
"""Independent check that the repo's KBase-dependent files agree with live KBase (cache.json from verify_kbase.py)
and with the independent truth (prep simulation map, prep R232 summary, HEAD 9e087c3).  Reads only; writes only
apply_scratch/verify_kbase/results_repo.tsv.
    nice -n 10 ~/Documents/py_venv/bin/python -I -B verify_repo.py
"""
import csv, hashlib, io, json, os, re, subprocess, sys, time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = os.path.dirname(os.path.dirname(HERE))
PREP = os.path.dirname(PLAN)
OUT = os.path.join(PLAN, "apply_scratch", "verify_kbase")
EK = os.path.expanduser("~/Documents/EmilyKin")
HEAD = "9e087c3"
R13 = "community_modeling/analysis/r13/annotations"
STAGE_START = time.mktime(time.strptime("2026-11-12 13:40:00", "%Y-%m-%d %H:%M:%S"))   # CST, kbase stage start
R = []


def check(name, ok, detail=""):
    R.append((name, "PASS" if ok else "FAIL", detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def note(name, detail):
    R.append((name, "NOTE", detail))
    print(f"[NOTE] {name}: {detail}")


def git(*a):
    return subprocess.run(["git", "-C", EK, *a], capture_output=True, check=True).stdout


def tsv_text(s):
    return list(csv.DictReader(io.StringIO(s), delimiter="\t"))


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    c = json.load(open(os.path.join(OUT, "cache.json")))
    latest = {int(o): i for o, i in c["latest"].items()}
    asm = c["asm_objid"]
    gbin = {int(o): b for o, b in c["gbin"].items()}
    obj_of_g = {b: o for o, b in gbin.items()}
    sim = {r["MAG"]: r for r in tsv(os.path.join(PREP, "results/iterativeid_simulation_r232.tsv"))}
    new_id = {b: r["simulated_r232_id"] for b, r in sim.items()}
    old_id = {b: r["committed_id"] for b, r in sim.items()}
    r232 = {r["user_genome"]: r["classification"] for r in tsv(os.path.join(PREP, "results/gtdbtk.r232.summary.tsv"))}
    inv = {int(r["objid"]): r for r in tsv(os.path.join(PLAN, "kbase/inventory_262319.tsv"))}
    bins = sorted(sim)
    live_ref = {b: f"262319/{asm[b]}/{latest[asm[b]][4]}" for b in bins}

    # ---------------------------------------------------------------- R1 kbase_mag_rename_log.json
    hb = git("show", f"{HEAD}:kbase_mag_rename_log.json")
    wb = open(os.path.join(EK, "kbase_mag_rename_log.json"), "rb").read()
    h, w = json.loads(hb), json.loads(wb)
    fmt_h = json.dumps(h, indent=1, sort_keys=True).encode() == hb
    fmt_w = json.dumps(w, indent=1, sort_keys=True).encode() == wb
    check("log.format", fmt_h and fmt_w and not wb.endswith(b"\n") and b"\r" not in wb,
          f"HEAD and working copy both == json.dumps(indent=1, sort_keys=True), no trailing newline, LF")
    bad = []
    for b in bins:
        e = w.get(b, {})
        if set(e) != {"name_r220", "new_name", "upa", "upa_r220"}:
            bad.append(f"{b} keys {sorted(e)}")
        elif not (e["name_r220"] == h[b]["new_name"] == old_id[b] and e["upa_r220"] == h[b]["upa"]
                  and e["new_name"] == new_id[b] == latest[asm[b]][1] and e["upa"] == live_ref[b]
                  and e["upa"].split("/")[:2] == h[b]["upa"].split("/")[:2]):
            bad.append(f"{b} {e}")
    check("log.fields", not bad and set(w) == set(bins),
          f"291/291: name_r220/upa_r220 == HEAD; new_name == sim R232 == live object name; upa == live latest "
          f"Assembly ref (v2), same objid" + (f"; bad {bad[:4]}" if bad else ""))
    note("log.upa_choice", "upa was moved to v2 (plan left 'update upa' vs 'add upa_r232' open); readers: fetch_rast.py "
         "keys on objid only (unaffected); check_fetch.py fixed to compare ws/objid; build_taxonomy.py and build_panel.py "
         "copy it into taxonomy.tsv kbase_assembly_ref / mag_panel.tsv kbase_ref (291-row cascade, no code reads those "
         "columns); validate_panel.py and rename_mags_kbase.py use keys / write only")

    # ---------------------------------------------------------------- R2 kbase_r232/ run records
    d = os.path.join(EK, "kbase_r232")
    files = sorted(os.listdir(d))
    check("rec.files", files == ["rename_mags_kbase_r232_log.tsv", "report_saved.json", "sets_saved.tsv",
                                 "writeback_gtdbtk_r232_log.tsv"], f"kbase_r232/ holds {files}")
    rl = tsv(os.path.join(d, "rename_mags_kbase_r232_log.tsv"))
    bad = []
    for r in rl:
        o = int(r["objid"])
        b = r["bin"]
        kind_ok = (r["kind"] == "Assembly" and asm.get(b) == o) or (r["kind"] == "Genome" and gbin.get(o) == b)
        if not (kind_ok and r["result"] == "renamed" and r["old_name"] == inv[o]["name"] and r["new_name"] == latest[o][1]
                and r["ver"] == inv[o]["ver"]):
            bad.append(r)
    ren_objs = {int(r["objid"]) for r in rl}
    exp_objs = {o for o in inv if inv[o]["name"] != latest[o][1]}
    check("rec.rename_log", len(rl) == 102 and not bad and ren_objs == exp_objs and sorted(int(r["order"]) for r in rl) == list(range(1, 103)),
          f"{len(rl)} rows, orders 1..102; objid/bin/kind/old name/new name/version match live KBase and the pre-run inventory; "
          f"objids == the 102 renamed objects" + (f"; bad {bad[:3]}" if bad else ""))
    wl = tsv(os.path.join(d, "writeback_gtdbtk_r232_log.tsv"))
    hlog = json.loads(hb)
    bad = [r for r in wl if not (r["old_ref"] == hlog[r["bin"]]["upa"] and r["new_ref"] == live_ref[r["bin"]]
                                 and r["new_lineage"] == r232[r["bin"]] and r["assembly_md5"] == c["classified"][r["bin"]][0]
                                 and r["name"] in (old_id[r["bin"]], new_id[r["bin"]]))]
    check("rec.writeback_log", len(wl) == 291 and len({r["bin"] for r in wl}) == 291 and not bad,
          f"{len(wl)} rows: old_ref == HEAD upa, new_ref == live latest, lineage == prep R232, md5 == recomputed classified MD5"
          + (f"; bad {bad[:2]}" if bad else ""))
    named_new = sum(1 for r in wl if r["name"] == new_id[r["bin"]])
    note("rec.writeback_log_names", f"name column holds the R232 ID for {named_new}/291 rows (writeback ran after the renames)")
    ss = tsv(os.path.join(d, "sets_saved.tsv"))
    check("rec.sets_saved", len(ss) == 1 and ss[0]["new_ref"] == f"262319/585/{latest[585][4]}" == "262319/585/3"
          and ss[0]["old_ref"] == "262319/585/2" and ss[0]["items"] == "291" and ss[0]["labels_changed"] == "51",
          f"{ss[0]['old_ref']} -> {ss[0]['new_ref']}, items {ss[0]['items']}, labels changed {ss[0]['labels_changed']}")
    rs = json.load(open(os.path.join(d, "report_saved.json")))
    check("rec.report_saved", rs["report"] == c["report"]["ref"] == "262319/892/1" and rs["name"] == c["report"]["name"]
          and rs["file_links"] == c["report"]["file_links"] and rs["objects_created"] == len(c["report"]["objects_created"]) == 292
          and rs["handles_resolved"] == 13, f"{rs['report']} {rs['name']}: file list and counts == live report")

    # ---------------------------------------------------------------- R3 r13 _manifest.tsv
    p = f"{R13}/rast/genomes/_manifest.tsv"
    hm = {r["bin"]: r for r in tsv_text(git("show", f"{HEAD}:{p}").decode())}
    wm = {r["bin"]: r for r in tsv(os.path.join(EK, p))}
    bad, other = [], Counter()
    for b in bins:
        H, W = hm[b], wm[b]
        go = int(W["genome_ref"].split("/")[1])
        if not (W["genome_ref"] == H["genome_ref"] and W["assembly_ref"] == H["assembly_ref"] and gbin.get(go) == b
                and int(W["assembly_ref"].split("/")[1]) == asm[b] and W["iterativeID"] == new_id[b]
                and W["genome_name"] == latest[go][1] == f"{W['genomeset']}_{new_id[b]}"):
            bad.append(b)
        for k in H:
            if k not in ("iterativeID", "genome_name") and H[k] != W.get(k):
                other[k] += 1
    check("repo.manifest", not bad and not other and set(wm) == set(bins) and list(next(iter(wm.values()))) == list(next(iter(hm.values()))),
          f"291/291: genome_ref/assembly_ref == HEAD and map to the bin through live KBase; iterativeID == sim; "
          f"genome_name == live Genome name; no other column changed" + (f"; bad {bad[:4]} other {dict(other)}" if bad or other else ""))
    gref_v = Counter(r["genome_ref"].split("/")[2] for r in wm.values())
    aref_v = Counter(r["assembly_ref"].split("/")[2] for r in wm.values())
    note("repo.manifest_versions", f"manifest genome_ref versions {dict(gref_v)}, assembly_ref versions {dict(aref_v)}")

    # ---------------------------------------------------------------- R4 r13 taxonomy.tsv / mag_panel.tsv KBase refs
    tx = {r["bin"]: r for r in tsv(os.path.join(EK, f"{R13}/taxonomy.tsv"))}
    htx = {r["bin"]: r for r in tsv_text(git("show", f"{HEAD}:{R13}/taxonomy.tsv").decode())}
    bad = [b for b in bins if not (tx[b]["kbase_assembly_ref"] == live_ref[b] and tx[b]["iterativeID"] == new_id[b]
                                   and tx[b]["kbase_rast_genome_ref"] == htx[b]["kbase_rast_genome_ref"]
                                   and (not tx[b]["kbase_rast_genome_ref"] or int(tx[b]["kbase_rast_genome_ref"].split("/")[1]) == obj_of_g[b])
                                   and tx[b]["gtdb_classification"] == r232[b])]
    check("repo.taxonomy_refs", not bad and len(tx) == 291, f"taxonomy.tsv 291/291: kbase_assembly_ref == live v2 ref, "
          f"kbase_rast_genome_ref == HEAD and is the bin's Genome, iterativeID == sim, gtdb_classification == prep R232"
          + (f"; bad {bad[:4]}" if bad else ""))
    mp = {r["bin"]: r for r in tsv(os.path.join(EK, f"{R13}/panel/mag_panel.tsv"))}
    bad = [b for b in bins if not (mp[b]["kbase_ref"] == live_ref[b] and mp[b]["iterativeID"] == new_id[b])]
    check("repo.mag_panel_refs", not bad and len(mp) == 291, f"mag_panel.tsv 291/291: kbase_ref == live v2 ref, iterativeID == sim"
          + (f"; bad {bad[:4]}" if bad else ""))

    # ---------------------------------------------------------------- R5 historical KBase records untouched (D5)
    hist = ["kbase_mag_upload_log.json", "kbase_mag_upload.out", f"{R13}/kbase_kofam_push", "community_modeling/plan_r13/kbase_rast",
            f"{R13}/rast/check_fetch.out", f"{R13}/rast/fetch_rast_full.log", f"{R13}/rast/genomes"]
    diff = git("diff", "--name-only", "HEAD", "--", *hist).decode().split()
    diff = [x for x in diff if x != f"{R13}/rast/genomes/_manifest.tsv"]
    newer = []
    for root in [f"{R13}/kbase_kofam_push", "community_modeling/plan_r13/kbase_rast", f"{R13}/rast"]:
        for dp, dn, fn in os.walk(os.path.join(EK, root)):
            for f in fn:
                q = os.path.join(dp, f)
                rel = os.path.relpath(q, EK)
                if os.path.getmtime(q) >= STAGE_START and rel not in (f"{R13}/rast/genomes/_manifest.tsv", f"{R13}/rast/check_fetch.py"):
                    newer.append(rel)
    for f in ("kbase_mag_upload_log.json", "kbase_mag_upload.out"):
        if os.path.getmtime(os.path.join(EK, f)) >= STAGE_START:
            newer.append(f)
    check("repo.d5_kbase_records", not diff and not newer,
          f"tracked KBase push/upload/fetch records unchanged vs HEAD; no file (tracked, untracked or ignored, incl. "
          f"payloads/ and __pycache__/) under kbase_kofam_push/, plan_r13/kbase_rast/, rast/ written since 13:40 CST "
          f"other than _manifest.tsv and check_fetch.py" + (f"; diff {diff[:5]} newer {newer[:8]}" if diff or newer else ""))

    # ---------------------------------------------------------------- R6 hygiene: token, syntax, stray outputs
    tok = open(os.path.expanduser("~/.kbase/token")).read().strip().encode()
    st = git("status", "--porcelain", "-z", "--untracked-files=all").decode().split("\0")
    paths = [x[3:] for x in st if x]
    scan = [os.path.join(EK, x) for x in paths if os.path.isfile(os.path.join(EK, x))]
    for dp, dn, fn in os.walk(os.path.join(PLAN, "apply_scratch")):
        scan += [os.path.join(dp, f) for f in fn]
    for dp, dn, fn in os.walk(os.path.join(PLAN, "apply_log")):
        scan += [os.path.join(dp, f) for f in fn]
    hits = []
    for q in scan:
        try:
            if os.path.getsize(q) < 300_000_000 and tok in open(q, "rb").read():
                hits.append(q)
        except OSError:
            pass
    check("hygiene.token", not hits, f"KBase token not present in {len(scan)} modified/untracked repo files, apply_scratch/ and "
          f"apply_log/" + (f"; FOUND in {hits}" if hits else ""))
    import ast
    syn = []
    for f in ("scripts/rename_mags_kbase_r232.py", "scripts/writeback_gtdbtk_r232_kbase.py", "scripts/rename_mags_kbase.py",
              f"{R13}/rast/check_fetch.py"):
        try:
            ast.parse(open(os.path.join(EK, f)).read())
        except SyntaxError as e:
            syn.append(f"{f}: {e}")
    check("hygiene.syntax", not syn, "the 2 new KBase scripts, rename_mags_kbase.py and check_fetch.py parse" + (f"; {syn}" if syn else ""))
    stray = [x for x in paths if re.search(r"report_stage|emilykin_r232_report|output\.zip|\.pyc$", x)]
    check("hygiene.no_stray_report_files", not stray, f"no report staging, zip or bytecode among git status paths" + (f": {stray}" if stray else ""))

    # name-based references to MAG objects in living repo text (cite bins/UPAs, not names)
    allnames = {n for b in bins for n in (old_id[b], new_id[b])}
    pat = re.compile(r"262319/([A-Za-z][\w.|-]*)")
    gg = subprocess.run(["git", "-C", EK, "grep", "-nIE", r"262319/[A-Za-z]", "--", ".", ":!*.json"], capture_output=True).stdout.decode()
    refs = Counter()
    for line in gg.splitlines():
        for m in pat.findall(line):
            nm = m.rstrip(".")
            if nm in allnames or any(nm.endswith("_" + x) for x in allnames):
                refs[line.split(":", 1)[0]] += 1
    note("repo.name_refs", f"tracked non-JSON files citing a MAG object by '262319/<name>': {dict(refs) if refs else 'none'}")

    with open(os.path.join(OUT, "results_repo.tsv"), "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["check", "status", "detail"])
        wr.writerows(R)
    print(Counter(s for _, s, _ in R))


if __name__ == "__main__":
    main()
