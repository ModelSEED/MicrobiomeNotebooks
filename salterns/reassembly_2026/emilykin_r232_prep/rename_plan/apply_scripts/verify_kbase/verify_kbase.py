#!/usr/bin/env python3
"""Independent READ-ONLY verifier of the KBase side of the EmilyKin GTDB R232 update (ws 262319 and neighbours).

Expected values are derived only from sources the kbase stage did not write:
  R232 ID per bin ........ prep results/iterativeid_simulation_r232.tsv
  R232 lineage per bin ... prep results/gtdbtk.r232.summary.tsv     R220 lineage: results/r220_vs_r232_per_mag.tsv
  bin -> Assembly objid .. git show 9e087c3:kbase_mag_rename_log.json (HEAD, before the run)
  Genome -> bin .......... KBase itself: Genome data.assembly_ref -> Assembly data.assembly_id
  pre-run state .......... rename_plan/kbase/inventory_262319.tsv, preflight_snapshot.json, no_action.tsv (built 2026-11-11)
  classified MD5 ......... recomputed here from ~/gtdbtk_kbase/emilykin_r232/batchfile.tsv FASTAs (KBase recipe)

Writes only to apply_scratch/verify_kbase/ (results.tsv, cache.json).
    nice -n 10 ~/Documents/py_venv/bin/python -I -B verify_kbase.py
"""
import csv, hashlib, json, os, subprocess, sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kbro import HS, rpc, ws, blob_get, upa  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = os.path.dirname(os.path.dirname(HERE))                       # rename_plan/
PREP = os.path.dirname(PLAN)                                         # emilykin_r232_prep/
OUT = os.path.join(PLAN, "apply_scratch", "verify_kbase")
EK = os.path.expanduser("~/Documents/EmilyKin")
HEAD = "9e087c3"
WSID = 262319
RUN_START, RUN_END = "2026-11-12T19:40:00", "2026-11-12T20:20:00"   # UTC window of the kbase stage writes
os.makedirs(OUT, exist_ok=True)

R = []                                                               # (check, status, detail)


def check(name, ok, detail=""):
    R.append((name, "PASS" if ok else "FAIL", detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return ok


def note(name, detail):
    R.append((name, "NOTE", detail))
    print(f"[NOTE] {name}: {detail}")


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kbase_md5(path):
    seqs, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif line:
            seqs[cur].append(line)
    md5s = sorted(hashlib.md5("".join(v).upper().encode()).hexdigest() for v in seqs.values())
    return hashlib.md5(",".join(md5s).encode()).hexdigest(), len(seqs)


def taxon_id(c):
    parts = [x for x in c.split(";") if "__" in x and x.split("__", 1)[1]]
    return parts[-1] if parts else ""


def get_objects(refs, included=None, chunk=50):
    out = []
    for i in range(0, len(refs), chunk):
        objs = [{"ref": r} if included is None else {"ref": r, "included": included} for r in refs[i:i + chunk]]
        out += ws("get_objects2", [{"objects": objs}])[0]["data"]
    return out


def main():
    # ------------------------------------------------------------ independent truth
    sim = {r["MAG"]: r for r in tsv(os.path.join(PREP, "results/iterativeid_simulation_r232.tsv"))}
    old_id = {b: r["committed_id"] for b, r in sim.items()}
    new_id = {b: r["simulated_r232_id"] for b, r in sim.items()}
    r232 = {r["user_genome"]: r for r in tsv(os.path.join(PREP, "results/gtdbtk.r232.summary.tsv"))}
    r220 = {r["mag"]: r["r220_classification"] for r in tsv(os.path.join(PREP, "results/r220_vs_r232_per_mag.tsv"))}
    head_log = json.loads(subprocess.check_output(["git", "-C", EK, "show", f"{HEAD}:kbase_mag_rename_log.json"]))
    inv = {int(r["objid"]): r for r in tsv(os.path.join(PLAN, "kbase/inventory_262319.tsv"))}
    snap = json.load(open(os.path.join(PLAN, "kbase/preflight_snapshot.json")))
    batch = {l.rstrip("\n").split("\t")[1]: l.split("\t")[0]
             for l in open(os.path.expanduser("~/gtdbtk_kbase/emilykin_r232/batchfile.tsv")) if l.strip()}
    bins = sorted(sim)
    changed = [b for b in bins if old_id[b] != new_id[b]]
    check("truth.sources", len(bins) == 291 and set(bins) == set(r232) == set(r220) == set(head_log) == set(batch)
          and len(changed) == 51 and len(set(new_id.values())) == 291,
          f"291 bins in sim/R232/R220/HEAD log/batchfile; {len(changed)} changed; R232 IDs unique")
    check("truth.head_log_names", all(head_log[b]["new_name"] == old_id[b] for b in bins),
          "HEAD rename log new_name == sim committed_id 291/291")
    asm_objid = {b: int(head_log[b]["upa"].split("/")[1]) for b in bins}
    bin_of_asm = {o: b for b, o in asm_objid.items()}
    check("truth.head_log_v1", all(head_log[b]["upa"] == f"{WSID}/{asm_objid[b]}/1" for b in bins),
          "HEAD rename log upa all 262319/<objid>/1")
    cmp_ = {}
    for b in bins:
        cmp_[b] = kbase_md5(batch[b])

    # ------------------------------------------------------------ 1. census of 262319
    wi = ws("get_workspace_info", [{"id": WSID}])[0]
    rows = []
    for lo in range(1, wi[4] + 1, 150):
        rows += ws("list_objects", [{"ids": [WSID], "minObjectID": lo, "maxObjectID": lo + 149, "showHidden": 1,
                                     "showDeleted": 1, "showAllVersions": 1, "includeMetadata": 1}])[0]
    deleted = ws("list_objects", [{"ids": [WSID], "showOnlyDeleted": 1, "showHidden": 1}])[0]
    visible = {i[0] for i in ws("list_objects", [{"ids": [WSID], "minObjectID": 880}])[0]}
    vers = defaultdict(dict)
    for i in rows:
        vers[i[0]][i[4]] = i
    latest = {o: v[max(v)] for o, v in vers.items()}
    check("census.max_objid", wi[4] == 892 and max(vers) == 892 and sorted(vers) == list(range(1, 893)),
          f"ws max objid {wi[4]}; objids 1..{max(vers)} contiguous ({len(vers)} objects); ws moddate {wi[3]}")
    check("census.no_deleted", deleted == [], f"{len(deleted)} deleted objects")
    newobj = sorted(set(vers) - set(inv))
    check("census.new_objects", newobj == [892] and latest[892][1] == "GTDBTk_R232_local_report_emilykin"
          and latest[892][2].startswith("KBaseReport.Report") and list(vers[892]) == [1],
          f"objects not in the 2026-11-11 inventory: {[(o, latest[o][1], latest[o][2]) for o in newobj]}")
    note("census.report_hidden", f"report 892 {'is hidden' if 892 not in visible else 'is visible'} in the data panel "
                                 f"(299 older reports in inventory hidden: "
                                 f"{sum(1 for r in inv.values() if r['type'].startswith('KBaseReport') and r['hidden'] == 'True')})")

    kinds = Counter()
    bad = []
    asm_ids = sorted(o for o, r in inv.items() if r["type"].startswith("KBaseGenomeAnnotations.Assembly"))
    gen_ids = sorted(o for o, r in inv.items() if r["type"].startswith("KBaseGenomes.Genome"))
    check("census.asm_objids", set(asm_ids) == set(asm_objid.values()),
          f"inventory Assembly objids ({len(asm_ids)}) == HEAD log objids")
    for o, r in sorted(inv.items()):
        v, iv = vers.get(o, {}), int(r["ver"])
        if not v:
            bad.append(f"{o} missing")
            continue
        old = v.get(iv)
        if o in asm_objid.values():
            exp = {1, 2}
            kinds["asm"] += 1
        elif o == 585:
            exp = {1, 2, 3}
            kinds["set585"] += 1
        else:
            exp = set(range(1, iv + 1))
            kinds[r["type"].split("-")[0]] += 1
        if set(v) != exp:
            bad.append(f"{o} {r['name']} versions {sorted(v)} != {sorted(exp)}")
        if old is None or old[8] != r["chsum"] or old[3] != r["saved"] or old[2] != r["type"]:
            bad.append(f"{o} {r['name']} v{iv} chsum/saved/type changed")
        if o not in asm_objid.values() and o != 585 and max(v) != iv:
            bad.append(f"{o} has a version newer than inventory v{iv}")
        names = {i[1] for i in v.values()}
        if len(names) != 1:
            bad.append(f"{o} versions carry different names {names}")
    check("census.versions_checksums", not bad,
          f"{dict(kinds)}; inventory version chsum/save date/type unchanged for all 891; "
          f"Assemblies exactly v1,v2; set 585 exactly v1-v3; every other object no new version"
          + (f"; problems: {bad[:10]}" if bad else ""))
    for o, s in snap["objects"].items():
        i = vers[int(o)].get(s["ver"])
        if i is None or i[8] != s["chsum"]:
            bad.append(f"snapshot {o}")
    check("census.snapshot", not [b for b in bad if b.startswith("snapshot")],
          f"preflight_snapshot.json: {len(snap['objects'])} objects, version+chsum still present unchanged")
    asm_new = [o for o in asm_objid.values() if 2 in vers[o]]
    tw = sorted(vers[o][2][3] for o in asm_new)
    check("census.asm_v2_window", len(asm_new) == 291 and RUN_START <= tw[0] and tw[-1] <= RUN_END,
          f"291 Assembly v2 saved {tw[0]} .. {tw[-1]}; set 585 v3 saved {vers[585][3][3]}; report {latest[892][3]}")
    names_all = Counter(i[1] for i in latest.values())
    check("census.unique_names", max(names_all.values()) == 1, f"{len(names_all)} distinct names over 892 objects")

    # ------------------------------------------------------------ 2. Assemblies v1 vs v2
    v1 = get_objects([f"{WSID}/{asm_objid[b]}/1" for b in bins])
    v2 = get_objects([f"{WSID}/{asm_objid[b]}/2" for b in bins])
    a_bad = Counter()
    a_ex = defaultdict(list)
    for b, x1, x2 in zip(bins, v1, v2):
        d1, d2, i1, i2 = x1["data"], x2["data"], x1["info"], x2["info"]
        lin = r232[b]["classification"]
        def bad_(k, msg):
            a_bad[k] += 1
            a_ex[k].append(f"{b}: {msg}")
        if d1.get("assembly_id") != b or d2.get("assembly_id") != b:
            bad_("assembly_id", d2.get("assembly_id"))
        if "std_lineages" in d1:
            bad_("v1_has_lineage", d1["std_lineages"])
        rest2 = {k: v for k, v in d2.items() if k != "std_lineages"}
        if rest2 != d1:
            bad_("v2_data_ne_v1", sorted(k for k in set(d1) | set(rest2) if d1.get(k) != rest2.get(k)))
        sl = d2.get("std_lineages") or {}
        if list(sl) != ["gtdb"]:
            bad_("std_lineages_keys", list(sl))
        g = sl.get("gtdb") or {}
        if g != {"lineage": lin, "source_ver": "232.0", "taxon_id": taxon_id(lin)}:
            bad_("lineage", f"{g} != {lin}")
        md5, n = cmp_[b]
        m1, m2 = i1[10] or {}, i2[10] or {}
        if not (d1.get("md5") == d2.get("md5") == m1.get("MD5") == m2.get("MD5") == md5):
            bad_("md5", f"{d1.get('md5')} {m2.get('MD5')} vs FASTA {md5}")
        if not (str(d2.get("num_contigs")) == m2.get("N Contigs") == m1.get("N Contigs") == str(n)):
            bad_("n_contigs", f"{d2.get('num_contigs')} vs {n}")
        if i1[2] != i2[2] or not i2[2].startswith("KBaseGenomeAnnotations.Assembly-5.1"):
            bad_("type", i2[2])
        exp_meta = dict(m1)
        exp_meta.update({"GTDB_lineage": lin, "GTDB_source_ver": "232.0", "GTDB_taxon_id": taxon_id(lin),
                         "GTDB_classified_by": "local GTDB-Tk 2.8.0, GTDB R232 package"})
        if m2 != exp_meta:
            bad_("metadata", {k: (m1.get(k), m2.get(k)) for k in set(m1) | set(m2) if exp_meta.get(k) != m2.get(k)})
        if i1[9] != i2[9] - 0 and False:
            pass
        pv = x2.get("provenance") or [{}]
        c = pv[0].get("custom") or {}
        if pv[0].get("input_ws_objects") != [f"{WSID}/{asm_objid[b]}/1"]:
            bad_("prov_input", pv[0].get("input_ws_objects"))
        if (c.get("bin"), c.get("iterativeID_r232"), c.get("iterativeID_r220"), c.get("previous_gtdb_r220")) != \
                (b, new_id[b], old_id[b], r220[b][:900]):
            bad_("prov_custom", (c.get("bin"), c.get("iterativeID_r232"), c.get("iterativeID_r220")))
        if i2[1] != new_id[b] or i1[1] != new_id[b]:
            bad_("name", f"v1 {i1[1]} v2 {i2[1]} != {new_id[b]}")
    for k in ["assembly_id", "v1_has_lineage", "v2_data_ne_v1", "std_lineages_keys", "lineage", "md5", "n_contigs",
              "type", "metadata", "prov_input", "prov_custom", "name"]:
        check(f"asm.{k}", a_bad[k] == 0, f"{291 - a_bad[k]}/291 ok" + (f"; e.g. {a_ex[k][:3]}" if a_bad[k] else ""))
    lc = sum(1 for b in bins if r232[b]["classification"] != r220[b])
    note("asm.lineage_changed", f"{lc} of 291 R232 lineages differ from R220 (plan: 130)")
    sizes1 = sum(x["info"][9] for x in v1)
    note("asm.size", f"v1 data size total {sizes1/1e6:.1f} MB; v2-v1 size delta per object "
                     f"{min(x2['info'][9]-x1['info'][9] for x1, x2 in zip(v1, v2))}..{max(x2['info'][9]-x1['info'][9] for x1, x2 in zip(v1, v2))} bytes")

    # ------------------------------------------------------------ 3. Genomes: bin via assembly_ref, names, events
    gen = get_objects([f"{WSID}/{o}/2" for o in gen_ids], included=["/assembly_ref", "/ontology_events", "/id"], chunk=25)
    gbin, g_bad = {}, []
    ev_ok = 0
    for o, x in zip(gen_ids, gen):
        ar = x["data"].get("assembly_ref", "")
        ws_, ao, av = (ar.split("/") + ["", "", ""])[:3]
        b = bin_of_asm.get(int(ao)) if ao.isdigit() else None
        if ws_ != str(WSID) or av != "1" or b is None:
            g_bad.append(f"{o} assembly_ref {ar}")
            continue
        gbin[o] = b
        ev = x["data"].get("ontology_events") or []
        if [e.get("id") for e in ev] == ["SSO", "KO"] and ev[1].get("method_version") == "r13-kofam-2":
            ev_ok += 1
        iname, nm = inv[o]["name"], x["info"][1]
        if not iname.endswith("_" + old_id[b]):
            g_bad.append(f"{o} inventory name {iname} does not end with _{old_id[b]}")
            continue
        exp = iname[: -len(old_id[b])] + new_id[b]
        if nm != exp:
            g_bad.append(f"{o} {b}: name {nm} != {exp}")
        if x["data"].get("id") != iname:
            g_bad.append(f"{o} data.id {x['data'].get('id')} != pre-run name {iname} (should be frozen)")
        if x["info"][8] != inv[o]["chsum"] or x["info"][4] != 2:
            g_bad.append(f"{o} v{x['info'][4]} chsum changed")
    check("genome.bins", len(gbin) == 291 and sorted(gbin.values()) == bins,
          f"291 Genomes map 1:1 to the 291 bins through assembly_ref (all pin Assembly v1)")
    check("genome.names", not g_bad, f"name == pre-run prefix + R232 ID for 291/291; data.id frozen; v2 chsum unchanged"
          + (f"; problems {g_bad[:8]}" if g_bad else ""))
    check("genome.ontology_events", ev_ok == 291, f"{ev_ok}/291 Genome v2 carry [SSO RAST, KO r13-kofam-2] events")
    gnew = {o: latest[o][1] for o in gen_ids}

    # ------------------------------------------------------------ 3b. renames by objid, swap chain, name resolution
    exp_asm_name = {asm_objid[b]: new_id[b] for b in bins}
    a_names_bad = [(o, latest[o][1], n) for o, n in exp_asm_name.items() if latest[o][1] != n]
    check("rename.assemblies", not a_names_bad, f"291/291 Assembly objids carry the sim R232 ID ({len(changed)} renamed)"
          + (f"; bad {a_names_bad[:5]}" if a_names_bad else ""))
    renamed_a = [o for o in asm_objid.values() if inv[o]["name"] != latest[o][1]]
    renamed_g = [o for o in gen_ids if inv[o]["name"] != latest[o][1]]
    others_renamed = [o for o in vers if o in inv and o not in asm_objid.values() and o not in gen_ids
                      and inv[o]["name"] != latest[o][1]]
    check("rename.counts", len(renamed_a) == 51 and len(renamed_g) == 51 and not others_renamed
          and {bin_of_asm[o] for o in renamed_a} == set(changed) == {gbin[o] for o in renamed_g},
          f"{len(renamed_a)} Assemblies + {len(renamed_g)} Genomes renamed, exactly the 51 changed bins; "
          f"no other object renamed ({others_renamed})")
    chain = {195: "Chitinophagaceae.2_m", 150: "Chitinophagaceae.1_m"}
    chain_ok = all(latest[o][1] == n for o, n in chain.items()) and bin_of_asm[195] == "coasm_bin.117" \
        and bin_of_asm[150] == "CAN_5_bin.149" and gnew[744].endswith("_Chitinophagaceae.2_m") \
        and gnew[785].endswith("_Chitinophagaceae.1_m") and gbin[744] == "coasm_bin.117" and gbin[785] == "CAN_5_bin.149"
    old_v = ws("get_object_info3", [{"objects": [{"ref": f"{WSID}/{o}/1"} for o in (150, 195, 744, 785)]}])[0]["infos"]
    check("rename.swap_chain", chain_ok and [i[1] for i in old_v] == [latest[o][1] for o in (150, 195, 744, 785)],
          f"195 {latest[195][1]} (coasm_bin.117), 150 {latest[150][1]} (CAN_5_bin.149), 744 {gnew[744]}, "
          f"785 {gnew[785]}; their v1 carry the same names")
    rl = tsv(os.path.join(EK, "kbase_r232/rename_mags_kbase_r232_log.tsv"))
    by_o = {int(r["objid"]): r for r in rl}
    order_ok = int(by_o[195]["order"]) < int(by_o[150]["order"]) and int(by_o[744]["order"]) < int(by_o[785]["order"]) \
        and by_o[195]["time"] <= by_o[150]["time"] and by_o[744]["time"] <= by_o[785]["time"]
    check("rename.chain_order_logged", order_ok,
          f"log: 195 order {by_o[195]['order']} {by_o[195]['time']} before 150 order {by_o[150]['order']} {by_o[150]['time']}; "
          f"744 order {by_o[744]['order']} before 785 order {by_o[785]['order']}")
    # name resolution: every R232 Assembly and Genome name resolves to its objid; retired names do not resolve
    refs = [f"{WSID}/{new_id[b]}" for b in bins] + [f"{WSID}/{gnew[o]}" for o in gen_ids]
    want = [asm_objid[b] for b in bins] + list(gen_ids)
    infos = []
    for i in range(0, len(refs), 200):
        infos += ws("get_object_info3", [{"objects": [{"ref": r} for r in refs[i:i + 200]], "ignoreErrors": 1}])[0]["infos"]
    res_bad = [(r, w, i and i[0]) for r, w, i in zip(refs, want, infos) if not i or i[0] != w]
    check("rename.resolve_new", not res_bad, f"{len(refs)}/582 current names resolve to the expected objid"
          + (f"; bad {res_bad[:5]}" if res_bad else ""))
    retired = sorted({old_id[b] for b in changed} - set(new_id.values()))
    retired_g = sorted({inv[o]["name"] for o in renamed_g} - set(gnew.values()))
    ri = ws("get_object_info3", [{"objects": [{"ref": f"{WSID}/{n}"} for n in retired + retired_g], "ignoreErrors": 1}])[0]["infos"]
    alive = [n for n, i in zip(retired + retired_g, ri) if i]
    check("rename.retired_dead", not alive, f"{len(retired)} retired Assembly names + {len(retired_g)} retired Genome names "
          f"no longer resolve in 262319 (reused old name: {sorted({old_id[b] for b in changed} & set(new_id.values()))})"
          + (f"; still alive {alive}" if alive else ""))
    sp = ws("get_object_info3", [{"objects": [{"ref": f"{WSID}/Simplicispira.1"}, {"ref": f"{WSID}/Simplicispira.1_m"}]}])[0]["infos"]
    check("rename.prefix_pair", (sp[0][0], sp[1][0]) == (169, 308)
          and (bin_of_asm[169], bin_of_asm[308]) == ("CAN_5_bin.40", "coasm_bin.249"),
          f"Simplicispira.1 -> {sp[0][0]} ({bin_of_asm[sp[0][0]]}), Simplicispira.1_m -> {sp[1][0]} ({bin_of_asm[sp[1][0]]})")

    # ------------------------------------------------------------ 4. AssemblySet 585
    s3, s2 = get_objects([f"{WSID}/585/3", f"{WSID}/585/2"])
    items = s3["data"]["items"]
    it_bad = []
    for it in items:
        w, o, v = it["ref"].split("/")
        o = int(o)
        if w != str(WSID) or v != "2" or o not in bin_of_asm:
            it_bad.append(it)
            continue
        b = bin_of_asm[o]
        if it["label"] != new_id[b] or it["label"] != latest[o][1]:
            it_bad.append(it)
    labels = [it["label"] for it in items]
    check("set585.items", len(items) == 291 and not it_bad and len({it["ref"] for it in items}) == 291
          and {int(it["ref"].split("/")[1]) for it in items} == set(asm_objid.values()),
          f"{len(items)} items; each ref 262319/<Assembly objid>/2 (the lineage-carrying latest), label == R232 ID "
          f"== object name; covers all 291" + (f"; bad {it_bad[:3]}" if it_bad else ""))
    check("set585.sorted_unique", labels == sorted(labels) and len(set(labels)) == 291,
          f"labels unique and Python-sorted (first {labels[0]}, last {labels[-1]})")
    old_items = {int(it["ref"].split("/")[1]): it for it in s2["data"]["items"]}
    lc = sum(1 for it in items if old_items[int(it["ref"].split("/")[1])]["label"] != it["label"])
    v2refs_v1 = all(it["ref"].endswith("/1") for it in s2["data"]["items"])
    check("set585.vs_v2", lc == 51 and v2refs_v1 and len(old_items) == 291,
          f"vs v2: {lc} labels changed; v2 items all pinned v1 (unchanged history)")
    desc = s3["data"].get("description", "")
    pv = (s3.get("provenance") or [{}])[0]
    meta = s3["info"][10] or {}
    check("set585.desc_meta_prov", "mag_iterativeID_r220_to_r232.json" in desc and "R232" in desc
          and pv.get("input_ws_objects") == [f"{WSID}/585/2"] and meta.get("naming", "").startswith("iterativeID (GTDB R232")
          and meta.get("n_assemblies") == "291" and s3["info"][1] == "emilykin_MAGs_291",
          f"name {s3['info'][1]}; meta naming={meta.get('naming')!r} n_assemblies={meta.get('n_assemblies')}; "
          f"prov input {pv.get('input_ws_objects')}; description names the R232 map")
    st = {(int(r["item_index_new"]), r["item_ref_new"], r["label_new"]) for r in tsv(os.path.join(PLAN, "kbase/sets.tsv"))
          if r["objid"] == "585"}
    live = {(k, it["ref"], it["label"]) for k, it in enumerate(items)}
    check("set585.vs_plan", st == live, f"plan sets.tsv (index_new, ref_new, label_new) == live v3 for {len(st & live)}/291")

    # ------------------------------------------------------------ 5. GenomeSets untouched
    gs = get_objects([f"{WSID}/891/1"])[0]
    el = gs["data"]["elements"]
    el_ok = {int(v["ref"].split("/")[1]) for v in el.values()} == set(gen_ids) and all(v["ref"].endswith("/2") for v in el.values())
    gsets = [o for o, r in inv.items() if r["type"].startswith("KBaseSearch.GenomeSet")]
    gs_unch = all(set(vers[o]) == {1} and vers[o][1][8] == inv[o]["chsum"] for o in gsets)
    check("genomeset.untouched", el_ok and gs_unch and len(gsets) == 8,
          f"8 GenomeSets v1 only, chsum unchanged; 891 elements = the 291 Genome v2 refs")

    # ------------------------------------------------------------ 6. report 892
    rep = get_objects([f"{WSID}/892/1"])[0]
    rd = rep["data"]
    oc = [x["ref"] for x in rd["objects_created"]]
    exp_oc = {f"{WSID}/{asm_objid[b]}/2" for b in bins} | {f"{WSID}/585/3"}
    check("report.objects_created", len(oc) == 292 and set(oc) == exp_oc,
          f"{len(oc)} refs == the 291 Assembly v2 + 585/3")
    links = rd["file_links"] + rd["html_links"]
    hs = rpc(HS, "AbstractHandle.hids_to_handles", [[l["handle"] for l in links]])[0]
    hmap = {h["hid"]: h for h in hs}
    lk_bad, blob = [], {}
    for l in links:
        h = hmap.get(l["handle"])
        node = l["URL"].rsplit("/", 1)[1]
        if not h or h.get("id") != node:
            lk_bad.append(f"{l['name']}: handle {l['handle']} -> {h and h.get('id')} != {node}")
            continue
        st_, js = blob_get(l["URL"])
        f = ((js or {}).get("data") or {}).get("file") or {}
        blob[l["name"]] = {"status": st_, "size": f.get("size"), "md5": (f.get("checksum") or {}).get("md5"),
                           "fname": f.get("name"), "node": node}
        if st_ != 200 or not f.get("size"):
            lk_bad.append(f"{l['name']}: GET {st_} size {f.get('size')}")
    check("report.handles_blobs", not lk_bad and len(hs) == len(links) == 13,
          f"{len(links)} links ({len(rd['file_links'])} files + {len(rd['html_links'])} html); all handles resolve to "
          f"their node; every node GET 200 with size > 0" + (f"; bad {lk_bad}" if lk_bad else ""))
    # compare blob md5s with the staged local copies
    stage = [d for d in os.listdir(os.path.join(PLAN, "apply_scratch/kbase")) if d.startswith("emilykin_r232_report_stage_")]
    md5_cmp = {}
    srcs = [os.path.expanduser("~/gtdbtk_kbase/emilykin_r232/gtdbtk_out")]
    if len(stage) == 1:
        srcs.insert(0, os.path.join(PLAN, "apply_scratch/kbase", stage[0]))
    for sd in srcs:
        for root, _, files in os.walk(sd):
            for fn in files:
                p = os.path.join(root, fn)
                if fn in md5_cmp or fn not in blob:
                    continue
                md5_cmp[fn] = hashlib.md5(open(p, "rb").read()).hexdigest()
    matched = {n: (blob[n]["md5"] == md5_cmp.get(n)) for n in blob if n in md5_cmp}
    note("report.blob_vs_stage", f"stage dir {stage}; blob md5 == staged file md5 for "
                                 f"{sum(matched.values())}/{len(matched)} same-named files; unmatched names "
                                 f"{sorted(set(blob) - set(md5_cmp))}")
    note("report.text", (rd.get("text_message") or "")[:300].replace("\n", " "))
    refs892 = ws("list_referencing_objects", [[{"ref": f"{WSID}/{asm_objid[b]}/2"} for b in bins]])[0]
    ext = [(r[6], r[0], r[1]) for lst in refs892 for r in lst if r[6] != WSID]
    who = Counter((r[0], r[4]) for lst in refs892 for r in lst)
    check("report.referencing_v2", not ext and set(who) == {(585, 3), (892, 1)} and all(c == 291 for c in who.values()),
          f"each Assembly v2 is referenced by exactly set 585/3 and report 892/1; outside 262319: {ext[:3]}")

    # ------------------------------------------------------------ 7. other workspaces
    na = tsv(os.path.join(PLAN, "kbase/no_action.tsv"))
    objs = [r for r in na if r["objid"].isdigit()]
    oi = ws("get_object_info3", [{"objects": [{"wsid": int(r["ws"]), "objid": int(r["objid"])} for r in objs],
                                  "ignoreErrors": 1}])[0]["infos"]
    o_bad = [(r["ws"], r["objid"], r["name"], i and (i[1], i[4], i[3])) for r, i in zip(objs, oi)
             if not i or i[1] != r["name"] or str(i[4]) != r["ver"] or i[2] != r["type"] or i[3] >= "2026-11-11"]
    check("otherws.objects", not o_bad, f"{len(objs)} ID-named objects in 261948/233212: name, version, type unchanged and "
          f"last saved before 2026-11-11" + (f"; changed {o_bad}" if o_bad else ""))
    for w in (261948, 233212, 275568, 221042):
        try:
            x = ws("get_workspace_info", [{"id": w}])[0]
            hit = RUN_START <= x[3] <= RUN_END
            msg = f"{x[1]} moddate {x[3]}, max objid {x[4]}" + (" (INSIDE the kbase-stage write window)" if hit else "")
            if hit:
                check(f"otherws.{w}", False, msg)
            else:
                note(f"otherws.{w}", msg)
        except RuntimeError as e:
            note(f"otherws.{w}", f"not readable: {str(e)[:120]}")

    # ------------------------------------------------------------ cache for the repo checks
    cache = {
        "latest": {o: list(i[:10]) for o, i in latest.items()},
        "asm_objid": asm_objid, "gbin": {str(o): b for o, b in gbin.items()}, "gnew": {str(o): n for o, n in gnew.items()},
        "set585_items": items, "report": {"name": rep["info"][1], "ref": upa(rep["info"]), "objects_created": oc,
                                          "file_links": [l["name"] for l in rd["file_links"]], "blob": blob},
        "classified": {b: list(cmp_[b]) for b in bins},
    }
    json.dump(cache, open(os.path.join(OUT, "cache.json"), "w"), indent=1, sort_keys=True)
    with open(os.path.join(OUT, "results_kbase.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["check", "status", "detail"])
        w.writerows(R)
    print(Counter(s for _, s, _ in R))


if __name__ == "__main__":
    main()
