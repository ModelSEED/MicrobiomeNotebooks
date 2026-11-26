#!/usr/bin/env python3
"""Build the KBase change plan for the EmilyKin GTDB R232 iterativeID update.  READ-ONLY.

Only read-only Workspace methods are callable (whitelist in rpc()); nothing is saved, renamed or deleted.
Token: ~/.kbase/token (never printed).  Standard library only.

Inputs
  ../../results/iterativeid_simulation_r232.tsv   bin -> committed_id (current) / simulated_r232_id (new)
  ../../results/gtdbtk.r232.summary.tsv           R232 classification per bin (GTDB-Tk 2.8.0)
  ../../results/r220_vs_r232_per_mag.tsv          R220 classification (reference only)
  ~/Documents/EmilyKin/kbase_mag_rename_log.json  bin -> {new_name, upa} of the 262319 Assemblies
  ~/Documents/EmilyKin/iterativeIDs.json          EmilyKin ASV iterativeID -> ASV md5 (context check)
  ~/gtdbtk_kbase/emilykin_r232/batchfile.tsv      the FASTAs GTDB-Tk classified (KBase-style MD5 check)

Outputs (this folder)
  inventory_262319.tsv      every object in ws 262319 (incl. hidden), context and planned action
  assemblies_lineage.tsv    (a) 291 Assemblies: new version with data.std_lineages.gtdb (R232)
  renames.tsv               (b) 51 Assemblies + 51 RAST Genomes: Workspace.rename_object, in `order`
  sets.tsv                  (c) AssemblySet 262319/585 items (ref + label) and the 8 GenomeSets
  genomes_taxonomy.tsv      (d) optional, deferred: Genome taxonomy / std_lineages.gtdb / taxon_assignments
  no_action.tsv             ID-like names outside 262319 (ASV / other-project context) and non-workspaces
  preflight_snapshot.json   names, versions and checksums the write scripts must re-assert before writing

    ~/Documents/py_venv/bin/python -I build_kbase_plan.py
"""
import csv, hashlib, json, os, re, socket, sys, time, urllib.error, urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PREP = os.path.dirname(os.path.dirname(HERE))
EK = os.path.expanduser("~/Documents/EmilyKin")
WS_URL = "https://kbase.us/services/ws"
WSID = 262319
SRC_VER = "232.0"
READ_ONLY = {"Workspace.get_workspace_info", "Workspace.list_objects", "Workspace.get_object_info3",
             "Workspace.get_objects2", "Workspace.list_referencing_objects"}
_T = None


def rpc(method, params, tries=4):
    global _T
    if method not in READ_ONLY:
        raise SystemExit(f"refusing non-read-only call {method}")
    if _T is None:
        _T = open(os.path.expanduser("~/.kbase/token")).read().strip()
    body = json.dumps({"version": "1.1", "id": "1", "method": method, "params": params}).encode()
    for k in range(tries):
        req = urllib.request.Request(WS_URL, data=body, headers={"Authorization": _T, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                res = json.load(r)
        except urllib.error.HTTPError as e:
            try:
                res = json.load(e)
            except Exception:
                res = None
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError):
            res = None
        if res is not None:
            if res.get("error"):
                raise RuntimeError(f"{method}: {str(res['error'].get('message') or res['error'])[:300]}")
            return res["result"]
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"{method}: no response")


def list_all(wsid, **kw):
    out, mn = [], 1
    while True:
        p = {"ids": [wsid], "minObjectID": mn, "includeMetadata": 1}
        p.update(kw)
        objs = rpc("Workspace.list_objects", [p])[0]
        if not objs:
            return out
        out += objs
        mn = max(o[0] for o in objs) + 1


def upa(i):
    return f"{i[6]}/{i[0]}/{i[4]}"


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(name, cols, rows):
    with open(os.path.join(HERE, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {name}: {len(rows)} rows")


def taxon_id(c):
    parts = [x for x in c.split(";") if "__" in x and x.split("__", 1)[1]]
    return parts[-1] if parts else ""


def kbase_md5(path):
    """KBase AssemblyUtil MD5: md5 of ','.join(sorted(md5(upper(contig seq))))."""
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


NAME_OK = re.compile(r"^[\w|.-]+$")
TOK = re.compile(r"[A-Za-z0-9_]+\.[0-9]+(?:_m)?(?![0-9])")


def main():
    sim = read_tsv(os.path.join(PREP, "results/iterativeid_simulation_r232.tsv"))
    r232 = {r["user_genome"]: r for r in read_tsv(os.path.join(PREP, "results/gtdbtk.r232.summary.tsv"))}
    r220 = {r["mag"]: r["r220_classification"] for r in read_tsv(os.path.join(PREP, "results/r220_vs_r232_per_mag.tsv"))}
    log = json.load(open(os.path.join(EK, "kbase_mag_rename_log.json")))
    asv_ids = json.load(open(os.path.join(EK, "iterativeIDs.json")))          # ASV id -> ASV md5
    asv_md5 = set(asv_ids.values())
    batch = {l.rstrip("\n").split("\t")[1]: l.split("\t")[0]
             for l in open(os.path.expanduser("~/gtdbtk_kbase/emilykin_r232/batchfile.tsv")) if l.strip()}
    old = {r["MAG"]: r["committed_id"] for r in sim}
    new = {r["MAG"]: r["simulated_r232_id"] for r in sim}
    case = {r["MAG"]: (r["case_r220"], r["case_r232"]) for r in sim}
    changed = [b for b in sorted(old) if old[b] != new[b]]
    assert len(old) == 291 and len(set(new.values())) == 291 and len(changed) == 51, (len(old), len(changed))
    assert all(log[b]["new_name"] == old[b] for b in old), "rename log != committed ids"
    universe = set(old.values()) | set(new.values()) | set(old)

    def ids_in(s):
        out = set()
        for t in TOK.findall(s):
            out |= {c for c in {t} | {t[i + 1:] for i, ch in enumerate(t) if ch == "_"} if c in universe}
        return out

    # ---------------------------------------------------------------- live inventory (read-only)
    wsi = rpc("Workspace.get_workspace_info", [{"id": WSID}])[0]
    objs = list_all(WSID, showHidden=1, showDeleted=1)
    visible = {o[0] for o in list_all(WSID)}
    by_id = {o[0]: o for o in objs}
    names = {o[1]: o for o in objs}
    A = {o[0]: o for o in objs if o[2].startswith("KBaseGenomeAnnotations.Assembly")}
    G = {o[0]: o for o in objs if o[2].startswith("KBaseGenomes.Genome")}
    print(f"ws {WSID} {wsi[1]} modified {wsi[3]} max objid {wsi[4]}; {len(objs)} objects "
          f"({len(visible)} visible): {Counter(o[2] for o in objs).most_common()}")

    # assemblies: identity by objid (rename log) + data.assembly_id == bin
    bin_objid = {b: int(log[b]["upa"].split("/")[1]) for b in old}
    asm_refs = [upa(A[bin_objid[b]]) for b in sorted(old)]
    asm_data = {}
    for i in range(0, len(asm_refs), 100):
        part = asm_refs[i:i + 100]
        res = rpc("Workspace.get_objects2", [{"objects": [{"ref": r, "included": ["/assembly_id", "/md5", "/std_lineages",
                                                                                   "/num_contigs", "/taxon_ref"]}
                                                          for r in part]}])[0]["data"]
        asm_data.update({r: d["data"] for r, d in zip(part, res)})
    # genomes: map to bin through assembly_ref
    gen_refs = [upa(g) for g in G.values()]
    gen_data = {}
    for i in range(0, len(gen_refs), 25):          # included paths still count whole-object size (1 GB cap)
        part = gen_refs[i:i + 25]
        res = rpc("Workspace.get_objects2", [{"objects": [{"ref": r, "included": [
            "/id", "/source_id", "/scientific_name", "/taxonomy", "/taxon_ref", "/taxon_assignments", "/std_lineages",
            "/assembly_ref", "/genetic_code", "/domain"]} for r in part]}])[0]["data"]
        gen_data.update({r: d["data"] for r, d in zip(part, res)})
    objid_bin = {v: b for b, v in bin_objid.items()}
    gen_bin = {}
    for gid, g in G.items():
        d = gen_data[upa(g)]
        b = objid_bin.get(int(d["assembly_ref"].split("/")[1]))
        assert b, f"genome {upa(g)} assembly_ref {d['assembly_ref']} not a MAG Assembly"
        gen_bin[b] = gid
    assert len(gen_bin) == 291

    # sets
    sets = {o[0]: o for o in objs if "Set" in o[2]}
    set_data = {}
    res = rpc("Workspace.get_objects2", [{"objects": [{"ref": upa(o)} for o in sets.values()]}])[0]["data"]
    for o, d in zip(sets.values(), res):
        set_data[o[0]] = d["data"]

    problems = []
    # ---------------------------------------------------------------- (a) assemblies_lineage.tsv
    a_rows = []
    for b in sorted(old):
        o = A[bin_objid[b]]
        d = asm_data[upa(o)]
        cls = r232[b]["classification"]
        local_md5, n_contigs = kbase_md5(batch[b])
        md5_ok = local_md5 == (o[10] or {}).get("MD5") == d.get("md5") and n_contigs == int((o[10] or {}).get("N Contigs", -1))
        id_ok = d.get("assembly_id") == b and o[1] == old[b]
        if not (md5_ok and id_ok and o[4] == 1 and not d.get("std_lineages")):
            problems.append(f"assembly {b} {upa(o)}: md5_ok={md5_ok} id_ok={id_ok} ver={o[4]} std_lineages={d.get('std_lineages')}")
        if not cls or cls.startswith("Unclassified"):
            problems.append(f"{b}: no R232 lineage")
        a_rows.append(dict(ws=WSID, objid=o[0], name=o[1], type=o[2], ver=o[4],
                           action="save_new_version:std_lineages.gtdb", new_name=new[b] if new[b] != old[b] else "",
                           expected_new_ver=o[4] + 1, bin=b, assembly_id=d.get("assembly_id"),
                           current_gtdb=json.dumps(d.get("std_lineages") or {}), current_GTDB_meta=(o[10] or {}).get("GTDB_lineage", ""),
                           kbase_md5=(o[10] or {}).get("MD5"), classified_fasta_kbase_md5=local_md5,
                           md5_match=md5_ok, n_contigs=n_contigs, r220_lineage=r220.get(b, ""),
                           r232_lineage=cls, source_ver=SRC_VER, taxon_id=taxon_id(cls),
                           r232_method=r232[b]["classification_method"],
                           lineage_changed=(r220.get(b, "") != cls)))
    write_tsv("assemblies_lineage.tsv",
              ["ws", "objid", "name", "type", "ver", "action", "new_name", "expected_new_ver", "bin", "assembly_id",
               "current_gtdb", "current_GTDB_meta", "kbase_md5", "classified_fasta_kbase_md5", "md5_match", "n_contigs",
               "r220_lineage", "r232_lineage", "source_ver", "taxon_id", "r232_method", "lineage_changed"], a_rows)

    # ---------------------------------------------------------------- (b) renames.tsv
    # order: a rename whose target name is currently held by another changing object goes after that object's rename
    holder = {old[b]: b for b in changed}
    order, seen = [], set()

    def visit(b, stack=()):
        if b in seen:
            return
        if b in stack:
            raise SystemExit(f"rename cycle through {b}: needs a temporary name")
        h = holder.get(new[b])
        if h and h != b:
            visit(h, stack + (b,))
        seen.add(b)
        order.append(b)
    for b in changed:
        visit(b)
    r_rows, final_names = [], {o[1] for o in objs}
    for kind in ("Assembly", "Genome"):
        for n, b in enumerate(order, 1):
            if kind == "Assembly":
                o = A[bin_objid[b]]
                cur, tgt, derived = o[1], new[b], ""
            else:
                o = G[gen_bin[b]]
                cur = o[1]
                if not cur.endswith("_" + old[b]):
                    problems.append(f"genome {upa(o)} name {cur} does not end with _{old[b]}")
                tgt = cur[: -len(old[b])] + new[b]
                derived = f"{WSID}/{bin_objid[b]} (assembly_ref {gen_data[upa(o)]['assembly_ref']})"
            h = holder.get(new[b])
            chain = ""
            if h and h != b:
                chain = f"target held now by {h} ({old[h]}); rename {h} first"
            elif any(holder.get(new[x]) == b for x in changed):
                dep = [x for x in changed if holder.get(new[x]) == b]
                chain = f"must precede {','.join(dep)} (frees this name)"
            occupied = names.get(tgt)
            coll = "free" if occupied is None else (
                f"held by objid {occupied[0]} ({occupied[1]}) which is renamed earlier in this list"
                if (kind == "Assembly" and occupied[0] in A and objid_bin.get(occupied[0]) in holder.values()) or
                   (kind == "Genome" and occupied[0] in G) else f"COLLISION with objid {occupied[0]} {occupied[2]}")
            if coll.startswith("COLLISION") or not NAME_OK.match(tgt) or tgt.isdigit():
                problems.append(f"rename {cur} -> {tgt}: {coll}")
            r_rows.append(dict(ws=WSID, objid=o[0], name=cur, type=o[2], ver=o[4], action="rename_object", new_name=tgt,
                               order=(n if kind == "Assembly" else 51 + n), kind=kind, bin=b, old_id=old[b], new_id=new[b],
                               case_r220=case[b][0], case_r232=case[b][1], derived_of=derived, chain=chain,
                               name_check=coll, hidden=o[0] not in visible))
    # simulate the renames in order against all 891 names
    cur_names = {o[0]: o[1] for o in objs}
    for r in sorted(r_rows, key=lambda r: r["order"]):
        if r["new_name"] in set(cur_names.values()):
            problems.append(f"simulated rename order fails at {r['name']} -> {r['new_name']}")
        cur_names[r["objid"]] = r["new_name"]
    assert len(set(cur_names.values())) == len(cur_names)
    write_tsv("renames.tsv", ["ws", "objid", "name", "type", "ver", "action", "new_name", "order", "kind", "bin", "old_id",
                              "new_id", "case_r220", "case_r232", "derived_of", "chain", "name_check", "hidden"],
              sorted(r_rows, key=lambda r: r["order"]))

    # ---------------------------------------------------------------- (c) sets.tsv
    s_rows = []
    for sid, so in sorted(sets.items()):
        d = set_data[sid]
        if so[2].startswith("KBaseSets.AssemblySet"):
            items = d["items"]
            new_items = []
            for j, it in enumerate(items):
                oid = int(it["ref"].split("/")[1])
                b = objid_bin[oid]
                assert it["label"] == old[b], (it, old[b])
                new_items.append((new[b], j, it, b, oid))
            for k, (nl, j, it, b, oid) in enumerate(sorted(new_items)):
                s_rows.append(dict(ws=WSID, objid=sid, name=so[1], type=so[2], ver=so[4],
                                   action="resave_set:item_ref+label" if nl != it["label"] else "resave_set:item_ref",
                                   new_name="", item_index_old=j, item_index_new=k, item_ref_old=it["ref"],
                                   item_ref_new=f"{WSID}/{oid}/{A[oid][4] + 1}", label_old=it["label"], label_new=nl,
                                   bin=b, note="set saved once, after (a) and (b): expected new set version "
                                               f"{so[4] + 1}; item_ref_new = latest Assembly version carrying the R232 "
                                               "lineage (expected v2; resolve at write time as salterns set_moves() does); "
                                               "items re-sorted by new label"))
        else:
            el = d.get("elements") or {}
            labels_are_refs = all(k == v.get("ref") for k, v in el.items())
            label_ids = sorted({x for k in el for x in ids_in(k)})
            if sid == 891:
                act, note = ("none_now; resave_if_(d)",
                             "KOfam GenomeSet pins Genome v2; re-point only if (d) saves Genome v3")
            else:
                act, note = "none", "RAST app output set, pinned to Genome v1 by design; keep as provenance"
            s_rows.append(dict(ws=WSID, objid=sid, name=so[1], type=so[2], ver=so[4], action=act, new_name="",
                               item_index_old=f"{len(el)} elements", item_ref_old="", item_ref_new="",
                               label_old=("element keys are the refs (no ID labels)" if labels_are_refs else
                                          f"ID-like keys: {label_ids}"), label_new="", bin="", note=note))
    write_tsv("sets.tsv", ["ws", "objid", "name", "type", "ver", "action", "new_name", "item_index_old", "item_index_new",
                           "item_ref_old", "item_ref_new", "label_old", "label_new", "bin", "note"], s_rows)

    # ---------------------------------------------------------------- (d) genomes_taxonomy.tsv (optional)
    g_rows = []
    for b in sorted(old):
        o = G[gen_bin[b]]
        d = gen_data[upa(o)]
        cls = r232[b]["classification"]
        names_only = "; ".join(x.split("__", 1)[1] for x in cls.split(";") if x.split("__", 1)[1])
        g_rows.append(dict(ws=WSID, objid=o[0], name=o[1], type=o[2], ver=o[4],
                           action="OPTIONAL_DEFERRED:save_new_version(taxonomy,std_lineages.gtdb,taxon_assignments.gtdb)",
                           new_name=next((r["new_name"] for r in r_rows if r["objid"] == o[0]), ""),
                           bin=b, size_bytes=o[9], data_id=d.get("id"), source_id=d.get("source_id"),
                           cur_scientific_name=d.get("scientific_name"), cur_taxonomy=d.get("taxonomy"),
                           cur_taxon_ref=d.get("taxon_ref"), cur_taxon_assignments=json.dumps(d.get("taxon_assignments") or {}),
                           cur_std_lineages=json.dumps(d.get("std_lineages") or {}),
                           proposed_taxonomy=names_only, proposed_std_lineages_gtdb=json.dumps(
                               {"lineage": cls, "source_ver": SRC_VER, "taxon_id": taxon_id(cls)}),
                           proposed_taxon_assignments_gtdb=taxon_id(cls),
                           keep="id, source_id, scientific_name, taxon_ref, features/cdss ids unchanged"))
    write_tsv("genomes_taxonomy.tsv", ["ws", "objid", "name", "type", "ver", "action", "new_name", "bin", "size_bytes",
                                       "data_id", "source_id", "cur_scientific_name", "cur_taxonomy", "cur_taxon_ref",
                                       "cur_taxon_assignments", "cur_std_lineages", "proposed_taxonomy",
                                       "proposed_std_lineages_gtdb", "proposed_taxon_assignments_gtdb", "keep"], g_rows)

    # ---------------------------------------------------------------- inventory_262319.tsv
    inv = []
    rn = {r["objid"]: r for r in r_rows}
    for oid, o in sorted(by_id.items()):
        t = o[2]
        if oid in A:
            b = objid_bin[oid]
            ctx, act = "MAG (Assembly; data.assembly_id = bin)", "new_version(a)" + ("; rename(b)" if oid in rn else "")
        elif oid in G:
            b = next(x for x, g in gen_bin.items() if g == oid)
            ctx, act = "MAG (RAST Genome of the Assembly)", ("rename(b)" if oid in rn else "none") + "; optional(d)"
        elif t.startswith("KBaseSets.AssemblySet"):
            b, ctx, act = "", "MAG set (labels = MAG iterativeIDs)", "resave(c): refs->v2, 51 labels"
        elif t.startswith("KBaseSearch.GenomeSet"):
            b, ctx, act = "", "MAG genome set (keys = refs)", "none" + ("; re-point if (d)" if oid == 891 else "")
        elif t.startswith("KBaseReport"):
            b, ctx, act = "", "app report (hidden; historical)", "none (immutable record of the run)"
        else:
            b, ctx, act = "", "narrative", "none (no cell references a MAG object name)"
        meta = o[10] or {}
        inv.append(dict(ws=WSID, objid=oid, name=o[1], type=t, ver=o[4], action=act,
                        new_name=rn[oid]["new_name"] if oid in rn else "", bin=b, context=ctx,
                        hidden=oid not in visible, saved=o[3], chsum=o[8], size=o[9],
                        GTDB_lineage_meta=meta.get("GTDB_lineage", ""), Taxonomy_meta=meta.get("Taxonomy", "")))
    write_tsv("inventory_262319.tsv", ["ws", "objid", "name", "type", "ver", "action", "new_name", "bin", "context",
                                       "hidden", "saved", "chsum", "size", "GTDB_lineage_meta", "Taxonomy_meta"], inv)

    # ---------------------------------------------------------------- no_action.tsv (other workspaces)
    na = []
    changed_ids = {old[b] for b in changed} | {new[b] for b in changed}
    for w in (275568, 261948, 233212, 221042):
        try:
            info = rpc("Workspace.get_workspace_info", [{"id": w}])[0]
        except RuntimeError as e:
            na.append(dict(ws=w, objid="", name="", type="", ver="", action="none", new_name="",
                           context=f"no such workspace ({e})", evidence=""))
            continue
        wo = list_all(w)
        hits = [o for o in wo if ids_in(o[1])]
        na.append(dict(ws=w, objid="", name=info[1], type="workspace", ver="", action="none", new_name="",
                       context=f"'{info[8].get('narrative_nice_name', '')}': {len(wo)} objects, {len(hits)} with "
                               f"MAG-ID-like names", evidence=json.dumps(Counter(o[2].split('-')[0] for o in wo).most_common(4))))
        for o in hits:
            m = o[10] or {}
            ev = []
            if m.get("asv"):
                ev.append(f"asv md5 {m['asv']} {'IS' if m['asv'] in asv_md5 else 'not'} in EmilyKin iterativeIDs.json")
            if m.get("Source"):
                ev.append(f"Source={m['Source']}")
            ids = ids_in(o[1])
            ns = {261948: "ASV context, other project (codiffusion bioreactor ASV synthetic genomes)",
                  233212: "ASV context (probably EmilyKin's 2025 ASV synthetic genomes: Source=Sludge, .ASV.mdl models; "
                          "the 262319 narrative's batch-import cell creates an ASV_genomes set)"}.get(w, "ASV context")
            na.append(dict(ws=w, objid=o[0], name=o[1], type=o[2], ver=o[4], action="none (never rename)",
                           new_name="", context=(ns + "; string equals a "
                                                 + ("CHANGED MAG id" if ids & changed_ids else "MAG id")
                                                 + f" {sorted(ids)}"), evidence="; ".join(ev)))
    na.append(dict(ws=278957, objid="", name="", type="", ver="", action="none", new_name="",
                   context="repo hits for '278957' are numbers (0.278957, gene coordinates), not a workspace", evidence=""))
    # references from outside 262319
    refs = [f"{WSID}/{o[0]}/{v}" for o in list(A.values()) + list(G.values()) for v in range(1, o[4] + 1)]
    outside = Counter()
    for i in range(0, len(refs), 100):
        for lst in rpc("Workspace.list_referencing_objects", [[{"ref": r} for r in refs[i:i + 100]]])[0]:
            for y in lst:
                if y[6] != WSID:
                    outside[(y[6], y[0], y[1], y[2])] += 1
    na.append(dict(ws="*", objid="", name="list_referencing_objects", type="", ver="", action="none", new_name="",
                   context=f"{len(refs)} Assembly/Genome versions in {WSID}: {len(outside)} referencing objects in other "
                           f"workspaces", evidence=json.dumps(sorted(outside)[:20])))
    write_tsv("no_action.tsv", ["ws", "objid", "name", "type", "ver", "action", "new_name", "context", "evidence"], na)

    # ---------------------------------------------------------------- preflight snapshot
    snap = {"built": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "ws": WSID, "ws_moddate": wsi[3], "ws_max_objid": wsi[4],
            "objects": {str(o[0]): {"name": o[1], "ver": o[4], "chsum": o[8], "type": o[2]} for o in objs
                        if o[0] in A or o[0] in G or o[0] in sets},
            "rename_order": [(r["objid"], r["name"], r["new_name"]) for r in sorted(r_rows, key=lambda r: r["order"])]}
    json.dump(snap, open(os.path.join(HERE, "preflight_snapshot.json"), "w"), indent=1)
    print(f"rename order (assemblies): {order}")
    print(f"problems: {len(problems)}")
    for p in problems:
        print("  ", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
