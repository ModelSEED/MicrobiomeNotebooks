#!/usr/bin/env python3
"""Bin the 11 re-assembled wetland metagenomes exactly like the other 770 (per Andrew, 2026-11-17: "start the rest of the
pipeline"), in KBase narrative 272312:

  per assembly: MetaBAT2 + MaxBin2 + CONCOCT (reads from ws 261700)  ->  DAS Tool consensus over the binners that succeeded

Settings, release hashes and output names are copied from the completed jobs of the original run (EE2 records of ws 272312;
every non-input parameter was constant across all 756-768 runs of each app):
  metabat.run_metabat        caf4e3c0…  min_contig_length 2500
  kb_maxbin.run_max_bin      3574f0b2…  thread 4, prob_threshold 0.8, markerset "40", min_contig_length 2500, plotmarker 0
  kb_concoct.run_kb_concoct  936c9bfe…  bowtie2_default, min_contig_length 2500, split 10000/0, kmer 4, vgmm 400/500, pca 90
  kb_das_tool.run_kb_das_tool c70fdbde… diamond, score 0.5, duplicate 0.6, megabin 0.5, bin evals + plots
Resumable and idempotent: a ledger (binning_ledger.json) records every submitted job; EE2 states decide what to do next.
A failed binner is resubmitted at most twice; DAS Tool runs once all three binners are terminal, over those that completed
(skipped if none did). One pass per call; --watch loops every 30 min until every assembly has a terminal DAS Tool job.
Dry run unless --apply. Usage: KB_AUTH_TOKEN=... python -I -B binning_driver.py [--apply] [--watch]
"""
import json, os, sys, time
sys.path.insert(0, "/home/freiburger/gtdbtk_kbase")
import kbase_mags as k

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER = os.path.join(HERE, "binning_ledger.json")
EV = os.path.join(HERE, "binning.events")
SRC_WS, SRC_NAME = 261700, "freiburgermsu:narrative_1780514279067"
DEST_WS, DEST_NAME = 272312, "freiburgermsu:narrative_1790959438780"
APPLY, WATCH = "--apply" in sys.argv, "--watch" in sys.argv
MAX_TRIES = 3
BINNERS = {
    "MetaBAT2": ("metabat.run_metabat", "metabat/run_metabat", "caf4e3c0c415602f3e11f1473734a7b682c6fdfb",
                 {"min_contig_length": 2500}),
    "MaxBin2": ("kb_maxbin.run_max_bin", "kb_maxbin/run_maxbin2", "3574f0b2cc343f4f530bd52b86c8cf8b6b1427f3",
                {"thread": 4, "prob_threshold": 0.8, "markerset": "40", "min_contig_length": 2500, "plotmarker": 0}),
    "CONCOCT": ("kb_concoct.run_kb_concoct", "kb_concoct/run_kb_concoct", "936c9bfed7e07dd1ddcb76ac4fa775b73d60afbb",
                {"read_mapping_tool": "bowtie2_default", "min_contig_length": 2500, "contig_split_size": 10000,
                 "contig_split_overlap": 0, "kmer_size": 4, "max_clusters_for_vgmm": 400, "max_iterations_for_vgmm": 500,
                 "total_percentage_pca": 90, "no_cov_normalization": 0, "no_total_coverage": 0}),
}
DASTOOL = ("kb_das_tool.run_kb_das_tool", "kb_das_tool/run_kb_das_tool", "c70fdbdef91d5b6614dfd387fa127c1492769a8f",
           {"search_engine": "diamond", "score_threshold": 0.5, "duplicate_penalty": 0.6, "megabin_penalty": 0.5,
            "write_bin_evals": 1, "create_plots": 1})
TERM = {"completed", "error", "terminated"}


def ev(msg):
    line = f"[{time.strftime('%F %T')}] {msg}"
    print(line, flush=True)
    open(EV, "a").write(line + "\n")


def assemblies():
    jobs = json.load(open(os.path.join(HERE, "jobs.json")))
    out = {}
    infos = k.rpc(k.WS, "Workspace.get_object_info3", [{"objects": [{"wsid": SRC_WS, "name": n} for n in jobs],
                                                         "ignoreErrors": 1}])[0]["infos"]
    for n, i in zip(jobs, infos):
        if i:
            acc = n.split("_")[0]
            out[acc] = {"asm_name": n, "asm_ref": k.upa(i), "reads_name": f"{acc}_reads",
                        "reads_ref": jobs[n]["read_library_ref"]}
    return out


def submit(method, app_id, ver, params, src):
    job = {"method": method, "app_id": app_id, "service_ver": ver, "wsid": DEST_WS, "params": [params],
           "source_ws_objects": src}
    if not APPLY:
        print(f"  [dry run] {app_id}: {json.dumps(params)[:200]}")
        return None
    return k.rpc(k.EE2, "execution_engine2.run_job", [job], tries=1)[0]


def states(job_ids):
    if not job_ids:
        return {}
    res = k.rpc(k.EE2, "execution_engine2.check_jobs", [{"job_ids": job_ids, "exclude_fields": ["job_input", "job_output"]}])[0]
    return {s["job_id"]: s for s in res["job_states"]}


def binned_ref(name):
    i = k.rpc(k.WS, "Workspace.get_object_info3", [{"objects": [{"wsid": DEST_WS, "name": name}], "ignoreErrors": 1}])[0]["infos"][0]
    return k.upa(i) if i else None


def one_pass():
    led = json.load(open(LEDGER)) if os.path.exists(LEDGER) else {}
    asm = assemblies()
    st = states([j for a in led.values() for t in a.values() for j in t])
    all_done = True
    for acc, a in sorted(asm.items()):
        L = led.setdefault(acc, {})
        tool_state = {}
        for tool, (method, app, ver, extra) in BINNERS.items():
            jobs = L.setdefault(tool, [])
            s = [st.get(j, {}).get("status") for j in jobs]
            if "completed" in s:
                tool_state[tool] = "completed"
            elif any(x not in TERM for x in s):
                tool_state[tool] = "active"
            elif len(jobs) >= MAX_TRIES:
                tool_state[tool] = "gave_up"
            else:
                params = {"workspace_name": DEST_NAME, "assembly_ref": a["asm_ref"],
                          "binned_contig_name": f"{a['asm_name']}.{tool}.BinnedContigs",
                          "reads_list": [f"{SRC_NAME}/{a['reads_name']}"], **extra}
                jid = submit(method, app, ver, params, [a["asm_ref"], a["reads_ref"]])
                if jid:
                    jobs.append(jid)
                    ev(f"SUBMIT {acc} {tool} attempt {len(jobs)} -> {jid}")
                tool_state[tool] = "active"
        d = L.setdefault("DASTool", [])
        ds = [st.get(j, {}).get("status") for j in d]
        if "completed" in ds:
            continue
        if any(x not in TERM for x in ds):
            all_done = False
            continue
        if any(v == "active" for v in tool_state.values()):
            all_done = False
            continue
        ok = [t for t, v in tool_state.items() if v == "completed"]
        if not ok:
            ev(f"SKIP {acc}: no binner completed; no DAS Tool")
            continue
        if len(d) >= MAX_TRIES:
            ev(f"GAVE_UP {acc}: DAS Tool failed {len(d)} times")
            continue
        names = [f"{DEST_NAME}/{a['asm_name']}.{t}.BinnedContigs" for t in ok]
        refs = [binned_ref(f"{a['asm_name']}.{t}.BinnedContigs") for t in ok]
        params = {"workspace_name": DEST_NAME, "assembly_ref": a["asm_ref"], "input_binned_contig_names": names,
                  "output_binned_contig_name": f"{a['asm_name']}.DASTool.BinnedContigs", **DASTOOL[3]}
        jid = submit(*DASTOOL[:3], params, [a["asm_ref"]] + [r for r in refs if r])
        if jid:
            d.append(jid)
            ev(f"SUBMIT {acc} DASTool attempt {len(d)} over {ok} -> {jid}")
        all_done = False
    if APPLY:
        json.dump(led, open(LEDGER, "w"), indent=1)
    summ = {acc: {t: st.get(L[t][-1], {}).get("status") if L.get(t) else None for t in list(BINNERS) + ["DASTool"]}
            for acc, L in sorted(led.items())}
    return all_done, summ


if __name__ == "__main__":
    ev(f"START binning_driver apply={APPLY} watch={WATCH}")
    while True:
        try:
            done, summ = one_pass()
            ev("STATUS " + json.dumps(summ))
            if done or not WATCH:
                if done:
                    ev("ALLDONE every assembly has a terminal DAS Tool outcome")
                break
        except Exception as e:
            ev(f"WARN {e!r:.300}")
        time.sleep(1800)
