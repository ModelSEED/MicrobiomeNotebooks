#!/usr/bin/env python
"""Assemble and bin the 133 GROWdb river metagenomes (171 SRA runs) in KBase, then refine, QC and classify the bins.

Adapted 2026-11-09 from salterns/reassembly_2026/pipeline.py (v2.2): the same apps, tags, settings and safeguards (write-
ahead ledger, recovery of lost submits, output validation, retry budgets, per-stage pilots, local DAS Tool 1.1.7, local
GTDB-Tk 2.8.0/R232 under a memory guard). GROW-specific (see README.md):
  * unit = sample (133). Multi-run samples (37) have their trimmed runs merged (read and base totals must equal the sums)
    before assembly; a sample starts as soon as every one of its runs is 'done' in the import ledger.
  * no co-assembly; binning coverage is per-sample (each sample's own library), as GROWdb did.
  * reads come from ws 275929 (imported by raw_metagenomes_sra/kbase_import.py; its ledger is only READ here); every
    output (trimmed/merged reads, assemblies, BinnedContigs, bins, CheckM) is written to ws 276048.
  * KBase capacity is shared with the salterns watcher, which has priority: caps extreme 3, bigmem 6, njs 6 while the
    salterns ledger has trim/merge/assembly/binning work pending or live; afterwards up to the salterns totals
    (12/20/15) minus the live salterns and import-watcher jobs (never below 3/6/6).
  * sample pilots: deadrun_0006 (13 Gbp, 1 run) and eastriver_0100 (67.7 Gbp, 2 runs, the largest). Other samples are
    trimmed/merged once both pilots passed that stage, and assembled/binned once both pilots finished binning.
  * GTDB-Tk (scope gtdbtk-grow-assembly) never starts while the salterns GTDB-Tk is pending or its scope is active.
  * local DAS Tool workers <= 2; local steps need >= 100 GB free in work/.

Stages (app, KBase tag, client group):
  import      read-only mirror of ../raw_metagenomes_sra/kbase_import_ledger.json (run -> PairedEndLibrary in ws 275929)
  trim        kb_trimmomatic/run_trimmomatic                 release  bigmem  per run; Trimmomatic 0.39 PE: TruSeq3-PE-2
                clip 2:30:10, LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50; output '<sample>.<SRR>.trim_paired'.
                If trimming fails for a run, its raw import is used instead.
  merge       kb_ReadsUtilities/KButil_Merge_MultipleReadsLibs_to_OneLibrary release njs  multi-run samples only;
                read_count and total_bases must equal the input sums (else the sample is not assembled)
  assemble    MEGAHIT/run_megahit                            release  bigmem  MEGAHIT 1.2.9 meta-sensitive for EVERY sample
                (user decision 2026-11-10; metaSPAdes is not used). After an intrinsic MEGAHIT failure or an assembly
                < 10 Mbp: min-count 2 + k 27..127, then the MEGAHIT defaults (the salterns co-assembly ladder)
  metabat     metabat/run_metabat                            beta     bigmem  \
  maxbin      kb_maxbin/run_maxbin2                          release  extreme  > per-sample coverage (the sample's own
  concoct     kb_concoct/run_kb_concoct                      release  extreme /  library, trimmed and merged)
  dastool     LOCAL DAS Tool 1.1.7 (~/Documents/dastool_env) via local/dastool_worker.py; saves '<sample>.DASTool.
                BinnedContigs' to ws 276048
  extract     MetagenomeUtils/extract_bins_as_assemblies     beta     njs     every refined bin -> Assembly (+ set)
  checkm      kb_Msuite/run_checkM_lineage_wf                beta     njs     CheckM 1.1.2 reduced tree; then
                local/collect_checkm_select.py picks bins >= 50% complete and < 10% contaminated
  gtdbtk      LOCAL GTDB-Tk 2.8.0 / R232 (local/run_gtdbtk_reassembly.sh) in a 170 GB systemd scope, started when
                MemAvailable >= 150 GB, once every sample's CheckM is finished and the salterns GTDB-Tk is not pending
The first unit of dastool, extract and checkm runs alone (pilot) before the rest fan out.

Do not delete work/: work/dastool/*/attempt*/<unit>.save_objects.json is read until extract has finished, and
work/checkm/*/fasta plus the batchfiles until GTDB-Tk has finished.

Usage (KBUtilLib from ~/Documents/KBUtilLib/src; token from ~/.kbase/token, never printed):
    python pipeline.py plan [--through STAGE]      # what would be submitted now; submits nothing
    python pipeline.py advance [--through STAGE]   # poll, record outputs, submit what is ready
    python pipeline.py status
    python pipeline.py watch --through STAGE --interval 1800   # see README.md for the detached systemd-run start
"""

import argparse
import csv
import fcntl
import json
import os
import signal
import subprocess
import sys
import time
import traceback
from collections import Counter, defaultdict
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
LOCAL = os.path.join(HERE, "local")
WORK = os.path.join(HERE, "work")                      # large local intermediates (git-ignored)
_KBUTILLIB_SRC = os.path.normpath(os.path.join(HERE, "..", "..", "..", "KBUtilLib", "src"))
if os.path.isdir(_KBUTILLIB_SRC) and _KBUTILLIB_SRC not in sys.path:
    sys.path.insert(0, _KBUTILLIB_SRC)

from kbutillib.kb_app_runner.nms import NMSSpecCache  # noqa: E402
from kbutillib.kb_app_runner.runner import AppRunner  # noqa: E402
from kbutillib.kb_genome_utils import KBGenomeUtils  # noqa: E402
from kbutillib.kb_job_utils import KBJobUtils  # noqa: E402

WS_ID = 276048                                         # outputs: "GROWdb metagenomes: assembly and binning (2026)"
WS_NAME = "freiburgermsu:narrative_1794250260895"
READS_WS_NAME = "freiburgermsu:narrative_1794152940209"   # ws 275929: the imported reads (read by ref only)
MANIFEST = os.path.join(HERE, "manifest.tsv")
LEDGER = os.path.join(HERE, "ledger.json")
IMPORT_LEDGER = os.path.normpath(os.path.join(HERE, "..", "raw_metagenomes_sra", "kbase_import_ledger.json"))  # READ ONLY
SALTERNS_LEDGER = os.path.normpath(os.path.join(HERE, "..", "..", "salterns", "reassembly_2026", "ledger.json"))  # READ ONLY
MAX_TRIES = 3
SPADES_TRIES_BEFORE_MEGAHIT = 2
HARD_CAP_ATTEMPTS = 6                                  # repeated cancels cannot loop forever
LIVE = {"created", "estimating", "queued", "running"}
STAGES = ["import", "trim", "merge", "assemble", "metabat", "maxbin", "concoct", "dastool", "extract", "checkm", "gtdbtk"]
APPS = {  # variant -> (app id, KBase tag)
    "trim": ("kb_trimmomatic/run_trimmomatic", "release"),
    "megahit": ("MEGAHIT/run_megahit", "release"),       # assemble (2026-11-10: MEGAHIT for every sample)
    "megahit_retry": ("MEGAHIT/run_megahit", "release"),
    "megahit_default": ("MEGAHIT/run_megahit", "release"),
    "merge": ("kb_ReadsUtilities/KButil_Merge_MultipleReadsLibs_to_OneLibrary", "release"),
    "metabat": ("metabat/run_metabat", "beta"),
    "maxbin": ("kb_maxbin/run_maxbin2", "release"),
    "concoct": ("kb_concoct/run_kb_concoct", "release"),
    "extract": ("MetagenomeUtils/extract_bins_as_assemblies", "beta"),
    "checkm": ("kb_Msuite/run_checkM_lineage_wf", "beta"),
}
POOL = {"import": "njs", "trim": "bigmem", "assemble": "extreme", "assemble_megahit": "bigmem", "merge": "njs",
        "megahit": "bigmem", "megahit_retry": "bigmem", "megahit_default": "bigmem",   # assemble since 2026-11-10
        "coassemble": "bigmem", "coassemble_retry": "bigmem", "coassemble_default": "bigmem", "metabat": "bigmem",
        "maxbin": "extreme", "concoct": "extreme", "extract": "njs", "checkm": "njs"}   # coassemble*: salterns' variants
POOL_CAP_BUSY = {"extreme": 3, "bigmem": 6, "njs": 6}  # while salterns has trim/merge/assembly/binning work left
POOL_CAP_TOTAL = {"extreme": 12, "bigmem": 20, "njs": 15}   # the salterns caps: shared total once salterns has cleared
SALTERNS_HEAVY = ("import", "trim", "assemble", "merge", "coassemble", "metabat", "maxbin", "concoct")
LOCAL_CAP = 2
MIN_FREE_GB = 100                                      # free space in work/ before a local step starts
BINNERS = ("metabat", "maxbin", "concoct")
BINNER_LABEL = {"metabat": "metabat2", "maxbin": "maxbin2", "concoct": "concoct"}
MIN_ASSEMBLY_BP = 10_000_000                           # "something broke" floor for a per-sample assembly
STALE_H = {"queued": 7 * 24, "running": 96}
PILOT = ("dastool", "extract", "checkm")
PILOT_SAMPLES = ("deadrun_2019_sw_WHONDRS-S19S_0006",   # small: 1 run, 13.0 Gbp
                 "eastriver_2019_sw_WHONDRS-S19S_0100")  # largest: 2 runs, 67.7 Gbp
OOM_MARKERS = ("return code: 12", "return code: 28", "Cannot allocate memory", "No space left on device")
TRANSIENT_MARKERS = ("URLError", "timed out", "TimeoutError", "no response from", "RemoteDisconnected", "ConnectionError",
                     "ConnectionResetError", "ConnectionRefusedError", "ConnectionAbortedError", "BrokenPipeError",
                     "Bad Gateway", "Service Unavailable", "Gateway Time-out", "Gateway Timeout",
                     "HTTP Error 5",                   # urllib 500-599, incl. Cloudflare's 52x in front of kbase.us
                     " Server Error: ",                # requests' raise_for_status on a 5xx
                     "ServerError: Unknown: 0.",       # KBUtilLib baseclient: a 5xx whose body is not JSON-RPC
                     "IncompleteRead", "ChunkedEncodingError", "SSLError", "Max retries exceeded", "Connection aborted",
                     "Temporary failure")
TRANSIENT_BUDGET_S = 24 * 3600                         # a run of transient local attempts gives up after a day
KBASE_TRANSIENT_MARKERS = ("Unexpected Job Error",)    # EE2-side kills/outages (incl. monitor_jobrunner_logs.py): user
KBASE_TRANSIENT_BUDGET_S = 72 * 3600                   # decision 2026-11-12; not counted, 72 h from the first such failure
BIG_SAMPLE_BP = 50_000_000_000                         # samples >= this wait for the eastriver pilot (2026-11-12)
COLLECT_TRANSIENT_MAX = 24                             # CheckM collections hit by KBase transport errors: ~12 h
LOCAL_STALE_H = 24
CHECK_FAILURES_MAX = 24                                # completion checks that raise: re-polled for ~12 h, then rejected
GTDBTK_SCOPE = "gtdbtk-grow-assembly"
SALTERNS_GTDBTK_SCOPE = "gtdbtk-salterns-reassembly"
SCOPE_RUNNING = ("active", "activating", "deactivating", "reloading", "refreshing")
SCOPE_ENDED = ("inactive", "failed")
GTDBTK_DIR = os.path.expanduser("~/gtdbtk_kbase/grow_assembly_2026")
GTDBTK_MIN_AVAIL_GB = 150
# Hybrid local/KBase mode (user decision 2026-11-23, BRIEF.md): these stages run locally (../hybrid/hw.py, the KBase
# apps' own tool versions and commands) and are saved back into ws 276048; queued KBase jobs are taken over.
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "hybrid")))
import glue as hybrid  # noqa: E402
HYBRID = ("assemble", "metabat", "maxbin", "concoct", "extract", "checkm")
DRIVER = "per-sample"


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(msg):
    print(f"[{now()}] {msg}", flush=True)


# ------------------------------------------------------------------------------------------------ units and naming
def load_manifest():
    with open(MANIFEST, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 171 and len({r["run"] for r in rows}) == 171 and len({r["sample"] for r in rows}) == 133, \
        "manifest must hold 171 unique runs from 133 samples"
    return rows


def units(rows):
    """The 133 sample units (pilots first) and the 171 run units (trim)."""
    by = defaultdict(list)
    for r in rows:
        by[r["sample"]].append(r)
    order = list(PILOT_SAMPLES) + sorted(s for s in by if s not in PILOT_SAMPLES)
    samples = [{"unit": s, "kind": "sample", "runs": by[s]} for s in order]
    runs = [{"unit": r["run"], "kind": "run", "sample": s, "row": r} for su in samples for r in su["runs"]
            for s in [su["unit"]]]
    return samples, runs


def name(stage, u):
    if stage == "import":
        return f"{u['sample']}.{u['unit']}.reads"
    if stage == "trim":
        return f"{u['sample']}.{u['unit']}.trim"                 # kb_trimmomatic saves '<this>_paired'
    s = u["unit"]
    return {"merge": f"{s}.runs{len(u['runs'])}.trim_paired.merged.reads",
            "assemble": f"{s}.MEGAHIT.Assembly",
            "metabat": f"{s}.MetaBAT2.BinnedContigs", "maxbin": f"{s}.MaxBin2.BinnedContigs",
            "concoct": f"{s}.CONCOCT.BinnedContigs", "dastool": f"{s}.DASTool.BinnedContigs",
            "extract": f"{s}.DASTool_bins.AssemblySet"}[stage]


# ------------------------------------------------------------------------------------------------ ledger
def load_ledger():
    if os.path.exists(LEDGER):
        with open(LEDGER) as fh:
            return json.load(fh)
    return {"workspace_id": WS_ID, "workspace_name": WS_NAME, "created_utc": now(), "tasks": {}}


def save_ledger(ledger):
    tmp = LEDGER + ".part"
    with open(tmp, "w") as fh:
        json.dump(ledger, fh, indent=1)
    os.replace(tmp, LEDGER)


def task(ledger, stage, unit):
    return ledger["tasks"].setdefault(f"{stage}:{unit}", {"stage": stage, "sample": unit, "state": "pending",
                                                          "attempts": []})


def get_task(ledger, stage, unit):
    return ledger["tasks"].get(f"{stage}:{unit}")


def output_of(ledger, stage, unit):
    t = get_task(ledger, stage, unit)
    return t.get("output_ref") if t and t["state"] == "done" else None


def terminal(ledger, stage, unit):
    t = get_task(ledger, stage, unit)
    return bool(t) and t["state"] in ("done", "empty", "failed")


def dead(ledger, stage, unit):
    """terminal without an output: everything downstream of it is 'empty'."""
    return terminal(ledger, stage, unit) and not output_of(ledger, stage, unit)


def counted(a):
    """Attempts that say something about the tool (a cancel or a KBase transport error does not)."""
    return a.get("status") != "terminated" and not a.get("transient")


def retire_submit_errors(t):
    if t.get("submit_errors"):
        t["submit_errors_history"] = t.get("submit_errors_history", []) + t.pop("submit_errors")


def in_flight(t):
    """A submit that raised or was interrupted: this unit's job may already exist in EE2 (a pilot in flight)."""
    return bool(t.get("submit_errors") or t.get("submitting")) and t["state"] not in ("done", "empty", "failed")


def run_reads(ledger, run):
    """A run's trimmed reads; its raw import if trimming ended without output; None while undecided."""
    ref = output_of(ledger, "trim", run)
    if ref:
        return ref
    if dead(ledger, "trim", run):
        return output_of(ledger, "import", run)
    return None


def run_settled(ledger, run):
    return terminal(ledger, "trim", run)


def lib_settled(ledger, u):
    """Is the sample's assembly library decided (single run: trimmed; multi-run: merge finished)?"""
    if len(u["runs"]) == 1:
        return run_settled(ledger, u["runs"][0]["run"])
    return terminal(ledger, "merge", u["unit"])


def lib_of(ledger, u):
    """The one PairedEndLibrary that is assembled and mapped for this sample; None if undecided or lost."""
    if len(u["runs"]) == 1:
        return run_reads(ledger, u["runs"][0]["run"])
    return output_of(ledger, "merge", u["unit"])


def sync_imports(ledger, rows):
    """Mirror the import watcher's ledger (read only; it is replaced atomically) into import:<SRR> tasks. A run the
    import watcher gave up on is NOT mirrored as failed: its sample waits, since a re-import would be wanted."""
    try:
        with open(IMPORT_LEDGER) as fh:
            runs = json.load(fh)["runs"]
    except (OSError, ValueError, KeyError) as e:
        log(f"import: could not read the import ledger ({type(e).__name__}: {str(e)[:150]}); no change this cycle")
        return
    failed = []
    for r in rows:
        x = runs.get(r["run"]) or {}
        t = task(ledger, "import", r["run"])
        if x.get("state") == "done" and x.get("ref") and t["state"] != "done":
            t.update(state="done", output_ref=x["ref"], mirrored_utc=now())
        elif x.get("state") == "failed":
            failed.append(r["run"])
    if failed and failed != ledger.get("import_failed_logged"):
        ledger["import_failed_logged"] = failed
        log(f"import: {len(failed)} runs failed in the import ledger; their samples wait (and so does GTDB-Tk): {failed}")


def salterns_load():
    """(busy, live jobs per pool, gtdbtk state) from the salterns ledger; it is only read, never locked."""
    try:
        with open(SALTERNS_LEDGER) as fh:
            tasks = json.load(fh)["tasks"]
    except (OSError, ValueError, KeyError):
        return True, Counter(), "unknown"               # unreadable: assume salterns still needs the capacity
    busy = any(t["stage"] in SALTERNS_HEAVY and t["state"] not in ("done", "empty", "failed") for t in tasks.values())
    live = Counter(POOL.get(a.get("variant", t["stage"]), "?") for t in tasks.values() for a in t["attempts"]
                   if a.get("status") in LIVE and not a.get("local"))
    g = tasks.get("gtdbtk:ALL")
    return busy, live, (g or {}).get("state", "pending")


def import_live_njs():
    try:
        with open(IMPORT_LEDGER) as fh:
            return sum(1 for x in json.load(fh)["runs"].values() if x.get("state") == "active")
    except (OSError, ValueError, KeyError):
        return 0


def pool_caps(ledger):
    """Our per-pool caps: reduced while salterns still has heavy work; afterwards the shared total minus the live
    salterns and import-watcher jobs, never below the reduced caps."""
    busy, live_s, _ = salterns_load()
    if busy:
        caps = dict(POOL_CAP_BUSY)
    else:
        imp = import_live_njs()
        caps = {p: max(POOL_CAP_BUSY[p], POOL_CAP_TOTAL[p] - live_s[p] - (imp if p == "njs" else 0))
                for p in POOL_CAP_BUSY}
    key = json.dumps([busy, caps], sort_keys=True)
    if ledger.get("caps_logged") != key:
        ledger["caps_logged"] = key
        log(f"pool caps {caps} (salterns {'busy' if busy else 'cleared'}; salterns live {dict(live_s)})")
    return caps


def free_gb(path=WORK):
    import shutil
    return shutil.disk_usage(path).free / 1e9


# ------------------------------------------------------------------------------------------------ KBase helpers
def retry_rpc(fn, *args, tries=4, **kw):
    for k in range(tries):
        try:
            return fn(*args, **kw)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(10 * (k + 1))


class KB:
    def __init__(self, dry_run=False):
        self.u = KBGenomeUtils(kb_version="prod")
        self.jobs = KBJobUtils(self.u, kb_version="prod")
        self.ws = self.u.ws_client()
        self.ee2 = self.jobs.ee2_client
        self.nms = NMSSpecCache()
        self.runner = None if dry_run else AppRunner(ws=self.u, kb_version="prod", job_store=self.jobs,
                                                     nms_cache=self.nms)

    def states(self, job_ids):
        out = {}
        for i in range(0, len(job_ids), 100):
            res = retry_rpc(self.ee2.check_jobs, {"job_ids": job_ids[i:i + 100], "exclude_fields": ["job_input"]})
            for s in (res.get("job_states", res) if isinstance(res, dict) else res):
                if isinstance(s, dict) and s.get("job_id"):
                    out[s["job_id"]] = s
        return out

    def info(self, name_or_ref):
        obj = {"ref": name_or_ref} if "/" in name_or_ref else {"wsid": WS_ID, "name": name_or_ref}
        infos = retry_rpc(self.ws.get_object_info3, {"objects": [obj], "ignoreErrors": 1, "includeMetadata": 1})["infos"]
        return infos[0] if infos else None

    def object_ref(self, name):
        i = self.info(name)
        return f"{i[6]}/{i[0]}/{i[4]}" if i else None

    def service_ver(self, app_id, tag):
        return self.nms.get(app_id, None if tag == "release" else tag).service_ver

    def submit(self, app_id, tag, params, meta):
        handle = self.runner.run_app(app_id, params, WS_ID, tag=None if tag == "release" else tag, meta=meta)
        return handle.job_id

    def find_job(self, app_id, needle, since_utc, exclude=()):
        """The newest EE2 job of this app whose params contain `needle` verbatim and that EE2 created after
        since_utc - 2 min, not already in `exclude`: a submit whose reply, or whose ledger write, was lost."""
        since_ms = int(datetime.fromisoformat(since_utc).timestamp() * 1000) - 120_000
        res = retry_rpc(self.ee2.check_workspace_jobs, {"workspace_id": WS_ID, "return_list": 1,
                                                        "exclude_fields": ["job_output"]})
        for j in sorted(res.get("job_states", []), key=lambda j: j.get("created") or 0, reverse=True):
            ji = j.get("job_input") or {}
            if (ji.get("app_id") == app_id and (j.get("created") or 0) >= since_ms and j.get("job_id") not in exclude
                    and j.get("status") in (LIVE | {"completed"}) and json.dumps(needle) in json.dumps(ji.get("params"))):
                return j
        return None


def recon_key(stage, params, out_name):
    """A string the submitted params contain verbatim, to find the EE2 job of a lost submit. Never None or '':
    json.dumps(None) == 'null' would match any null parameter (e.g. MEGAHIT's "min_count": null)."""
    p0 = params[0] if isinstance(params, list) and params and isinstance(params[0], dict) else {}
    key = {"trim": p0.get("output_reads_name"),            # the app saves '<output_reads_name>_paired'
           "extract": p0.get("assembly_set_name"),          # present for one bin and for many
           "checkm": p0.get("input_ref")}.get(stage, out_name)
    return key or None


def err_full(s):
    e = s.get("error")
    return str((e.get("message") if isinstance(e, dict) else None) or s.get("errormsg") or "")


def clip(msg, n=2000):
    """Head and tail of a long message (multi-library job logs end with the informative lines)."""
    return msg if len(msg) <= n else msg[:n // 2] + " ... " + msg[-(n // 2):]


def worker_alive(pid, workdir):
    """Is a local worker still running? Our own exited child is reaped here, since a SIGKILLed worker would
    otherwise stay a zombie whose /proc entry looks alive; a worker started by an earlier watcher must still be
    running with this attempt's workdir on its command line (a reused PID does not count)."""
    if not pid or pid <= 0:
        return False
    try:
        return os.waitpid(pid, os.WNOHANG)[0] == 0
    except ChildProcessError:
        pass
    try:
        with open(f"/proc/{pid}/stat") as fh:
            if fh.read().rsplit(")", 1)[1].split()[0] in ("Z", "X"):
                return False
        with open(f"/proc/{pid}/cmdline", "rb") as fh:
            return workdir.encode() in fh.read()
    except (OSError, IndexError):
        return False


def scope_state(scope):
    """`systemctl --user is-active`: 'active', 'inactive', 'failed', ...; '' when the user bus is unreachable."""
    return (subprocess.run(["systemctl", "--user", "is-active", scope], capture_output=True, text=True).stdout
            or "").strip()


def mem_available_gb():
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1048576
    return 0.0


# ------------------------------------------------------------------------------------------------ parameters
def params_for(variant, u, ledger, kb, rows):
    """(service-shape params, output object name), None while inputs are pending, or 'skip' if they ended empty."""
    s = u["unit"]
    if variant == "trim":                                # u is a run; its import is in ws 275929 (by absolute ref)
        reads = output_of(ledger, "import", s)
        if not reads:
            return None                                  # imports are never 'dead' here: the sample waits
        base = name("trim", u)
        return [{"input_ws": READS_WS_NAME, "output_ws": WS_NAME, "input_reads_ref": reads, "output_reads_name": base,
                 "translate_to_phred33": 1,
                 "adapter_clip": {"adapterFa": "TruSeq3-PE-2.fa", "seed_mismatches": 2, "palindrome_clip_threshold": 30,
                                  "simple_clip_threshold": 10},
                 "sliding_window": {"sliding_window_size": 4, "sliding_window_min_quality": 15},
                 "leading_min_quality": 3, "trailing_min_quality": 3, "crop_length": 0, "head_crop_length": 0,
                 "min_length": 50}], base + "_paired"
    if variant == "merge":                               # multi-run samples only
        runs = [r["run"] for r in u["runs"]]
        if not all(run_settled(ledger, r) for r in runs):
            return None
        refs = [run_reads(ledger, r) for r in runs]
        if not all(refs):
            return "skip"
        assert len(set(refs)) == len(refs), f"{s}: duplicate read refs {refs}"
        kinds = ["Trimmomatic-paired" if output_of(ledger, "trim", r) else "raw (trim failed)" for r in runs]
        return [{"workspace_name": WS_NAME, "input_refs": refs, "output_name": name("merge", u),
                 "desc": f"{s}: {len(refs)} SRA runs ({', '.join(f'{r} {k}' for r, k in zip(runs, kinds))}) merged "
                         f"for a per-sample assembly"}], name("merge", u)
    if variant.startswith("megahit"):
        if not lib_settled(ledger, u):
            return None
        reads = lib_of(ledger, u)
        if not reads:
            return "skip"
        p = {"workspace_name": WS_NAME, "read_library_ref": reads, "output_contigset_name": name("assemble", u),
             "megahit_parameter_preset": "meta-sensitive", "min_count": None, "k_min": None, "k_max": None,
             "k_step": None, "k_list": [], "min_contig_length": 1000}
        if variant == "megahit_retry":
            p.update(megahit_parameter_preset=None, min_count=2, k_list=[27, 37, 47, 57, 67, 77, 87, 97, 107, 117, 127])
        elif variant == "megahit_default":
            p.update(megahit_parameter_preset=None)
        return [p], name("assemble", u)
    asm = output_of(ledger, "assemble", s)
    if variant in BINNERS:
        if not asm:
            return "skip" if dead(ledger, "assemble", s) else None
        reads = lib_of(ledger, u)                        # per-sample coverage: the library that was assembled
        if not reads:
            return "skip"
        base = {"workspace_name": WS_NAME, "assembly_ref": asm, "binned_contig_name": name(variant, u),
                "reads_list": [reads]}
        if variant == "metabat":
            return [{**base, "min_contig_length": 2500}], name(variant, u)
        if variant == "maxbin":
            return [{**base, "thread": 32, "prob_threshold": 0.8, "markerset": "40", "min_contig_length": 1000,
                     "plotmarker": 0}], name(variant, u)
        return [{**base, "read_mapping_tool": "bowtie2_default", "min_contig_length": 2500, "contig_split_size": 10000,
                 "contig_split_overlap": 0, "kmer_size": 4, "max_clusters_for_vgmm": 400, "max_iterations_for_vgmm": 500,
                 "total_percentage_pca": 90, "no_cov_normalization": 0, "no_total_coverage": 0}], name(variant, u)
    if variant == "dastool":
        if not asm:
            return "skip" if dead(ledger, "assemble", s) else None
        if not all(terminal(ledger, b, s) for b in BINNERS):
            return None
        binners = [f"{BINNER_LABEL[b]}={output_of(ledger, b, s)}" for b in BINNERS if output_of(ledger, b, s)]
        if not binners:
            return "skip"
        return {"assembly_ref": asm, "binners": binners}, name("dastool", u)
    if variant == "extract":
        das = output_of(ledger, "dastool", s)
        if not das:
            return "skip" if dead(ledger, "dastool", s) else None
        pred = json.load(open(get_task(ledger, "dastool", s)["predicted"]))["predicted_assemblies"]
        bids = [p["bid"] for p in pred]
        out = name("extract", u) if len(bids) > 1 else f"{bids[0]}__{s}"   # one bin: no AssemblySet is saved
        return [{"workspace_name": WS_NAME, "binned_contig_obj_ref": das, "extracted_assemblies": ",".join(bids),
                 "assembly_suffix": f"__{s}", "assembly_set_name": name("extract", u)}], out
    if variant == "checkm":
        target = output_of(ledger, "extract", s)
        if not target:
            return "skip" if dead(ledger, "extract", s) else None
        return [{"workspace_name": WS_NAME, "input_ref": target, "reduced_tree": "1", "save_plots_dir": "0",
                 "threads": "4"}], None
    raise ValueError(variant)


def trailing_transient(att):
    """The attempts at the end of the list that failed with a KBase transport error (local workers only)."""
    run = []
    for a in reversed(att):
        if not a.get("transient"):
            break
        run.append(a)
    return run[::-1]


def choose_variant(stage, t):
    """Variant for the next attempt, or None if the task should give up."""
    att = t["attempts"]
    if sum(1 for a in att if not a.get("transient")) >= HARD_CAP_ATTEMPTS:
        return None
    tr = trailing_transient(att)                       # an outage gets a time budget, not a count
    if tr and tr[0].get("local") and \
            time.time() - datetime.fromisoformat(tr[0]["submitted_utc"]).timestamp() > TRANSIENT_BUDGET_S:
        return None
    if tr and not tr[0].get("local") and time.time() - datetime.fromisoformat(       # an EE2 job may queue for hours:
            tr[0].get("finished_utc") or tr[0]["submitted_utc"]).timestamp() > KBASE_TRANSIENT_BUDGET_S:  # from its end
        return None
    if stage == "assemble":                            # the salterns co-assembly ladder
        intrinsic = [a for a in att if counted(a) and a.get("intrinsic")]
        other = [a for a in att if counted(a) and not a.get("intrinsic")]
        if len(other) >= MAX_TRIES:
            return None
        return ("megahit", "megahit_retry", "megahit_default", None)[min(len(intrinsic), 3)]
    return stage if sum(map(counted, att)) < MAX_TRIES else None


# ------------------------------------------------------------------------------------------------ completion checks
def validate_output(kb, stage, ref, params):
    """(ok, reason) for an output object, fresh or adopted: assembly size floor; merged reads/bases == input sums."""
    if stage == "assemble":
        size = int((kb.info(ref)[10] or {}).get("Size") or 0)
        floor = MIN_ASSEMBLY_BP
        if size < floor:
            return False, f"assembly {ref} is only {size} bp (< {floor})"
    if stage == "merge":
        ins = params[0]["input_refs"]
        infos = retry_rpc(kb.ws.get_object_info3, {"objects": [{"ref": r} for r in [ref] + ins], "includeMetadata": 1})["infos"]
        m = [i[10] or {} for i in infos]
        got = (int(m[0]["read_count"]), int(m[0]["total_bases"]))
        exp = (sum(int(x["read_count"]) for x in m[1:]), sum(int(x["total_bases"]) for x in m[1:]))
        if got != exp or got[0] % 2:
            return False, f"merged reads/bases {got} != input sums {exp}"
    return True, None


def on_completed(kb, ledger, key, t, a, s):
    """Validate a finished KBase job; return (ok, output_ref, reason)."""
    stage, unit = t["stage"], t["sample"]
    result = ((s.get("job_output") or {}).get("result") or [{}])
    result = result[0] if isinstance(result, list) and result else {}
    if stage == "extract":
        refs = result.get("assembly_ref_list") or []
        target = result.get("assembly_set_ref") or (refs[0] if len(refs) == 1 else None)
        pred = json.load(open(get_task(ledger, "dastool", unit)["predicted"]))["predicted_assemblies"]
        if not target or len(refs) != len(pred):
            return False, None, f"extract made {len(refs)} assemblies for {len(pred)} bins"
        infos = retry_rpc(kb.ws.get_object_info3, {"objects": [{"ref": r} for r in refs], "includeMetadata": 1})["infos"]
        want = {p.get("assembly_name") or f"{p['bid']}__{unit}": p for p in pred}   # dastool_to_binnedcontigs manifest
        bad, md5_other = [], []
        for i in infos:
            m, p = i[10] or {}, want.get(i[1])
            if (not p or ("dna_size" in p and int(m.get("Size") or -1) != p["dna_size"])
                    or ("num_contigs" in p and int(m.get("N Contigs") or -1) != p["num_contigs"])):
                bad.append(i[1])
            elif p.get("md5") and m.get("MD5") != p["md5"]:
                md5_other.append(i[1])
        if bad:
            return False, None, f"extracted assemblies differ from the DAS Tool bins (name, Size or N Contigs): {bad[:3]}"
        if md5_other:                    # logged only: CheckM and GTDB-Tk both use the KBase copy (collect_checkm_select)
            log(f"{key}: {len(md5_other)} extracted assemblies have an MD5 other than predicted: {md5_other[:3]}")
        return True, target, None
    if stage == "checkm":
        rep = result.get("report_ref")
        return (True, rep, None) if rep else (False, None, "CheckM job returned no report_ref")
    out_name = a.get("output_name")
    ref = kb.object_ref(out_name) if out_name else None
    if ref in ledger.get("superseded_refs", ()):       # only the metaSPAdes-era version exists: nothing new was saved
        ref = None
    if out_name and not ref:
        return None, None, f"completed without {out_name}"
    ok, reason = validate_output(kb, stage, ref, a["params"])
    return (True, ref, None) if ok else (False, ref, reason)


def collect_checkm(ledger, unit, report_ref, target_ref):
    """Download CheckM tables, pick >= 50 / < 10 bins, write their FASTA and a GTDB-Tk batchfile (local)."""
    out = os.path.join(WORK, "checkm", unit)
    das = get_task(ledger, "dastool", unit)
    cmd = [sys.executable, os.path.join(LOCAL, "collect_checkm_select.py"), "--report", report_ref, "--set", target_ref,
           "--tag", unit, "--out", out]
    if das and das.get("bins_dir"):
        cmd += ["--local-fasta-dir", das["bins_dir"]]
    rc = subprocess.run(cmd, stdout=open(os.path.join(WORK, "checkm", f"{unit}.collect.log"), "w"),
                        stderr=subprocess.STDOUT).returncode
    return rc == 0, os.path.join(out, f"gtdbtk_batchfile_{unit}.tsv")


# ------------------------------------------------------------------------------------------------ polling
def poll(kb, ledger, dry_run=False):
    for t in ledger["tasks"].values():                 # self-heal attempts stranded by the v1 poll ordering
        if t["state"] == "active" and t["attempts"] and t["attempts"][-1].get("status") not in LIVE \
                and t["attempts"][-1].get("job_id") and not t["attempts"][-1].get("local"):
            t["attempts"][-1]["status"] = "running"
    live = {a["job_id"]: (k, a) for k, t in ledger["tasks"].items() for a in t["attempts"]
            if a.get("status") in LIVE and a.get("job_id") and not a.get("local")}
    st = kb.states(list(live)) if live else {}
    for jid, (key, a) in live.items():
        s = st.get(jid)
        if not s:
            continue
        t, status = ledger["tasks"][key], s["status"]
        if status in LIVE:
            a["status"] = status
            t["state"] = "active"
            since = s.get("running") or s.get("queued") or s.get("created")
            age_h = (time.time() * 1000 - since) / 3.6e6 if since else 0
            if age_h > STALE_H.get(status, STALE_H["queued"]) and not a.get("stale_logged"):
                a["stale_logged"] = now()
                log(f"{key}: job {jid} {status} for {age_h:.0f} h - check its log; cancel it if silent (it is retried)")
            continue
        if status == "completed":
            try:
                ok, ref, reason = on_completed(kb, ledger, key, t, a, s)
            except Exception as e:                     # lookup failed: leave the attempt live, re-poll next cycle
                a["check_failures"] = a.get("check_failures", 0) + 1
                if a["check_failures"] < CHECK_FAILURES_MAX:
                    log(f"{key}: completion check failed ({type(e).__name__}: {str(e)[:200]}); re-polled next cycle")
                    continue
                a.update(status="error", finished_utc=now(), error=f"completion check failed {a['check_failures']} "
                                                                  f"times: {type(e).__name__}: {str(e)[:300]}")
                t["state"] = "retry"                   # never re-poll forever; a plain retry, not a rejection
                log(f"{key}: {a['error'][:300]} - retried")
                continue
            a["status"], a["finished_utc"] = status, now()
            if ok:
                t.update(state="done", output_ref=ref, job_id=jid)
                log(f"{key}: done -> {ref}")
                if t["stage"] == "checkm":
                    t["collect_pending"] = True
            elif ok is None:
                t["state"] = "empty"
                log(f"{key}: {reason}")
            else:
                a.update(status="error", error=reason, rejected=True, rejected_ref=ref)
                if t["stage"] == "assemble":
                    a["intrinsic"] = True            # walk the MEGAHIT ladder (min-count 2 + k-list, then defaults)
                t["state"] = "failed" if t["stage"] == "merge" else "retry"   # the same merge cannot change the counts
                log(f"{key}: rejected: {reason}")
        else:
            full = err_full(s)
            a["status"], a["finished_utc"], a["error"] = status, now(), clip(full)
            if status == "error" and any(m in full[:2000] for m in KBASE_TRANSIENT_MARKERS):
                a["transient"] = True                  # a KBase-side kill or outage: not counted, 72-h budget
                t["state"] = "retry"
                log(f"{key}: transient KBase error (attempt {len(t['attempts'])}, not counted): {a['error'][:160]}")
                continue
            if t["stage"] == "assemble" and ("Error running MEGAHIT" in full[:2000] or "return code" in full[:2000]
                                             or any(m in full[:2000] for m in OOM_MARKERS)):
                a["intrinsic"] = True                  # next rung: min-count 2 (less memory), then defaults
            if t["stage"] == "maxbin" and status == "error" and "dataset cannot be binned" in full:
                t.update(state="empty", reason="MaxBin2: too few marker genes in this assembly to bin it")
                log(f"{key}: empty - {t['reason']}")   # deterministic; DAS Tool goes on with the other binners
                continue
            t["state"] = "retry"
            log(f"{key}: {status} (attempt {len(t['attempts'])}): {a['error'][:200]}")
    for key, t in ledger["tasks"].items():             # local DAS Tool workers (gtdbtk_stage polls its own attempt)
        if t["stage"] != "dastool":
            continue
        for a in t["attempts"]:
            if not (a.get("local") and a.get("status") in LIVE and a.get("workdir")):
                continue
            sp = os.path.join(a["workdir"], "status.json")
            try:
                stj = json.load(open(sp)) if os.path.exists(sp) else {"state": "running"}
            except ValueError:                         # status.json is replaced atomically, so only if hand-edited
                stj = {"state": "running"}
            pid = a.get("pid") or stj.get("pid")       # pid None: the ledger write after the launch was lost
            if stj["state"] == "running" and worker_alive(pid, a["workdir"]):
                age_h = (time.time() - datetime.fromisoformat(a["submitted_utc"]).timestamp()) / 3600
                if age_h > LOCAL_STALE_H and not a.get("stale_logged"):
                    a["stale_logged"] = now()
                    log(f"{key}: local worker {pid} still running after {age_h:.0f} h - see {a['workdir']}/worker.log")
                continue
            a["finished_utc"] = now()
            if stj["state"] == "done":
                a["status"] = "completed"
                t.update(state="done", output_ref=stj["output_ref"], predicted=stj["predicted"],
                         bins_dir=stj["bins_dir"], n_bins=stj["n_bins"])
                log(f"{key}: done -> {stj['output_ref']} ({stj['n_bins']} bins)")
            elif stj["state"] == "empty":
                a["status"] = "completed"
                t.update(state="empty", reason=stj.get("reason"))
                log(f"{key}: empty - {stj.get('reason')}")
            else:
                a["status"], a["error"] = "error", stj.get("error") or "worker exited without a status"
                if any(m in a["error"] for m in TRANSIENT_MARKERS):
                    a["transient"] = True              # a KBase outage says nothing about the unit: not counted
                t["state"] = "retry"
                log(f"{key}: local worker failed{' (transient)' if a.get('transient') else ''}: {a['error'][:200]}")
    if not dry_run:
        hybrid.poll(sys.modules[__name__], kb, ledger)  # local hybrid steps (assemble, binners, extract, CheckM)
    if dry_run:                                        # plan: never collect alongside the watcher
        return
    for key, t in ledger["tasks"].items():             # CheckM collection (local, after the job is done)
        if t.get("collect_pending"):
            if free_gb() < MIN_FREE_GB:
                log(f"{key}: CheckM collection waits: less than {MIN_FREE_GB} GB free in {WORK}")
                continue
            target = get_task(ledger, "extract", t["sample"])["output_ref"]
            raised = False
            try:
                ok, batch = collect_checkm(ledger, t["sample"], t["output_ref"], target)
            except Exception as e:
                ok, batch, raised = False, None, True
                log(f"{key}: CheckM collection raised {type(e).__name__}: {str(e)[:200]}")
            if ok:
                t.update(collect_pending=False, batchfile=batch)
                log(f"{key}: CheckM collected -> {batch}")
            else:
                last = ""
                if not raised:                         # collect.log was rewritten by this try: its final line is
                    try:                               # the exception (earlier lines may be recovered retries)
                        with open(os.path.join(WORK, "checkm", f"{t['sample']}.collect.log"), errors="replace") as fh:
                            lines = [x for x in fh.read()[-4000:].splitlines() if x.strip()]
                        last = lines[-1] if lines else ""
                    except OSError:
                        pass
                if any(m in last for m in TRANSIENT_MARKERS) and t.get("collect_transient", 0) < COLLECT_TRANSIENT_MAX:
                    t["collect_transient"] = t.get("collect_transient", 0) + 1     # a KBase outage: not counted
                    log(f"{key}: CheckM collection hit a KBase transport error ({last[:150]}); retried next cycle")
                    continue
                t["collect_failures"] = t.get("collect_failures", 0) + 1
                if t["collect_failures"] >= 3:
                    t.update(collect_pending=False, collect_failed=True)
                    log(f"{key}: CheckM collection failed 3 times - unit left out of GTDB-Tk (see work/checkm/*.collect.log)")


# ------------------------------------------------------------------------------------------------ advance
def guard_alive(pid):
    try:
        with open(f"/proc/{pid}/cmdline", "rb") as fh:   # a zombie or a reused PID reads empty / other
            return bool(pid) and b"mem_guard.sh" in fh.read()
    except (OSError, TypeError):
        return False


def launch_guard(scope):
    """mem_guard.sh for a GTDB-Tk scope: FREEZE/THAW it when the shared machine runs out of memory (logs in GTDBTK_DIR)."""
    return subprocess.Popen(["bash", os.path.join(LOCAL, "mem_guard.sh")], env=dict(os.environ, UNIT=scope, GUARD_DIR=GTDBTK_DIR),
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True).pid


def launch_local_dastool(u, prm, wd):
    os.makedirs(wd, exist_ok=True)
    cmd = [sys.executable, os.path.join(LOCAL, "dastool_worker.py"), "--unit", u["unit"], "--assembly-ref",
           prm["assembly_ref"], "--workdir", wd] + [x for b in prm["binners"] for x in ("--binner", b)]
    p = subprocess.Popen(cmd, stdout=open(os.path.join(wd, "stdout.log"), "w"), stderr=subprocess.STDOUT,
                         start_new_session=True)
    return p.pid


def gtdbtk_stage(ledger, dry_run, planned, samples):
    t = task(ledger, "gtdbtk", "ALL")
    ck = [f"checkm:{u['unit']}" for u in samples]       # every sample, not just those that have reached CheckM
    if t["state"] in ("done", "failed"):
        return
    if not all(k in ledger["tasks"] and terminal(ledger, "checkm", k.split(":", 1)[1])
               and not ledger["tasks"][k].get("collect_pending") for k in ck):
        return
    scope = f"{GTDBTK_SCOPE}.scope"
    out = os.path.join(GTDBTK_DIR, "gtdbtk_out")
    code_file = out + ".exit_code"
    if t["state"] == "active":
        if dry_run:
            return
        running = scope_state(scope)
        if running in SCOPE_RUNNING:
            a = t["attempts"][-1]
            if not guard_alive(a.get("guard_pid")):     # it died with its scope, or the run was adopted: the new guard
                a["guard_pid"] = launch_guard(scope)    # adopts a freeze left behind and decides when to thaw
                log(f"gtdbtk: mem_guard for {scope} was not running - restarted (pid {a['guard_pid']})")
            return
        if running not in SCOPE_ENDED:                  # never conclude a run we cannot see
            log(f"gtdbtk: 'systemctl --user is-active' gave {running!r} (user bus unreachable?); not concluding")
            return
        a = t["attempts"][-1]
        result = subprocess.run(["systemctl", "--user", "show", "-p", "Result", "--value", scope], capture_output=True,
                                text=True).stdout.strip()    # 'oom-kill' if the scope was stopped by an OOM kill
        rc = open(code_file).read().strip() if os.path.exists(code_file) else "missing"
        a.update(status="completed" if rc == "0" else "error", finished_utc=now(), exit_code=rc, scope_state=running,
                 scope_result=result)
        if guard_alive(a.get("guard_pid")):             # only if the PID is still that mem_guard (no PID reuse)
            try:
                os.kill(a["guard_pid"], signal.SIGTERM)
            except OSError:
                pass
        if rc == "0":
            dest = os.path.join(HERE, "gtdbtk")
            os.makedirs(dest, exist_ok=True)
            for f in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
                if os.path.exists(os.path.join(out, f)):
                    subprocess.run(["cp", os.path.join(out, f), dest])
            t.update(state="done", output_ref=dest)
            log(f"gtdbtk: done -> {dest}")
        else:
            t["state"] = "retry" if sum(map(counted, t["attempts"])) < 2 else "failed"
            log(f"gtdbtk: classify_wf exit {rc} (scope {running}, Result={result}) -> {t['state']}")
        return
    batches = [ledger["tasks"][k].get("batchfile") for k in ck if ledger["tasks"][k]["state"] == "done"]
    lines = []
    for b in batches:
        if b and os.path.exists(b):
            lines += [x for x in open(b).read().splitlines() if x.strip()]
    nf = sum(1 for k in ck if ledger["tasks"][k].get("collect_failed"))
    if not lines:
        state = "failed" if nf else "empty"
        if t["state"] != state:
            t["state"] = state
            log(f"gtdbtk: nothing to classify ({nf} units failed CheckM collection) -> {state}")
        return
    planned.append(("gtdbtk", "ALL", "local", "GTDB-Tk 2.8.0/R232", "local"))
    if dry_run:
        return
    st = scope_state(scope)
    _, _, salt_g = salterns_load()                       # salterns has priority: never alongside its GTDB-Tk
    salt_scope = scope_state(f"{SALTERNS_GTDBTK_SCOPE}.scope")
    if st not in SCOPE_RUNNING and (salt_g not in ("done", "failed", "empty") or salt_scope in SCOPE_RUNNING):
        if t.get("salterns_wait_logged") != salt_g:
            t["salterns_wait_logged"] = salt_g
            log(f"gtdbtk: {len(lines)} genomes ready; waiting for the salterns GTDB-Tk (ledger {salt_g}, scope "
                f"{salt_scope or '?'})")
        return
    if st not in SCOPE_RUNNING and free_gb(GTDBTK_DIR if os.path.isdir(GTDBTK_DIR) else WORK) < MIN_FREE_GB:
        log(f"gtdbtk: less than {MIN_FREE_GB} GB free; not launching")
        return
    if st in SCOPE_RUNNING:                              # launched by a cycle whose ledger write was lost: adopt it
        t["attempts"].append({"local": True, "status": "running", "submitted_utc": now(), "out_dir": out, "adopted": True})
        t["state"] = "active"
        log(f"gtdbtk: adopted the running scope {scope}")
        return
    if st not in SCOPE_ENDED:                            # systemd-run would fail too
        log(f"gtdbtk: 'systemctl --user is-active' gave {st!r} (user bus unreachable?); not launching")
        return
    avail = mem_available_gb()
    if avail < GTDBTK_MIN_AVAIL_GB:
        last = t.get("waiting_logged")
        if not last or (datetime.now(timezone.utc) - datetime.fromisoformat(last)).total_seconds() >= 6 * 3600:
            t["waiting_logged"] = now()
            log(f"gtdbtk: {len(lines)} genomes ready; waiting for MemAvailable >= {GTDBTK_MIN_AVAIL_GB} GB "
                f"(now {avail:.0f} GB)")
        return
    os.makedirs(GTDBTK_DIR, exist_ok=True)
    batch = os.path.join(GTDBTK_DIR, "batchfile.tsv")
    with open(batch, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    subprocess.run(["systemctl", "--user", "reset-failed", scope], capture_output=True)  # an OOM-killed scope stays loaded
    if os.path.exists(code_file):
        os.remove(code_file)
    pplacer = "1" if any(counted(x) for x in t["attempts"]) else "8"   # a retry: least pplacer memory
    cmd = ["systemd-run", "--user", "--scope", f"--unit={GTDBTK_SCOPE}", "-p", "MemoryMax=170G", "-p",
           "MemorySwapMax=0", "bash", os.path.join(LOCAL, "run_gtdbtk_reassembly.sh"), batch, out]
    a = {"local": True, "pid": None, "guard_pid": None, "batchfile": batch, "n_genomes": len(lines),
         "pplacer_cpus": int(pplacer), "units_failed_collection": nf, "status": "running", "submitted_utc": now(),
         "out_dir": out}
    t["attempts"].append(a)
    t["state"] = "active"
    save_ledger(ledger)              # write-ahead: a lost write can never make a later cycle rerun a finished run
    p = subprocess.Popen(cmd, stdout=open(os.path.join(GTDBTK_DIR, "gtdbtk_run.log"), "a"), stderr=subprocess.STDOUT,
                         start_new_session=True, env=dict(os.environ, PPLACER_CPUS=pplacer))
    a.update(pid=p.pid, guard_pid=launch_guard(scope))
    log(f"gtdbtk: launched classify_wf on {len(lines)} genomes in scope {scope} (MemAvailable {avail:.0f} GB, "
        f"pplacer_cpus {pplacer}; {nf} units left out after failed CheckM collection)")


SMALL_GATED = ("assemble", "metabat", "maxbin", "concoct")


def pilot_gates(ledger, samples):
    """Sample pilots: which stages may go beyond the two pilot samples, and pilot tasks that ended badly."""
    pil = [u for u in samples if u["unit"] in PILOT_SAMPLES]
    runs = [r["run"] for u in pil for r in u["runs"]]
    multi = [u["unit"] for u in pil if len(u["runs"]) > 1]
    trim_ok = all(output_of(ledger, "trim", r) for r in runs)
    merge_ok = trim_ok and all(output_of(ledger, "merge", s) for s in multi)
    bin_ok = all(output_of(ledger, "assemble", u["unit"]) and all(terminal(ledger, b, u["unit"]) for b in BINNERS)
                 and any(output_of(ledger, b, u["unit"]) for b in BINNERS) for u in pil)
    bad = [k for k, t in ledger["tasks"].items() if t["state"] in ("failed", "empty") and t["stage"] in
           ("trim", "merge", "assemble") + BINNERS and (t["sample"] in PILOT_SAMPLES or t["sample"] in runs)]
    small = [u for u in pil if u["unit"] == PILOT_SAMPLES[0]]   # deadrun: the full chain proven (2026-11-12)
    small_ok = all(output_of(ledger, "assemble", u["unit"]) and all(terminal(ledger, b, u["unit"]) for b in BINNERS)
                   and any(output_of(ledger, b, u["unit"]) for b in BINNERS) for u in small) and bool(small)
    gates = {"trim": trim_ok, "merge": merge_ok, "assemble": bin_ok, "metabat": bin_ok, "maxbin": bin_ok,
             "concoct": bin_ok}
    gates.update({"small_" + st: small_ok for st in SMALL_GATED})   # samples < 50 Gbp need only deadrun
    return gates, bad


def advance(kb, ledger, rows, through, dry_run=False):
    poll(kb, ledger, dry_run)
    if not dry_run:
        save_ledger(ledger)
    samples, run_units = units(rows)
    stages = STAGES[:STAGES.index(through) + 1]
    sync_imports(ledger, rows)
    caps = pool_caps(ledger)
    gates, pilot_bad = pilot_gates(ledger, samples)
    if pilot_bad and not all(gates.values()) and ledger.get("pilot_bad_logged") != pilot_bad and not dry_run:
        ledger["pilot_bad_logged"] = pilot_bad
        log(f"sample pilot: {pilot_bad} ended without output - fan-out held where it depends on them; inspect")
    live = Counter(POOL[a.get("variant", t["stage"])] for t in ledger["tasks"].values() for a in t["attempts"]
                   if a.get("status") in LIVE and not a.get("local"))
    local_live = sum(1 for t in ledger["tasks"].values() for a in t["attempts"] if a.get("local") and a.get("status") in LIVE
                     and t["stage"] == "dastool")
    planned = []
    order = list(reversed(stages))                       # downstream first: finishing a sample beats starting one
    if "trim" in order and sum(1 for t in ledger["tasks"].values() if t["stage"] == "trim"
                               and t["state"] in ("pending", "retry")) < 10:
        order.remove("trim")                             # a few (re)trims left: they go first, or ~120 queued MEGAHIT
        order.insert(0, "trim")                          # jobs would hold bigmem and starve the samples they block
    if "assemble" in order:                              # hybrid: the one local MEGAHIT slot is filled first (the long
        order.remove("assemble")                         # pole); the binners and extracts fill the rest of the gate
        order.insert(1 if order and order[0] == "trim" else 0, "assemble")
    cyc = hybrid.Cycle(sys.modules[__name__], ledger, DRIVER)
    cyc.extra = (rows,)
    cyc.extract_name = lambda u: name("dastool", u)
    for stage in order:
        if stage == "import":                            # (168 trims would otherwise fill bigmem before MetaBAT2)
            continue                                     # mirrored above; never submitted by this driver
        if stage == "gtdbtk":
            gtdbtk_stage(ledger, dry_run, planned, samples)
            continue
        if stage == "trim":
            todo = run_units
        elif stage == "merge":
            todo = [u for u in samples if len(u["runs"]) > 1]
        else:
            todo = samples
        if stage in gates and not gates[stage]:          # sample pilot: only the two pilot samples for now
            for u in todo:
                task(ledger, stage, u["unit"])
            small_ok = stage in SMALL_GATED and gates.get("small_" + stage)
            todo = [u for u in todo if (u.get("sample") or u["unit"]) in PILOT_SAMPLES
                    or (small_ok and sum(int(r["bases"]) for r in u["runs"]) < BIG_SAMPLE_BP)]
        # pilot: one unit at a time until one succeeds (CheckM: until its tables were also collected); a pilot
        # that fails, or whose CheckM collection fails, holds the stage for inspection
        tasks_ = [x for x in ledger["tasks"].values() if x["stage"] == stage]
        sts = Counter(x["state"] for x in tasks_)
        ok_done = sum(1 for x in tasks_ if x["state"] == "done" and (stage != "checkm" or x.get("batchfile")))
        bad = sts["failed"] + sum(1 for x in tasks_ if x.get("collect_failed"))
        fan_out = stage not in PILOT or ok_done > 0
        allow_new = fan_out or (sts["active"] == 0 and sts["retry"] == 0 and bad == 0
                                and not any(x.get("collect_pending") for x in tasks_)
                                and not any(in_flight(x) for x in tasks_))
        if stage in PILOT and not fan_out and bad and not dry_run:
            log(f"{stage}: pilot failed - stage held; inspect the failed task before fanning out")
        for u in todo:
            s = u["unit"]
            t = task(ledger, stage, s)
            if t["state"] in ("done", "empty", "failed"):
                continue
            if t["state"] == "active" and not (stage in HYBRID and hybrid.queued_kbase(t)):
                continue                                 # (a queued KBase job of a hybrid stage may be taken over)
            if not allow_new and not t["attempts"] and not in_flight(t):   # an in-flight pilot is reconciled below
                continue
            variant = choose_variant(stage, t)
            if variant is None:
                t["state"] = "failed"
                log(f"{stage}:{s}: giving up after {len(t['attempts'])} attempts")
                if stage in PILOT and not fan_out:      # a failed pilot: no second unit in this cycle
                    allow_new = False
                continue
            got = params_for(variant, u, ledger, kb, rows)
            if got is None:
                continue
            if got == "skip":
                t.update(state="empty", reason="upstream ended without output")
                log(f"{stage}:{s}: skipped - upstream ended without output")
                continue
            params, out_name = got
            ref = kb.object_ref(out_name) if out_name and stage not in ("dastool", "extract") else None
            if ref and ref in ledger.get("superseded_refs", ()):
                ref = None                               # a metaSPAdes-era version: never adopted (2026-11-10)
            if ref and not any(x.get("rejected_ref") == ref for x in t["attempts"]):
                try:                                     # an existing output is checked like a fresh one
                    ok, reason = validate_output(kb, stage, ref, params)
                except Exception as e:
                    t["adopt_check_failures"] = t.get("adopt_check_failures", 0) + 1
                    if t["adopt_check_failures"] < CHECK_FAILURES_MAX:
                        log(f"{stage}:{s}: check of existing {out_name} failed ({type(e).__name__}: {str(e)[:200]}); "
                            f"next cycle")
                        continue
                    t.update(state="failed", reason=f"check of existing {out_name} raised {t['adopt_check_failures']} "
                                                    f"times: {type(e).__name__}: {str(e)[:300]}")
                    log(f"{stage}:{s}: {t['reason']} - failed")   # downstream goes 'empty'; barriers still open
                    continue
                t.pop("adopt_check_failures", None)
                if ok:
                    t.update(state="done", output_ref=ref, adopted=now())
                    log(f"{stage}:{s}: adopted existing {out_name} -> {ref}")
                    continue
                log(f"{stage}:{s}: existing {out_name} ({ref}) not adopted: {reason}")
            if stage == "dastool":
                tr = trailing_transient(t["attempts"])  # back off during a KBase outage: 0.5, 1, 2, 4, 4 ... h
                if tr and tr[-1].get("finished_utc") and time.time() - datetime.fromisoformat(
                        tr[-1]["finished_utc"]).timestamp() < min(1800 * 2 ** (len(tr) - 1), 4 * 3600) - 120:
                    continue
                if local_live >= LOCAL_CAP:
                    continue
                if free_gb() < MIN_FREE_GB:
                    if not dry_run:
                        log(f"dastool:{s}: waits - less than {MIN_FREE_GB} GB free in {WORK}")
                    continue
                planned.append((stage, s, variant, "local DAS Tool 1.1.7", "local"))
                if not fan_out:
                    allow_new = False
                local_live += 1
                if dry_run:
                    continue
                stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                wd = os.path.join(WORK, "dastool", s, f"attempt{len(t['attempts']) + 1}-{stamp}")   # never reused
                a = {"local": True, "pid": None, "workdir": wd, "variant": "dastool", "params": params,
                     "status": "running", "output_name": out_name, "submitted_utc": now()}
                t["attempts"].append(a)
                t["state"] = "active"
                save_ledger(ledger)                     # write-ahead: a full disk fails here, before the worker starts
                a["pid"] = launch_local_dastool(u, params, wd)
                save_ledger(ledger)                     # if this write is lost, poll() takes the pid from status.json
                log(f"dastool:{s}: launched local worker (pid {a['pid']})")
                continue
            if stage in HYBRID:                          # never a new KBase job: local (or a queued one taken over)
                if not fan_out and not dry_run:
                    allow_new = False
                cyc.try_run(kb, stage, u, t, variant, params, out_name, s, dry_run)
                continue
            app_id, tag = APPS[variant]
            pool = POOL[variant]
            outstanding = t.get("submit_errors", []) + ([t["submitting"]] if t.get("submitting") else [])
            if outstanding and not dry_run:             # a submit that raised, or whose ledger write was lost:
                if not fan_out:                         # this unit is the pilot, possibly already running in EE2
                    allow_new = False
                needle = recon_key(stage, params, out_name)   # did EE2 accept it after all?
                since = min(e.get("started_utc") or e["utc"] for e in outstanding)
                try:
                    orphan = (kb.find_job(app_id, needle, since, exclude={x.get("job_id") for x in t["attempts"]})
                              if needle else None)
                except Exception as e:
                    log(f"{stage}:{s}: orphan lookup failed ({type(e).__name__}: {str(e)[:200]}); next cycle")
                    continue
                if orphan:
                    t["attempts"].append({"job_id": orphan["job_id"], "variant": variant, "app_id": app_id, "tag": tag,
                                          "params": params, "output_name": out_name, "submitted_utc": now(),
                                          "status": "queued", "reconciled": True})
                    t.pop("submitting", None)
                    retire_submit_errors(t)
                    t["state"] = "active"
                    live[pool] += 1
                    save_ledger(ledger)
                    log(f"{stage}:{s}: adopted orphan job {orphan['job_id']}")
                    continue
                if t.pop("submitting", None) is not None:
                    log(f"{stage}:{s}: an interrupted submit left no EE2 job; submitting afresh")
            if live[pool] >= caps[pool]:
                continue
            planned.append((stage, s, variant, app_id, tag))
            if not fan_out:
                allow_new = False
            live[pool] += 1
            if dry_run:
                continue
            t0 = now()
            t["submitting"] = {"utc": t0, "started_utc": t0, "variant": variant}
            save_ledger(ledger)                         # write-ahead: a full disk fails here, before anything is submitted
            try:
                ver = kb.service_ver(app_id, tag)
                jid = kb.submit(app_id, tag, params, meta={"stage": stage, "unit": s})
            except Exception as e:
                log(f"{stage}:{s}: submit raised {type(e).__name__}: {str(e)[:300]}")
                t.pop("submitting", None)
                t.setdefault("submit_errors", []).append({"utc": now(), "started_utc": t0, "variant": variant,
                                                          "error": str(e)[:1000]})
                save_ledger(ledger)
                continue
            t["attempts"].append({"job_id": jid, "variant": variant, "app_id": app_id, "tag": tag, "service_ver": ver,
                                  "params": params, "output_name": out_name, "submitted_utc": t0, "status": "queued"})
            t.pop("submitting", None)
            retire_submit_errors(t)
            t["state"] = "active"
            save_ledger(ledger)
            log(f"{stage}:{s}: submitted {app_id}@{tag} ({ver[:8]}) -> {jid}")
    if not dry_run:
        cyc.finish()
        needs = [lib_of(ledger, u) for u in samples if lib_of(ledger, u) and any(
            not terminal(ledger, st, u["unit"]) for st in ("assemble",) + BINNERS)]
        hybrid.write_needs(DRIVER, needs)
        hybrid.G.reads_cache_cleanup(hybrid.READS_CACHE, needs, log)
        save_ledger(ledger)
    return planned


def status(ledger):
    by = defaultdict(Counter)
    for t in ledger["tasks"].values():
        by[t["stage"]][t["state"]] += 1
    for stage in STAGES:
        if by.get(stage):
            print(f"{stage:10s} " + "  ".join(f"{k}={v}" for k, v in sorted(by[stage].items())), flush=True)
    live = [(k, a) for k, t in ledger["tasks"].items() for a in t["attempts"] if a.get("status") in LIVE]
    pools = Counter("local" if a.get("local") else POOL.get(a.get("variant", ledger["tasks"][k]["stage"]), "?")
                    for k, a in live)
    print(f"live: {len(live)} {dict(pools)} | statuses {dict(Counter(a['status'] for _, a in live))}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=("plan", "advance", "status", "watch"))
    ap.add_argument("--through", default="import", choices=STAGES, help="last stage to run (default import)")
    ap.add_argument("--interval", type=int, default=1800, help="watch: seconds between cycles")
    args = ap.parse_args()
    rows = load_manifest()
    os.makedirs(os.path.join(WORK, "checkm"), exist_ok=True)
    if args.command == "status":
        return status(load_ledger())
    if args.command == "plan":
        ledger = load_ledger()
        kb = KB(dry_run=True)
        planned = advance(kb, ledger, rows, args.through, dry_run=True)
        print(f"would start {len(planned)} now:", dict(Counter((p[0], p[3], p[4]) for p in planned)))
        return None
    lock = open(LEDGER + ".lock", "w")                          # one writer at a time
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit("another pipeline.py advance/watch holds ledger.json.lock; stop it first")
    kb = KB()
    if args.command == "advance":
        ledger = load_ledger()
        planned = advance(kb, ledger, rows, args.through)
        status(ledger)
        return planned
    while True:                                                  # watch: reload the ledger every cycle
        try:
            ledger = load_ledger()
            advance(kb, ledger, rows, args.through)
            status(ledger)
        except Exception:
            log("cycle failed:\n" + traceback.format_exc())
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
