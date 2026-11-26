"""Shared local-resource gate for the two GROW drivers' hybrid workers (user decision 2026-11-23).

One flocked JSON (grow/hybrid/slots.json) lists every local hybrid worker of either driver with the threads and
memory it was given. A driver may start a worker only if, with it, the hybrid workers stay within
  TOTAL_THREADS 48 / TOTAL_MEM_GB 160  (BRIEF limits),
MemAvailable stays >= request + MEM_MARGIN_GB, the machine has request/2 idle cores by load average (the other
users' jobs, e.g. codiffusion, keep theirs), the work disk keeps >= MIN_FREE_GB + the step's disk need free, and no
GTDB-Tk scope is active.
Shares (user decision 2026-11-24, replacing the strict per-sample priority of 2026-11-23): on FM the sitecov driver has
a reserved 12 threads / 40 GB and the per-sample driver the other 36 threads / 120 GB. A driver may always start a step
that keeps it within its own share (if the totals allow; running steps are never preempted, so a share that the other
driver borrowed comes back as its steps finish). It may go beyond its share (borrow the other's idle capacity) only
while the other driver reported no ready-but-waiting steps in its last cycle (< 2 h ago) - so borrowing stops, and the
capacity is yielded back, as soon as the owner has work. TF is shared the same way, half each (10 threads / 25 GB),
using the drivers' TF-eligible waiting counts.
The DAS Tool workers stay under each driver's own LOCAL_CAP and are counted here too (by their declared size).
"""
import fcntl
import json
import os
import shutil
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SLOTS = os.path.join(HERE, "slots.json")
TOTAL_THREADS = 48
TOTAL_MEM_GB = 160
MEM_MARGIN_GB = 30
MIN_FREE_GB = 100
GTDBTK_SCOPES = ("gtdbtk-grow-assembly", "gtdbtk-grow-sitecov", "gtdbtk-salterns-reassembly")
PRIORITY_DRIVER = "per-sample"
WAITING_TTL_S = 2 * 3600
TF_THREADS, TF_MEM_GB = 20, 50
SHARE = {"per-sample": (36, 120), "sitecov": (12, 40)}          # FM threads / GB (user decision 2026-11-24)
TF_SHARE = {"per-sample": (10, 25), "sitecov": (10, 25)}        # TF: half each, borrowable when the other is idle
POPLAR_THREADS = 32                                              # poplar overflow (user decision 2026-11-25)
POPLAR_MIN_FREE_GB = 200


def _alive(pid, workdir):
    try:
        with open(f"/proc/{pid}/stat") as fh:
            if fh.read().rsplit(")", 1)[1].split()[0] in ("Z", "X"):
                return False
        with open(f"/proc/{pid}/cmdline", "rb") as fh:
            return workdir.encode() in fh.read()
    except (OSError, IndexError):
        return False


def mem_available_gb():
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1048576
    return 0.0


def gtdbtk_running():
    for s in GTDBTK_SCOPES:
        st = subprocess.run(["systemctl", "--user", "is-active", s], capture_output=True, text=True).stdout.strip()
        if st in ("active", "activating", "deactivating", "reloading", "refreshing"):
            return s
    return None


class Gate:
    """with Gate() as g: g.can(...); g.add(...)  -- holds the flock for the whole block."""

    def __enter__(self):
        self.fh = open(SLOTS + ".lock", "w")
        fcntl.flock(self.fh, fcntl.LOCK_EX)
        try:
            self.d = json.load(open(SLOTS))
        except (OSError, ValueError):
            self.d = {"workers": {}, "waiting": {}}
        self.d.setdefault("workers", {})
        self.d.setdefault("waiting", {})
        for wd, w in list(self.d["workers"].items()):            # drop finished workers
            if w.get("pid") and not _alive(w["pid"], wd.split("#", 1)[0]):
                del self.d["workers"][wd]
            elif not w.get("pid") and time.time() - w.get("t", 0) > 600:
                del self.d["workers"][wd]
        return self

    def __exit__(self, *exc):
        tmp = SLOTS + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.d, f, indent=1)
        os.replace(tmp, SLOTS)
        fcntl.flock(self.fh, fcntl.LOCK_UN)
        self.fh.close()

    def used(self, host="FM", driver=None):
        ws = [w for w in self.d["workers"].values() if w.get("host", "FM") == host
              and (driver is None or w.get("driver") == driver)]
        return sum(w["threads"] for w in ws), sum(w["mem_gb"] for w in ws)

    def other_waiting(self, driver, key="n"):
        """Ready-but-waiting steps the other driver reported in its last cycle (0 if stale)."""
        n = 0
        for d, w in self.d["waiting"].items():
            if d != driver and time.time() - w.get("t", 0) < WAITING_TTL_S:
                n += w.get(key, 0) or 0
        return n

    def share_ok(self, driver, threads, mem_gb, host="FM"):
        """(ok, reason): within the driver's own share, or borrowing while the other driver has nothing waiting."""
        st, sm = (SHARE if host == "FM" else TF_SHARE).get(driver, (0, 0))
        t, m = self.used(host, driver)
        if t + threads <= st and m + mem_gb <= sm:
            return True, None
        ow = self.other_waiting(driver, "n" if host == "FM" else "n_tf")
        if ow:
            return False, (f"{host} share of {driver} {t}+{threads}/{st} threads, {m}+{mem_gb}/{sm} GB; "
                           f"other driver has {ow} steps waiting")
        return True, None

    def can_tf(self, driver, threads, mem_gb):
        """TF (user plan 2026-11-23): <= 20 threads / 50 GB for our jobs, >= 100 GB free, and back off when the
        user's own work loads it (load average leaves fewer than threads idle cores)."""
        ok, why = self.share_ok(driver, threads, mem_gb, "TF")
        if not ok:
            return False, why
        t, m = self.used("TF")
        if t + threads > TF_THREADS or m + mem_gb > TF_MEM_GB:
            return False, f"TF budget {t}+{threads} threads / {m}+{mem_gb} GB"
        try:
            out = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", "TF",
                                  "cat /proc/loadavg; nproc; grep MemAvailable /proc/meminfo; df -k --output=avail ~ | tail -1"],
                                 capture_output=True, text=True, timeout=40).stdout.split("\n")
            load, ncpu = float(out[0].split()[0]), int(out[1])
            avail_gb, free_gb = int(out[2].split()[1]) / 1048576, int(out[3]) / 1048576
        except Exception as e:  # noqa: BLE001
            return False, f"TF unreachable ({type(e).__name__})"
        if ncpu - load < threads:
            return False, f"TF busy (load {load:.1f} on {ncpu})"
        if avail_gb < mem_gb + 8:
            return False, f"TF MemAvailable {avail_gb:.0f} GB"
        if free_gb < 100:
            return False, f"TF disk {free_gb:.0f} GB free"
        return True, None

    def can(self, driver, threads, mem_gb, disk_gb=0, disk_path=HERE, share=True):
        """(ok, reason). share=False for the small FM footprint of a step that computes on poplar."""
        ok, why = self.share_ok(driver, threads, mem_gb, "FM") if share else (True, None)
        if not ok:
            return False, why
        g = gtdbtk_running()
        if g:
            return False, f"GTDB-Tk scope {g} is active"
        t, m = self.used()
        if t + threads > TOTAL_THREADS:
            return False, f"threads {t}+{threads} > {TOTAL_THREADS}"
        if m + mem_gb > TOTAL_MEM_GB:
            return False, f"memory {m}+{mem_gb} GB > {TOTAL_MEM_GB}"
        avail = mem_available_gb()
        if avail < mem_gb + MEM_MARGIN_GB:
            return False, f"MemAvailable {avail:.0f} GB < {mem_gb}+{MEM_MARGIN_GB}"
        idle = os.cpu_count() - os.getloadavg()[0]
        if share and idle < threads / 2:                   # (a poplar step's FM side mostly copies and polls)
            return False, f"only {idle:.0f} idle cores by load average"
        free = shutil.disk_usage(disk_path).free / 1e9
        if free < MIN_FREE_GB + disk_gb:
            return False, f"{free:.0f} GB free < {MIN_FREE_GB}+{disk_gb}"
        return True, None

    def can_poplar(self, driver, threads, mem_gb, disk_gb):
        """poplar (shared ANL host, MFA ControlMaster): <= 32 threads of ours, /scratch keeps >= 200 GB free after the
        step's need, MemAvailable >= the step's memory + 100 GB, enough idle cores; per-sample keeps priority here."""
        if driver != PRIORITY_DRIVER:
            w = self.d["waiting"].get(PRIORITY_DRIVER)
            if w and w.get("n_poplar") and time.time() - w.get("t", 0) < WAITING_TTL_S:
                return False, "per-sample priority on poplar"
        t, _ = self.used("poplar")
        if t + threads > POPLAR_THREADS:
            return False, f"poplar budget {t}+{threads} > {POPLAR_THREADS} threads"
        try:
            out = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "poplar",
                                  "cat /proc/loadavg; nproc; grep MemAvailable /proc/meminfo; "
                                  "df -k --output=avail /scratch | tail -1"],
                                 capture_output=True, text=True, timeout=60).stdout.split("\n")
            load, ncpu = float(out[0].split()[0]), int(out[1])
            avail_gb, free_gb = int(out[2].split()[1]) / 1048576, int(out[3]) / 1048576
        except Exception as e:  # noqa: BLE001
            return False, f"poplar unreachable ({type(e).__name__}; MFA ControlMaster down?)"
        if ncpu - load < threads / 2:                     # same rule as FM: half the step's threads idle
            return False, f"poplar busy (load {load:.0f} on {ncpu})"
        if avail_gb < mem_gb + 100:
            return False, f"poplar MemAvailable {avail_gb:.0f} GB"
        if free_gb < POPLAR_MIN_FREE_GB + disk_gb:
            return False, f"poplar /scratch {free_gb:.0f} GB free < {POPLAR_MIN_FREE_GB}+{disk_gb:.0f}"
        return True, None

    def add(self, workdir, driver, stage, unit, threads, mem_gb, pid=None, host="FM"):
        self.d["workers"][workdir] = {"driver": driver, "stage": stage, "unit": unit, "threads": threads,
                                      "mem_gb": mem_gb, "pid": pid, "t": time.time(), "host": host}

    def set_pid(self, workdir, pid):
        if workdir in self.d["workers"]:
            self.d["workers"][workdir]["pid"] = pid

    def waiting(self, driver, n, n_tf=0, n_poplar=0):
        self.d["waiting"][driver] = {"n": n, "n_tf": n_tf, "n_poplar": n_poplar, "t": time.time()}


def reads_cache_cleanup(cache, needed_refs, log=print):
    """Delete cached read files whose library ref no driver still needs (needs_*.json in this dir) and that no live
    worker is using (a worker records the refs it reads in its spec)."""
    needed = set(needed_refs)
    for f in os.listdir(HERE):
        if f.startswith("needs_") and f.endswith(".json"):
            try:
                needed |= set(json.load(open(os.path.join(HERE, f))))
            except (OSError, ValueError):
                return                                          # unreadable: delete nothing this time
    with Gate() as g:
        for wd in g.d["workers"]:
            try:
                sp = json.load(open(os.path.join(wd, "spec.json")))
                needed |= set(sp.get("reads_refs", []))
            except (OSError, ValueError):
                pass
    if not os.path.isdir(cache):
        return
    for f in os.listdir(cache):
        if f.endswith(".part"):
            continue
        ref = "/".join(f.split(".", 1)[0].split("_")[:3])
        if ref not in needed:
            for p in (os.path.join(cache, f),):
                try:
                    os.remove(p)
                    log(f"reads cache: removed {f} (no longer needed)")
                except OSError:
                    pass
