#!/usr/bin/env python3
"""GROW hybrid: detached job runner on poplar (overflow host; user decisions 2026-11-23 / 2026-11-25).

    ssh poplar 'python3 /scratch/afreiburger/grow_hybrid/bin/poplar_runner.py <rundir>'  < token-on-stdin

<rundir>/job.json: {"downloads": [{"node": shock id, "md5": ..., "path": abs path}], "cmd": bash command,
                    "outputs": [abs paths to md5 when done], "min_free_gb": floor while running,
                    "free_path": the file system that floor is for (default <rundir>),
                    "moves": [[src, dst], ...] after success (before the md5s), "rm_after": [dirs removed after them]}
(since 2026-11-25, 3rd decision: the run dir with reads and outputs is on /scratch1; MEGAHIT's work dir is on local
/scratch when it has room - on success the outputs move to /scratch1 and the work dir is removed at once.)

- The KBase token is read from stdin only when there is something to download. It lives only in this process's
  memory (never in a file, never in the child's environment). The process detaches (double fork) right after reading
  stdin, so the ssh call returns at once.
- STATE holds one word: running | done | failed. state.json has the details (rc, md5s, reason, times).
- Standard library only (poplar's /usr/bin/python3); metag-hybrid is not touched.
"""
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import urllib.request

BLOB = "https://kbase.us/services/shock-api"


def write_state(rd, word, **kw):
    st = {"state": word, "t": time.strftime("%Y-%m-%dT%H:%M:%S%z"), **kw}
    with open(os.path.join(rd, "state.json.tmp"), "w") as f:
        json.dump(st, f, indent=1)
    os.replace(os.path.join(rd, "state.json.tmp"), os.path.join(rd, "state.json"))
    with open(os.path.join(rd, "STATE.tmp"), "w") as f:
        f.write(word + "\n")
    os.replace(os.path.join(rd, "STATE.tmp"), os.path.join(rd, "STATE"))


def log(rd, msg):
    with open(os.path.join(rd, "runner.log"), "a") as f:
        f.write(time.strftime("[%Y-%m-%d %H:%M:%S] ") + msg + "\n")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def download(rd, token, d, tries=5):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(f"{BLOB}/node/{d['node']}?download_raw",
                                         headers={"Authorization": "OAuth " + token})
            h = hashlib.md5()
            with urllib.request.urlopen(req, timeout=600) as r, open(d["path"] + ".part", "wb") as out:
                for b in iter(lambda: r.read(8 << 20), b""):
                    out.write(b)
                    h.update(b)
            if d.get("md5") and h.hexdigest() != d["md5"]:
                last = f"md5 {h.hexdigest()} != {d['md5']}"
                continue
            os.replace(d["path"] + ".part", d["path"])
            log(rd, f"downloaded node {d['node']} -> {d['path']} (md5 ok)")
            return
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
            time.sleep(min(600, 15 * 2 ** i))
    raise RuntimeError(f"download of node {d['node']} failed: {last}")


def free_gb(path):
    return shutil.disk_usage(path).free / 1e9


def main():
    rd = os.path.abspath(sys.argv[1])
    job = json.load(open(os.path.join(rd, "job.json")))
    token = sys.stdin.readline().strip() if job.get("downloads") else None
    if job.get("downloads") and not token:
        write_state(rd, "failed", reason="no token on stdin")
        return
    # detach: the ssh session returns now; the work continues without a terminal
    if os.fork() > 0:
        return
    os.setsid()
    if os.fork() > 0:
        os._exit(0)
    devnull = os.open(os.devnull, os.O_RDWR)
    for fd in (0, 1, 2):
        os.dup2(devnull, fd)
    write_state(rd, "running", pid=os.getpid(), host=os.uname().nodename)
    try:
        for d in job.get("downloads", []):
            if not (os.path.exists(d["path"]) and d.get("md5") and md5(d["path"]) == d["md5"]):
                download(rd, token, d)
        token = None                                               # not needed past the downloads
        env = {k: v for k, v in os.environ.items() if "TOKEN" not in k.upper()}
        log(rd, "RUN: " + job["cmd"])
        t0 = time.time()
        p = subprocess.Popen(["/bin/bash", "-o", "pipefail", "-c", job["cmd"]], cwd=rd, env=env,
                             stdout=open(os.path.join(rd, "cmd.log"), "ab"), stderr=subprocess.STDOUT,
                             start_new_session=True)
        floor = job.get("min_free_gb", 150)
        fpath = job.get("free_path", rd)
        reason = None
        while p.poll() is None:
            time.sleep(30)
            if free_gb(fpath) < floor:                             # keep /scratch usable for everyone else
                reason = f"{fpath} free {free_gb(fpath):.0f} GB < {floor} GB: stopped"
                os.killpg(p.pid, signal.SIGTERM)
                time.sleep(30)
                if p.poll() is None:
                    os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                break
        rc = p.returncode
        log(rd, f"exit {rc} after {time.time() - t0:.0f} s" + (f" ({reason})" if reason else ""))
        if rc != 0 or reason:
            write_state(rd, "failed", rc=rc, reason=reason or f"command exit {rc}")
            return
        for src, dst in job.get("moves", []):
            if os.path.exists(src):
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.move(src, dst)
        for d in job.get("rm_after", []):
            shutil.rmtree(d, ignore_errors=True)
            log(rd, f"removed work dir {d}")
        sums = {o: md5(o) for o in job.get("outputs", []) if os.path.exists(o)}
        missing = [o for o in job.get("outputs", []) if o not in sums]
        if missing:
            write_state(rd, "failed", rc=rc, reason=f"missing outputs {missing}")
            return
        write_state(rd, "done", rc=rc, md5=sums, seconds=round(time.time() - t0))
    except Exception as e:  # noqa: BLE001
        write_state(rd, "failed", reason=f"{type(e).__name__}: {e}"[:1500])


if __name__ == "__main__":
    main()
