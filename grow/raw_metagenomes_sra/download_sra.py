"""Download the raw metagenomes behind the 3,825 GROWdb MAGs (KBase ws 145226) from ENA as fastq.gz.

runs.tsv lists 171 SRA runs = every metagenome the MAGs were binned from (133 samples):
  - 65 NASQAN + 2 Yukon samples, BioProject PRJNA946291 (GROWdb / Colorado State)
  - 66 WHONDRS-2019 surface-water samples, JGI per-sample BioProjects
    (the 11 WHONDRS runs CSU re-deposited in PRJNA946291 are exact duplicates and are excluded)

usage: python download_sra.py runs.tsv OUTDIR [workers]
Files land in OUTDIR/<sample>/<run>_{1,2}.fastq.gz, each checked against ENA's size + md5
(a <file>.md5 is written on success). Re-runnable and resumable (HTTP Range); finished files are skipped.
"""
import csv, hashlib, os, shutil, sys, time, concurrent.futures as cf
import requests

runs_tsv, OUT = sys.argv[1], sys.argv[2]
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 6
MARGIN = 20e9                      # never fill the disk past this much free space
CHUNK = 8 << 20

def log(*a):
    print(time.strftime('%F %T'), *a, flush=True)

tasks = []
for r in csv.DictReader(open(runs_tsv), delimiter='\t'):
    for url, size, md5 in zip(r['fastq_ftp'].split(';'), r['fastq_bytes'].split(';'), r['fastq_md5'].split(';')):
        tasks.append(dict(sample=r['sample'], run=r['Run'], url='https://' + url, size=int(size), md5=md5,
                          dst=os.path.join(OUT, r['sample'], url.rsplit('/', 1)[1])))
TOTAL = sum(t['size'] for t in tasks)

def fetch(t):
    dst, part, size = t['dst'], t['dst'] + '.part', t['size']
    if os.path.exists(dst) and os.path.getsize(dst) == size and os.path.exists(dst + '.md5'):
        return t, 'skip'
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    for attempt in range(1, 8):
        try:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            if have > size:
                os.remove(part); have = 0
            while shutil.disk_usage(OUT).free < (size - have) + MARGIN:
                log(f'WAIT disk nearly full, holding {os.path.basename(dst)}'); time.sleep(600)
            h = hashlib.md5()
            if have:
                with open(part, 'rb') as fh:
                    for b in iter(lambda: fh.read(CHUNK), b''):
                        h.update(b)
            t0, n0 = time.time(), have
            hdr = {'Range': f'bytes={have}-'} if have else {}
            with requests.get(t['url'], headers=hdr, stream=True, timeout=(60, 300)) as r:
                if have and r.status_code != 206:          # server ignored the range: start over
                    have, h = 0, hashlib.md5()
                r.raise_for_status()
                with open(part, 'ab' if have else 'wb') as fh:
                    for b in r.iter_content(CHUNK):
                        fh.write(b); h.update(b); have += len(b)
            got = h.hexdigest()
            if have != size or got != t['md5']:
                os.remove(part)
                raise IOError(f'mismatch: {have}/{size} bytes, md5 {got} vs {t["md5"]}')
            os.replace(part, dst)
            with open(dst + '.md5', 'w') as fh:
                fh.write(f'{got}  {os.path.basename(dst)}\n')
            log(f'OK   {t["sample"]}/{os.path.basename(dst)} {size/1e9:.2f} GB {(have-n0)/1e6/max(1e-6,time.time()-t0):.0f} MB/s')
            return t, 'ok'
        except Exception as e:
            log(f'RETRY {attempt} {os.path.basename(dst)}: {str(e)[:200]}'); time.sleep(min(600, 30 * attempt))
    return t, 'FAILED'

log(f'start: {len(tasks)} files, {len({t["run"] for t in tasks})} runs, {len({t["sample"] for t in tasks})} samples, '
    f'{TOTAL/1e12:.3f} TB -> {OUT} ({WORKERS} workers, {shutil.disk_usage(OUT).free/1e12:.2f} TB free)')
done_bytes, status = 0, {}
with cf.ThreadPoolExecutor(WORKERS) as ex:
    for i, (t, st) in enumerate(ex.map(fetch, tasks), 1):
        status[t['dst']] = st
        if st != 'FAILED':
            done_bytes += t['size']
        if i % 10 == 0 or i == len(tasks):
            log(f'PROGRESS {i}/{len(tasks)} files, {done_bytes/1e12:.3f}/{TOTAL/1e12:.3f} TB')
with open(os.path.join(OUT, 'manifest.tsv'), 'w') as fh:
    fh.write('sample\trun\tfile\tbytes\tmd5\tstatus\n')
    for t in tasks:
        fh.write(f"{t['sample']}\t{t['run']}\t{os.path.relpath(t['dst'], OUT)}\t{t['size']}\t{t['md5']}\t{status[t['dst']]}\n")
bad = [os.path.basename(k) for k, v in status.items() if v == 'FAILED']
log(f'done: {len(tasks)-len(bad)}/{len(tasks)} files ok; failed: {bad}')
