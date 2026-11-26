"""Review check on a running chunk's partial tblout/domtblout (complete KOs only).
1. tblout best-1-domain score == max domtblout domain score (for pairs with a domtblout row)
2. domain score - full score excess distribution; lowest thresholds of domain-type KOs
3. significance under DRAM's domtblout rule == significance under postprocess.py's tblout rule
usage: check_partial.py <W> <chunk>"""
import sys, os
W, c = sys.argv[1:3]
info = {}
with open(os.path.join(W, "db/hmm_info.tsv")) as fh:
    fh.readline()
    for l in fh:
        f = l.rstrip("\n").split("\t")
        info[f[0]] = (f[2], f[3])
def rows(p, n):
    with open(p) as fh:
        lines = fh.readlines()
    out = []
    for l in lines[:-1]:            # drop the last (possibly partial) line
        if l.startswith("#"):
            continue
        out.append(l.split(None, n))
    return out
tb = rows(os.path.join(W, f"hmm/chunk_{c}.tblout"), 18)
dt = rows(os.path.join(W, f"hmm/chunk_{c}.domtblout"), 22)
kos_t = list(dict.fromkeys(r[2] for r in tb)); kos_d = list(dict.fromkeys(r[3] for r in dt))
done = set(kos_t[:-1]) & set(kos_d[:-1])       # KOs complete in both files
tbd = {(r[0], r[2]): (float(r[5]), float(r[8])) for r in tb if r[2] in done}
dmax, drule = {}, {}
for r in dt:
    ko = r[3]
    if ko not in done:
        continue
    k = (r[0], ko); ds = float(r[13]); fs = float(r[7])
    dmax[k] = max(dmax.get(k, -1e9), ds)
    thr, st = info[ko]
    if st == "domain":
        s = ds
    elif st == "full":
        s = fs
    else:
        continue
    drule[k] = drule.get(k, False) or (s > float(thr))
print("KOs complete", len(done), "tblout pairs", len(tbd), "pairs with dom rows", len(dmax))
mism = [(k, tbd[k][1], v) for k, v in dmax.items() if abs(tbd[k][1] - v) > 0.051]
print("best1dom != max domtbl domain score:", len(mism), mism[:5])
exc = sorted(tbd[k][1] - tbd[k][0] for k in tbd)
print("domain-full excess: max %.1f; n>0: %d; n>3 bits: %d" % (exc[-1], sum(e > 0 for e in exc), sum(e > 3 for e in exc)))
# rule equivalence: postprocess (tblout) vs DRAM (domtblout rows)
pp = {}
for k, (fs, ds) in tbd.items():
    thr, st = info[k[1]]
    if st in ("full", "domain"):
        pp[k] = (fs if st == "full" else ds) > float(thr)
diff = [k for k in pp if pp[k] != drule.get(k, False)]
print("significant pairs: postprocess", sum(pp.values()), "DRAM-domtblout", sum(drule.values()), "differ", len(diff), diff[:5])
dthr = sorted(float(t) for t, s in info.values() if s == "domain")
print("lowest domain-type thresholds", dthr[:8], "n domain-type < 25:", sum(t < 25 for t in dthr))
ftr = sorted(float(t) for t, s in info.values() if s == "full")
print("lowest full-type thresholds", ftr[:5])
