#!/usr/bin/env python
"""Step 11b (local, internet): KEGG organism code -> genome assembly accession (KEGG REST get/gn:<code>, field
DATA_SOURCE 'Assembly:GCF_/GCA_...'), so KEGG references can be given their GTDB R232 genus
(gtdbtk_data/release232/taxonomy/gtdb_taxonomy.tsv, matched on the accession number without version, GCF<->GCA).
Usage: python -I 11b_kegg_org_assembly.py <org_codes.txt> <out.tsv>"""
import sys, time, re, os
from concurrent.futures import ThreadPoolExecutor
import requests
codes = [l.strip() for l in open(sys.argv[1]) if l.strip()]
out = sys.argv[2]
done = set()
if os.path.exists(out):
    done = {l.split("\t")[0] for l in open(out)}
codes = [c for c in codes if c not in done]
S = requests.Session()


def get(url, tries=8):
    for i in range(tries):
        try:
            r = S.get(url, timeout=120)
            if r.status_code == 200:
                return r.text
            if r.status_code == 404:
                return ""
        except requests.RequestException:
            pass
        time.sleep(3 * (i + 1))
    return None


def fetch(batch):
    txt = get("https://rest.kegg.jp/get/" + "+".join("gn:" + c for c in batch))
    time.sleep(0.3)
    res = {}
    if not txt:
        return batch, res
    for ent in txt.split("///"):
        m = re.search(r"^ENTRY\s+(T\d+)", ent, re.M)
        if not m:
            continue
        oc = re.search(r"^ORG_CODE\s+(\S+)", ent, re.M)
        code = oc.group(1) if oc else ""
        acc = re.search(r"Assembly:\s*(GC[AF]_\d+\.\d+)", ent)
        tx = re.search(r"^TAXONOMY\s+TAX:(\d+)", ent, re.M)
        org = re.search(r"^NAME\s+(.+)$", ent, re.M)
        res[code] = (m.group(1), acc.group(1) if acc else "", tx.group(1) if tx else "", org.group(1).strip() if org else "")
    return batch, res


batches = [codes[i:i + 10] for i in range(0, len(codes), 10)]
with open(out, "a") as fo, ThreadPoolExecutor(3) as ex:
    for batch, res in ex.map(fetch, batches):
        for c in batch:
            t = res.get(c)
            fo.write(f"{c}\t" + ("\t".join(t) if t else "\t\t\t") + "\n")
        fo.flush()
print("done")
