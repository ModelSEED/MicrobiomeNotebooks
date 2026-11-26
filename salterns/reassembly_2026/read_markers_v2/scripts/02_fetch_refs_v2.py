#!/usr/bin/env python
"""Step 2 (v2; runs on poplar, which has internet): external protein references for the v2 panel.

  * KOs of the v1 panel (source == v1): the v1 raw files (read_markers/refs_raw/kegg/<KO>.faa, link/<KO>.txt) are
    copied unchanged, so every v1 reference stays byte-identical (needed for the v1 -> v2 regression check).
  * New KOs: KEGG REST link/genes/<KO> + get/<10 genes>/aaseq (as v1). Prokaryote genes; if more than kegg_cap:
      targets -> every archaeal gene + a random sample of kegg_cap bacterial genes (the v1 marker rule, so archaeal
                 homologs of shared methanogenesis/C1 genes are all present),
      decoys  -> a random sample of kegg_cap prokaryote genes (the v1 decoy rule).
    Seed per KO: random.Random(int(KO[1:])) (v1 used one shared Random(1), whose draws depend on the panel order).
  * Eukaryote genes of every counted target KO (euk_cap random genes, same per-KO seed + 1) -> kegg_euk/<KO>.faa;
    labelled '<gene>_euk' (decoy_reported) in 04_build_db_v2.py, so that reads of algal / fungal homologs (the waters
    are Dunaliella-rich) are not counted as prokaryotic copies.
  * UniProtKB: the v1 Pfam/name sets (copied), plus reviewed (Swiss-Prot) Bacteria + Archaea members of the Pfam
    families of the generic targets (labelled by KOfam on poplar; members that pass no panel KO become
    '<family>_other' decoys, as v1's SDIMO_other / GRE_other / Rieske_other).
Usage: python -I 02_fetch_refs_v2.py <ko_panel_v2.tsv> <v1 refs_raw> <outdir> [<KO subset file> <tag>]
  (subset: only these KOs, so two instances on different hosts can share the work, as in v1; env NO_UNIPROT=1 skips
  the UniProt part. Selection is seeded per KO, so it does not depend on which instance fetches a KO.)
"""
import os, sys, time, random, re, csv, shutil, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import requests

panel_p, v1raw, out = sys.argv[1:4]
SUBSET = set(open(sys.argv[4]).read().split()) if len(sys.argv) > 4 else None
TAG = sys.argv[5] if len(sys.argv) > 5 else "all"
for d in ("kegg/link", "kegg_euk", "uniprot"):
    os.makedirs(f"{out}/{d}", exist_ok=True)
S = requests.Session()
LOG = open(f"{out}/fetch_log.tsv", "a")


def get(url, tries=8):
    for i in range(tries):
        try:
            r = S.get(url, timeout=180)
            if r.status_code == 200:
                return r.text
            if r.status_code == 404:
                return ""
        except requests.RequestException:
            pass
        time.sleep(4 * (i + 1))
    raise RuntimeError(url)


br_p = f"{out}/kegg/br08601.txt"
if not os.path.exists(br_p):
    shutil.copy(f"{v1raw}/kegg/br08601.txt", br_p)     # the v1 copy (2026-11-24), so organism taxonomy is identical
tax = {}
A = B = C = D = ""
for line in open(br_p):
    if line.startswith("A"):
        A = re.sub(r"\s*\(\d+\)$", "", line[1:].strip())
    elif line.startswith("B"):
        B = re.sub(r"\s*\(\d+\)$", "", line[1:].strip()); C = D = ""
    elif line.startswith("C"):
        C = re.sub(r"\s*\(\d+\)$", "", line[1:].strip()); D = ""
    elif line.startswith("D"):
        D = re.sub(r"\s*\(\d+\)$", "", line[1:].strip())
    elif line.startswith("E"):
        f = line[1:].strip().split(None, 1)
        if len(f) == 2:
            tax[f[0]] = (A, B, C, D, f[1])

panel = list(csv.DictReader(open(panel_p), delimiter="\t"))
genes_meta, todo, todo_euk = [], {}, {}
for row in panel:
    ko, cap, ecap = row["ko"], int(row["kegg_cap"]), int(row["euk_cap"])
    if SUBSET is not None and ko not in SUBSET:
        continue
    lp = f"{out}/kegg/link/{ko}.txt"
    if not os.path.exists(lp):
        if row["source"] == "v1" and os.path.exists(f"{v1raw}/kegg/link/{ko}.txt"):
            shutil.copy(f"{v1raw}/kegg/link/{ko}.txt", lp)
        else:
            open(lp, "w").write(get(f"https://rest.kegg.jp/link/genes/{ko}")); time.sleep(0.4)
    genes = [l.split("\t")[1].strip() for l in open(lp) if "\t" in l]
    prok = sorted(g for g in genes if tax.get(g.split(":")[0], ("",))[0] == "Prokaryotes")
    euk = sorted(g for g in genes if tax.get(g.split(":")[0], ("",))[0] == "Eukaryotes")
    unk = [g for g in genes if g.split(":")[0] not in tax]
    if row["source"] == "v1":
        fp = f"{out}/kegg/{ko}.faa"
        if not os.path.exists(fp + ".ok"):
            if os.path.exists(f"{v1raw}/kegg/{ko}.faa"):
                shutil.copy(f"{v1raw}/kegg/{ko}.faa", fp); open(fp + ".ok", "w").write("ok (copied from v1)\n")
            else:   # v1 had no KEGG gene for this KO (e.g. K12905: link/genes is empty)
                open(fp, "w").close(); open(fp + ".ok", "w").write("ok (v1 had no KEGG genes)\n")
        sel = None
    else:
        rng = random.Random(int(ko[1:]))
        sel = prok
        if cap and len(sel) > cap:
            if row["role"] == "target":
                arc = [g for g in sel if tax[g.split(":")[0]][1] == "Archaea"]
                bac = [g for g in sel if tax[g.split(":")[0]][1] != "Archaea"]
                sel = sorted(arc + (rng.sample(bac, cap) if len(bac) > cap else bac))
            else:
                sel = sorted(rng.sample(sel, cap))
        todo[ko] = sel
    seleuk = []
    if row["role"] == "target" and ecap > 0 and euk:
        rng2 = random.Random(int(ko[1:]) + 1)
        seleuk = sorted(rng2.sample(euk, ecap)) if len(euk) > ecap else euk
        todo_euk[ko] = seleuk
    LOG.write(f"kegg_link\t{ko}\tall={len(genes)}\tprok={len(prok)}\teuk={len(euk)}\tunknown_org={len(unk)}\t"
              f"selected_prok={'v1_copy' if sel is None else len(sel)}\tselected_euk={len(seleuk)}\n")
    for g in (sel or []) + seleuk:
        t = tax[g.split(":")[0]]
        genes_meta.append((g, ko, t[0], t[1], t[2], t[3], t[4]))
LOG.flush()
print("link done; prok genes to fetch", sum(len(v) for v in todo.values()), "euk", sum(len(v) for v in todo_euk.values()), flush=True)


def fetch_batch(batch):
    txt = get("https://rest.kegg.jp/get/" + "+".join(batch) + "/aaseq")
    time.sleep(0.35)
    return txt


def fetch_set(gl, fp, tag, ko):
    if os.path.exists(fp + ".ok"):
        return
    batches = [gl[i:i + 10] for i in range(0, len(gl), 10)]
    with ThreadPoolExecutor(3) as ex:
        res = list(ex.map(fetch_batch, batches))
    seqs = {}
    for txt in res:
        cur = None
        for l in txt.splitlines():
            if l.startswith(">"):
                cur = l[1:].split()[0]; seqs[cur] = []
            elif cur:
                seqs[cur].append(l.strip())
    with open(fp, "w") as fh:
        for g in gl:
            if g in seqs:
                fh.write(f">{g}\n{''.join(seqs[g])}\n")
    LOG.write(f"{tag}\t{ko}\trequested={len(gl)}\treceived={sum(1 for g in gl if g in seqs)}\n"); LOG.flush()
    open(fp + ".ok", "w").write("ok\n")
    print(tag, ko, len(gl), sum(1 for g in gl if g in seqs), flush=True)


for ko, gl in todo.items():
    fetch_set(gl, f"{out}/kegg/{ko}.faa", "kegg_aaseq", ko)
for ko, gl in todo_euk.items():
    fetch_set(gl, f"{out}/kegg_euk/{ko}.faa", "kegg_aaseq_euk", ko)
with open(f"{out}/kegg/kegg_genes_v2_new.{TAG}.tsv", "w") as fh:
    fh.write("gene\tko\tA\tB\tC\tD\torganism\n")
    for r in genes_meta:
        fh.write("\t".join(r) + "\n")

if os.environ.get("NO_UNIPROT"):
    print("done (no UniProt)", flush=True); sys.exit(0)
# ---- UniProt: v1 sets (copied) + reviewed members of the generic targets' Pfam families (Bacteria + Archaea)
for f in os.listdir(f"{v1raw}/uniprot"):
    if not os.path.exists(f"{out}/uniprot/{f}"):
        shutil.copy(f"{v1raw}/uniprot/{f}", f"{out}/uniprot/{f}")
NEWFAMS = {"PF00763": "folD", "PF02882": "folD", "PF00551": "purU", "PF01268": "fhs", "PF00464": "glyA",
           "PF03313": "sdaA", "PF03315": "sdaA", "PF00389": "fdh", "PF02826": "fdh", "PF00909": "amtB",
           "PF00756": "fghA", "PF08240": "frmA", "PF00107": "frmA", "PF00296": "ssuD", "PF00871": "ackA",
           "PF01515": "pta", "PF00501": "acs"}
for pf in NEWFAMS:
    name = f"rev_{pf}"
    fp = f"{out}/uniprot/{name}.tsv"
    if os.path.exists(fp + ".ok"):
        continue
    q = f"(xref:pfam-{pf}) AND (reviewed:true) AND ((taxonomy_id:2) OR (taxonomy_id:2157))"
    url = ("https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,organism_name,organism_id,lineage,length,sequence&query="
           + urllib.parse.quote(q))
    txt = get(url)
    open(fp, "w").write(txt)
    n = max(0, txt.count("\n") - 1)
    LOG.write(f"uniprot\t{name}\t{q}\trows={n}\n"); LOG.flush()
    open(fp + ".ok", "w").write("ok\n")
    print(name, n, flush=True)
print("done", flush=True)
