#!/usr/bin/env python
"""Step 2 (local, needs internet): external protein references for the KO panel.

Sources (fetched 2026-11-24):
  A. KEGG GENES members of every panel KO (meta/ko_panel.tsv), via KEGG REST link/genes/<KO> + get/<10 genes>/aaseq.
     Organism taxonomy from KEGG br08601 (A = Eukaryotes/Prokaryotes, B = Bacteria/Archaea/..., C, D, organism name).
     Prokaryote genes are kept for every KO. KOs with kegg_cap > 0 are randomly subsampled to that many prokaryote genes
     (random.Random(1)). For the marker KOs, eukaryote genes (cap 1000 per KO) are fetched too and later become a
     'euk_marker' decoy class (cytosolic/organellar ribosomal proteins of eukaryotes must not count as prokaryote GE).
     NOTE: UniProt no longer exposes KO cross-references (query xref:ko-K00399 returns 0 in release 2026_04) and maps
     only ~55 % of KEGG genes, so KEGG itself is the KO-defined source.
  B. UniProtKB (release 2026_04) members of Pfam families around selected targets, Bacteria + Archaea only, for breadth
     and taxonomy (labelled later by KOfam on poplar, never by the family name alone).
Writes <out>/kegg/<KO>.faa, <out>/kegg/kegg_genes.tsv, <out>/uniprot/<fam>.tsv (accession, organism, taxid, lineage,
length, sequence) and <out>/fetch_log.tsv.
Usage: python -I 02_fetch_refs.py <ko_panel.tsv> <outdir>
"""
import os, sys, time, random, re, csv, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import requests

panel_p, out = sys.argv[1], sys.argv[2]
os.makedirs(f"{out}/kegg/link", exist_ok=True)
os.makedirs(f"{out}/uniprot", exist_ok=True)
S = requests.Session()
LOG = open(f"{out}/fetch_log.tsv", "a")


def get(url, tries=6):
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
    raise RuntimeError(url)


# ---- KEGG organism taxonomy (br08601)
br_p = f"{out}/kegg/br08601.txt"
if not os.path.exists(br_p):
    open(br_p, "w").write(get("https://rest.kegg.jp/get/br:br08601"))
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
rng = random.Random(1)
genes_meta = []
for row in panel:
    ko, cap = row["ko"], int(row["kegg_cap"])
    lp = f"{out}/kegg/link/{ko}.txt"
    if not os.path.exists(lp):
        open(lp, "w").write(get(f"https://rest.kegg.jp/link/genes/{ko}")); time.sleep(0.4)
    genes = [l.split("\t")[1].strip() for l in open(lp) if "\t" in l]
    prok = [g for g in genes if tax.get(g.split(":")[0], ("",))[0] == "Prokaryotes"]
    euk = [g for g in genes if tax.get(g.split(":")[0], ("",))[0] == "Eukaryotes"]
    unk = [g for g in genes if g.split(":")[0] not in tax]
    sel = sorted(prok)
    if cap and len(sel) > cap:
        if row["role"] == "marker":   # markers: every archaeal gene + a random sample of bacterial genes
            arc = [g for g in sel if tax[g.split(":")[0]][1] == "Archaea"]
            bac = [g for g in sel if tax[g.split(":")[0]][1] != "Archaea"]
            sel = sorted(arc + (rng.sample(bac, cap) if len(bac) > cap else bac))
        else:
            sel = sorted(rng.sample(sel, cap))
    seleuk = []
    if row["role"] == "marker":
        seleuk = sorted(euk)
        if len(seleuk) > 1000:
            seleuk = sorted(rng.sample(seleuk, 1000))
    LOG.write(f"kegg_link\t{ko}\tall={len(genes)}\tprok={len(prok)}\teuk={len(euk)}\tunknown_org={len(unk)}\tselected_prok={len(sel)}\tselected_euk={len(seleuk)}\n")
    for g in sel + seleuk:
        t = tax[g.split(":")[0]]
        genes_meta.append((g, ko, t[0], t[1], t[2], t[3], t[4]))
LOG.flush()

# ---- KEGG aaseq, 10 genes per request, 2 workers
todo = {}
for g, ko, *_ in genes_meta:
    todo.setdefault(ko, []).append(g)


def fetch_batch(batch):
    txt = get("https://rest.kegg.jp/get/" + "+".join(batch) + "/aaseq")
    time.sleep(0.35)
    return txt


for ko, gl in todo.items():
    fp = f"{out}/kegg/{ko}.faa"
    if os.path.exists(fp) and os.path.exists(fp + ".ok"):
        continue
    batches = [gl[i:i + 10] for i in range(0, len(gl), 10)]
    with ThreadPoolExecutor(2) as ex:
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
    LOG.write(f"kegg_aaseq\t{ko}\trequested={len(gl)}\treceived={sum(1 for g in gl if g in seqs)}\n"); LOG.flush()
    open(fp + ".ok", "w").write("ok\n")
    print(ko, len(gl), sum(1 for g in gl if g in seqs), flush=True)

with open(f"{out}/kegg/kegg_genes.tsv", "w") as fh:
    fh.write("gene\tko\tA\tB\tC\tD\torganism\n")
    for r in genes_meta:
        fh.write("\t".join(r) + "\n")

# ---- UniProt Pfam families (Bacteria + Archaea), with lineage
FAMS = {"PF02249": "mcrA", "PF02745": "mcrA", "PF02241": "mcrB", "PF02783": "mcrB", "PF02240": "mcrG",
        "PF06007": "phnJ", "PF02461": "CuMMO_AMO", "PF02332": "SDIMO", "PF07355": "GRDB", "PF02901": "GRE_PFL",
        "PF00848": "Rieske_alpha", "PF12176": "mtaB", "PF06253": "MttB", "PF09505": "mtbB", "PF05369": "mtmB"}
QUERIES = {f"pfam_{pf}": f"(xref:pfam-{pf}) AND ((taxonomy_id:2) OR (taxonomy_id:2157))" for pf in FAMS}
QUERIES["name_mpnS"] = '(protein_name:"methylphosphonate synthase") OR (gene:mpnS)'
QUERIES["name_hepD"] = '(protein_name:"hydroxyethylphosphonate dioxygenase") OR (gene:hepD)'
for name, q in QUERIES.items():
    fp = f"{out}/uniprot/{name}.tsv"
    if os.path.exists(fp + ".ok"):
        continue
    url = ("https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,organism_name,organism_id,lineage,length,sequence&query="
           + urllib.parse.quote(q))
    txt = get(url)
    open(fp, "w").write(txt)
    n = max(0, txt.count("\n") - 1)
    LOG.write(f"uniprot\t{name}\t{q}\trows={n}\n"); LOG.flush()
    open(fp + ".ok", "w").write("ok\n")
    print(name, n, flush=True)
print("done")
