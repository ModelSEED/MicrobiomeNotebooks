import sys, re, hashlib
from collections import defaultdict
from pathlib import Path
import pandas as pd
SP = Path(sys.argv[1]); MARK = sys.argv[2]; KOL = sys.argv[3]
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
markers = pd.read_csv(MARK, sep="\t", index_col=0)
kol = pd.read_csv(KOL, sep="\t", index_col=0)
thr = pd.to_numeric(kol.threshold, errors="coerce"); st = kol.score_type
def rf(p):
    s, n = {}, None
    for line in open(p):
        if line.startswith(">"): n = line[1:].split()[0]; s[n] = []
        else: s[n].append(line.strip())
    return {k: "".join(v) for k, v in s.items()}
rows = []
for line in open(SP / "markers.domtblout"):
    if line.startswith("#"): continue
    f = line.split()
    rows.append((f[0], int(f[2]), f[3], int(f[5]), float(f[6]), float(f[7]), float(f[13]), int(f[15]), int(f[16])))
h = pd.DataFrame(rows, columns=["gene", "glen", "ko", "qlen", "full_e", "full", "dom", "hf", "ht"])
h["score"] = h.full.where(h.ko.map(st) == "full", h.dom); h["thr"] = h.ko.map(thr)
h["sig"] = h.ko.map(st).isin(["full", "domain"]) & (h.score > h.thr)
g = h.groupby(["gene", "ko"]).agg(full_e=("full_e", "min"), score=("score", "max"), sig=("sig", "any"), thr=("thr", "first"),
                                  hf=("hf", "min"), ht=("ht", "max"), qlen=("qlen", "first"), glen=("glen", "first")).reset_index()
g["cov"] = ((g.ht - g.hf + 1) / g.qlen).round(2)
sig = g[g.sig].sort_values(["gene", "full_e"]).drop_duplicates("gene").set_index("gene")
anyb = g.sort_values(["gene", "full_e"]).drop_duplicates("gene").set_index("gene")
g["method"] = g.gene.str.split("|").str[0]; g["label"] = g.gene.str.split("|").str[1]
labs = sorted({p.name.split(".")[0] for p in (SP / "bakta_tmp").glob("*.cds.faa")})
# contig map contig_i -> original id by sequence identity
cmap = {}
for l in labs:
    orig = {hashlib.md5(v.upper().encode()).hexdigest(): k for k, v in rf(SP / f"../fasta/{l}.fna").items()}
    for k, v in rf(SP / f"bakta_tmp/{l}.sequences.fna").items():
        cmap[(l, k)] = orig[hashlib.md5(v.upper().encode()).hexdigest()]
calls = sig.reset_index()
calls["method"] = calls.gene.str.split("|").str[0]; calls["label"] = calls.gene.str.split("|").str[1]
calls = calls[calls.label.isin(labs)]
calls["mk"] = calls.ko.map(markers.gene)
t = calls.groupby(["mk", "label", "method"]).size().unstack("method", fill_value=0)
t["PB"] = t.get("prodigal", 0).astype(str) + "/" + t.get("bakta", 0).astype(str)
w = t.PB.unstack("label").fillna("0/0")
diff = w[[any(x.split("/")[0] != x.split("/")[1] for x in r) for _, r in w.iterrows()]]
print("marker calls Prodigal/Bakta-CDS (68 marker profiles; Bakta = proteins sent to PSC search, i.e. without the UPS exact-match CDS)")
print(w.to_string()); print("\nrows that differ:"); print(diff.to_string())
# grd / Sec detail
print("\n=== grd genes (K21579/K21578/K10672/K10670) best hit per gene, both methods ===")
x = g[g.ko.isin(["K21579", "K21578", "K10672", "K10670"]) & g.label.isin(labs) & (g.score > 80)].copy()
x["call"] = x.gene.map(sig.ko).fillna("-").map(lambda k: markers.gene.get(k, k))
x["gene_id"] = x.gene.map(lambda s: s.split("|", 2)[2])
print(x.sort_values(["label", "method", "gene", "ko"])[["label", "method", "gene_id", "glen", "ko", "score", "thr", "cov", "call"]].to_string(index=False))
seqs = rf(SP / "all_markers_input.faa")
for k, v in seqs.items():
    if "U" in v:
        print("\nU-containing Bakta protein:", k, len(v), "aa; U at", v.find("U") + 1)
        print(g[g.gene == k][["ko", "score", "thr", "cov", "sig"]].to_string(index=False))
