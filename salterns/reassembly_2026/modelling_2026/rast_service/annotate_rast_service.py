#!/usr/bin/env python
"""RAST functions for the 57 modelling genomes through the public RAST annotation service (Andrew, 2026-11-20:
"RastClient, cancel KBase"). ModelSEEDpy RastClient -> https://tutorial.theseed.org/services/genome_annotation,
GenomeAnnotation.run_pipeline with stages annotate_proteins_kmer_v2, annotate_proteins_kmer_v1,
annotate_proteins_similarity (RAST's function-assignment stages) on the Prodigal 2.6.3 (-p meta -g 11) proteins that
the KOfam tracks use (sources.tsv), so RAST and KOfam annotate the same genes.
Writes functions/<mag_id>.tsv (protein_id, function; '' = no RAST call) and summary.tsv. Resumable.
"""
import glob, os, time
import pandas as pd
from modelseedpy.core.rast_client import RastClient

HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "functions"), exist_ok=True)
src = pd.read_csv(os.path.join(HERE, "sources.tsv"), sep="\t")
rast = RastClient()
rows = []
for r in src.itertuples():
    faa = os.path.join(HERE, "faa_src", r.poplar_faa.lstrip("/"))
    out = os.path.join(HERE, "functions", f"{r.mag_id}.tsv")
    seqs, name, cur = {}, None, []
    for l in open(faa):
        if l.startswith(">"):
            if name: seqs[name] = "".join(cur)
            name, cur = l[1:].split()[0], []
        else:
            cur.append(l.strip().rstrip("*"))
    seqs[name] = "".join(cur)
    if not os.path.exists(out):
        for k in range(4):
            try:
                t = time.time(); res = rast.annotate_protein_sequences(seqs); dt = time.time() - t; break
            except Exception as e:
                print(r.mag_id, "retry", k, type(e).__name__, str(e)[:120]); time.sleep(20 * (k + 1))
        else:
            raise SystemExit(f"{r.mag_id}: RAST service failed 4 times")
        pd.DataFrame({"protein_id": list(seqs), "function": [res.get(p, "") for p in seqs]}).to_csv(out, sep="\t", index=False)
    else:
        dt = float("nan")
    f = pd.read_csv(out, sep="\t", keep_default_na=False)
    called = f.function != ""
    hyp = f.function.str.lower().str.contains("hypothetical")
    rows.append({"mag_id": r.mag_id, "proteins": len(f), "with_rast_function": int(called.sum()),
                 "non_hypothetical": int((called & ~hyp).sum()), "seconds": round(dt, 1)})
    print(rows[-1])
s = pd.DataFrame(rows)
s["frac_non_hypothetical"] = (s.non_hypothetical / s.proteins).round(3)
s.to_csv(os.path.join(HERE, "summary.tsv"), sep="\t", index=False)
print(s.describe().round(3).to_string())
