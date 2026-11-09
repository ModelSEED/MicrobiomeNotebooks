"""B9. Non-canonical homoserine kinase check with TIGRFAM 15.0 HMMs (extracted to verify_bio/tigr_extract, data only):
TIGR00191 thrB (GHMP), TIGR00938 thrB_alt (type II), TIGR02535 'proposed homoserine kinase' (members adjacent to Thr
genes in genomes lacking any known HSK), TIGR01488 HAD-IB (ThrH/SerB-like). Trusted cutoffs (TC) from the HMM files.
For every TIGR02535 hit: gene neighbourhood (+-4 CDS) to test adjacency to threonine-pathway genes."""
import sys, json
from pathlib import Path
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_bio import *
import pyhmmer

TD = VB / "tigr_extract"
alpha = pyhmmer.easel.Alphabet.amino()
hmms = []
for f in ["TIGR00191.HMM", "TIGR00938.HMM", "TIGR02535.HMM", "TIGR01488.HMM"]:
    with pyhmmer.plan7.HMMFile(str(TD / f)) as hf:
        hmms.append(hf.read())
seqs, feats_by_mag = [], {}
for mid in model_ids():
    g = rast(mid)
    feats_by_mag[iid(mid)] = g["features"]
    for f in g["features"]:
        p = f.get("protein_translation")
        if p:
            seqs.append(pyhmmer.easel.TextSequence(name=f"{iid(mid)}|{short(f['id'])}".encode(), sequence=p.replace("*", "")).digitize(alpha))
hits = {}
for top in pyhmmer.hmmsearch(hmms, seqs, cpus=8, E=1e-5, Z=len(seqs), bit_cutoffs=None):
    q = top.query.name.decode() if isinstance(top.query.name, bytes) else top.query.name
    tc = top.query.cutoffs.trusted[0] if top.query.cutoffs.trusted_available() else None
    for h in top:
        n = h.name.decode() if isinstance(h.name, bytes) else h.name
        mag, gene = n.split("|")
        hits.setdefault(q, []).append((mag, gene, round(h.score, 1), (h.score >= tc) if tc else None, tc))
mags = sorted(feats_by_mag, key=lambda x: (x not in OBLIGATE, x))
print("MAG".ljust(22) + "  thrB(TIGR00191)  thrB_alt(TIGR00938)  hypHSK(TIGR02535)  HAD-IB(TIGR01488)   [n hits >= TC; best score]")
for m in mags:
    cells = []
    for q in ["TIGR00191", "TIGR00938", "TIGR02535", "TIGR01488"]:
        hs = [x for x in hits.get(q, []) if x[0] == m]
        npass = sum(1 for x in hs if x[3])
        best = max([x[2] for x in hs], default=0)
        cells.append(f"{npass} [{best:6.1f}]".rjust(19))
    print(("*" if m in OBLIGATE else " ") + m.ljust(21) + " ".join(cells))
print("\nTIGR02535 hits >= TC with gene neighbourhood (+-4):")
for mag, gene, sc, ok, tc in sorted(hits.get("TIGR02535", []), key=lambda x: (x[0] not in OBLIGATE, x[0])):
    if not ok:
        continue
    feats = feats_by_mag[mag]
    idx = {short(f["id"]): i for i, f in enumerate(feats)}
    i = idx[gene]
    print(f"  {mag} {gene} score {sc} (TC {tc}) RAST: {'; '.join(feats[i].get('functions') or [])}")
    for j in range(max(0, i - 4), min(len(feats), i + 5)):
        f = feats[j]
        loc = f["location"][0]
        print(f"      {'>>' if j == i else '  '} {short(f['id']):10s} {loc[0][:22]:22s} {loc[1]:>8} {loc[2]} {'; '.join(f.get('functions') or [])[:90]}")
json.dump(hits, open(VB / "b9_tigr_hsk.json", "w"), indent=0)
