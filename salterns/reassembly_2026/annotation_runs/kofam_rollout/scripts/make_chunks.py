#!/usr/bin/env python3
"""Pool the per-MAG Prodigal proteins into N_CHUNKS hmmsearch chunks of WHOLE species clusters, de-duplicated
within each chunk.

Why: 28.4 % of the 3,911,146 proteins are exact duplicates of another MAG's protein (the same strain binned from
a per-sample assembly and a site co-assembly). Identical sequences get identical hmmsearch results (scores do not
depend on the database; E-values do only through -Z, which is fixed to the total protein count), so each chunk
searches one representative per distinct sequence and the postprocessor copies its hits to the other genes
(chunks/chunk_XX.members.tsv: member gene -> representative gene). Duplicates are mostly within a species
cluster (skani 95 % ANI, comparison/mags_old_vs_new.tsv), so whole clusters are kept in one chunk; every chunk
still holds whole MAGs, so per-MAG results are complete as soon as its chunk finishes.

- The species clusters of the control bins (first argument file) are forced into chunk 01, so the first chunk
  can be checked against annotation_choice/kofam_control_test/marker_calls_by_mag.tsv.
- Remaining clusters: longest-processing-time greedy on de-duplicated residues.

kofam_rollout change (2026-11-19): the manifest needs only a 'name' column. Without a 'species_cluster' column
every MAG is its own group (LPT over whole MAGs), which is how the 935 rollout MAGs were chunked (no species
clusters existed when the track was set up). The controls file may be empty (production run: no forced chunk).
Still exact: identical proteins are expanded per chunk, so the grouping only changes how much is saved.
stdlib only.  usage: make_chunks.py <controls.txt> <mag_manifest.tsv> <genes_dir> <out_dir> [n_chunks=40]
"""
import hashlib
import os
import sys
from collections import defaultdict

controls_f, manifest_f, genes_dir, out_dir = sys.argv[1:5]
n_chunks = int(sys.argv[5]) if len(sys.argv) > 5 else 40
controls = [l.strip() for l in open(controls_f) if l.strip()]
cluster = {}
with open(manifest_f) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        cluster[r["name"]] = r.get("species_cluster") or r["name"]     # no cluster column -> one group per MAG
mags = sorted(f[:-5] for f in os.listdir(genes_dir) if f.endswith(".done") and f != "ALL_DONE")
assert set(mags) == set(cluster), (len(mags), len(cluster))


def read_faa(m):
    out, name, seq = [], None, []
    with open(os.path.join(genes_dir, m + ".faa")) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(seq)))
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name:
        out.append((name, "".join(seq)))
    return out


# per-MAG protein counts and residues; per-cluster de-duplicated residues; global uniqueness
mag_stats, clu_hashes, glob = {}, defaultdict(dict), set()
for m in mags:
    prots = read_faa(m)
    mag_stats[m] = (len(prots), sum(len(s) for _, s in prots))
    d = clu_hashes[cluster[m]]
    for _, s in prots:
        h = hashlib.md5(s.encode()).digest()
        d[h] = len(s)
        glob.add(h)
clu_res = {c: sum(d.values()) for c, d in clu_hashes.items()}
del clu_hashes

ctrl_clusters = sorted({cluster[c] for c in controls})
load = [0] * n_chunks
assign_c = {}
for c in ctrl_clusters:
    assign_c[c] = 0
    load[0] += clu_res[c]
for c in sorted((c for c in clu_res if c not in assign_c), key=lambda c: (-clu_res[c], c)):
    i = min(range(n_chunks), key=lambda k: load[k])
    assign_c[c] = i
    load[i] += clu_res[c]

os.makedirs(out_dir, exist_ok=True)
tot_n = sum(v[0] for v in mag_stats.values())
tot_res = sum(v[1] for v in mag_stats.values())
searched_n = searched_res = 0
chunk_lines = []
for i in range(n_chunks):
    cid = f"{i + 1:02d}"
    cm = [m for m in mags if assign_c[cluster[m]] == i]
    seen = {}
    n_all = n_uni = r_uni = 0
    with open(os.path.join(out_dir, f"chunk_{cid}.faa.part"), "w") as fo, \
            open(os.path.join(out_dir, f"chunk_{cid}.members.tsv.part"), "w") as fmem:
        fmem.write("member_gene\trepresentative_gene\n")
        for m in cm:
            for name, s in read_faa(m):
                n_all += 1
                h = hashlib.md5(s.encode()).digest()
                rep = seen.get(h)
                if rep is None:
                    seen[h] = name
                    n_uni += 1
                    r_uni += len(s)
                    fo.write(f">{name}\n")
                    for k in range(0, len(s), 80):
                        fo.write(s[k:k + 80] + "\n")
                else:
                    fmem.write(f"{name}\t{rep}\n")
    for suf in (".faa", ".members.tsv"):
        os.replace(os.path.join(out_dir, f"chunk_{cid}{suf}.part"), os.path.join(out_dir, f"chunk_{cid}{suf}"))
    searched_n += n_uni
    searched_res += r_uni
    chunk_lines.append(f"chunk_{cid}\tmags\t{len(cm)}\tclusters\t{sum(1 for c in assign_c if assign_c[c] == i)}\t"
                       f"proteins\t{n_all}\tsearched\t{n_uni}\tresidues_searched\t{r_uni}\n")

with open(os.path.join(out_dir, "chunk_assignment.tsv"), "w") as fo:
    fo.write("mag\tchunk\tn_proteins\tresidues\tspecies_cluster\tforced_control_cluster\n")
    for m in mags:
        c = cluster[m]
        fo.write(f"{m}\t{assign_c[c] + 1:02d}\t{mag_stats[m][0]}\t{mag_stats[m][1]}\t{c}\t{int(c in ctrl_clusters)}\n")
with open(os.path.join(out_dir, "chunk_stats.txt"), "w") as fo:
    fo.write(f"mags\t{len(mags)}\nproteins\t{tot_n}\nresidues\t{tot_res}\n"
             f"unique_protein_sequences_global\t{len(glob)}\tfraction_duplicate\t{1 - len(glob) / tot_n:.4f}\n"
             f"searched_proteins_within_chunk_dedup\t{searched_n}\tfraction_saved\t{1 - searched_n / tot_n:.4f}\n"
             f"searched_residues\t{searched_res}\tfraction_saved\t{1 - searched_res / tot_res:.4f}\n"
             f"chunks\t{n_chunks}\tdedup_cluster_residues_min\t{min(load)}\tmax\t{max(load)}\n"
             f"control_clusters_in_chunk_01\t{len(ctrl_clusters)}\n")
    fo.writelines(chunk_lines)
print(open(os.path.join(out_dir, "chunk_stats.txt")).read())
