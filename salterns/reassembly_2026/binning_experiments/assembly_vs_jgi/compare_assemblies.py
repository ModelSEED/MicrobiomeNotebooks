#!/usr/bin/env python
"""Per-sample contiguity of the three assemblies of each salterns metagenome:
  jgi2  = (ii)  JGI IMG-team 'SPAdes version 2' reassembly, 2018 (Ga0209xxx; IMG 3300026094-251): exact lengths of
          all contigs >= 500 bp from the CoverM BAM headers (lengths_jgi_v2.tsv.gz; 33 samples)
  orig  = (iii) the non-public Tringe-lab Tadpole + SPAdes per-sample assembly the original 450 MAGs were binned from:
          contigs >= 2,500 bp reconstructed from the MAG contig NODE ranks (mag_source_reconstruction.tsv; 39 samples)
  new   = 2026 KBase metaSPAdes 3.15.3 scaffolds >= 1,000 bp (lengths_2026.tsv.gz; 39 samples)
The columns of provenance/assembly_vs_bins.tsv that need contigs < 1,000 bp (n_contigs_ge500, asm_Mbp_ge500,
median_len, pct_bp_ge1000, pct_bp_ge2500, N50_ge500) cannot be computed for 'new' (kb_SPAdes min_contig_length 1000)
and are left empty there; for 'orig' only contigs >= 2,500 bp are known.
Also joins: binned Mbp per sample (original 450 / 310 MAGs; 2026 selected MAGs), read and run statistics
(spades_run_stats.tsv), JGI read recruitment (assembly_vs_bins.tsv).
Writes assembly_vs_jgi.tsv (one row per sample) and assembly_vs_jgi_summary.tsv (by sample type / site), and
coassembly_2026.tsv (the 4 MEGAHIT sediment co-assemblies, for context).
Usage: ~/Documents/py_venv/bin/python -I compare_assemblies.py
"""
import collections, csv, gzip, os, re, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SALTERNS = os.path.normpath(os.path.join(ROOT, ".."))


def load_lengths(path):
    d = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            u, _, l = line.rstrip("\n").split("\t")
            d[u].append(int(l))
    for u in d:
        d[u].sort(reverse=True)
    return d


def nx(lens, f=0.5):
    tot, acc = sum(lens), 0
    for l in lens:
        acc += l
        if acc >= f * tot:
            return l
    return 0


def stats(lens, prefix, low):
    """lens sorted descending; low = smallest length present (500 for jgi2, 1000 for new)."""
    o = {}
    ge = lambda t: [l for l in lens if l >= t]  # noqa: E731
    if low <= 500:
        g5 = ge(500)
        o[f"{prefix}_n_contigs_ge500"] = len(g5)
        o[f"{prefix}_asm_Mbp_ge500"] = round(sum(g5) / 1e6, 1)
        o[f"{prefix}_median_len"] = int(st.median(g5))
        o[f"{prefix}_N50_ge500"] = nx(g5)
    for t in (1000, 2500, 5000, 10000, 20000):
        g = ge(t)
        o[f"{prefix}_n_ge{t}"] = len(g)
        o[f"{prefix}_Mbp_ge{t}"] = round(sum(g) / 1e6, 1)
        if t in (1000, 2500):
            o[f"{prefix}_N50_ge{t}"] = nx(g)
    if low <= 500:
        tot = sum(ge(500))
        o[f"{prefix}_pct_bp_ge1000"] = round(100 * sum(ge(1000)) / tot, 1)
        o[f"{prefix}_pct_bp_ge2500"] = round(100 * sum(ge(2500)) / tot, 1)
    o[f"{prefix}_pct_of_ge1000_in_ge2500"] = round(100 * sum(ge(2500)) / sum(ge(1000)), 1)
    o[f"{prefix}_max_len"] = lens[0]
    return o


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main():
    new = load_lengths(os.path.join(HERE, "lengths_2026.tsv.gz"))
    jgi = load_lengths(os.path.join(HERE, "lengths_jgi_v2.tsv.gz"))
    rec = {r["sample"]: r for r in tsv(os.path.join(HERE, "mag_source_reconstruction.tsv"))}
    runs = {r["sample"]: r for r in tsv(os.path.join(HERE, "spades_run_stats.tsv"))}
    avb = {r["sample"][:-3]: r for r in tsv(os.path.join(ROOT, "provenance", "assembly_vs_bins.tsv"))}
    man = {r["sample"]: r for r in tsv(os.path.join(ROOT, "manifest.tsv"))}
    portal = {r["sample"][:-3] if r["sample"].endswith("_MG") else r["sample"]: r
              for r in tsv(os.path.join(ROOT, "provenance", "jgi_portal.tsv"))}

    # binned Mbp per sample: original MAGs (450 and the 310 dereplicated) and 2026 selected MAGs
    orig_bin = collections.defaultdict(lambda: [0, 0, 0, 0])   # bp450, n450, bp310, n310
    for r in tsv(os.path.join(SALTERNS, "MAG_contigs", "index.tsv")):
        m = re.match(r"Salt_Pond_MetaG_?(.+?)_MG_DASTool", r["mag"])
        s = m.group(1).replace("R2_restored_DShore", "R2_restored_D_shore")
        orig_bin[s][0] += int(r["dna_size"])
        orig_bin[s][1] += 1
        if r["population"] == "derep310":
            orig_bin[s][2] += int(r["dna_size"])
            orig_bin[s][3] += 1
    new_bin = collections.defaultdict(lambda: [0, 0])
    for u in new:
        p = os.path.join(ROOT, "work", "checkm", u, f"selected_{u}.tsv")
        if os.path.exists(p):
            for r in tsv(p):
                if r["pass_50_10"] == "1":
                    new_bin[u][0] += int(r["size_bp"])
                    new_bin[u][1] += 1

    rows = []
    for s in sorted(man):
        row = {"sample": s, "type": man[s]["type"], "site": "R2_restored" if s.startswith("R2_restored") else man[s]["site"],
               "raw_Gbp": round(int(man[s]["bases"]) / 1e9, 2)}
        row["jgi2_Ga"] = portal.get(s, {}).get("sp2_Ga", "") or ""
        row["jgi2_in_BAM"] = s in jgi
        row["jgi_pct_reads_to_jgi2_asm"] = avb.get(s, {}).get("pct_reads_to_assembly", "")
        if s in jgi:
            row.update(stats(jgi[s], "jgi2", 500))
        row.update(stats(new[s], "new", 1000))
        r = rec[s]
        row.update({"orig_n_ge2500_est": int(r["n_ge_est"]), "orig_Mbp_ge2500_lower": round(int(r["bp_lower"]) / 1e6, 1),
                    "orig_Mbp_ge2500_est": round(float(r["bp_est"]) / 1e6, 1), "orig_N50_ge2500_est": int(r["N50_est"]),
                    "orig_Mbp_ge5000_est": round(float(r["bp_ge5000_est"]) / 1e6, 1),
                    "orig_Mbp_ge10000_est": round(float(r["bp_ge10000_est"]) / 1e6, 1),
                    "orig_Mbp_ge20000_est": round(float(r.get("bp_ge20000_est") or 0) / 1e6, 1),
                    "orig_rank_points": int(r["n_points"])})
        # ratios on the binner-relevant basis (contigs >= 2,500 bp)
        row["new_over_orig_Mbp_ge2500"] = round(row["new_Mbp_ge2500"] / row["orig_Mbp_ge2500_est"], 3)
        row["new_over_orig_N50_ge2500"] = round(row["new_N50_ge2500"] / row["orig_N50_ge2500_est"], 3)
        row["new_over_orig_Mbp_ge10000"] = round(row["new_Mbp_ge10000"] / row["orig_Mbp_ge10000_est"], 3)
        if s in jgi:
            row["new_over_jgi2_Mbp_ge2500"] = round(row["new_Mbp_ge2500"] / row["jgi2_Mbp_ge2500"], 3)
            row["new_over_jgi2_N50_ge2500"] = round(row["new_N50_ge2500"] / row["jgi2_N50_ge2500"], 3)
            row["new_over_jgi2_Mbp_ge1000"] = round(row["new_Mbp_ge1000"] / row["jgi2_Mbp_ge1000"], 3)
            row["orig_over_jgi2_Mbp_ge2500"] = round(row["orig_Mbp_ge2500_est"] / row["jgi2_Mbp_ge2500"], 3)
        ob = orig_bin[s]
        row.update({"orig_MAGs_450": ob[1], "orig_binned_Mbp_450": round(ob[0] / 1e6, 1), "orig_MAGs_310": ob[3],
                    "orig_binned_Mbp_310": round(ob[2] / 1e6, 1), "new_MAGs": new_bin[s][1],
                    "new_binned_Mbp": round(new_bin[s][0] / 1e6, 1)})
        row["orig_pct_ge2500_binned_450"] = round(100 * ob[0] / 1e6 / row["orig_Mbp_ge2500_est"], 1)
        row["new_pct_ge2500_binned"] = round(100 * new_bin[s][0] / 1e6 / row["new_Mbp_ge2500"], 1)  # MaxBin2 bins
        # also hold contigs 1-2.5 kb, so this can exceed the binned share of >= 2.5 kb contigs
        ru = runs[s]
        row.update({"trim_pct_bases_kept": float(ru["pct_bases_kept_by_trim"]),
                    "trim_Gbp": round(int(ru["trim_bases"]) / 1e9, 2),
                    "spades_insert_size": float(ru["insert_size"]),
                    "spades_pct_pairs_on_long_edges_K127": float(ru["pct_pairs_on_long_edges_last"])})
        rows.append(row)

    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(HERE, "assembly_vs_jgi.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)

    # ---- summaries
    def summ(group, rs):
        o = {"group": group, "n_samples": len(rs)}
        for k in ("new_Mbp_ge2500", "orig_Mbp_ge2500_est", "new_N50_ge2500", "orig_N50_ge2500_est",
                  "new_over_orig_Mbp_ge2500", "new_over_orig_N50_ge2500", "new_over_orig_Mbp_ge10000"):
            o["median_" + k] = round(st.median(r[k] for r in rs), 3)
        o["sum_new_Mbp_ge2500"] = round(sum(r["new_Mbp_ge2500"] for r in rs), 1)
        o["sum_orig_Mbp_ge2500_est"] = round(sum(r["orig_Mbp_ge2500_est"] for r in rs), 1)
        o["sum_new_Mbp_ge10000"] = round(sum(r["new_Mbp_ge10000"] for r in rs), 1)
        o["sum_orig_Mbp_ge10000_est"] = round(sum(r["orig_Mbp_ge10000_est"] for r in rs), 1)
        j = [r for r in rs if r["jgi2_in_BAM"]]
        o["n_with_jgi2"] = len(j)
        if j:
            for k in ("new_over_jgi2_Mbp_ge2500", "new_over_jgi2_N50_ge2500", "new_over_jgi2_Mbp_ge1000",
                      "orig_over_jgi2_Mbp_ge2500"):
                o["median_" + k] = round(st.median(r[k] for r in j), 3)
            o["sum_jgi2_Mbp_ge2500"] = round(sum(r["jgi2_Mbp_ge2500"] for r in j), 1)
            o["sum_new_Mbp_ge2500_same_samples"] = round(sum(r["new_Mbp_ge2500"] for r in j), 1)
            o["sum_orig_Mbp_ge2500_same_samples"] = round(sum(r["orig_Mbp_ge2500_est"] for r in j), 1)
        o["n_new_below_orig_Mbp_ge2500"] = sum(r["new_over_orig_Mbp_ge2500"] < 1 for r in rs)
        o["n_new_below_jgi2_Mbp_ge2500"] = sum(r.get("new_over_jgi2_Mbp_ge2500", 9) < 1 for r in j)
        return o

    groups = [("all", rows), ("Sediment", [r for r in rows if r["type"] == "Sediment"]),
              ("Water", [r for r in rows if r["type"] == "Water"])]
    for site in ("R1", "R2", "R2A", "SF2", "R2_restored"):
        groups.append((f"site {site}", [r for r in rows if r["site"] == site]))
    summ_rows = [summ(g, rs) for g, rs in groups]
    keys = list(dict.fromkeys(k for r in summ_rows for k in r))
    with open(os.path.join(HERE, "assembly_vs_jgi_summary.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(summ_rows)

    co = []
    for u in sorted(x for x in new if x.endswith("_sediment")):
        o = {"unit": u, **stats(new[u], "new", 1000)}
        cores = [r for r in rows if r["site"] == u.split("_")[0] and r["type"] == "Sediment"]
        o["sum_6_core_new_Mbp_ge2500"] = round(sum(r["new_Mbp_ge2500"] for r in cores), 1)
        o["sum_6_core_orig_Mbp_ge2500_est"] = round(sum(r["orig_Mbp_ge2500_est"] for r in cores), 1)
        co.append(o)
    with open(os.path.join(HERE, "coassembly_2026.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(co[0]), delimiter="\t")
        w.writeheader()
        w.writerows(co)
    print("wrote", len(rows), "samples")


if __name__ == "__main__":
    main()
