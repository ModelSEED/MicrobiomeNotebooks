#!/usr/bin/env python
"""How GTDB R220 names map to R226 and R232, genome by genome (desk part of the EmilyKin R232 re-run).

Inputs (downloaded, md5-checked against GTDB's MD5SUM.txt, in gtdb_downloads/rNNN/):
  {ar53,bac120}_taxonomy_rNNN.tsv.gz  every genome in the release with its GTDB lineage
  sp_clusters_rNNN.tsv                species representative -> clustered genomes
Genomes are matched across releases by the numeric part of the accession (GCA/GCF and the version ignored:
GTDB sometimes swaps GCA<->GCF or bumps the version between releases).

Outputs (tables/):
  release_stats.tsv                   genomes, species clusters and distinct taxa per rank, per release
  rep_fate_<a>_<b>.tsv                every species representative of release a: its name in b, whether it is still a
                                      representative, which b cluster holds it, the deepest rank whose name changed
  rep_fate_summary_<a>_<b>.tsv        counts of the above
  taxon_map_<a>_<b>.tsv               every genus/family/order/class/phylum of a: the b name(s) its genomes carry
                                      (dominant name + share; unchanged / renamed / split / merged)
Usage: python -I gtdb_release_diff.py [--downloads DIR] [--out DIR]
"""
import argparse, collections, csv, gzip, os, re, sys

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RELEASES = ["220", "226", "232"]
HERE = os.path.dirname(os.path.abspath(__file__))


def key(acc):
    """RS_GCF_000016525.1 / GB_GCA_000016525.2 / GCA_000016525.1 -> 000016525"""
    m = re.search(r"GC[AF]_(\d+)", acc)
    return m.group(1) if m else acc


def read_taxonomy(d, r):
    tax = {}
    for dom in ("ar53", "bac120"):
        with gzip.open(os.path.join(d, f"r{r}", f"{dom}_taxonomy_r{r}.tsv.gz"), "rt") as fh:
            for line in fh:
                gid, lin = line.rstrip("\n").split("\t")
                tax[key(gid)] = (gid, lin.split(";"))
    return tax


def read_clusters(d, r):
    """genome key -> representative key; and the set of representative keys"""
    member_of, reps = {}, set()
    with open(os.path.join(d, f"r{r}", f"sp_clusters_r{r}.tsv")) as fh:
        next(fh)
        for line in fh:
            row = line.rstrip("\n").split("\t")
            rep = key(row[0])
            reps.add(rep)
            member_of[rep] = rep
            for g in row[9].split(","):
                if g and g != "N/A":
                    member_of[key(g)] = rep
    return member_of, reps


def strip_suffix(name):
    """g__Halomonas_A -> g__Halomonas (GTDB alphabetic suffixes mark polyphyletic splits of one name)"""
    return re.sub(r"_[A-Z]+$", "", name) if not name.startswith("s__") else \
        " ".join(re.sub(r"_[A-Z]+$", "", t) for t in name.split(" "))


def first_change(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return RANKS[i], ("suffix_only" if strip_suffix(x) == strip_suffix(y) else "name")
    return "none", ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--downloads", default=os.path.join(HERE, "gtdb_downloads"))
    ap.add_argument("--out", default=os.path.join(HERE, "tables"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tax = {r: read_taxonomy(a.downloads, r) for r in RELEASES}
    clu = {r: read_clusters(a.downloads, r) for r in RELEASES}

    with open(os.path.join(a.out, "release_stats.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["release", "domain", "genomes", "species_clusters"] + [f"n_{x}" for x in RANKS[1:6]])
        for r in RELEASES:
            for dom in ("d__Bacteria", "d__Archaea"):
                rows = [lin for _, lin in tax[r].values() if lin[0] == dom]
                reps = sum(1 for k in clu[r][1] if k in tax[r] and tax[r][k][1][0] == dom)
                w.writerow([f"R{r}", dom[3:], len(rows), reps] + [len({lin[i] for lin in rows}) for i in range(1, 6)])

    for ra, rb in (("220", "226"), ("226", "232"), ("220", "232")):
        ta, tb = tax[ra], tax[rb]
        mem_b, reps_b = clu[rb]
        reps_a = clu[ra][1]
        cnt = collections.Counter()
        with open(os.path.join(a.out, f"rep_fate_r{ra}_r{rb}.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["rep_key", f"gid_r{ra}", f"lineage_r{ra}", f"gid_r{rb}", f"lineage_r{rb}", f"still_rep_r{rb}",
                        f"cluster_rep_r{rb}", f"cluster_species_r{rb}", "fate", "species_name_changed",
                        "first_changed_rank", "change_kind"])
            for k in sorted(reps_a):
                if k not in ta:
                    continue
                gid_a, lin_a = ta[k]
                if k not in tb:
                    # the representative genome left GTDB (QC / NCBI suppression)
                    w.writerow([k, gid_a, ";".join(lin_a), "", "", "", "", "", "dropped", "", "", ""])
                    cnt["dropped"] += 1
                    continue
                gid_b, lin_b = tb[k]
                rep_b = mem_b.get(k, "")
                still = k in reps_b
                fate = "still_rep" if still else ("absorbed_into_other_cluster" if rep_b else "not_clustered")
                sp_changed = lin_a[6] != lin_b[6]
                rank, kind = first_change(lin_a, lin_b)
                w.writerow([k, gid_a, ";".join(lin_a), gid_b, ";".join(lin_b), still, tb[rep_b][0] if rep_b in tb else "",
                            tb[rep_b][1][6] if rep_b in tb else "", fate, sp_changed, rank, kind])
                cnt[fate] += 1
                cnt[f"species_name_changed={sp_changed}"] += 1
                cnt[f"first_changed_rank={rank}{'/' + kind if kind else ''}"] += 1
                cnt[f"domain={lin_a[0]}|species_name_changed={sp_changed}"] += 1
        with open(os.path.join(a.out, f"rep_fate_summary_r{ra}_r{rb}.tsv"), "w") as fh:
            fh.write("metric\tcount\n")
            for m, c in sorted(cnt.items()):
                fh.write(f"{m}\t{c}\n")

        # taxon-level mapping (genomes present in both releases)
        with open(os.path.join(a.out, f"taxon_map_r{ra}_r{rb}.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["rank", f"taxon_r{ra}", "n_genomes_both", f"dominant_r{rb}", "dominant_share", f"n_r{rb}_names",
                        f"all_r{rb}(name:n)", "status", "sources_of_dominant"])
            for i in range(1, 6):
                fwd = collections.defaultdict(collections.Counter)
                for k, (_, lin_a) in ta.items():
                    if k in tb:
                        fwd[lin_a[i]][tb[k][1][i]] += 1
                inv = collections.defaultdict(set)   # b name -> a names that send >= 10 % of their genomes to it
                for t, c in fwd.items():
                    n = sum(c.values())
                    for nb, x in c.items():
                        if x / n >= 0.10:
                            inv[nb].add(t)
                for t in sorted(fwd):
                    c = fwd[t]
                    n = sum(c.values())
                    dom, dn = c.most_common(1)[0]
                    share = dn / n
                    if share < 0.90:
                        status = "split"
                    elif dom == t:
                        status = "unchanged"
                    elif strip_suffix(dom) == strip_suffix(t):
                        status = "suffix_changed"
                    else:
                        status = "renamed"
                    if status != "split" and len(inv[dom]) > 1:
                        status += "+merged"
                    w.writerow([RANKS[i], t, n, dom, f"{share:.3f}", len(c),
                                ";".join(f"{x}:{y}" for x, y in c.most_common(8)), status,
                                ";".join(sorted(inv[dom])) if len(inv[dom]) > 1 else ""])
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
