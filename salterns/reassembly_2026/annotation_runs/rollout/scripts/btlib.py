#!/usr/bin/env python3
"""Small stdlib helpers for the salterns binning test (driver.py task scripts call these with python -I).

Subcommands:
  filter-fasta IN MINLEN OUT          keep sequences >= MINLEN (length with N's); print JSON stats of input and output
  bbduk-counts DIR                    JSON of Input/Result reads and bases from DIR/bbduk_{adapt,filter,qtrim}.log
  merge-depth ASM OUT T1 [T2 ...]     merge per-BAM jgi depth tables exactly like KBase's MetaBAT2 app
                                      (combineDepthTablesUtil.py: zero-fill, totalAvgDepth = sum, assembly order)
  abund IN OUT                        per-library MaxBin2 abundance file (contig <tab> mean depth) from a 1-BAM jgi table
  c2b-metabat2 BINDIR OUT             contig2bin from MetaBAT2 bin.N.fa (unbinned/lowDepth/tooShort excluded)
  c2b-maxbin2 OUTDIR OUT              contig2bin from MaxBin2 <prefix>.NNN.fasta (noclass/tooshort excluded)
  c2b-concoct CSV OUT                 contig2bin from CONCOCT clustering_merged.csv
  das-inputs TSV...                   print "<comma tsvs> <comma labels>" of the non-empty contig2bin files
  fasta-stats DIR [EXT]               TSV of genome, size_bp, n_contigs, n50 for DIR/*.EXT
"""
import json
import os
import re
import sys


def read_fasta_lengths(path):
    """[(id, length)] in file order; id = first word of the header."""
    out, cur, n = [], None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out.append((cur, n))
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if cur is not None:
        out.append((cur, n))
    return out


def n50(lengths):
    ls = sorted(lengths, reverse=True)
    tot, acc = sum(ls), 0
    for x in ls:
        acc += x
        if 2 * acc >= tot:
            return x
    return 0


def stats(lengths):
    return {"n": len(lengths), "bp": sum(lengths), "n50": n50(lengths), "max": max(lengths) if lengths else 0}


def filter_fasta(inp, minlen, out):
    minlen = int(minlen)
    all_l, kept_l = [], []
    keep = False
    with open(inp) as fh, open(out, "w") as fo:
        cur_len, cur_lines, cur_hdr = 0, [], None

        def flush():
            if cur_hdr is None:
                return
            all_l.append(cur_len)
            if cur_len >= minlen:
                kept_l.append(cur_len)
                fo.write(cur_hdr)
                fo.writelines(cur_lines)

        for line in fh:
            if line.startswith(">"):
                flush()
                cur_hdr, cur_lines, cur_len = line, [], 0
            else:
                cur_lines.append(line)
                cur_len += len(line.strip())
        flush()
    print(json.dumps({"input": inp, "min_len": minlen, "all": stats(all_l), "kept": stats(kept_l)}, indent=1))


def bbduk_counts(d):
    res = {}
    for step in ("adapt", "filter", "qtrim"):
        p = os.path.join(d, f"bbduk_{step}.log")
        if not os.path.exists(p):
            continue
        r = {}
        for line in open(p, errors="replace"):
            m = re.match(r"^(Input|Result):\s+(\d+) reads.*?\s(\d+) bases", line)
            if m:
                k = m.group(1).lower()
                r[f"{k}_reads"], r[f"{k}_bases"] = int(m.group(2)), int(m.group(3))
        res[step] = r
    print(json.dumps(res, indent=1))


def read_depth_table(path):
    """KBase createDictFromDepthfile: {contig: [len, total, depth, var], 'header': [...]}"""
    d = {}
    with open(path) as f:
        header = f.readline().rstrip().split("\t")
        for line in f:
            vals = line.rstrip().split("\t")
            d[vals[0].split()[0]] = vals[1:]
    d["header"] = header
    return d


def merge_depth(asm, out, tables):
    data = [read_depth_table(t) for t in tables]
    contig_list = [c for c, _ in read_fasta_lengths(asm)]
    header = ["contigName", "contigLen", "totalAvgDepth"]
    main = {}
    for dd in data:
        for c in dd:
            if c == "header":
                continue
            main[c] = [dd[c][0], 0]
    for dd in data:
        header += [dd["header"][3], dd["header"][4]]
        for c in main:
            if c in dd:
                main[c][1] = float(main[c][1]) + float(dd[c][2])
                main[c] += [dd[c][2], dd[c][3]]
            else:
                main[c] += [float(0), float(0)]
    with open(out, "w") as f:
        f.write("\t".join(header) + "\n")
        for c in contig_list:
            if c in main:
                f.write("{}\t{}\n".format(c, "\t".join(str(i) for i in main[c])))
    print(f"merged {len(tables)} tables, {sum(1 for c in contig_list if c in main)} contigs -> {out}")


def abund(inp, out):
    n = 0
    with open(inp) as f, open(out, "w") as fo:
        f.readline()
        for line in f:
            v = line.rstrip().split("\t")
            fo.write(f"{v[0]}\t{v[3]}\n")
            n += 1
    print(f"{n} contigs -> {out}")


def write_c2b(pairs, out):
    seen = {}
    with open(out, "w") as fo:
        for contig, b in pairs:
            if contig in seen:
                raise SystemExit(f"contig {contig} in two bins ({seen[contig]}, {b})")
            seen[contig] = b
            fo.write(f"{contig}\t{b}\n")
    print(f"{len(seen)} contigs in {len(set(seen.values()))} bins -> {out}")


def c2b_metabat2(bindir, out):
    pairs = []
    for fn in sorted(os.listdir(bindir)):
        m = re.match(r"^bin\.(\d+)\.fa$", fn)
        if not m:
            continue
        for c, _ in read_fasta_lengths(os.path.join(bindir, fn)):
            pairs.append((c, f"metabat2.bin.{int(m.group(1)):03d}"))
    write_c2b(pairs, out)


def c2b_maxbin2(outdir, out):
    pairs = []
    for fn in sorted(os.listdir(outdir)):
        m = re.match(r"^.+\.(\d+)\.fasta$", fn)
        if not m:
            continue
        for c, _ in read_fasta_lengths(os.path.join(outdir, fn)):
            pairs.append((c, f"maxbin2.Bin.{m.group(1)}"))
    write_c2b(pairs, out)


def c2b_concoct(csv, out):
    pairs = []
    with open(csv) as f:
        f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            c, k = line.split(",")
            pairs.append((c, f"concoct.bin.{int(k):03d}"))
    write_c2b(pairs, out)


def das_inputs(tsvs):
    keep = [t for t in tsvs if os.path.exists(t) and os.path.getsize(t) > 0]
    labels = [os.path.basename(os.path.dirname(os.path.abspath(t))) for t in keep]
    print(",".join(keep), ",".join(labels))


def fasta_stats(d, ext="fa"):
    print("genome\tsize_bp\tn_contigs\tn50")
    for fn in sorted(os.listdir(d)):
        if fn.endswith("." + ext):
            ls = [x for _, x in read_fasta_lengths(os.path.join(d, fn))]
            print(f"{fn[:-len(ext) - 1]}\t{sum(ls)}\t{len(ls)}\t{n50(ls)}")


def main():
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "filter-fasta":
        filter_fasta(*args)
    elif cmd == "bbduk-counts":
        bbduk_counts(*args)
    elif cmd == "merge-depth":
        merge_depth(args[0], args[1], args[2:])
    elif cmd == "abund":
        abund(*args)
    elif cmd == "c2b-metabat2":
        c2b_metabat2(*args)
    elif cmd == "c2b-maxbin2":
        c2b_maxbin2(*args)
    elif cmd == "c2b-concoct":
        c2b_concoct(*args)
    elif cmd == "das-inputs":
        das_inputs(args)
    elif cmd == "fasta-stats":
        fasta_stats(*args)
    else:
        raise SystemExit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
