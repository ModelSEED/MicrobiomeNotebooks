#!/usr/bin/env python3
"""Timing subsets: 300 random KOfam profiles (seed 1) and the first 8 of them. usage: make_test_hmms.py <kofam_db> <db_dir>"""
import os, random, sys
kdb, db = sys.argv[1:3]
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(db, "hmm_info.tsv"))][1:]
random.seed(1)
pick = sorted(random.sample(rows, 300), key=lambda r: r[0])
for name, sub in (("random300", pick), ("random8", pick[:8])):
    with open(os.path.join(db, f"{name}.hmm"), "w") as fo:
        for r in sub:
            fo.write(open(os.path.join(kdb, "profiles", r[0] + ".hmm")).read())
    tot = sum(int(r[1]) for r in rows)
    print(name, len(sub), "sum_LENG", sum(int(r[1]) for r in sub), "of", tot)
