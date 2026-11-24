#!/usr/bin/env python3
"""Concatenate the KOfam profiles (one release) into one HMM file for hmmsearch, as DRAM does, and check it.

Reads (read-only) the shared release dir  <kofam_db>/profiles/K*.hmm  and  <kofam_db>/ko_list.
Writes into <out_dir>:
  kofam_all.hmm          every profile, sorted by K number (DRAM searches all profiles; '-' KOs are searched, never called)
  markers.hmm            the marker KOs of annotation_choice/marker_ko_table.tsv (test runs)
  hmm_info.tsv           ko, LENG, threshold, score_type, profile_index (1-based order in kofam_all.hmm), definition
  db_check.txt           profile/ko_list consistency report
stdlib only.  usage: build_db.py <kofam_db> <marker_ko_table.tsv> <out_dir>
"""
import os
import sys

kdb, markers_f, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
ko_list = {}
with open(os.path.join(kdb, "ko_list")) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        ko_list[f[0]] = dict(zip(hdr, f))
markers = []
with open(markers_f) as fh:
    fh.readline()
    for line in fh:
        markers.append(line.split("\t")[0])

files = sorted(f for f in os.listdir(os.path.join(kdb, "profiles")) if f.endswith(".hmm"))
info, names_seen = [], set()
with open(os.path.join(out, "kofam_all.hmm.part"), "w") as fa, open(os.path.join(out, "markers.hmm"), "w") as fm:
    for i, f in enumerate(files, 1):
        txt = open(os.path.join(kdb, "profiles", f)).read()
        name = leng = None
        for line in txt.splitlines()[:30]:
            if line.startswith("NAME"):
                name = line.split()[1]
            elif line.startswith("LENG"):
                leng = int(line.split()[1])
        assert name == f[:-4], (f, name)
        assert txt.count("\n//") == 1 and txt.rstrip().endswith("//"), f   # exactly one profile per file
        names_seen.add(name)
        fa.write(txt if txt.endswith("\n") else txt + "\n")
        if name in markers:
            fm.write(txt if txt.endswith("\n") else txt + "\n")
        k = ko_list.get(name, {})
        info.append((name, leng, k.get("threshold", "NA"), k.get("score_type", "NA"), i, k.get("definition", "")))
os.replace(os.path.join(out, "kofam_all.hmm.part"), os.path.join(out, "kofam_all.hmm"))
with open(os.path.join(out, "hmm_info.tsv"), "w") as fo:
    fo.write("ko\tleng\tthreshold\tscore_type\tprofile_index\tdefinition\n")
    for r in info:
        fo.write("\t".join(map(str, r)) + "\n")
no_prof = sorted(set(ko_list) - names_seen)
no_row = sorted(names_seen - set(ko_list))
thr = [float(v["threshold"]) for v in ko_list.values() if v["threshold"] != "-"]
with open(os.path.join(out, "db_check.txt"), "w") as fo:
    fo.write(f"profiles\t{len(files)}\nko_list_rows\t{len(ko_list)}\n"
             f"profiles_without_ko_list_row\t{len(no_row)}\t{','.join(no_row[:50])}\n"
             f"ko_list_rows_without_profile\t{len(no_prof)}\t{','.join(no_prof[:50])}\n"
             f"ko_list_dash_threshold\t{sum(1 for v in ko_list.values() if v['threshold'] == '-')}\n"
             f"ko_list_min_numeric_threshold\t{min(thr)}\n"
             f"sum_LENG\t{sum(r[1] for r in info)}\n"
             f"markers_in_table\t{len(markers)}\tmarkers_with_profile\t{sum(m in names_seen for m in markers)}\t"
             f"missing\t{','.join(m for m in markers if m not in names_seen)}\n")
print(open(os.path.join(out, "db_check.txt")).read())
