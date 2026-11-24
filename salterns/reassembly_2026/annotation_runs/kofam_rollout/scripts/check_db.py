#!/usr/bin/env python3
"""Check the reused KOfam 2026-10-26 HMM DB (kofam/db/, built by build_db.py on 2026-11-11) against the release
files whose md5 is documented (kofam_db/md5sums.txt: ko_list.gz 697af5b3..., profiles.tar.gz bc68ad2a...).

 1. index kofam_all.hmm: profiles in order, byte ranges, LENG; no bytes outside a profile; NAMEs sorted
 2. stream profiles.tar.gz: every K*.hmm member is byte-identical to its range in kofam_all.hmm (build_db adds a
    final newline if a file lacks one); the member set = the profile set (27,981)
 3. markers.hmm = the marker profiles (meta/marker_ko_table.tsv) in K order, byte-identical
 4. hmm_info.tsv = regenerated from kofam_all.hmm (LENG, order) + ko_list.gz (threshold, score_type, definition)
 5. the bakta_test profile chunks (used by sec_full_db.sh) concatenate to kofam_all.hmm byte-for-byte
stdlib only. usage: check_db.py <kofam_db> <db_dir> <marker_ko_table.tsv> <bakta_chunk_dir>
"""
import glob
import gzip
import hashlib
import os
import sys
import tarfile

kdb, db, mk_tbl, bch = sys.argv[1:5]
allp = os.path.join(db, "kofam_all.hmm")
out = []

# 1. index
idx, order, leng = {}, [], {}
md5_all = hashlib.md5()
with open(allp, "rb") as fh:
    pos = start = 0
    name = None
    for line in fh:
        md5_all.update(line)
        if line.startswith(b"NAME"):
            name = line.split()[1].decode()
        elif line.startswith(b"LENG"):
            leng[name] = int(line.split()[1])
        pos += len(line)
        if line.startswith(b"//"):
            assert name and name not in idx, name
            idx[name] = (start, pos)
            order.append(name)
            start, name = pos, None
    total = pos
out.append(f"kofam_all.hmm\tbytes\t{total}\tmd5\t{md5_all.hexdigest()}\tprofiles\t{len(order)}")
out.append(f"kofam_all.hmm\ttrailing_bytes_after_last_profile\t{total - start}")
out.append(f"kofam_all.hmm\tnames_sorted\t{order == sorted(order)}")

# 2. tarball vs ranges
n_mem = n_eq = 0
bad = []
with tarfile.open(os.path.join(kdb, "profiles.tar.gz"), "r:gz") as tf, open(allp, "rb") as fh:
    members = set()
    for m in tf:
        if not (m.isfile() and m.name.endswith(".hmm")):
            continue
        ko = os.path.basename(m.name)[:-4]
        members.add(ko)
        n_mem += 1
        data = tf.extractfile(m).read()
        if not data.endswith(b"\n"):
            data += b"\n"
        if ko not in idx:
            bad.append(ko)
            continue
        s, e = idx[ko]
        fh.seek(s)
        if fh.read(e - s) == data:
            n_eq += 1
        else:
            bad.append(ko)
out.append(f"tarball\thmm_members\t{n_mem}\tbyte_identical_to_kofam_all_range\t{n_eq}\tdiffering\t{len(bad)}\t{','.join(bad[:20])}")
out.append(f"tarball\tmember_set_equals_profile_set\t{members == set(order)}")

# 3. markers.hmm
markers = [l.split("\t")[0] for l in open(mk_tbl).read().splitlines()[1:] if l.strip()]
h = hashlib.md5()
with open(allp, "rb") as fh:
    for ko in order:
        if ko in markers:
            s, e = idx[ko]
            fh.seek(s)
            h.update(fh.read(e - s))
md5_mk = hashlib.md5(open(os.path.join(db, "markers.hmm"), "rb").read()).hexdigest()
out.append(f"markers.hmm\tmd5\t{md5_mk}\texpected_from_kofam_all\t{h.hexdigest()}\tequal\t{md5_mk == h.hexdigest()}"
           f"\tmarkers\t{len(markers)}\twith_profile\t{sum(m in idx for m in markers)}")

# 4. hmm_info.tsv
ko_list = {}
with gzip.open(os.path.join(kdb, "ko_list.gz"), "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        ko_list[f[0]] = dict(zip(hdr, f))
lines = ["ko\tleng\tthreshold\tscore_type\tprofile_index\tdefinition\n"]
for i, ko in enumerate(order, 1):
    k = ko_list.get(ko, {})
    lines.append("\t".join(map(str, (ko, leng[ko], k.get("threshold", "NA"), k.get("score_type", "NA"), i,
                                      k.get("definition", "")))) + "\n")
exp = "".join(lines)
got = open(os.path.join(db, "hmm_info.tsv")).read()
out.append(f"hmm_info.tsv\tregenerated_from_ko_list.gz_equal\t{exp == got}\trows\t{len(lines) - 1}"
           f"\tsum_LENG\t{sum(leng.values())}")

# 5. bakta_test profile chunks
hb = hashlib.md5()
chs = sorted(glob.glob(os.path.join(bch, "chunk_*.hmm")))
for p in chs:
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            hb.update(blk)
out.append(f"bakta_test_chunks\tfiles\t{len(chs)}\tconcat_md5\t{hb.hexdigest()}\tequal_kofam_all\t{hb.hexdigest() == md5_all.hexdigest()}")
print("\n".join(out))
