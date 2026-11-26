#!/usr/bin/env bash
# Build the smoke-test inputs under binning_test/smoke/ (small, real data): ~30 Mbp of the longest R1_A_D1 Arm C
# scaffolds, 4 M R1_A_D1 + 2 M SF2_C_D1 cleaned read pairs (pilot qc/clean), 6 OLD bins and 4 originals of R1_A_D1.
# SMOKE_B never gets an assembly (tests the SF2 timeout path) and has empty OLD/original sets (tests the EMPTY path).
set -euo pipefail
B=/scratch1/afreiburger/salterns_2026/binning_test
P=/scratch1/afreiburger/salterns_2026/tadpole_pilot
PY=/scratch1/afreiburger/salterns_2026/envs/bt_map/bin/python
S=$B/smoke; mkdir -p $S/input $S/inputs/old/SMOKE_A $S/inputs/orig/SMOKE_A $S/inputs/old/SMOKE_B $S/inputs/orig/SMOKE_B $S/inputs/kbase_checkm
[ -s $S/input/smoke_asm.fa ] || $PY -I - $P/R1_A_D1/armC/spades/scaffolds.fasta $S/input/smoke_asm.fa 30000000 <<'PYEOF'
import sys
inp, out, target = sys.argv[1], sys.argv[2], int(sys.argv[3])
seqs, name, buf = [], None, []
for line in open(inp):
    if line.startswith(">"):
        if name: seqs.append((name, "".join(buf)))
        name, buf = line[1:].strip(), []
    else:
        buf.append(line.strip())
if name: seqs.append((name, "".join(buf)))
seqs.sort(key=lambda x: -len(x[1]))
tot = 0
with open(out, "w") as fo:
    for n, s in seqs:
        if tot >= target: break
        fo.write(f">{n}\n" + "\n".join(s[i:i+80] for i in range(0, len(s), 80)) + "\n"); tot += len(s)
print("smoke asm bp", tot)
PYEOF
for x in "SMK1 R1_A_D1 16000000" "SMK2 SF2_C_D1 8000000"; do
  set -- $x
  mkdir -p $S/input/pilotlike/$1/qc $S/input/pilotlike/$1/logs
  # seqkit head stops after N records (zcat | head under pipefail dies of SIGPIPE)
  for m in 1 2; do $B/../envs/bt_map/bin/seqkit head -n $(( $3 / 4 )) -j 4 -o $S/input/pilotlike/$1/qc/clean_$m.fq.gz $P/$2/qc/clean_$m.fq.gz; done
  cp $P/$2/logs/bbduk_adapt.log $P/$2/logs/bbduk_filter.log $S/input/pilotlike/$1/logs/
done
for b in $(ls $B/inputs/old/R1_A_D1 | head -6); do cp $B/inputs/old/R1_A_D1/$b $S/inputs/old/SMOKE_A/; done
head -1 $B/inputs/orig/R1_A_D1/originals.tsv > $S/inputs/orig/SMOKE_A/originals.tsv
grep -E "metabat\.(45|47|51)\.contigs__|concoct_out\.3\.contigs\.fa" $B/inputs/orig/R1_A_D1/originals.tsv >> $S/inputs/orig/SMOKE_A/originals.tsv
tail -n +2 $S/inputs/orig/SMOKE_A/originals.tsv | cut -f3 | while read f; do cp $B/inputs/orig/R1_A_D1/$f $S/inputs/orig/SMOKE_A/; done
head -1 $B/inputs/orig/R1_A_D1/originals.tsv > $S/inputs/orig/SMOKE_B/originals.tsv
sed 's/__R1_A_D1/__SMOKE_A/' $B/inputs/kbase_checkm/checkm_R1_A_D1.tsv > $S/inputs/kbase_checkm/checkm_SMOKE_A.tsv
touch $S/input/asm.done
ls -la $S/input $S/input/pilotlike/* $S/inputs/*/*
echo SMOKE_PREP_OK
