#!/usr/bin/env bash
# Download (segmented, resumable) and verify every raw file DRAM 1.5.0 prepare_databases needs (--skip_uniref).
# Output: $RAW/<files>, $RAW/MANIFEST.tsv (url, last-modified, bytes, md5), $RAW/FETCH_OK when all checks pass.
set -euo pipefail
RAW=${RAW:-/scratch1/afreiburger/salterns_2026/dram_db/raw}
NSEG=${NSEG:-8}
mkdir -p "$RAW"; cd "$RAW"
log(){ echo "[$(date '+%F %T')] $*"; }

hdr(){ curl -sIL -m 60 "$1" | tr -d '\r'; }
seg_dl(){ # url out [nseg]
  local url=$1 out=$2 n=${3:-$NSEG} size chunk i s e need
  size=$(hdr "$url" | grep -i '^content-length' | tail -1 | tr -dc 0-9)
  [ -n "$size" ] || { log "no size for $url"; return 1; }
  if [ -s "$out" ] && [ "$(stat -c %s "$out")" = "$size" ]; then log "have $out ($size)"; return 0; fi
  [ "$size" -lt 50000000 ] && n=1
  mkdir -p "$out.parts"; chunk=$(( (size + n - 1) / n ))
  for i in $(seq 0 $((n-1))); do
    s=$((i*chunk)); e=$((s+chunk-1)); [ $e -ge $size ] && e=$((size-1)); need=$((e-s+1))
    ( p="$out.parts/$i"
      for try in $(seq 1 30); do
        have=0; [ -f "$p" ] && have=$(stat -c %s "$p")
        [ "$have" -ge "$need" ] && break
        curl -sS --fail -L -m 7200 -r $((s+have))-$e "$url" >> "$p" || { log "retry $try seg $i $out"; sleep 20; }
      done ) &
  done
  wait
  for i in $(seq 0 $((n-1))); do
    s=$((i*chunk)); e=$((s+chunk-1)); [ $e -ge $size ] && e=$((size-1)); need=$((e-s+1))
    [ "$(stat -c %s "$out.parts/$i")" = "$need" ] || { log "segment $i of $out incomplete"; return 1; }
  done
  for i in $(seq 0 $((n-1))); do cat "$out.parts/$i"; done > "$out.tmp" && mv "$out.tmp" "$out" && rm -r "$out.parts"
  [ "$(stat -c %s "$out")" = "$size" ] || { log "size mismatch $out"; return 1; }
  log "downloaded $out ($size)"
}
record(){ # url file
  printf '%s\t%s\t%s\t%s\t%s\n' "$2" "$1" "$(hdr "$1" | grep -i '^last-modified' | tail -1 | cut -d' ' -f2-)" \
    "$(stat -c %s "$2")" "$(md5sum "$2" | cut -d' ' -f1)" >> MANIFEST.tsv.tmp
}

KOFAM=https://www.genome.jp/ftp/db/kofam
PFAM=http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
DBCAN=https://pro.unl.edu/dbCAN2/download/Databases/V11      # bcb.unl.edu (DRAM's URL) now redirects to the dbCAN2 home page
MEROPS=https://ftp.ebi.ac.uk/pub/databases/merops/current_release
VIRAL=https://ftp.ncbi.nlm.nih.gov/refseq/release/viral
VOG=https://fileshare.lisc.univie.ac.at/vog/latest           # DRAM's fileshare.csb.univie.ac.at redirects here

declare -A FILES=(
  [kofam_profiles.tar.gz]=$KOFAM/profiles.tar.gz
  [kofam_ko_list.gz]=$KOFAM/ko_list.gz
  [Pfam-A.full.gz]=$PFAM/Pfam-A.full.gz
  [Pfam-A.hmm.dat.gz]=$PFAM/Pfam-A.hmm.dat.gz
  [Pfam.version.gz]=$PFAM/Pfam.version.gz
  [dbCAN-HMMdb-V11.txt]=$DBCAN/dbCAN-HMMdb-V11.txt
  [CAZyDB.08062022.fam-activities.txt]=$DBCAN/CAZyDB.08062022.fam-activities.txt
  [CAZyDB.08062022.fam.subfam.ec.txt]=$DBCAN/CAZyDB.08062022.fam.subfam.ec.txt
  [pepunit.lib]=$MEROPS/pepunit.lib
  [viral.1.protein.faa.gz]=$VIRAL/viral.1.protein.faa.gz
  [vog.hmm.tar.gz]=$VOG/vog.hmm.tar.gz
  [vog.annotations.tsv.gz]=$VOG/vog.annotations.tsv.gz
  [vog_release.txt]=$VOG/release.txt
)
# small ones first, Pfam last (40 GB)
ORDER="kofam_ko_list.gz kofam_profiles.tar.gz Pfam-A.hmm.dat.gz Pfam.version.gz dbCAN-HMMdb-V11.txt CAZyDB.08062022.fam-activities.txt CAZyDB.08062022.fam.subfam.ec.txt pepunit.lib viral.1.protein.faa.gz vog.hmm.tar.gz vog.annotations.tsv.gz vog_release.txt Pfam-A.full.gz"
rm -f MANIFEST.tsv.tmp
for f in $ORDER; do log "fetch $f"; seg_dl "${FILES[$f]}" "$f"; done
curl -sS --fail -m 120 -o pfam_md5_checksums "$PFAM/md5_checksums"
curl -sS --fail -m 120 -o vog.hmm.tar.gz.md5 "$VOG/vog.hmm.tar.gz.md5"
curl -sS --fail -m 120 -o vog.annotations.tsv.gz.md5 "$VOG/vog.annotations.tsv.gz.md5"

log "verify"
grep -E ' (Pfam-A.full.gz|Pfam-A.hmm.dat.gz)$' pfam_md5_checksums | md5sum -c -
for f in vog.hmm.tar.gz vog.annotations.tsv.gz; do [ "$(md5sum $f | cut -d' ' -f1)" = "$(grep -o '^[0-9a-f]\{32\}' $f.md5)" ] || { log "md5 FAIL $f"; exit 1; }; done
# KOfam ko_list must be the 2026-10-26 release used in annotation_choice/ (sources/kegg/kofam_ko_list_2026-10-26.gz)
[ "$(md5sum kofam_ko_list.gz | cut -d' ' -f1)" = 697af5b3b22a919306a1e60f99cec1a6 ] || { log "ko_list is not the 2026-10-26 release"; exit 1; }
gzip -t kofam_profiles.tar.gz viral.1.protein.faa.gz Pfam-A.hmm.dat.gz
tar -tzf kofam_profiles.tar.gz | grep -c '\.hmm$' | sed 's/^/kofam profiles: /'
tar -tzf vog.hmm.tar.gz | grep -c 'VOG.*\.hmm$' | sed 's/^/VOG hmms: /'
head -1 dbCAN-HMMdb-V11.txt | grep -q '^HMMER3' || { log "dbCAN file is not HMMER3"; exit 1; }
grep -c '^NAME' dbCAN-HMMdb-V11.txt | sed 's/^/dbCAN models: /'
head -c 1 pepunit.lib | grep -q '>' || { log "pepunit.lib not FASTA"; exit 1; }
zcat Pfam.version.gz
for f in $ORDER; do record "${FILES[$f]}" "$f"; done
mv MANIFEST.tsv.tmp MANIFEST.tsv
touch FETCH_OK
log "FETCH_OK"
