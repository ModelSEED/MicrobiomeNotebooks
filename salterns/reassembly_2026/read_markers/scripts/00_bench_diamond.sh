DM=/scratch1/afreiburger/salterns_2026/envs/bt_dastool/bin/diamond
export TMPDIR=/scratch1/afreiburger/salterns_2026/read_markers/tmp
mode=$1; tag=$2; shift 2
/usr/bin/time -v -o bench_$tag.time nice -n 10 $DM blastx -d prov/prov -q q8M.fq.gz -o bench_$tag.tsv.gz --compress 1 $mode -k 5 -e 1e-3 --threads 16 "$@" -t $TMPDIR --quiet --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
echo done > bench_$tag.done
