#!/bin/bash
set -euo pipefail
if [ $# -lt 3 ]; then
  echo "usage: subset_bams.sh <bam_list.txt> <contig:start-end> <out_dir> [threads] [reference.fa for CRAM]" >&2
  exit 1
fi
LIST=$1
REGION=$2
OUT=$3
THREADS=${4:-4}
REF=${5:-}
SAMTOOLS=${SAMTOOLS:-samtools}
mkdir -p "$OUT"
one() {
  in=$1
  name=$(basename "$in"); name=${name%.*}
  out="$OUT/$name.bam"
  [ -s "$out.bai" ] && return 0
  if [ -n "$REF" ]; then
    "$SAMTOOLS" view -b -T "$REF" -o "$out" "$in" "$REGION"
  else
    "$SAMTOOLS" view -b -o "$out" "$in" "$REGION"
  fi
  "$SAMTOOLS" index "$out"
}
export -f one
export OUT REGION REF SAMTOOLS
grep -v '^#' "$LIST" | awk 'NF {print $NF}' | xargs -P "$THREADS" -I{} bash -c 'one "$@"' _ {}
echo "subset BAMs in $OUT: $(ls "$OUT"/*.bam | wc -l)"
