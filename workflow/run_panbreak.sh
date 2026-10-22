#!/bin/bash
set -euo pipefail
if [ $# -ne 1 ]; then
  echo "usage: run_panbreak.sh <config.sh>" >&2
  exit 1
fi
source "$1"
PANBREAK=${PANBREAK:-panbreak}
THREADS=${THREADS:-8}
ENGINE=${ENGINE:-docker}
FLANK=${FLANK:-20000}
OUT=${OUT:-panbreak_out}
mkdir -p "$OUT"

if [ -n "${GENE:-}" ]; then
  $PANBREAK regions -a "$ASSEMBLIES" -r "$REF" -g "$GENE" -l "$LOCUS" -f "$FLANK" -t "$THREADS" -o "$OUT/01_regions"
else
  $PANBREAK regions -a "$ASSEMBLIES" -r "$REF" -L "$REF_LOCUS" -l "$LOCUS" -f "$FLANK" -t "$THREADS" -o "$OUT/01_regions"
fi
$PANBREAK graph -i "$OUT/01_regions/locus.fa.gz" -r "$REF" -e "$ENGINE" -t "$THREADS" ${PGGB_ARGS:+--pggb-args "$PGGB_ARGS"} \
  ${TRIM_UP:+--regions "$OUT/01_regions/regions.tsv" --trim-up "$TRIM_UP" --trim-down "${TRIM_DOWN:-$TRIM_UP}"} -o "$OUT/02_graph"
$PANBREAK bubbles -g "$OUT/02_graph/locus.gfa" -r "$REF" --regions "$OUT/01_regions/regions.tsv" -m "${MIN_SV:-50}" -o "$OUT/03_bubbles"
$PANBREAK annotate -b "$OUT/03_bubbles" -r "$REF" --regions "$OUT/01_regions/regions.tsv" ${REF_GFF:+--gff "$REF_GFF"} \
  ${TE_LIB:+--te-lib "$TE_LIB"} -t "$THREADS"
$PANBREAK kmers -b "$OUT/03_bubbles" --fasta "$OUT/01_regions/locus.fa.gz" ${BACKGROUND:+--background $BACKGROUND} -t "$THREADS" -o "$OUT/04_kmers"
$PANBREAK genotype -k "$OUT/04_kmers" -b "$OUT/03_bubbles" -B "$BAMS" ${BAM_REF:+--bam-ref "$BAM_REF"} ${BAM_REGION:+--bam-region "$BAM_REGION"} \
  --ploidy "${PLOIDY:-2}" -t "$THREADS" -o "$OUT/05_genotypes"
$PANBREAK popstats -G "$OUT/05_genotypes/genotypes.tsv" ${SAMPLES:+-s "$SAMPLES" -c "${GROUP_COL:-population}"} -b "$OUT/03_bubbles" \
  --ploidy "${PLOIDY:-2}" -o "$OUT/06_popstats"
$PANBREAK plot -b "$OUT/03_bubbles" -r "$REF" --gfa "$OUT/02_graph/locus.gfa" --regions "$OUT/01_regions/regions.tsv" ${REF_GFF:+--gff "$REF_GFF"} \
  --genotypes "$OUT/05_genotypes" --popstats "$OUT/06_popstats" ${SAMPLES:+-s "$SAMPLES" -c "${GROUP_COL:-population}"} -o "$OUT/07_plots"
