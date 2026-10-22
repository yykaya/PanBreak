#!/bin/bash
#SBATCH --job-name=panbreak_gt
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/panbreak_gt_%A_%a.log
set -euo pipefail
KMERS=${KMERS:-results/04_kmers}
BUBBLES=${BUBBLES:-results/03_bubbles}
BAM_LIST=${BAM_LIST:-bams.txt}
OUT=${OUT:-results/05_genotypes}
CHUNKS=${SLURM_ARRAY_TASK_MAX:-20}
TASK=${SLURM_ARRAY_TASK_ID:-1}
mkdir -p logs "$OUT/chunks"
chunk="$OUT/chunks/list_$TASK.txt"
awk -v n="$CHUNKS" -v t="$TASK" 'NF && !/^#/ && (NR - 1) % n == t - 1' "$BAM_LIST" > "$chunk"
panbreak genotype -k "$KMERS" -b "$BUBBLES" -B "$chunk" ${BAM_REF:+--bam-ref "$BAM_REF"} ${BAM_REGION:+--bam-region "$BAM_REGION"} \
  -t "${SLURM_CPUS_PER_TASK:-8}" -o "$OUT/chunks/chunk_$TASK"
