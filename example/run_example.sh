#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
PANBREAK=${PANBREAK:-panbreak}
THREADS=${THREADS:-8}
ENGINE=${ENGINE:-docker}
R=results

[ -s data/samples.tsv ] || python3 make_example.py -t "$THREADS"

$PANBREAK regions -a data/assemblies.tsv -r ATG1 -g geneX -l x -f 15000 -o $R/01_regions
$PANBREAK graph -i $R/01_regions/locus.fa.gz -r ATG1 -e "$ENGINE" -t "$THREADS" -o $R/02_graph
$PANBREAK bubbles -g $R/02_graph/locus.gfa -r ATG1 --regions $R/01_regions/regions.tsv -o $R/03_bubbles
$PANBREAK annotate -b $R/03_bubbles -r ATG1 --gff data/assemblies/ATG1.gff3 --regions $R/01_regions/regions.tsv \
  --te-lib data/te_library.fa -t "$THREADS"
$PANBREAK kmers -b $R/03_bubbles --fasta $R/01_regions/locus.fa.gz -o $R/04_kmers
$PANBREAK genotype -k $R/04_kmers -b $R/03_bubbles -B data/bams --bam-ref data/assemblies/ATG1.fa \
  --bam-region chr1:20000-80000 -t "$THREADS" -o $R/05_genotypes
$PANBREAK popstats -G $R/05_genotypes/genotypes.tsv -s data/samples.tsv -b $R/03_bubbles -o $R/06_popstats
$PANBREAK plot -b $R/03_bubbles -r ATG1 --gff data/assemblies/ATG1.gff3 --gfa $R/02_graph/locus.gfa --regions $R/01_regions/regions.tsv \
  --genotypes $R/05_genotypes --popstats $R/06_popstats -s data/samples.tsv -o $R/07_plots
python3 evaluate.py
