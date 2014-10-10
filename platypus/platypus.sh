#!/bin/bash
set -e
set -x
set -o pipefail

[ "$input"     != "" ] && input=$(./swe fetch $input)
[ "$gatk_jar"  != "" ] && gatk_jar=$(./swe fetch $gatk_jar)
[ "$interval"  != "" ] && interval_file=$(./swe fetch $interval)
[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)

cpu_cores=32

samtools index $input

platypus callVariants \
    --refFile $gatk_data/hg19/ucsc.hg19.fasta \
    --bamFiles $input \
    -o raw.vcf

./swe emit file raw.vcf
