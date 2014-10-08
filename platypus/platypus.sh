#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input | ./swe fetch -)
interval_file=$(./swe get interval|./swe fetch -)
interval=$(cat $interval_file)
gatk_data=$(./swe get GATK_DATA)
cpu_cores=32

samtools index $input

platypus callVariants \
    --refFile $gatk_data/hg19/ucsc.hg19.fasta \
    --bamFiles $input \
    -o raw.vcf

./swe emit file raw.vcf
