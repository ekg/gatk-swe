#!/bin/bash
set -e
set -x
set -o pipefail

[ "$input"     != "" ] && input=$(./swe fetch $input)
[ "$interval"  != "" ] && interval_file=$(./swe fetch $interval)
[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)
[ "$candidates" != "" ] && candidates=$(./swe fetch $candidates)

interval=$(cat $interval_file)
tabix -p vcf $candidates

<$input glia -Rru -w 1000 -S 100 -Q 100 -G 4 -f $gatk_data/hg19/ucsc.hg19.fasta -v $candidates \
    | freebayes -f $gatk_data/hg19/ucsc.hg19.fasta \
        --stdin --variant-input $candidates --only-use-input-alleles > raw.vcf

./swe emit file raw.vcf
