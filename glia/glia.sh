#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input | ./swe fetch -)
interval_file=$(./swe get interval|./swe fetch -)
interval=$(cat $interval_file)
gatk_data=$(./swe get GATK_DATA)
candidates=$(./swe get candidates | ./swe fetch -)

tabix -p vcf $candidates

<$input glia -Rru -w 1000 -S 100 -Q 100 -G 4 -f $gatk_data/hg19/ucsc.hg19.fasta -v $candidates \
    | freebayes -f $gatk_data/hg19/ucsc.hg19.fasta \
        --stdin --variant-input $candidates --only-use-input-alleles > raw.vcf

./swe emit file raw.vcf
