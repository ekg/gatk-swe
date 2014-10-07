#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input | ./swe fetch -)
interval_file=$(./swe get interval|./swe fetch -)
interval=$(cat $interval_file)
gatk_data=$(./swe get GATK_DATA)
cpu_cores=32

freebayes -f $gatk_data/hg19/ucsc.hg19.fasta $input > raw.vcf

./swe emit file raw.vcf
