#!/bin/bash
set -e
set -x
set -o pipefail

[ "$input"     != "" ] && input=$(./swe fetch $input)
[ "$interval"  != "" ] && interval_file=$(./swe fetch $interval)
[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)

interval=$(cat $interval_file)
cpu_cores=32

freebayes -f $gatk_data/hg19/ucsc.hg19.fasta $input > raw.vcf

./swe emit file raw.vcf
