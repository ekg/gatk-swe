#!/bin/bash
set -e
set -x
set -o pipefail

splits=$(./swe get splits)

input1=$(./swe get input1 | ./swe fetch -)
input2=$(./swe get input2 | ./swe fetch -)


./split_input_fastq.pl --input paired[$input1,$input2] --splits $splits
# split_input_fastq.pl will output 1.fastq.gz 2.fastq.gz .... N.fastq.gz  - interleaved FASTQ files

for i in $(seq 1 $splits)
do
    ./swe emit file $i.fastq.gz
done
