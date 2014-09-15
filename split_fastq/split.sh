#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input)
splits=$(./swe get splits)

if [[ $input =~ paired\[(.*),(.*)\] ]] ;
    then
	 file1=$(./swe fetch ${BASH_REMATCH[1]})
	 file2=$(./swe fetch ${BASH_REMATCH[2]})
    else
	echo Not implemented
	false
fi


./split_input_fastq.pl --input paired[$file1,$file2] --splits $splits

for i in $(seq 1 $splits)
do
    ./swe emit file $i.fastq.gz
done
