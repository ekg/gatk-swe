#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input)
chr=$(./swe get chr)
cpu_cores=32


# compute number of input files
words=( $input )

if [ "${#words[@]}" == "1" ]
then
    #single input file, nothing to merge
    cp $(./swe fetch $input) $chr.bam
    cp $(./swe fetch $input.bai ) $chr.bam.bai
else
    # multiple input files

    local_bams=""
    for bam in $input
    do
	local_bam=$(./swe fetch $bam)
	local_bams="$local_bams $local_bam"
    done

    samtools merge -@ $cpu_cores $chr.bam $local_bams
    samtools index $chr.bam
fi

    ./swe emit file $chr.bam
    ./swe emit file $chr.bam.bai

exit 0
