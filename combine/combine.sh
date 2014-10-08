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

    # run merge and streaming duplicate marking
    samtools merge -@ $cpu_cores $chr.dups.bam $local_bams
    samtools index $chr.dups.bam
    java -Xmx2g -jar ./MarkDuplicates.jar INPUT=$chr.dups.bam OUTPUT=$chr.bam REMOVE_DUPLICATES=true METRICS_FILE=duplication.metrics
    samtools index $chr.bam

fi

    ./swe emit file $chr.bam
    ./swe emit file $chr.bam.bai

exit 0
