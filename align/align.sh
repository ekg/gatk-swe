#!/bin/bash
set -e
set -x
set -o pipefail

[ "$input" != "" ] # input must be defined
input=$(./swe fetch $input)

[ "$sample_id" != "" ] # sample_id must be defined

[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)


cpu_cores=32
GROUP_ID="@RG\tID:1\tPL:ILLUMINA\tPU:pu\tLB:group1\tSM:$sample_id"

#samtools sort is limited to 62M per thread, so that it won't take more than 2GB total
bwa mem -M -p -t $cpu_cores -R "$GROUP_ID" $gatk_data/hg19/ucsc.hg19.fasta $input \
    | samtools view -@ $cpu_cores -1 -bt   $gatk_data/hg19/ucsc.hg19.fasta.fai - \
    | samtools sort -@ $cpu_cores -m $[2000/$cpu_cores]M -l 0 - raw

samtools index raw.bam


#split alignments by chromosomes. Emit separate bam files for each contig:  chr1.bam chr2.bam ...
chr_list=$(samtools idxstats raw.bam| cut -f 1 |grep chr)
for chr in $chr_list
do
	samtools view -@ $cpu_cores -F 4 -b raw.bam  $chr > $chr.bam
	samtools index $chr.bam 
	# run emits in parallel
	{ ./swe emit file $chr.bam     || touch emit.failed & }
	{ ./swe emit file $chr.bam.bai || touch emit.failed & }
done

wait
[ ! -e emit.failed ]

#	samtools view -@ $cpu_cores -f 4 -b raw.bam  > unaligned.bam
#	samtools index unaligned.bam
#	./swe emit file unaligned.bam
#	./swe emit file unaligned.bam.bai

exit 0
