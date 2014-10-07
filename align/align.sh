#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input | ./swe fetch -)
sample=$(./swe get sample_id)
gatk_data=$(./swe get GATK_DATA)
cpu_cores=32
GROUP_ID="@RG\tID:1\tPL:ILLUMINA\tPU:pu\tLB:group1\tSM:Sample_XXX"

sort_dir=$$.sambamba_sort_dir
mkdir $sort_dir

#samtools sort is limited to 62M per thread, so that it won't take more than 2GB total
bwa mem -M -p -t $cpu_cores -R "$GROUP_ID" $gatk_data/hg19/ucsc.hg19.fasta $input \
    | sambamba view -S -f bam -l 0 /dev/stdin \
    | sambamba sort -t $cpu_cores -m $[2000/$cpu_cores]M --tmpdir=$sort_dir -o raw.bam /dev/stdin

    #| samtools view -@ $cpu_cores -1 -bt   $gatk_data/hg19/ucsc.hg19.fasta.fai - \
    #| samtools sort -@ $cpu_cores -m $[2000/$cpu_cores]M -l 0 - raw

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
