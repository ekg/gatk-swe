#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input | ./swe fetch -)
splits=$(./swe get splits)
chr=$(./swe get chr)
gatk_data=$(./swe get GATK_DATA)
cpu_cores=32

samtools index $input #fix it 
./equally_spaced_intervals.pl --bam $input --blocks $splits --chr $chr > es_intervals.txt

cat es_intervals.txt | parallel -j $cpu_cores "./advanced_splitter.pl --reference $gatk_data/hg19/ucsc.hg19.fasta --bam $input --region {} > {}.breakpoints"

[ ! -e breakpoints.lst ] || rm breakpoints.lst

for interval in `cat es_intervals.txt`
do
    cat $interval.breakpoints >> breakpoints.lst 
done

./breakpoints2intervals.pl --bam $input --breakpoints breakpoints.lst --splits=$splits > interval.lst 

for i in $(seq 1 $splits)
do
    interval=$(head -n $i interval.lst|tail -n 1)
    echo $interval > $i.interval
    samtools view -@ $cpu_cores -b $input $interval > $i.bam 
    samtools index $i.bam 
    { ./swe emit file $i.interval || touch emit.failed & }
    { ./swe emit file $i.bam      || touch emit.failed & }
    { ./swe emit file $i.bam.bai  || touch emit.failed & }
done
wait


[ ! -e emit.failed ]
