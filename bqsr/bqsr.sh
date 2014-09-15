#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input)
bai=$(./swe fetch $input.bai)
input=$(./swe fetch $input)


chr=$(./swe get chr)
gatk_data=$(./swe get GATK_DATA)
gatk_jar=$(./swe get gatk_jar)
[ -e $gatk_jar ]
cpu_cores=32


samtools index $input #fix me

tabix -h $gatk_data/dbsnp_137.hg19.vcf.gz $chr:1-300000000 > interval.dbsnp_137.hg19.vcf
java -Xmx6g -jar $gatk_jar \
    -T BaseRecalibrator \
    -I $input \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
    -knownSites interval.dbsnp_137.hg19.vcf \
    --covariate ContextCovariate  \
    --covariate ReadGroupCovariate  \
    --covariate QualityScoreCovariate \
    --covariate CycleCovariate  \
    -o bqsr.grp  \
    --disable_indel_quals \
    -nct $cpu_cores 

./swe emit file bqsr.grp
exit 0

