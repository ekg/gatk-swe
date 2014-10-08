#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input | ./swe fetch -)
gatk_jar=$(./swe get gatk_jar | ./swe fetch -)
#bqsr=$(./swe get bqsr | ./swe fetch -)
interval_file=$(./swe get interval|./swe fetch -)
interval=$(cat $interval_file)
gatk_data=$(./swe get GATK_DATA)
cpu_cores=32

tabix -h $gatk_data/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz $interval > VCF1.vcf
tabix -h $gatk_data/1000G_phase1.indels.hg19.vcf.gz $interval                  > VCF2.vcf
tabix -h $gatk_data/dbsnp_137.hg19.vcf.gz $interval         > interval.dbsnp_137.hg19.vcf

samtools index $input

java -Xmx2g -jar $gatk_jar \
    -I $input \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
    -T RealignerTargetCreator \
    -o realigned.intervals \
    --known VCF1.vcf \
    --known VCF2.vcf \
    -L $interval 

java -Xmx2g -jar $gatk_jar \
    -I $input \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
    -T IndelRealigner \
    -targetIntervals realigned.intervals \
    -o realigned.bam \
    -known VCF1.vcf \
    -known VCF2.vcf \
    --consensusDeterminationModel KNOWNS_ONLY \
    -L $interval  \
    -LOD 0.4  \
    --downsample_to_coverage 250 \
    -compress 0

# drop BQSR

#java -Xmx2g -jar $gatk_jar \
#    -R $gatk_data/hg19/ucsc.hg19.fasta \
#    -I realigned.bam  \
#    -T PrintReads  \
#    -o recalibrated.bam  \
#    --disable_indel_quals  \
#    -BQSR $bqsr \
#    -compress 0

# note that for GATK-UG we use the realigned bam, whereas for GATK-HC we use the raw input

if [[ $gatk_jar =~ GenomeAnalysisTKLite ]] ;
then
    java -Xmx6g -jar $gatk_jar \
        -R $gatk_data/hg19/ucsc.hg19.fasta \
        -T UnifiedGenotyper \
        -I realigned.bam \
        --dbsnp interval.dbsnp_137.hg19.vcf \
        -o raw.vcf \
        -glm BOTH \
        -L $interval \
        -stand_call_conf 30.0 \
        -stand_emit_conf 30.0 \
        -nct $cpu_cores \
        -rf BadCigar
else 
    #full featured GATK, run HaplotypeCaller
    java -Xmx6g -jar $gatk_jar \
        -R $gatk_data/hg19/ucsc.hg19.fasta \
        -T HaplotypeCaller \
        -I $input \
        --dbsnp interval.dbsnp_137.hg19.vcf \
        -stand_call_conf 30.0 \
        -stand_emit_conf 30.0 \
        -o raw.hc.vcf \
        -L $interval \
        -nct $cpu_cores \
        -rf BadCigar

    java -Xmx6g -jar $gatk_jar \
        -R $gatk_data/hg19/ucsc.hg19.fasta \
        -T VariantAnnotator \
        -I $input \
        --variant raw.hc.vcf \
        -A MappingQualityZero \
        --dbsnp interval.dbsnp_137.hg19.vcf \
        -o raw.vcf \
        -L $interval \
        -nt $cpu_cores \
        -rf BadCigar

fi

./swe emit file raw.vcf
