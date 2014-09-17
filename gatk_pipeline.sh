#!/bin/bash
set -e -x
set -o pipefail

. ./project.settings
#INPUT_FASTQ=paired[s3://gapp-west/test@clusterk.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_1.fastq.gz,s3://gapp-west/test@clusterk.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_2.fastq.gz]
GATK_JAR=s3://gapp-west/test@clusterk.com/GenomeAnalysisTK.jar

export K_GATK_DATA=/tmp/gatk-reference
export K_ANALYSIS=$K_ANALYSIS
export GATK_JAR=./bin/GenomeAnalysisTKLite.jar
export PATH=$PATH:./bin
export SHELL=/bin/bash

CHROMOSOMES="chr22 chr21 chr20 chr19 chr18 chr17 chr16 chr15 chr14 chr13 chr12 chr11 chr10 chr9 chr8 chr7 chr6 chr5 chr4 chr3 chr2 chr1 chrY chrX"
#CHROMOSOMES="chr22"
# chr21 chr20 chr19 chr18 chr17 chr16 chr15 chr14 chr13 chr12 chr11 chr10 chr9 chr8 chr7 chr6 chr5 chr4 chr3 chr2 chr1 chrY chrX"

BQSR_CHR=chr22

# split input files into chunks


#  SWE_ENGINE clusterk or local
export SWE_ENGINE=local



if [ "$SWE_ENGINE" = "clusterk" ]
then
    export SWE_QUEUE=default
    export SWE_S3_STORAGE=s3://gapp-west/swe-test
    export SWE_KSUB_EXTRA_PARAMS=" -u bin.tar.gz -t ANALYSIS=$K_ANALYSIS -c auto:ANALYSIS,STAGE -e auto:ANALYSIS,STAGE -m 10000 -dm 20000 -de auto:ANALYSIS,STAGE -th auto:ANALYSIS,STAGE  "
else
    export SWE_DEV_MODE=1 # cache successfull task IDs
    export SWE_KSUB_EXTRA_PARAMS=" -u bin.tar.gz "

fi

[ -e bin.tar.gz ] || tar czvf bin.tar.gz ./bin

NAME_PREFIX="Sample:$K_ANALYSIS";

#####################################################################################################
# Process input files:
# 1. Split each input file into chunks of 1GB
# 2. Launch alignment job per each chunk
# 3. Each alignment will return one sorted file per chromosome, which will later be combined.


#list of alignment job IDs
align_job_ids=""
for input in $INPUT_FASTQ
do
	#check input string format
	#compute number of alignment splits. one split per GB of input data.
	if [[ $input =~ paired\[(.*),(.*)\] ]] ;
	then
		 file1=${BASH_REMATCH[1]}
		 file2=${BASH_REMATCH[2]}
		 file1_size=$(es3 ls $file1 | head -n 1| cut -f 2)

	elif [[ $input =~ local\[(.*),(.*)\] ]] ;
	then
		file1=$(./swe store ${BASH_REMATCH[1]})
		file2=$(./swe store ${BASH_REMATCH[2]})
		file1_size=$(stat -c%s ${BASH_REMATCH[1]})
	
	else
		echo Not implemented
		false
	fi

	splits=$[$file1_size/1000000000+1]

	echo Processing $input, will split it in $splits chunks

	# input.sh:  accepts input fastq, splits in N chunks, saves them as N.fastq.gz
	split_job_id=$(./swe submit \
						 -u split_fastq/split.sh \
						 -u split_fastq/split_input_fastq.pl \
						 -t NAME=$NAME_PREFIX:interleave -t STAGE=interleave\
						 -isplits=$splits \
						 -iinput1=$file1 \
						 -iinput2=$file2 \
						 --wrap="bash split.sh")

	#align each split, reference them via $split_job_id

	for split in $(seq 1 $splits)
	do
		# align.sh: accepts input split file, produces alignment, split by chromosome
		align_job_id=$(./swe submit \
							 -u align/align.sh \
							 -d $split_job_id \
							 -t NAME=$NAME_PREFIX:align -t STAGE=align\
							 -isample_id=SAMPLE \
							 -iinput=$split_job_id:$split.fastq.gz \
							 --wrap="bash align.sh")
		align_job_ids="$align_job_ids $align_job_id"

	done

done
#exit 1

########### GATK  stage
# 1. combine chromosome files from each alignment job into single BAM file per chromosome
# 2. Use one of the chromosomes for base quality recalibration (chr22)
# 3. For each chromosome, run split_chr to find safe spot where chromosome can be broken
# 4. For each resulting interval run GATK HaplotypeCaller and the rest of Best Practices GATK pipeline
#
#

# number of splits per chromosomes
gatk_splits=10
#combine per chromosomes

for chr in $CHROMOSOMES
do
	#get chromosome size and compute optimal number of splits
	
	chr_size=$(grep "^$chr	" $K_GATK_DATA/hg19/ucsc.hg19.fasta.fai |cut -f 2)
	gatk_splits=$[$chr_size/10000000]
	[ "$gatk_splits" != "0" ]

	#create comma separated list of alignment jobs
	align_job_list=$(echo $align_job_ids |tr " " ",")

	#create list of input alignment files for current chromosome
	input_array=""
	for align_job in $align_job_ids
	do
		input_array="$input_array $align_job:$chr.bam"
	done

	#submit comine jobs, and pass list of alignent files as input, and all alignment jobs as prerequsite
	#combine.sh: accepts a list of aligned bam files, produced combined file for a given chromsome
	combine_job_id=$(./swe submit \
							-d $align_job_list \
							-ichr=$chr \
							-iinput="$input_array" \
							-u combine/combine.sh \
							-t NAME=$NAME_PREFIX:combine:$chr -t STAGE=combine \
							--wrap="bash combine.sh" )

	#for one of the chromosomes submit BQSR job to produce Base Quality Recalibration files
	if [ "$chr" == "$BQSR_CHR" ]
	then
		#bqsr.sh: runs Base Quality recalibration on chr22
		bqsr_job_id=$(./swe submit \
						     -d $combine_job_id \
						     -iinput=$combine_job_id:$chr.bam \
						     -t NAME=$NAME_PREFIX:bqsr:$chr -t STAGE=bqsr \
						     -u bqsr/bqsr.sh \
						     -ichr=$chr \
						     -igatk_jar=$GATK_JAR \
						     --wrap="bash bqsr.sh")

	fi
	

	# for each combined chromosomes, run split_chr.sh to obtain list of 
	#chr_split.sh: finds genomic locations where it is safe to split a chromosomes
	#               returns list of bam files:  split.$chr.$split_id.bam
	chr_split_id=$(./swe submit \
						  -d $combine_job_id \
						  -iinput=$combine_job_id:$chr.bam \
						  -isplits=$gatk_splits \
						  -ichr=$chr \
						  -t NAME=$NAME_PREFIX:chr_split:$chr -t STAGE=chr_split \
						  -u split_chr/split_chr.sh \
						  -u split_chr/advanced_splitter.pl \
						  -u split_chr/breakpoints2intervals.pl \
						  -u split_chr/equally_spaced_intervals.pl \
						  --wrap="bash split_chr.sh")

	#for each chromosome split, start a GATK analysis job
	for split_id in $(seq 1 $gatk_splits)
	do

		#gatk.sh runs gatk on a sub-interval and applies BQSR
		gatk_job_id=$(./swe submit \
							-d $chr_split_id,$bqsr_job_id \
							-iinput=$chr_split_id:$split_id.bam \
							-iinterval=$chr_split_id:$split_id.interval \
							-ibqsr=$bqsr_job_id:bqsr.grp \
							-t NAME=$NAME_PREFIX:gatk:$chr:$split_id -t STAGE=gatk \
							-u gatk/gatk.sh \
							-igatk_jar=$GATK_JAR \
							-u gatk/MarkDuplicates.jar \
							--wrap="bash gatk.sh")
	[ "$swe_wait" == "" ] || ./swe wait $gatk_job_id

		gatk_job_ids="$gatk_job_ids $gatk_job_id"
		
	done
done

#exit 1
# collect all GATK job ids, and create comma separated list for dependendencies (-d)
gatk_job_list=$(echo $gatk_job_ids |tr " " ",")

input_array=""
for gatk_job in $gatk_job_ids
do
    input_array="$input_array $gatk_job:raw.vcf"
done

#submit a job that combines all sub-region vcf files, into one sorted VCF
#combine_vcf.sh concatenate sub-region VCFs into final VCF
combine_vcf_job_id=$(./swe submit \
			    -d $gatk_job_list \
			    -iinput="$input_array" \
			    -t NAME=$NAME_PREFIX:combine_vcf -t STAGE=combine_vcf \
			    -u combine_vcf/vcf-sort \
			    -u combine_vcf/combine_vcf.sh \
			    --wrap="bash combine_vcf.sh" )

#submit a job that run 
#run variant quality recalibration
vqsr_job_id=$(./swe submit \
		    -d $combine_vcf_job_id \
		    -u vqsr/vqsr.sh \
		    -igatk_jar=$GATK_JAR \
		    -iinput=$combine_vcf_job_id:raw.vcf.gz \
		    --wrap=" bash vqsr.sh")

#wait for all jobs to finish
./swe wait $vqsr_job_id


exit 0


