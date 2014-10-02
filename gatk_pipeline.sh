#!/bin/bash
set -e -x
set -o pipefail

. ./general_settings.sh

. $1

#INPUT_FASTQ=paired[s3://gapp-west/test@clusterk.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_1.fastq.gz,s3://gapp-west/test@clusterk.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_2.fastq.gz]

export K_GATK_DATA=/tmp/gatk-reference
export K_ANALYSIS=$K_ANALYSIS
export PATH=$PATH:./bin
export SHELL=/bin/bash


[ "$SWE_ENGINE"  != "" ] || export SWE_ENGINE=clusterk




#chromosomes should start with BQSR_CHROMOSOME
BQSR_CHR=chr22
[ "$CHROMOSOMES" != "" ] || export CHROMOSOMES="chr22 chr21 chr20 chr19 chr18 chr17 chr16 chr15 chr14 chr13 chr12 chr11 chr10 chr9 chr8 chr7 chr6 chr5 chr4 chr3 chr2 chr1 chrY chrX"


#export SWE_DEV_MODE=1 # cache successfull task IDs


GENOME_FAI=./bin/ucsc.hg19.fasta.fai
[ -e $GENOME_FAI ] 



if [ "$SWE_ENGINE" = "clusterk" ]
then
    export SWE_QUEUE=default
    export SWE_S3_STORAGE=s3://gapp-east-temp
    export common_params=" -u bin.tar.gz -t ANALYSIS=$K_ANALYSIS -c auto:ANALYSIS,STAGE -e auto:ANALYSIS,STAGE -m 15000 -dm 40000 -de auto:ANALYSIS,STAGE -th auto:ANALYSIS,STAGE  "
    export SWE_DEV_MODE=1 # cache successfull task IDs
else
    export SWE_KSUB_EXTRA_PARAMS=" -u bin.tar.gz "

fi


common_params="$common_params -j 1"
#if a GATK 3.0 jar is specified - use it. Otherwise fall back to GATK_Lita
#upload GATK_JAR
if [ "$GATK_JAR" != "" ]
then
GATK_JAR=$(./swe store $GATK_JAR)
else
    # if HTTP link for 3.0 is not provided fall back on free GATK_Lite
GATK_JAR=$(./swe store ./GenomeAnalysisTKLite.jar)
fi


[ "$GATK_JAR" != "" ] #GATK_JAR must be defined
[ -e bin.tar.gz ] || tar czvf bin.tar.gz ./bin

NAME_PREFIX="$NAME:$K_ANALYSIS";

#####################################################################################################
# Process input files:
# 1. Split each input file into chunks of 1GB
# 2. Launch alignment job per each chunk
# 3. Each alignment will return one sorted file per chromosome, which will later be combined.
#
# supported input formats:   paired[s3_path_for_read1,s3_path_for_read2]
#                            local[local_path_for_read1,local_path_for_read2]
#
#

#approximate size of input splits in bytes. 500MB 
input_split_size=500000000

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

	splits=$[$file1_size/$input_split_size+1]

	echo Processing $input, will split it in $splits chunks

	# input.sh:  accepts input fastq, splits in N chunks, saves them as N.fastq.gz
	split_job_id=$(./swe submit $common_params \
						 -u split_fastq/split.sh \
						 -u split_fastq/split_input_fastq.pl \
						 -t NAME=$NAME_PREFIX:interleave -t STAGE=interleave\
						 -isplits=$splits \
						 -iinput1=$file1 \
						 -iinput2=$file2 \
						 -c 16 \
						 --wrap="bash split.sh")
	
	#align each split, reference them via $split_job_id

	for split in $(seq 1 $splits)
	do
		# align.sh: accepts input split file, produces alignment, split by chromosome
		align_job_id=$(./swe submit $common_params  \
							 -u align/align.sh \
							 -d $split_job_id \
							 -t NAME=$NAME_PREFIX:align -t STAGE=align\
							 -isample_id=SAMPLE \
							 -iinput=$split_job_id:$split.fastq.gz \
							 --wrap="bash align.sh")
		align_job_ids="$align_job_ids $align_job_id"
	done

done


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
	
	chr_size=$(grep "^$chr	" $GENOME_FAI |cut -f 2)
	gatk_splits=$[$chr_size/5000000]
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
	# output: $combine_job_id:$chr.bam
	combine_job_id=$(./swe submit  $common_params \
							-d $align_job_list \
							-ichr=$chr \
							-iinput="$input_array" \
							-u combine/combine.sh \
							-c 8 \
							-t NAME=$NAME_PREFIX:combine:$chr -t STAGE=combine \
							--wrap="bash combine.sh" )

	#for one of the chromosomes submit BQSR job to produce Base Quality Recalibration files
	if [ "$chr" == "$BQSR_CHR" ]
	then
		#bqsr.sh: runs Base Quality recalibration on chr22
		# output is $bqsr_job_id:bqsr.grp
		bqsr_job_id=$(./swe submit  $common_params \
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
	chr_split_id=$(./swe submit  $common_params \
						  -d $combine_job_id \
						  -iinput=$combine_job_id:$chr.bam \
						  -isplits=$gatk_splits \
						  -ichr=$chr \
						  -t NAME=$NAME_PREFIX:chr_split:$chr -t STAGE=chr_split \
						  -u split_chr/split_chr.sh \
						  -c 8 \
						  -u split_chr/advanced_splitter.pl \
						  -u split_chr/breakpoints2intervals.pl \
						  -u split_chr/equally_spaced_intervals.pl \
						  --wrap="bash split_chr.sh")

	#for each chromosome split, start a GATK analysis job
	for split_id in $(seq 1 $gatk_splits)
	do

		#gatk.sh runs gatk on a sub-interval and applies BQSR
		# output: $gatk_job_id:raw.vcf
		gatk_job_id=$(./swe submit  $common_params \
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

# collect all GATK job ids, and create comma separated list for dependendencies (-d)
gatk_job_list=$(echo $gatk_job_ids |tr " " ",")

input_array=""
for gatk_job in $gatk_job_ids
do
    input_array="$input_array $gatk_job:raw.vcf"
done

#submit a job that combines all sub-region vcf files, into one sorted VCF
#combine_vcf.sh concatenate sub-region VCFs into final VCF
combine_vcf_job_id=$(./swe submit  $common_params \
			    -d $gatk_job_list \
			    -iinput="$input_array" \
			    -t NAME=$NAME_PREFIX:combine_vcf -t STAGE=combine_vcf \
			    -u combine_vcf/vcf-sort \
			    -u combine_vcf/combine_vcf.sh \
			    --wrap="bash combine_vcf.sh" )

#submit a job that run 
#run variant quality recalibration
# output: $vqsr_job_id:recalibrated.filtered.vcf.gz
vqsr_job_id=$(./swe submit  $common_params \
		    -d $combine_vcf_job_id \
		    -u vqsr/vqsr.sh \
		    -igatk_jar=$GATK_JAR \
		    -t NAME=$NAME_PREFIX:vqsr -t STAGE=vqsr \
		    -iinput=$combine_vcf_job_id:raw.vcf.gz \
		    --wrap=" bash vqsr.sh")

#wait for all jobs to finish
./swe wait  $vqsr_job_id
./swe fetch $vqsr_job_id:recalibrated.filtered.vcf.gz

exit 0


