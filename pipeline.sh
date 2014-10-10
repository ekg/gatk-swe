#!/bin/bash
set -e -x
set -o pipefail

. ./general_settings.sh

. $1

#INPUT_FASTQ=paired[s3://gapp-west/test@clusterk.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_1.fastq.gz,s3://gapp-west/test@clusterk.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_2.fastq.gz]

export GATK_REFERENCE=s3://gapp-east/gatk-reference/
export ANALYSIS=$K_ANALYSIS
export PATH=$PATH:./bin
export SHELL=/bin/bash


[ "$SWE_ENGINE"  != "" ] || export SWE_ENGINE=clusterk




[ "$CHROMOSOMES" != "" ] || export CHROMOSOMES="chr22 chr21 chr20 chr19 chr18 chr17 chr16 chr15 chr14 chr13 chr12 chr11 chr10 chr9 chr8 chr7 chr6 chr5 chr4 chr3 chr2 chr1 chrY chrX"


#export SWE_DEV_MODE=1 # cache successfull task IDs


GENOME_FAI=./bin/ucsc.hg19.fasta.fai
[ -e $GENOME_FAI ] 


export SWE_S3_STORAGE=s3://gapp-east-temp
export common_params=" -p 10 -u ./swe -u pre_script.sh -u bin.tar.gz -t ANALYSIS=$K_ANALYSIS -c auto:ANALYSIS,STAGE -e auto:ANALYSIS,STAGE -m 15000 -dm 40000 -de auto:ANALYSIS,STAGE -th auto:ANALYSIS,STAGE  "


# create a new job
job_id=$(kjob add -n $NAME-`date "+%Y-%m-%d-%H:%M:%S"`)
common_params="$common_params -j $job_id"


#common_params="$common_params -j 1"
#if a GATK 3.0 jar is specified - use it. Otherwise fall back to GATK_Lita
#upload GATK_JAR
if [ "$GATK_FULL_JAR" != "" ]
then
    GATK_FULL_JAR=$(./swe store ./GenomeAnalysisTK.jar)
fi
if [ "$GATK_LITE_JAR" != "" ]
then
    # if HTTP link for 3.0 is not provided fall back on free GATK_Lite
    GATK_LITE_JAR=$(./swe store ./GenomeAnalysisTKLite.jar)
fi


[ "$GATK_LITE_JAR" != "" ] #GATK_JAR must be defined
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

# also, you can just use a number of line
#input_split_lines=8000000

#list of alignment job IDs
align_job_ids=""
for input in $INPUT_FASTQ
do
    echo $input
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
	split_job_id=$(ksub $common_params \
						 -u split_fastq/split.sh \
						 -u split_fastq/split_input_fastq.pl \
						 -t NAME=$NAME_PREFIX:split -t STAGE=split \
						 -v splits=$splits \
						 -v input1=$file1 \
						 -v input2=$file2 \
						 -c 16 \
						 --wrap="bash split.sh")

    # count the splits we've made
    #./swe wait $split_job_id
    #cat $(./swe fetch $split_job_id:splits.txt)

	#align each split, reference them via $split_job_id

	for split in $(seq 1 $splits)
	do
		# align.sh: accepts input split file, produces alignment, split by chromosome
		align_job_id=$(ksub $common_params  \
							 -u align/align.sh \
							 -d $split_job_id \
							 -t NAME=$NAME_PREFIX:align -t STAGE=align\
							 -v sample_id=SAMPLE \
							 -v input=$split_job_id:$split.fastq.gz \
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
region_splits=10
#combine per chromosomes


for chr in $CHROMOSOMES
do
	#get chromosome size and compute optimal number of splits
	
	chr_size=$(grep "^$chr	" $GENOME_FAI |cut -f 2)
	region_splits=$[$chr_size/5000000]
	[ "$region_splits" != "0" ]

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
 	combine_job_id=$(ksub  $common_params \
		-d $align_job_list \
		-v chr=$chr \
		-v input="$input_array" \
		-u combine/combine.sh \
        -u MarkDuplicates.jar \
		-c 8 \
		-t NAME=$NAME_PREFIX:combine:$chr -t STAGE=combine \
		--wrap="bash combine.sh" )

	#for one of the chromosomes submit BQSR job to produce Base Quality Recalibration files
	#if [ "$chr" == "$BQSR_CHR" ]
    #if [ false ]q
    #then
		#bqsr.sh: runs Base Quality recalibration on chr22
		# output is $bqsr_job_id:bqsr.grp
	#	bqsr_job_id=$(ksub  $common_params \
	#					     -d $combine_job_id \
	#					     -v input=$combine_job_id:$chr.bam \
	#					     -t NAME=$NAME_PREFIX:bqsr:$chr -t STAGE=bqsr \
	#					     -u bqsr/bqsr.sh \
	#					     -v chr=$chr \
	#					     -v gatk_jar=$GATK_JAR \
	#					     --wrap="bash bqsr.sh")
#
#	fi
	
	# for each combined chromosomes, run split_chr.sh to obtain list of 
	#chr_split.sh: finds genomic locations where it is safe to split a chromosomes
	#               returns list of bam files:  split.$chr.$split_id.bam
	chr_split_id=$(ksub  $common_params \
		-d $combine_job_id \
		-v input=$combine_job_id:$chr.bam \
		-v splits=$region_splits \
		-v chr=$chr \
		-t NAME=$NAME_PREFIX:chr_split:$chr -t STAGE=chr_split \
		-u split_chr/split_chr.sh \
		-c 8 \
		-u split_chr/advanced_splitter.pl \
		-u split_chr/breakpoints2intervals.pl \
		-u split_chr/equally_spaced_intervals.pl \
		--wrap="bash split_chr.sh")

	#for each chromosome split, start an analysis job from each caller
    
	for split_id in $(seq 1 $region_splits)
	do

		#gatk.sh runs gatk on a sub-v nterval and applies BQSR
		# output: $gatk_job_id:raw.vcf
		#					-d $chr_split_id,$bqsr_job_id \
		
        gatk_ug_job_id=$(ksub  $common_params \
            -d $chr_split_id \
		    -v input=$chr_split_id:$split_id.bam \
		    -v interval=$chr_split_id:$split_id.interval \
		    -t NAME=$NAME_PREFIX:gatk-ug:$chr:$split_id -t STAGE=gatk-ug \
	        -u gatk-ug/gatk-ug.sh \
		    -v gatk_jar=$GATK_LITE_JAR \
		    --wrap="bash gatk-ug.sh")

        gatk_hc_job_id=$(ksub  $common_params \
            -d $chr_split_id \
		    -v input=$chr_split_id:$split_id.bam \
		    -v interval=$chr_split_id:$split_id.interval \
		    -t NAME=$NAME_PREFIX:gatk-hc:$chr:$split_id -t STAGE=gatk-hc \
	        -u gatk-hc/gatk-hc.sh \
		    -v gatk_jar=$GATK_FULL_JAR \
		    --wrap="bash gatk-hc.sh")

        freebayes_job_id=$(ksub $common_params \
            -d $chr_split_id \
            -v input=$chr_split_id:$split_id.bam \
            -v interval=$chr_split_id:$split_id.interval \
            -t NAME=$NAME_PREFIX:freebayes:$chr:$split_id -t STAGE=freebayes \
            -u freebayes/freebayes.sh \
            --wrap="bash freebayes.sh")

        platypus_job_id=$(ksub $common_params \
            -d $chr_split_id \
            -v input=$chr_split_id:$split_id.bam \
            -v interval=$chr_split_id:$split_id.interval \
            -t NAME=$NAME_PREFIX:platypus:$chr:$split_id -t STAGE=platypus \
            -u platypus/platypus.sh \
            --wrap="bash platypus.sh")


	    #[ "$swe_wait" == "" ] || ./swe wait $caller_job_id

		caller_job_ids="$caller_job_ids $gatk_ug_job_id $gatk_hc_job_id $freebayes_job_id $platypus_job_id"
		
	done
done

# collect all caller job ids, and create comma separated list for dependendencies (-d)
caller_job_list=$(echo $caller_job_ids |tr " " ",")

input_array=""
for caller_job in $caller_job_ids
do
    input_array="$input_array $caller_job:raw.vcf"
done

#submit a job that combines all sub-region vcf files, into one sorted VCF
#combine_vcf.sh concatenate sub-region VCFs into final VCF
combine_vcf_job_id=$(ksub  $common_params \
			    -d $caller_job_list \
			    -v input="$input_array" \
			    -t NAME=$NAME_PREFIX:combine_vcf -t STAGE=combine_vcf \
			    -u combine_vcf/combine_vcf.sh \
			    --wrap="bash combine_vcf.sh" )


for chr in $CHROMOSOMES
do
	
	chr_size=$(grep "^$chr	" $GENOME_FAI |cut -f 2)
	region_splits=$[$chr_size/5000000]
	[ "$region_splits" != "0" ]

	#for each chromosome split, start an analysis job from each caller
    
	for split_id in $(seq 1 $region_splits)
	do

		#gatk.sh runs gatk on a sub-v nterval and applies BQSR
		# output: $gatk_job_id:raw.vcf
		#					-d $chr_split_id,$bqsr_job_id \
		
        glia_job_id=$(ksub $common_params \
            -d $combine_vcf_job_id \
            -v input=$chr_split_id:$split_id.bam \
            -v interval=$chr_split_id:$split_id.interval \
            -v candidates=$combine_vcf_job_id:raw.vcf.gz \
            -t NAME=$NAME_PREFIX:glia:$chr:$split_id -t STAGE=glia \
            -u glia/glia.sh \
            --wrap="bash glia.sh")

		gl_job_ids="$gl_job_ids $glia_job_id"
		
	done
done

gl_job_list=$(echo $gl_job_ids |tr " " ",")

input_array=""
for gl_job in $gl_job_ids
do
    input_array="$input_array $gl_job:raw.vcf"
done

combine_final_vcf_job_id=$(ksub  $common_params \
			    -d $gl_job_list \
			    -v input="$input_array" \
			    -t NAME=$NAME_PREFIX:combine_vcf -t STAGE=combine_vcf \
			    -u combine_vcf/combine_vcf.sh \
			    --wrap="bash combine_vcf.sh" )

#submit a job that run 
#run variant quality recalibration
# output: $vqsr_job_id:recalibrated.filtered.vcf.gz
#vqsr_job_id=$(ksub  $common_params \
#		    -d $combine_vcf_job_id \
#		    -u vqsr/vqsr.sh \
#		    -v gatk_jar=$GATK_JAR \
#		    -t NAME=$NAME_PREFIX:vqsr -t STAGE=vqsr \
#		    -v input=$combine_vcf_job_id:raw.vcf.gz \
#		    --wrap=" bash vqsr.sh")

#wait for all jobs to finish
kwait $combine_final_vcf_job_id
./swe fetch $combine_final_vcf_job_id:raw.vcf.gz

exit 0


