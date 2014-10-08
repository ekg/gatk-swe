#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input)
gatk_data=$(./swe get GATK_DATA)

[ ! -e raw.vcf.tmp ] || rm raw.vcf.tmp

run_in_parallel=1

# cat all input VCFs into raw.vcf.tmp
#if [ "$run_in_parallel" != "" ]
#then
    # save input files into input.XXXX.vcf in parallel
#    echo $input | parallel -j 10 'mv $(./swe fetch {}) input.$$.vcf '
#    cat input.*.vcf > raw.vcf.tmp
#else
    for vcf in $input
    do
	    cat $(./swe fetch $vcf) >> raw.vcf.tmp
    done
#fi


#obtain header from the first VCF in the list
set -- $input
first_vcf=$1

cat $(./swe fetch $first_vcf) |grep -P '^\#' >header.vcf


( cat header.vcf
  cat raw.vcf.tmp | grep -vP "^\#" ) \
      | vcfallelicprimitives --keep-info --keep-geno \
      | vt normalize -r $gatk_data/hg19/ucsc.hg19.fasta - \
      | vcf-sort -c \
      | vcfuniq > raw.vcf

bgzip -f raw.vcf 
tabix -p vcf raw.vcf.gz


./swe emit file raw.vcf.gz
./swe emit file raw.vcf.gz.tbi

exit 0

