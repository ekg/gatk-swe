#!/bin/bash
set -e
set -x
set -o pipefail

[ "$input" != "" ]
[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)

[ ! -e raw.vcf.tmp ] || rm raw.vcf.tmp

run_in_parallel=1

# cat all input VCFs into raw.vcf.tmp
#if [ "$run_in_parallel" != "" ]
#then
    # save input files into input.XXXX.vcf in parallel
echo $input | tr ' ' '\n' | parallel -j 16 'cat $(./swe fetch {})' >raw.vcf.tmp
#    cat input.*.vcf > raw.vcf.tmp
#else
#    for vcf in $input
#    do
#	    cat $(./swe fetch $vcf) >> raw.vcf.tmp
#    done
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

