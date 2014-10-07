#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input)

[ ! -e raw.vcf.tmp ] || rm raw.vcf.tmp

run_in_parallel=1

# cat all input VCFs into raw.vcf.tmp
if [ "$run_in_parallel" != "" ]
then
    # save input files into input.XXXX.vcf in parallel
    parallel -j 10 -i bash -c "mv \`./swe fetch {}\` input.\$\$.vcf " -- $input

    cat input.*.vcf > raw.vcf.tmp
else
    for vcf in $input
    do
	cat $(./swe fetch $vcf) >> raw.vcf.tmp
    done
fi


#obtain header from the first VCF in the list
set -- $input
first_vcf=$1

cat $(./swe fetch $first_vcf) |grep -P '^\#' >header.vcf


( cat header.vcf
  cat raw.vcf.tmp | grep -vP "^\#" | ./vcf-sort -c ) | vcfuniq > raw.vcf

bgzip -f raw.vcf 
tabix -p vcf raw.vcf.gz


./swe emit file raw.vcf.gz
./swe emit file raw.vcf.gz.tbi

exit 0

