#!/bin/bash
set -e
set -x
set -o pipefail

input=$(./swe get input)

[ ! -e raw.vcf.tmp ] || rm raw.vcf.tmp

# replace with parallel
for vcf in $input
do
    local_vcf=$(./swe fetch $vcf)

    cat $local_vcf |grep -P '^\#' >header.vcf
    cat $local_vcf  >> raw.vcf.tmp

done

cat header.vcf > raw.vcf
cat raw.vcf.tmp | grep -vP "^\#" | ./vcf-sort -c >> raw.vcf

bgzip -f raw.vcf 
tabix -p vcf raw.vcf.gz

./swe emit file raw.vcf.gz
./swe emit file raw.vcf.gz.tbi

exit 0

