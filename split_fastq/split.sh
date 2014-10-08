#!/bin/bash
set -e
set -x
set -o pipefail

splits=$(./swe get splits)

input1=$(./swe get input1 | ./swe fetch -)
input2=$(./swe get input2 | ./swe fetch -)


./split_input_fastq.pl --input paired[$input1,$input2] --splits $splits
# split_input_fastq.pl will output 1.fastq.gz 2.fastq.gz .... N.fastq.gz  - interleaved FASTQ files

for i in $(seq 1 $splits)
do
    ./swe emit file $i.fastq.gz
done

#seq 1 $splits | parallel ./swe emit file {}.fastq.gz


# attempt at dynamically determining splits

#pigz -c $input1 | split -l $splitlines - $input1.
#pigz -c $input2 | split -l $splitlines - $input2.

#rm $input1 $input2
#ls $input1.* $input2.* | parallel -j 16 'gzip {}'

# hack due to problems with swe creating "/./" in s3 URLs
#cd $(dirname $input1) && ../swe emit file *gz && cd -

# save the file list as a workaround for s3 eventual consistency
#ls $(dirname $input1)/*gz >splits.txt && ./swe emit file splits.txt
