#!/usr/bin/perl
set -e -x 
#Download reference files and unpack it to /tmp/gatk-reference


if [ ! -e /tmp/gatk-reference ]
then
wget -O /tmp/gatk-reference.tar.gz http://gapp-west.s3.amazonaws.com/gatk-reference.tar.gz

cd /tmp && tar xvf gatk-reference.tar.gz

echo Reference files downloaded successfully


fi


# download sample exome

if [ ! -e /tmp/gcat_set_025_2.fastq.gz ]
then
wget	-O /tmp/gcat_set_025_1.fastq.gz http://gapp-west.s3.amazonaws.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_1.fastq.gz
wget	-O /tmp/gcat_set_025_2.fastq.gz http://gapp-west.s3.amazonaws.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_2.fastq.gz

fi