#!/usr/bin/perl
# GATK pipeline execution script. Runs certain stages of the pipeline.

use Getopt::Long;
use strict;
use File::Basename;

my $cmd = $ARGV[0];

#GetOptions ('sample-id=i', "status-file=s");

#initial sanity check

print "gapp.pl(): ".join(" ",@ARGV)."\n";
#check that script is started from /tmp

$ENV{PATH}="./code/bin/:$ENV{PATH}";

#fetch code
print "Fetching the code:...\n";
unless (-e "code")
{
	die "Codebase not found. aborting";
}

validate_prerequisites();

#run("[ -e code ] || ln -s /home/dmitry/gatk code");

sub get_database_file
{
	die "GATK_DATA bucket is not defined!" unless $ENV{GATK_DATA};
	foreach my $file (@_)
	{
		run ("download-cached ".$ENV{GATK_DATA}."/$file") unless -e basename($file);
	}
}

sub get_gatk_jar {
	my $gatk_jar = shift;

	if ($gatk_jar eq "LITE")
	{
		run("ln -s code/bin/GenomeAnalysisTKLite.jar gatk.jar") unless -e "gatk.jar";
	}
	else
	{
		run ("es3 test $gatk_jar");
		run ("download-cached $gatk_jar");
		run ("mv ".basename($gatk_jar)." gatk.jar");
	}
}


my $cpu_cores=get_cpu_cores();




if ( $cmd eq "interleave")
{
	my $input;
	my $output;
	my $prefix;
	my $pairs;
	GetOptions ( 'input=s'=>\$input, "output=s"=>\$output,"prefix=s"=>\$prefix,"pairs=i"=>\$pairs);
	die "--input must be defined " unless $input;
	die "--output must be defined " unless $output;
	die "--prefix must be defined " unless $prefix;
	die "--pairs must be defined" unless $pairs;
	if ($input =~/^region/)
	{	#todo fix!
	    run("touch $prefix.parts.list");
	    exit(0);
	}
	status ("interleave");

	#print "Sync input file to local directory:\n";
	run("code/1-interleave/split_input_fastq.pl --input $input --prefix $prefix --pairs $pairs");

	# for bam files we output reads in Stabq format so that alignment jobs could take care of
	if ($input =~/^bam/)
	{
		run("ls $prefix"."*.stabq.gz > $prefix.parts.list");
	}
	else
	{
		run("ls $prefix"."*.fastq.gz > $prefix.parts.list");
	}

	status ("sync");

	system("env |grep PATH");
	run("es3 sync $prefix"."* $output/ ");
	print "interleave done\n";
	exit(0);
}

if ( $cmd eq "align")
{
	my $input;
	my $output;
	my $sample_id;
	my $mode;
	my $alignment_algorithm;
	GetOptions ( 	'input=s' =>\$input, 
					"output=s"=>\$output,
					"sample_id=s"=>\$sample_id,
					"mode=s"=>\$mode,
					"alignment_algorithm=s" =>\$alignment_algorithm);
	$alignment_algorithm="BWA-mem" unless $alignment_algorithm;
	die "Incorrect alignment algorithm " unless $alignment_algorithm=~/(BWA-mem|novoalign)/;

	die "--input must be defined " 		unless $input=~/\S+\.(stabq|fastq).gz/;
	die "--output must be defined " 	unless $output;
	die "--sample_id must be defined " 	unless $sample_id;
	die "--mode is incorrect" 			unless $mode=~/^(align_raw|align)$/;
	print "Sync input file to local directory:\n";
	status ("sync");

	run ("es3 sync $input ./");

	if ($input =~/\.stabq.gz/) #if file is in stabq format, shuffle it, and convert to intelreaved fastq
	{
		run(	"zcat ".basename($input). 
				"| sort --parallel $cpu_cores -S 6G ".
				"| code/2-align/stabq2interleaved.pl ".
				"> input.fastq");
		$input="input.fastq";
	}
	#download necessary database files:
	
	my $GROUP_ID='@RG\tID:1\tPL:ILLUMINA\tPU:pu\tLB:group1\tSM:Sample_'.$sample_id;
	
	status ("align");

	if ($alignment_algorithm eq "BWA-mem")
	{
		get_database_file(	"hg19/ucsc.hg19.fasta.amb", 
							"hg19/ucsc.hg19.fasta.ann", 
							"hg19/ucsc.hg19.fasta.bwt",
							"hg19/ucsc.hg19.fasta.pac",
							"hg19/ucsc.hg19.fasta.sa",
							"hg19/ucsc.hg19.fasta",
							"hg19/ucsc.hg19.fasta.fai");

		

		#run alignment 
		run( "bwa mem -M -p -t $cpu_cores -R '$GROUP_ID' ucsc.hg19.fasta ".basename($input).
		 	"| samtools view -\@ $cpu_cores -1 -bt ucsc.hg19.fasta.fai - ".
		 	"| samtools sort -\@ $cpu_cores -l 0 - raw");

	}
	else
	{
		die unless $alignment_algorithm eq "novoalign";
		get_database_file(	
							"hg19/13-3.index",
							"hg19/ucsc.hg19.fasta",
							"hg19/ucsc.hg19.fasta.fai");
		
		run( ( $input =~/\.gz$/ ? "zcat ":"cat ").basename($input)." | code/2-align/interleaved2fastq.pl");
		run("cp code/2-align/novoalign.lic ./");
		run("novoalign -d 13-3.index  -c $cpu_cores -f input_1.fastq input_2.fastq -o SAM '$GROUP_ID' > input.sam");
		run("  samtools view -\@ $cpu_cores -1 -bt ucsc.hg19.fasta.fai input.sam ".
			"| samtools sort -\@ $cpu_cores -l 0 - raw");
		run("rm input_1.fastq input_2.fastq input.sam");

	}

	status ("split by chromosomes");

	run("samtools index raw.bam");

	#align and split by chromosome
	if ($mode eq "align")
	{	
		# get list of chromosomes
		open(FL,"samtools idxstats raw.bam| cut -f 1 |grep chr |") || die "can't get list of chromosomes";
		my @chrs=<FL>; chomp(@chrs);
		die "no chromosomes were found " unless scalar(@chrs);

		#spilt alignment file by chromosomes
		run("parallel -j 3 -i bash -c \"samtools view -\@ $cpu_cores -F 4 -b raw.bam  {}:1-300000000 > {}.bam\" -- ".join(" ",@chrs));
		run("parallel -j 3 -i samtools index {}.bam -- ".join(" ",@chrs));

		run("samtools view -\@ $cpu_cores -f 4 -b raw.bam  > unaligned.bam ");
		run("samtools index unaligned.bam ");
		
		run("rm raw.bam*");
	}
	else
	{
		# just save raw.bam, don't split
	}

	#run("rm *.fastq*"); 
	run("es3 sync *.bam *.bam.bai $output");
	exit(0);
}

if ( $cmd eq "combine")
{
	my $input;
	my $output;
	my $fastq_list;
	my $chr;
	my $rename; #if rename is defined, final bam file will be renamed 
	GetOptions ( 	"input=s" 		=>\$input, 
					"output=s"	  	=>\$output,
					"fastq_list=s"	=>\$fastq_list,
					"chr=s"			=>\$chr,
					"rename=s"		=>\$rename);
	die "--input must be defined " unless $input;
	die "--output must be defined " unless $output;
	die "--fastq_list must be defined " unless $fastq_list;
	die "--chr must be defined" unless $chr;

	status ("download chunks");

	#get list of chunks
	run("es3 sync $fastq_list ./");
	open(FL,basename($fastq_list)) || die "can't open $fastq_list";
	my @chunks=<FL>;
	chomp(@chunks);

	die "no chunks were given as input " unless scalar (@chunks);

	# download chunks for given chromosomes to bams/
	run("[ -e bams ] || mkdir bams");
	run ("es3 sync -I \"*/$chr.bam*\" $input/ ./bams/");

	foreach (@chunks)
	{
		die "$chr is not found in chunk $_ " unless -e "bams/$_/$chr.bam";
	}

	status ("merge");

	#combine chunks into one file for given chromosome
	if (scalar (@chunks) == 1) 	{
		run ("mv bams/".$chunks[0]."/$chr.bam $chr.bam");
	}
	else	{
		run ("samtools merge -\@ $cpu_cores $chr.bam ".join(" ",map {"bams/$_/$chr.bam"} @chunks));
	}

	run("samtools index $chr.bam");


	status ("upload");

	if ($rename)
	{
		die  "--rename must end in .bam" unless $rename =~/\S+\.bam/;
		run( "mv $chr.bam $rename");
		run( "mv $chr.bam.bai $rename.bai");
		run("es3 sync $rename $rename.bai $output");

	}
	else
	{
		run("es3 sync $chr.bam $chr.bam.bai $output");
	}
	run ("rm -Rf bams/");

	exit(0);

}

if ( $cmd eq "bqsr")
{
	my $input;
	my $output;
	my $gatk_jar;
	my $chr;
	GetOptions ( 	"input=s" 		=>\$input, 
					"output=s"	  	=>\$output,
					"gatk-jar=s"	=>\$gatk_jar,
					"chr=s"			=>\$chr);
	die "--input must be defined " unless $input;
	die "--output must be defined " unless $output;
	die "--gatk_jar must be defined " unless $gatk_jar;
	die "--chr must be defined" unless $chr;

	get_gatk_jar ($gatk_jar);
	get_database_file(	"dbsnp_137.hg19.vcf.gz",
						"dbsnp_137.hg19.vcf.gz.tbi",
						"hg19/ucsc.hg19.fasta",
						"hg19/ucsc.hg19.dict",
						"hg19/ucsc.hg19.fasta.fai");

	#create reduced representation of the dbSNP for the chromosome in question
	run("tabix -h dbsnp_137.hg19.vcf.gz $chr:1-300000000 > interval.dbsnp_137.hg19.vcf");
	run("es3 sync $input ./");
	run("mv ".basename($input)." input.bam");
	run("mv ".basename($input).".bai input.bam.bai");
	
	my @cmd =();
	status ("VariantRecalibrator");

	push (@cmd, "java -Xmx6g -jar gatk.jar");
	push (@cmd, "-T BaseRecalibrator ");
	push (@cmd, "-I input.bam ");
	push (@cmd, "-R ucsc.hg19.fasta ");
	push (@cmd, "-knownSites interval.dbsnp_137.hg19.vcf ");
	push (@cmd, "--covariate ContextCovariate ");
	push (@cmd, "--covariate ReadGroupCovariate ");
	push (@cmd, "--covariate QualityScoreCovariate ");
	push (@cmd, "--covariate CycleCovariate ");
	push (@cmd, "-o bqsr.grp ");
	push (@cmd, "--disable_indel_quals ");
	push (@cmd, "-nct $cpu_cores ");

	run(join(" ",@cmd));
	run("es3 sync bqsr.grp $output");
	exit(0);
}

if ( $cmd eq "split")
{
	my $input;
	my $output;
	my $min_interval;
	my $chr;
	GetOptions ( 	"input=s" 		=>\$input, 
					"output=s"	  	=>\$output,
					"min-interval=i"=>\$min_interval,
					"chr=s"			=>\$chr);
	die "--input must be defined " unless $input;
	die "--output must be defined " unless $output;
	die "--min-interval must be defined " unless $min_interval;
	die "--chr must be defined" unless $chr;

	get_database_file(	"hg19/ucsc.hg19.fasta",
						"hg19/ucsc.hg19.fasta.fai");

	status ("download input");

	run ("es3 sync $input ./");
	run ("es3 sync $input.bai ./");

	status ("compute splits");

	#generate equally spaced coarse intervals 1/3 of the size if $interval
	run("code/4-split/equally_spaced_intervals.pl --bam ".basename($input)." --block-length=".int ($min_interval/3)." > 1mb.intervals");

	open(FL,"1mb.intervals")||die;
	my @intervals=<FL>; chomp (@intervals);	
	#within each coarsely spaced interval, generate a list of safe breakpoints
	
	validate_interval(@intervals);
	die "No coarse intervals!!" unless scalar(@intervals);


	# for each coarse interval indentify a few potential safe breakpoints
	run ("parallel -j $cpu_cores -i bash -c \"".
		"samtools view -1 -b ".basename($input)." {} > {}.bam ".
		"&& code/4-split/advanced_splitter.pl --reference ucsc.hg19.fasta --bam {}.bam > {}.breakpoints".
		"&& rm {}.bam ".
		"\" -- ".join (" ",@intervals));

	#combine these breakpoints and sort them	
	run ("parallel -j $cpu_cores -i cat {}.breakpoints -- ".join (" ",@intervals)." |sort -n >breakpoints.lst");
	run ("rm *.breakpoints");
	#construct intervals of at at leatn --min-interval
	run("code/4-split/breakpoints2intervals.pl --bam ".basename($input).
		" --breakpoints breakpoints.lst --interval-size $min_interval >interval.lst");

	# assign NEW intervals to @intervals
	@intervals=();
	open(FL,"interval.lst")||die;
	my @intervals=<FL>; chomp (@intervals);	
	validate_interval(@intervals);

	status ("create subset bams");

	#construct one bam file per interval
	run ("parallel -j 5 -i bash -c \"".
		"samtools view -\@ $cpu_cores -b ".basename($input)." {} > {}.bam && samtools index {}.bam ".
		"\" -- ".join (" ",@intervals));

	status ("upload");

	run ("rm ".basename($input)." ".basename($input).".bai");
	run ("es3 sync *.bam* interval.lst $output");

	exit(0);
}	


if ( $cmd eq "gatk")
{
	my $input;
	my $output;
	my $gatk_jar;
	my $bqsr;
	my $interval;

	GetOptions ( 	"input=s" 		=>\$input, 
					"output=s"	  	=>\$output,
					"gatk-jar=s"	=>\$gatk_jar,
					"interval=s"	=>\$interval,
					"bqsr=s"		=>\$bqsr);
	die "--input must be defined " unless $input;
	die "--output must be defined " unless $output;
	die "--gatk_jar must be defined " unless $gatk_jar;
	die "--interval must be defined" unless $interval;
	die "--bqsr  must be defined" unless $bqsr;

	validate_interval($interval);
	get_gatk_jar ($gatk_jar);

	get_database_file(	"dbsnp_137.hg19.vcf.gz",
						"dbsnp_137.hg19.vcf.gz.tbi",
						"1000G_phase1.indels.hg19.vcf.gz",
						"1000G_phase1.indels.hg19.vcf.gz.tbi",
						"Mills_and_1000G_gold_standard.indels.hg19.vcf.gz",
						"Mills_and_1000G_gold_standard.indels.hg19.vcf.gz.tbi",
						"hg19/ucsc.hg19.fasta",
						"hg19/ucsc.hg19.dict",
						"hg19/ucsc.hg19.fasta.fai");

	run ("tabix -h Mills_and_1000G_gold_standard.indels.hg19.vcf.gz $interval > VCF1.vcf");
	run ("tabix -h 1000G_phase1.indels.hg19.vcf.gz $interval                  > VCF2.vcf");
	run ("tabix -h dbsnp_137.hg19.vcf.gz $interval         > interval.dbsnp_137.hg19.vcf");

	die "Incorrect input format ($input) " unless $input =~/\S+\.bam$/;
	die "Incorrect bqsr format ($bqsr) " unless $bqsr =~/\S+\.grp$/;

	run("es3 sync $input ./");
	run("es3 sync $bqsr  ./");

	#run deduplicaion 
	status ("MarkDuplicates");

	run ("java -Xmx2g -jar code/5-gatk/MarkDuplicates.jar INPUT=".basename($input).
		" OUTPUT=deduplicated.bam REMOVE_DUPLICATES=true METRICS_FILE=duplication.metrics");
	run ("samtools index deduplicated.bam");

	my @cmd=();
	push(@cmd, "java -Xmx2g -jar gatk.jar");
	push(@cmd, "-I deduplicated.bam");
	push(@cmd, "-R ucsc.hg19.fasta");
	push(@cmd, "-T RealignerTargetCreator");
	push(@cmd, "-o realigned.intervals");
	push(@cmd, "--known VCF1.vcf");
	push(@cmd, "--known VCF2.vcf");
	push(@cmd, "-L $interval");
	run(join(" ",@cmd));

	status ("IndelRealigner");

	@cmd=();
	push(@cmd, "java -Xmx2g -jar gatk.jar");
	push(@cmd, "-I deduplicated.bam");
	push(@cmd, "-R ucsc.hg19.fasta");
	push(@cmd, "-T IndelRealigner");
	push(@cmd, "-targetIntervals realigned.intervals");
	push(@cmd, "-o realigned.bam ");
	push(@cmd, "-known VCF1.vcf");
	push(@cmd, "-known VCF2.vcf");
	push(@cmd, "--consensusDeterminationModel KNOWNS_ONLY");
	push(@cmd, "-L $interval ");
	push(@cmd, "-LOD 0.4 ");
	push(@cmd, "--downsample_to_coverage 250 ");
	push(@cmd, "-compress 0");
	run(join(" ",@cmd));

	status ("Apply BQSR");

	@cmd=();
	push(@cmd, "java -Xmx2g -jar gatk.jar ");
	push(@cmd, "-R ucsc.hg19.fasta ");
	push(@cmd, "-I realigned.bam ");
	push(@cmd, "-T PrintReads ");
	push(@cmd, "-o recalibrated.bam ");
	push(@cmd, "--disable_indel_quals ");
	push(@cmd, "-BQSR ".basename($bqsr));
	push(@cmd, "-compress 0");
	run(join(" ",@cmd));


	status ("Variant Calling");

	if ($gatk_jar eq "LITE") #run Unified Genotyper for GATK-Lite
	{

		@cmd=();
		push(@cmd, "java -Xmx6g -jar gatk.jar ");
	 	push(@cmd, "-R ucsc.hg19.fasta ");
	 	push(@cmd, "-T UnifiedGenotyper ");
	 	push(@cmd, "-I recalibrated.bam ");
	 	push(@cmd, "--dbsnp interval.dbsnp_137.hg19.vcf ");
	 	push(@cmd, "-o raw.vcf ");
	 	push(@cmd, "-glm BOTH ");
	 	push(@cmd, "-L $interval ");
	 	push(@cmd, "-stand_call_conf 30.0 ");
	 	push(@cmd, "-stand_emit_conf 30.0 ");
	 	push(@cmd, "-nct $cpu_cores ");
	  	push(@cmd, "-rf BadCigar");
		run(join(" ",@cmd));
	}
	else #full featured GATK, run HaplotypeCaller
	{
		@cmd=();
		push(@cmd,"java -Xmx6g -jar gatk.jar ");
	 	push(@cmd,"-R ucsc.hg19.fasta ");
	 	push(@cmd,"-T HaplotypeCaller ");
	  	push(@cmd,"-I recalibrated.bam ");
	   	push(@cmd,"--dbsnp interval.dbsnp_137.hg19.vcf ");
	   	push(@cmd,"-stand_call_conf 30.0 ");
	   	push(@cmd,"-stand_emit_conf 30.0 ");
	   	push(@cmd,"-o raw.hc.vcf ");
	   	push(@cmd,"-L $interval ");
	   	push(@cmd,"-nct $cpu_cores ");
	   	push(@cmd,"-rf BadCigar");
		run(join(" ",@cmd));

		@cmd=();
		push(@cmd,"java -Xmx6g -jar gatk.jar ");
	   	push(@cmd,"-R ucsc.hg19.fasta ");
	   	push(@cmd,"-T VariantAnnotator ");
	   	push(@cmd,"-I recalibrated.bam ");
	   	push(@cmd,"--variant raw.hc.vcf ");
	   	push(@cmd,"-A MappingQualityZero ");
	   	push(@cmd,"--dbsnp interval.dbsnp_137.hg19.vcf ");
	   	push(@cmd,"-o raw.vcf ");
	   	push(@cmd,"-L $interval ");
	   	push(@cmd,"-nt $cpu_cores ");
	   	push(@cmd,"-rf BadCigar");
		run(join(" ",@cmd));
	}
	status("output");

	run ("es3 sync raw.vcf $output");

	exit(0);
}


if ( $cmd eq "freebayes")
{
	my $input;
	my $output;
	my $gatk_jar;
	my $bqsr;
	my $interval;

	GetOptions ( 	"input=s" 		=>\$input, 
					"output=s"	  	=>\$output,
					"gatk-jar=s"	=>\$gatk_jar,
					"interval=s"	=>\$interval,
					"bqsr=s"		=>\$bqsr);
	die "--input must be defined " unless $input;
	die "--output must be defined " unless $output;
	die "--gatk_jar must be defined " unless $gatk_jar;
	die "--interval must be defined" unless $interval;
	die "--bqsr  must be defined" unless $bqsr;

	validate_interval($interval);
	get_gatk_jar ($gatk_jar);

	get_database_file(	"dbsnp_137.hg19.vcf.gz",
						"dbsnp_137.hg19.vcf.gz.tbi",
						"1000G_phase1.indels.hg19.vcf.gz",
						"1000G_phase1.indels.hg19.vcf.gz.tbi",
						"Mills_and_1000G_gold_standard.indels.hg19.vcf.gz",
						"Mills_and_1000G_gold_standard.indels.hg19.vcf.gz.tbi",
						"hg19/ucsc.hg19.fasta",
						"hg19/ucsc.hg19.dict",
						"hg19/ucsc.hg19.fasta.fai");

	run ("tabix -h Mills_and_1000G_gold_standard.indels.hg19.vcf.gz $interval > VCF1.vcf");
	run ("tabix -h 1000G_phase1.indels.hg19.vcf.gz $interval                  > VCF2.vcf");
	run ("tabix -h dbsnp_137.hg19.vcf.gz $interval         > interval.dbsnp_137.hg19.vcf");

	die "Incorrect input format ($input) " unless $input =~/\S+\.bam$/;
	die "Incorrect bqsr format ($bqsr) " unless $bqsr =~/\S+\.grp$/;

	run("es3 sync $input ./");
	run("es3 sync $bqsr  ./");

	#run deduplicaion 
	status ("MarkDuplicates");

	run ("java -Xmx2g -jar code/5-gatk/MarkDuplicates.jar INPUT=".basename($input).
		" OUTPUT=deduplicated.bam REMOVE_DUPLICATES=true METRICS_FILE=duplication.metrics");
	run ("samtools index deduplicated.bam");

	my @cmd=();
	push(@cmd, "java -Xmx2g -jar gatk.jar");
	push(@cmd, "-I deduplicated.bam");
	push(@cmd, "-R ucsc.hg19.fasta");
	push(@cmd, "-T RealignerTargetCreator");
	push(@cmd, "-o realigned.intervals");
	push(@cmd, "--known VCF1.vcf");
	push(@cmd, "--known VCF2.vcf");
	push(@cmd, "-L $interval");
	run(join(" ",@cmd));

	status ("IndelRealigner");

	@cmd=();
	push(@cmd, "java -Xmx2g -jar gatk.jar");
	push(@cmd, "-I deduplicated.bam");
	push(@cmd, "-R ucsc.hg19.fasta");
	push(@cmd, "-T IndelRealigner");
	push(@cmd, "-targetIntervals realigned.intervals");
	push(@cmd, "-o realigned.bam ");
	push(@cmd, "-known VCF1.vcf");
	push(@cmd, "-known VCF2.vcf");
	push(@cmd, "--consensusDeterminationModel KNOWNS_ONLY");
	push(@cmd, "-L $interval ");
	push(@cmd, "-LOD 0.4 ");
	push(@cmd, "--downsample_to_coverage 250 ");
	push(@cmd, "-compress 0");
	run(join(" ",@cmd));

	status ("Apply BQSR");

	@cmd=();
	push(@cmd, "java -Xmx2g -jar gatk.jar ");
	push(@cmd, "-R ucsc.hg19.fasta ");
	push(@cmd, "-I realigned.bam ");
	push(@cmd, "-T PrintReads ");
	push(@cmd, "-o recalibrated.bam ");
	push(@cmd, "--disable_indel_quals ");
	push(@cmd, "-BQSR ".basename($bqsr));
	push(@cmd, "-compress 0");
	run(join(" ",@cmd));


	status ("Variant Calling");

	run ("freebayes -f ucsc.hg19.fasta recalibrated.bam > fb.vcf ");

	status ("Variant Annotations");

		@cmd=();
		push(@cmd,"java -Xmx6g -jar gatk.jar ");
	   	push(@cmd,"-R ucsc.hg19.fasta ");
	   	push(@cmd,"-T VariantAnnotator ");
	   	push(@cmd,"-I recalibrated.bam ");
	   	push(@cmd,"--variant fb.vcf ");
	   	push(@cmd,"-A MappingQualityZero ");
	   	push(@cmd,"-A QualByDepth ");
 		push(@cmd,"-A MappingQualityRankSumTest ");
	   	push(@cmd,"-A FisherStrand ");
	   	push(@cmd,"-A ReadPosRankSumTest ");
	   	push(@cmd,"--dbsnp interval.dbsnp_137.hg19.vcf ");
	   	push(@cmd,"-o raw.vcf ");
	   	push(@cmd,"-L $interval ");
	   	push(@cmd,"-nt $cpu_cores ");
	   	push(@cmd,"-rf BadCigar");
		run(join(" ",@cmd));

	status("output");

	run ("es3 sync raw.vcf $output");

	exit(0);
}


if ( $cmd eq "calibrate")
{
	my $input;
	my $output;
	my $gatk_jar;
	my $exome;
	my $bed_filter;
	my $rename; #will change name to $rename if specified
	GetOptions ( 	"input=s" 		=>\$input, 
					"output=s"	=>\$output,
					"gatk-jar=s"	=>\$gatk_jar,
					"exome"		=>\$exome,
					"rename=s"	=>\$rename,
					"bed-filter=s"	=>\$bed_filter);

	die "--input must be defined " unless $input;
	die "--output must be defined " unless $output;
	die "--gatk_jar must be defined " unless $gatk_jar;

	run("es3 test $bed_filter");
	get_gatk_jar ($gatk_jar);

	get_database_file(	"dbsnp_137.hg19.vcf.gz",
						"dbsnp_137.hg19.vcf.gz.tbi",
						"Mills_and_1000G_gold_standard.indels.hg19.vcf.gz",
						"Mills_and_1000G_gold_standard.indels.hg19.vcf.gz.tbi",
						"1000G_omni2.5.hg19.vcf.gz",
						"1000G_omni2.5.hg19.vcf.gz.tbi",
						"hapmap_3.3.hg19.vcf.gz",
						"hapmap_3.3.hg19.vcf.gz.tbi",
						"hg19/ucsc.hg19.fasta",
						"hg19/ucsc.hg19.dict",
						"hg19/ucsc.hg19.fasta.fai");

	status ("download vcf");

	run ("es3 sync $input ./");
	run ("zcat ".basename($input)." >unsorted.vcf");
	run ("vcfsorter.pl ucsc.hg19.dict unsorted.vcf > input_pre_noise.vcf");
	
	if (0) { #somewhat controversial feature
		run ("code/6-calibrate/add_noise.pl QD MQRankSum FS MQ0 DP ReadPosRankSum < input_pre_noise.vcf >input.vcf ");
	}
	else{
		run ("mv input_pre_noise.vcf input.vcf");
	}


	############################### SNPS ###########################################
	status ("SNP recalibration");

	print "SNP recalibration:\n";
	run("code/6-calibrate/select_snps.pl < input.vcf > snps.raw.vcf");
	my $variants = `grep -vPc "\^\#" < snps.raw.vcf`; chomp($variants);
	print "raw.vcf: Total $variants SNPS found\n";
	die "Not enough variants to run recalibration" if $variants < 10000;

	my $min_variants;
	if ($exome)	{
		$min_variants = 3000;
	} 
	else {	
		$min_variants = int($variants/50);	
	}

	foreach my $extra_options ("" , " --maxGaussians 4 "," --maxGaussians 2 ")
	{
		next if -e "snps.recal";
		my @cmd=();
		push(@cmd," java -Xmx7g -jar gatk.jar ");
	    push(@cmd,"-T VariantRecalibrator ");
	    push(@cmd,"-R ucsc.hg19.fasta ");
	    push(@cmd,"-input snps.raw.vcf ");
	    if ($gatk_jar eq "LITE") {
	    	 push(@cmd,"--minNumBadVariants $min_variants ");
	    }
	    else {
	    	 push(@cmd,"--numBadVariants $min_variants ");
	    }
	    push(@cmd,"-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg19.vcf.gz ");
	    push(@cmd,"-resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg19.vcf.gz ");
	    push(@cmd,"-resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 dbsnp_137.hg19.vcf.gz ");
	    if ($exome) { #do not use depth-based annotations for exome
	    	push(@cmd,"-an QD -an MQRankSum -an FS ");
	    }
	    else {
	    	push(@cmd,"-an MQ0 -an QD -an MQRankSum -an FS -an DP");
	    }
	    push(@cmd,"-mode SNP ");
	    push(@cmd,"-recalFile snps.recal.tmp ");
	    push(@cmd,$extra_options );
	    push(@cmd,"-tranchesFile snps.tranches ");
	    push(@cmd,"-rscriptFile snps.plots.R");
	    push(@cmd,"-nt $cpu_cores");
		
		if (run(join(" ",@cmd),1))
		{
			run ("mv snps.recal.tmp snps.recal");
		}
	}

	my @cmd=();
	push(@cmd,"java  -Xmx7g  -jar gatk.jar ");
	push(@cmd,"-T ApplyRecalibration ");
	push(@cmd,"-R ucsc.hg19.fasta ");
	push(@cmd,"-input snps.raw.vcf");
	push(@cmd,"--ts_filter_level 99.0 ");
	push(@cmd,"-tranchesFile snps.tranches ");
	push(@cmd,"-recalFile snps.recal ");
	push(@cmd,"-mode SNP ");
	push(@cmd,"-o snps.recalibrated.filtered.vcf ");
	run(join(" ",@cmd));

	status ("Indel recalibration");

	############################################ Indels
	print "INDEL recalibration:\n";
	run("code/6-calibrate/select_indels.pl < input.vcf > indels.raw.vcf");
	my $variants = `grep -vPc "\^\#" < indels.raw.vcf`; chomp($variants);
	print "raw.vcf: Total $variants indels found\n";
	die "Not enough variants to run recalibration" if $variants < 3000;

	my $min_variants;
	if ($exome)	{
		$min_variants = 3000;
	} 
	else {	
		$min_variants = int($variants/50);	
	}

	foreach my $extra_options ("" , " --maxGaussians 4 "," --maxGaussians 2 ")
	{
		next if -e "indels.recal";
		my @cmd=();
		push(@cmd,"java -Xmx7g  -jar gatk.jar ");
	   	push(@cmd,"-T VariantRecalibrator ");
	   	push(@cmd,"-R ucsc.hg19.fasta ");
	   	push(@cmd,"-input indels.raw.vcf ");
	    if ($gatk_jar eq "LITE") {
	    	 push(@cmd,"--minNumBadVariants $min_variants ");
	    }
	    else {
	    	 push(@cmd,"--numBadVariants $min_variants ");
	    }
	   	push(@cmd,"-resource:mills,VCF,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg19.vcf.gz ");
	   	push(@cmd,"-resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 dbsnp_137.hg19.vcf.gz ");
	    if ($exome) { #do not use depth-based annotations for exome
	    	push(@cmd,"-an QD -an MQRankSum -an FS -an ReadPosRankSum");
	    }
	    else {
	    	push(@cmd,"-an QD -an MQRankSum -an FS -an DP -an ReadPosRankSum");
	    }
	   	push(@cmd,"-mode INDEL ");
	   	push(@cmd,"-recalFile indels.recal.tmp ");
	   	push(@cmd, $extra_options );
	   	push(@cmd,"-tranchesFile indels.tranches ");
	   	push(@cmd,"-rscriptFile indels.plots.R ");
	   	push(@cmd,"-nt $cpu_cores");

		if (run(join(" ",@cmd),1))
		{
			run ("mv indels.recal.tmp indels.recal");
		}
	}
	
	status ("Apply recalibration");

	my @cmd=();
	push(@cmd,"java -Xmx7g -jar gatk.jar ");
	push(@cmd,"   -T ApplyRecalibration ");
	push(@cmd,"   -R ucsc.hg19.fasta ");
	push(@cmd,"   -input indels.raw.vcf  ");
	push(@cmd,"   --ts_filter_level 95.0 ");
	push(@cmd,"   -tranchesFile indels.tranches ");
	push(@cmd,"   -recalFile    indels.recal");
	push(@cmd,"   -mode INDEL ");
	push(@cmd,"   -o indels.recalibrated.filtered.vcf");
	run(join(" ",@cmd));

	################################### Combine variants 

	status ("combine variants");

	my @cmd=();
	push(@cmd,"java -Xmx7g -jar gatk.jar ");
	push(@cmd,"-T CombineVariants  ");
	push(@cmd,"-R ucsc.hg19.fasta  ");
	push(@cmd,"--variant snps.recalibrated.filtered.vcf  ");
	push(@cmd,"--variant indels.recalibrated.filtered.vcf  ");
	push(@cmd,"-o merged.recalibrated.filtered.vcf ");
	run(join(" ",@cmd));


	#apply BED region filter that defined capture regions
	if ($bed_filter)
	{
		run("es3 sync $bed_filter ./");
		run("rm merged.recalibrated.filtered.vcf.idx");
		run("mv merged.recalibrated.filtered.vcf pre-filtering.vcf");
		run("code/6-calibrate/apply_bed_filter.pl ".basename($bed_filter)."< pre-filtering.vcf > merged.recalibrated.filtered.vcf");
	}

	run("bgzip merged.recalibrated.filtered.vcf");
	run("tabix -p vcf merged.recalibrated.filtered.vcf.gz");

	if ($rename)
	{
		die  "--rename must end in .vcf.gz" unless $rename =~/\S+\.vcf.gz/;
		run( "mv merged.recalibrated.filtered.vcf.gz $rename");
		run( "mv merged.recalibrated.filtered.vcf.gz.tbi $rename.tbi");
		run("es3 sync $rename* $output");
	}
	else
	{
		run("es3 sync merged.recalibrated.filtered.vcf* $output");
	}
	exit(0);
}

die "Unsupported command '$cmd'";


















################################################# subrotines

sub validate_interval
{
	foreach (@_)
	{
		die "Incorrect interval $_" unless /^(\S+):\d+\-\d+$/;
	}

}
sub validate_prerequisites 
{
	print "Validating validating prerequisites:";
	foreach (	"pigz --help", 
				"zcat --help",
				"basename --help", 
				"realpath --help",
				"parallel -i touch {}.test -- 1 2 3",
				"dirname --help",
				"which samtools",
				"rm 1.test 2.test 3.test" )
	{
		run ("$_") if system("$_ >/dev/null 2>/dev/null");
	}
	
	print " OK\n";

}

#get number of available cpu cores
sub get_cpu_cores
{
	my $cores;
	open(FL,"cat /proc/cpuinfo|")||die;
	while (<FL>) {$cores ++ if /^(P|p)rocessor/;}
	die "Can't measure number of available CPU cores" unless $cores;
	return $cores;

}

# send job status to the scheduler
sub status 
{
	my $msg=shift;
	print STDERR "\nreporter:status:$msg\n"

}

#Run shell command. If ignore errors enabled, will return success if command succeeded
sub run
{
	my ($cmd, $ignore_errors)=@_;
	my $time = localtime;
	print "run(): $cmd\n";

	my $ret=system("bash","-o","pipefail","-c",$cmd);
	if ($ret != 0)
	{
		print "run():Command ($cmd) failed with exit code $ret ....";
		if ($ignore_errors)
		{
			print "Error ignored !!!\n";
		}
		else
		{
			print "\n";
			exit(1);
		}
	}
	return ($ret == 0);
}

