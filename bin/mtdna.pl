#!/usr/bin/env perl
################################################################################
# File Name: mtdna.pl
# Author: samson-xu
# mail: xy_xu@foxmail.com
# Created Time: Mon 22 Jun 2020 02:57:53 PM CST
################################################################################

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use POSIX;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use ConfigParse;
use SampleStat;
use WriteShell;
use MTdna;

# File or tool path check
my $config = path_check("$Bin/config.txt");

# Global variable
my ($help, $stat, $fastqc_help, $fastp_help, $backtrack_help, $mem_help, %mtdna_shell, $main_shell);
my $project = strftime("%Y%m%d-%H%M%S",localtime());
my $workDir = $ENV{'PWD'};
my $ref = $config->{'chrM'};
my $step = 'cfbm';
my $thread = '35';
my $fastqc_arg = '';
my $fastp_arg = "--detect_adapter_for_pe -q 15 -u 40 -n 5 -l 50 -w $thread -d 3";
my $align_way = 'mem';
my $align_arg = '';
my $run = 'no'; 

# Guide
my $guide_separator = "=" x 150;
my $indent = " " x 15;
my $parameter_separator = "*" x 70;
my $guide=<<INFO;
$guide_separator
=
=$indent$indent$indent NAME: Mitochondria SNP/InDel NGS Analysis Pipeline(mtdna) 
=$indent$indent$indent AUTHOR: xuxiangyang(xy_xu\@foxmail.com) 
=$indent$indent$indent VERSION: v0.1	2020/08/31                            
=                                          
$guide_separator


FUNCTIONS
$indent c. Quality control and check(FastQC) of input data(FASTQ).
$indent f. Adapter cut and low quality sequence filter of fastq(Fastp).
$indent b. Fastq Alignment and quality control of sequence alignment results.
$indent m. Mitochondrial gene mutation detection.

PARAMETER
$indent $0 [options] sample.lst

$parameter_separator Basic $parameter_separator 
$indent --help                        Print this guide information 
$indent --project <str>               Project name, default "$project"
$indent --workDir <str>               Work directory, default "$workDir"
$indent --ref <str>                   Reference genome absolute path, default "$ref"
$indent --step <str>                  Set step for run, default "$step"
$indent --thread <i>                  Set the number of threads for the program to run, default "$thread"
$indent --run <str>                   whether run pipeline, yes or no, default "$run"
$indent --stat                        Whether stat sample information, default not stat
$parameter_separator Filter $parameter_separator 
$indent --fastqc_help                 Print fastqc help information
$indent --fastqc_arg                  Fastqc argument setting, default "$fastqc_arg"
$indent --fastp_help                  Print fastp help information
$indent --fastp_arg                   Fastp argument setting, default "$fastp_arg"
$parameter_separator Align $parameter_separator 
$indent --align_way                   Select align algorithm, 'backtrack', 'mem', default "$align_way"
$indent --backtrack_help              Print BWA-backtrack help information
$indent --mem_help                    Print BWA-mem help information
$indent --align_arg                   Align argument setting, this has to correspond to the align_way, default "$align_arg"

NOTE
$indent 1. Fastq quality system should be phred 33
$indent 2. If input is fastq, sample.lst format like: SampleId    fq1    fq2, if input is bam, sample.lst format like: SampleId    bam
$guide_separator

INFO

# Parameter
GetOptions(
	"h|help" => \$help,
	"project=s" => \$project,
	"workDir=s" => \$workDir,
	"ref=s" => \$ref,
	"step=s" => \$step,
	"thread=i" => \$thread,
	"run=s" => \$run,
	"stat" => \$stat,
	"fastqc_help" => \$fastqc_help,
	"fastqc_arg=s" => \$fastqc_arg,
	"fastp_help" => \$fastp_help,
	"fastp_arg=s" => \$fastp_arg,
	"align_way=s" => \$align_way, 
	"backtrack_help" => \$backtrack_help,
	"mem_help" => \$mem_help,
);

if (@ARGV == 0) {
	die `$config->{fastqc} -h` . "\n" if (defined $fastqc_help);
	die `$config->{fastp}` . "\n" if (defined $fastp_help);
	die `$config->{bwa} aln` . "\n" if (defined $backtrack_help);
	die `$config->{bwa} mem` . "\n" if (defined $mem_help);
}
die $guide if (@ARGV == 0 || defined $help);

# Main
my $projectDir = "$workDir/.$project";
my $sample_file = shift;
system("mkdir -p $projectDir") == 0 || die $!;
system("cp $sample_file $projectDir") == 0 || die $!;
# load sample info
my %sampleInfo;
my $fastq_label = `grep 'gz\$' $sample_file`;
open SF, $sample_file or die $!;
while (<SF>) {
	next if (/#/);
	chomp;
	my @arr = split /\s+/;
	if ($fastq_label) {
		@{$sampleInfo{$arr[0]}{'fastq'}} = ($arr[1], $arr[2]);
	} else {
		$sampleInfo{$arr[0]}{'align'} = $arr[1];
	}
}
close SF;
#print Dumper \%sampleInfo;
my $sample_total = keys %sampleInfo;

my ($step1_shell);

if ($fastq_label) {
	foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
		# Fastq quality control
		my $fastq = $sampleInfo{$sampleId}{'fastq'};
		if ($step =~ /c/) {
			my $fastqcDir = "$projectDir/$sampleId/00.fastqc";
			my $fastqc_shell = "$config->{fastqc} -t $thread -o $fastqcDir $fastqc_arg $fastq->[0] $fastq->[1]\n";
			$fastqc_shell .= "rm $fastqcDir/*.zip";
			write_shell($fastqc_shell, "$fastqcDir/$sampleId.fastqc.sh");
			$step1_shell = "sh $fastqcDir/$sampleId.fastqc.sh >$fastqcDir/$sampleId.fastqc.sh.o 2>$fastqcDir/$sampleId.fastqc.sh.e";
		}
		# Fastq filter
		if ($step =~ /f/) {
			my $filterDir = "$projectDir/$sampleId/01.filter";
			my $filter_shell = "$config->{fastp} -i $fastq->[0] -o $filterDir/$sampleId.clean.1.fq.gz -I $fastq->[1] -O $filterDir/$sampleId.clean.2.fq.gz $fastp_arg -j $filterDir/$sampleId.fastq.json -h $filterDir/$sampleId.fastq.html -R '$sampleId fastq report'\n";
			write_shell($filter_shell, "$filterDir/$sampleId.filter.sh");
			$mtdna_shell{$sampleId} .= $step1_shell . " &\n"if ($step =~ /c/);
			$mtdna_shell{$sampleId} .= "sh $filterDir/$sampleId.filter.sh >$filterDir/$sampleId.filter.sh.o 2>$filterDir/$sampleId.filter.sh.e &\n";
			$mtdna_shell{$sampleId} .= "\nwait\n\n";
			@{$sampleInfo{$sampleId}{'clean'}} = ("$filterDir/$sampleId.clean.1.fq.gz", "$filterDir/$sampleId.clean.2.fq.gz");
		}
		# Alignment
		if ($step =~ /b/) {
			my $alignDir = "$projectDir/$sampleId/02.align";
			@{$sampleInfo{$sampleId}{'clean'}} = ($fastq->[0], $fastq->[1]) unless ($sampleInfo{$sampleId}{'clean'}); 
			my $fq1 = $sampleInfo{$sampleId}{'clean'}->[0];
			my $fq2 = $sampleInfo{$sampleId}{'clean'}->[1];
			my $align_shell = "# BWA for alignment\n";
			if ($align_way eq 'backtrack') {
				$align_shell .= "$config->{'bwa'} aln $thread $ref $fq1 > $alignDir/$sampleId.aln_sa1.sai &\n";
				$align_shell .= "$config->{'bwa'} aln $thread $ref $fq2 > $alignDir/$sampleId.aln_sa2.sai &\n";
				$align_shell .= "wait\n";
				$align_shell .= "$config->{'bwa'} sampe -r '\@RG\\tID:$sampleId\\tSM:$sampleId\\tPL:ILLUMINA' $align_arg $ref $alignDir/$sampleId.aln_sa1.sai $alignDir/$sampleId.aln_sa2.sai $fq1 $fq2 | $config->{'samtools'} view -1 --threads $thread - > $alignDir/$sampleId.bam\n\n";
			} else {
				$align_shell .= "$config->{'bwa'} mem -K 100000000 -t $thread -Y -R '\@RG\\tID:$sampleId\\tSM:$sampleId\\tPL:ILLUMINA' $align_arg $ref $fq1 $fq2 | $config->{'samtools'} view -1 --threads $thread - > $alignDir/$sampleId.bam\n";
			}
			$align_shell .= "$config->{'samtools'} sort -T $alignDir/$sampleId --threads $thread -o $alignDir/$sampleId.final.bam $alignDir/$sampleId.bam\n";
			$align_shell .= "$config->{'samtools'} index -@ $thread $alignDir/$sampleId.final.bam $alignDir/$sampleId.final.bai\n";
			$align_shell .= "rm $alignDir/$sampleId.bam\n";
			write_shell($align_shell, "$alignDir/$sampleId.align.sh");
			$mtdna_shell{$sampleId} .= "sh $alignDir/$sampleId.align.sh >$alignDir/$sampleId.align.sh.o 2>$alignDir/$sampleId.align.sh.e\n";
			$sampleInfo{$sampleId}{'align'} = "$alignDir/$sampleId.final.bam"; 
		}	
	}
}

print STDERR "Because your input is alignment result, the program will automatically skip the first few steps!\n" if (!$fastq_label and $step =~ /c|f|b/);	
foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
	$sampleInfo{$sampleId}{'align'} = "$projectDir/$sampleId/02.align/$sampleId.final.bam" unless ($sampleInfo{$sampleId}{'align'});
	# Mitochondrial gene mutation detection
	if ($step =~ /m/) {
		my $mtdnaDir = "$projectDir/$sampleId/03.mtdna";
		mtdna($sampleInfo{$sampleId}{'align'}, $config->{'gatk'}, $config->{'bwa'}, $thread, $config->{'annot'}, $config->{'lib'}, $config->{'db'}, $mtdnaDir);
		$mtdna_shell{$sampleId} .= "bash $mtdnaDir/$sampleId.mtdna.sh >$mtdnaDir/$sampleId.mtdna.sh.o 2>$mtdnaDir/$sampleId.mtdna.sh.e\n";
		$sampleInfo{$sampleId}{'mtdna'} = ""; 
	}
}

$main_shell = "# Run mtdna pipeline for all samples\n";

if ($step =~ /c|f|b|m/) {
	foreach my $sampleId (sort {$a cmp $b} keys %sampleInfo) {
		write_shell($mtdna_shell{$sampleId}, "$projectDir/$sampleId/$sampleId.sh");
		$main_shell .= "sh $projectDir/$sampleId/$sampleId.sh >$projectDir/$sampleId/$sampleId.sh.o 2>$projectDir/$sampleId/$sampleId.sh.e\n";
	}
}

$main_shell .= "cp $projectDir/*/03.mtdna/*.xls $workDir\n";
$main_shell .= "rm -rf $projectDir $workDir/haplogrep.log\n";

write_shell($main_shell, "$projectDir/main.sh");

stat_log($sample_total, $Bin) if (defined $stat);

system("nohup sh $projectDir/main.sh >$projectDir/main.sh.o 2>$projectDir/main.sh.e &") == 0 || die $! if ($run =~ m/y/i);
