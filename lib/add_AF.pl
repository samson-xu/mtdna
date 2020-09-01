#!/usr/bin/env perl
use strict;
use warnings;

my ($vcf, $multianno) = @ARGV;
my %hash;
my $line_number = 1;
$hash{$line_number} = "AF";

open VCF, $vcf or die $!;
while (<VCF>) {
	chomp;
	$line_number++;
	next if(/^#/);
	next if(/^\s*$/);
	my @arr = split /\t/;
	my @arr2 = split /:/, $arr[-1];
	my @arr3 = split /,/, $arr2[2];
	$hash{$line_number} = "$arr3[0]";
}
close VCF;

my $line_count = 0;
open ANNO, $multianno or die $!;
while (<ANNO>) {
	chomp;
	$line_count++;
	next if(/^#/);	
	next if(/^\s*$/);
	my @arr = split /\t/;
	my @new_arr = (@arr[0..4], $hash{$line_count}, @arr[5..$#arr]);	
	my $new_line = join "\t", @new_arr;
	print "$new_line\n" if($hash{$line_count} >= 0.001 or $hash{$line_count} eq "AF");
}
close ANNO;
