#! /usr/bin/env perl
# This program aims for calculating pairwise nucleotide diversity between two populations (Dxy) in sliding window scanning the whole genome
# input format: vcftools 012 format

use strict;
use warnings;
use Data::Dumper;
use v5.10;
use Getopt::Long;
use List::Util qw(max min);

if(!@ARGV){
	die "Dxy_WGS.pl -012 vcftools.012 -pop sample.list -w window_size -out out_prefix\n";
	}

# variables
my $input;
my $pop;
my $w = 10000;
my $out;

# Options
GetOptions('012=s' => \$input,
		'pop=s' => \$pop,
		'w=i' => \$w,
		'out=s' => \$out,
		) or die;
		
# STOUT parameters
say "Input parameters:";
say "Dxy_Df_WGS_new.pl -012 $input -pop $pop -w $w -out $out";

# read genotype file and split into each scaffold
my %pos = read_pos($input);


unlink("$out.Dxy");
open OUT, '>', "$out.Dxy" or die;
say OUT "#CHR\tBIN_START\tN_SNP\tDxy";

my $header;
foreach my $scaffold (sort keys %pos){
	# remove the scaffold smaller than window size
	next if ( max(keys %{$pos{$scaffold}}) < $w );
	
	open TEMP, '>', "$scaffold.GT" or die "Cannot open $scaffold.GT file!\n";
	say TEMP $header;
	my @index;
	for my $snp (sort {$a<=>$b} keys %{$pos{$scaffold}}){

		say TEMP $pos{$scaffold}{$snp};
	}
	close TEMP;

	# execute Dxy_Df_new.pl for each scaffold
	system("perl /home/fanhan/private/scripts/Dxy.pl -012 $scaffold.GT -pop $pop -chr $scaffold -w $w -out $scaffold && cat $scaffold.Dxy >> $out.Dxy && rm $scaffold*");
}

close OUT;


sub read_pos{
	my $file = shift @_;
	open IN, '<', $input or die "No genotype file is found!\n";
	my %hash;
	$header=<IN>;
	chomp $header;
	while(<IN>){
		chomp;
		my @row = split(/\t| /);
		$hash{$row[0]}{$row[1]} = $_;
	}
	close IN;

	return %hash;
}
