#! /usr/bin/env perl
# This program aims for calculating number of fixed SNPs (Df) in sliding window scanning the whole genome
# input format: vcftools 012 format

use strict;
use warnings;
use Data::Dumper;
use v5.10;
use Getopt::Long;
use List::Util qw(max min);

if(!@ARGV){
	die "Df_WGS.pl -012 vcftools.012 -pop popuation.list -w window_size -out out_prefix\n";
	}

# variables
my $input;
my $pop;
my $w = 10000;
my $out;
my $repeat;

# Options
GetOptions('012=s' => \$input,
		'pop=s' => \$pop,
		'w=i' => \$w,
		'out=s' => \$out,
		) or die;
		
# STOUT parameters
say "Input parameters:";
say "Df.pl -012 $input -pop $pop -w $w  -out $out";

# read 012.pos file and split into each scaffold
my %pos = read_pos($input);

unlink("$out.Df");
open OUT, '>', "$out.Df" or die;
say OUT "#CHR\tBIN_START\tBIN_END\tN_SNP\tTotal_Df";


foreach my $scaffold (sort keys %pos){
	# remove the scaffold smaller than window size
	next if ( max(keys %{$pos{$scaffold}}) < $w );
	
	open TEMP, '>', "$scaffold.012.pos" or die;
	my @index;
	for my $snp (sort {$a<=>$b} keys %{$pos{$scaffold}}){
		say TEMP $scaffold, "\t", $snp;
		push @index, $pos{$scaffold}{$snp};
	}
	close TEMP;


	my $min = min(@index)+2;
	my $max = max(@index)+2;
	
	system("cut -f 1,$min\-$max $input > $scaffold.012");

	# execute Dxy_Df.pl for each scaffold
	system("cp $input.indv $scaffold.012.indv && perl ./Df.pl -012 $scaffold.012 -pop $pop -chr $scaffold  -w $w -out $scaffold && cat $scaffold.Df >> $out.Df && rm $scaffold*");
}

close OUT;


sub read_pos{
	my $file = shift @_;
	open IN, '<', "$input.pos" or die "No pos file is found!\n";
	my %hash;
	my $n = 0;
	while(<IN>){
		chomp;
		my @row = split(/\t| /);
		$hash{$row[0]}{$row[1]} = $n;
		$n++;
	}
	close IN;

	return %hash;
}
