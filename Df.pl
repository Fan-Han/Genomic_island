#! /usr/bin/env perl
# This program aims for calculating number of fixed SNPs (Df) in sliding window scanning for one chromosome
# input format: vcftools 012 format


use strict;
use warnings;
use Data::Dumper;
use v5.10;
use Getopt::Long;
use List::Util qw(max min);

if(!@ARGV){
	die "Df.pl -012 vcftools.012 -pop popuation.list -chr scaffold -w window_size -out out_prefix\n";
	}

# variables
my $input;
my $pop;
my $w = 10000;
my $out;
my $chr;
my $repeat;

# Options
GetOptions('012=s' => \$input,
		'pop=s' => \$pop,
		'w=i' => \$w,
		'out=s' => \$out,
		'chr=s' => \$chr,
		'repeat=s' => \$repeat,
		) or die;
		
# STOUT parameters
say "Input parameters:";
say "Df.pl -012 $input -pop $pop -w $w -chr $chr -out $out";

##### main functions #####
# read -012 and write into a hash
my %SNP = read_input($input);
my %pos = read_pos($input);
# read pop x and y
my ($popx, $popy) = read_pop($pop);
# window scanning
open OUT, '>', "$out.Df" or die;
window_scan();
close OUT;

##### subroutines #####
# read -012 (file.012, file.012.indv, file.012.pos)
sub read_input{
	my $name = shift @_;
	# read file.012.indv
	my @indv;
	open IN, '<', "$name.indv" or die "No indv file is found!\n";
	while(<IN>){
		chomp;
		push @indv, $_;
	}
	close IN;
	
	# read file.012
	my %hash;
	open IN, '<', $name or die "No input file is found!\n";
	while(<IN>){
		chomp;
		my @row = split(/\t| /, $_, 2);
		my $ID = shift @indv;
		$row[1] =~ s/\-1/N/g;
		$row[1] =~ s/\t| //g;
		$hash{$ID} = $row[1];
	}
	close IN;
	
	return %hash;	
	
}

sub read_pos{
	my $file = shift @_;
	open IN, '<', "$input.pos" or die "No pos file is found!\n";
	my %hash;
	my $n=0;
	while(<IN>){
		chomp;
		my @row = split(/\t| /);
		$hash{$row[1]} = $n;
		$n++;
	}
	close IN;
	
	return %hash;
}

# read population list $pop
##SAMPLE	POPULATION
sub read_pop{
	my $list = shift @_;
	my (@x, @y);
	open IN, '<', $list or die "No pop file is found!\n";
	my %hash;
	while(<IN>){
		chomp;
		my @row = split(/\t| /);
		push @{$hash{$row[1]}}, $row[0]; 
	}
	close IN;
	
	my $keyx = (keys %hash)[0];
	my $keyy = (keys %hash)[1];

	@x = @{$hash{$keyx}};
	@y = @{$hash{$keyy}};
	
	return (\@x, \@y);
}


# window scanning for a scaffold/chromosome
sub window_scan{
	for(my $i = 1; $i <= max(keys %pos); $i += $w ){
		my ($Df, $rep_count) = (0,0);

		my @pos_index;

		foreach my $p (sort {$a<=>$b} keys %pos){
			if ( $p >= $i && $p < ($i + $w) ){
				push @pos_index, $pos{$p};
			}
		}
		
		# Df per site in this window
		foreach (@pos_index){
			my ($x_pop_seq, $y_pop_seq);
			foreach my $sx (@$popx){
				$x_pop_seq .= substr($SNP{$sx}, $_, 1);
			}
			foreach my $sy (@$popy){
				$y_pop_seq .= substr($SNP{$sy}, $_, 1);
			}
			#say $x_pop_seq, "\t", $y_pop_seq;
			my $Df_1 =  Df_boolen($x_pop_seq, $y_pop_seq);
			$Df += $Df_1;
		}

		
		say  OUT $chr, "\t", $i, "\t", $i+$w-1, "\t", scalar(@pos_index), "\t", $Df;


	}
}

# calculation
sub segregating{
	my ($X, $Y) = @_;
	my $count = 0;

	for(0..(length($X)-1)){
		my $a = substr($X, $_, 1);
		my $b = substr($Y, $_, 1);
		if($a ne "N" && $b ne "N" && $a ne $b){
			$count += abs($a - $b);
		}

	}
	return $count;
}

sub Df_boolen{
	my ($X, $Y) = @_;
	my $count=0;

	if( ($X eq (0 x length($X)) && $Y eq (2 x length($Y))) || ($Y eq (0 x length($Y)) && $X eq (2 x length($X))) ){
		$count=1;
	}
	return $count;

}
