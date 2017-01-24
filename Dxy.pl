#! /usr/bin/env perl
# This program aims for calculating pairwise nucleotide diversity between two populations (Dxy) in sliding window scanning for one chromosome
# input format: vcftools genotype format


use strict;
use warnings;
use Data::Dumper;
use v5.10;
use Getopt::Long;
use List::Util qw(max min);

if(!@ARGV){
	die "Dxy.pl -012 vcftools.GT -pop sample.list -chr scaffold -w window_size -out out_prefix\n";
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
		) or die;
		
# STOUT parameters
say "Input parameters:";
say "Dxy_Df.pl -012 $input -pop $pop -w $w -chr $chr -out $out";

##### main functions #####
# read sample.list
my ($popx, $popy) = read_pop($pop);
# read genotype 
my %GTX;
my %GTY;
my @pos;
read_GT($input);
#window scan
open OUT, '>', "$out.Dxy";
#say OUT "CHR\tBIN_START\tN_SNP\tDxy";
my $start_point=0;
for(my $ind=0; $ind < max(@pos); $ind+=$w){
	# number of SNPs in the window
	my $N_SNP=0;
	foreach (@pos){
		if($_>=($ind+1) && $_ <= ($ind+$w-1)){
			$N_SNP++;
		}
	}


	next if($N_SNP==0);

	#extract the corresponding genotype and calculate
	my $Dxy=0;
	foreach my $kX (keys %GTX){
		my $fragX = substr($GTX{$kX}, $start_point, $N_SNP);
		foreach my $kY (keys %GTY){
			my $fragY = substr($GTY{$kY}, $start_point, $N_SNP);
			$Dxy += cal_Dxy($fragX, $fragY);
		

		}
	}

	$Dxy = $Dxy/(2*scalar(@$popx)*2*scalar(@$popy)*$w);
	say OUT $chr, "\t", $ind, "\t", $N_SNP, "\t", $Dxy;
	$start_point += $N_SNP;

}


##### subroutines #####
close OUT;

# read sample.list into hash
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

# read genotype file
sub read_GT{
	my $in = shift @_;
	open IN, '<', $input or die "No 012 file $input is found!\n";
	my $in_head = <IN>;
	chomp $in_head;
	my @in_head = split(/\t| /, $in_head);

	my (@indexX, @indexY);
	my $i=0;
	while($i<scalar(@in_head)){
		if($in_head[$i] ~~ @$popx){
			push @indexX, $i;
		}elsif($in_head[$i] ~~ @$popy){
			push @indexY, $i;
		}
		$i++;
	}
	
	while(<IN>){
		chomp;
		my @row = split(/\t| /, $_);
		push @pos, $row[1];

		foreach my $j (@indexX){
			$GTX{$j."1"}.=substr($row[$j], 0, 1);
			$GTX{$j."2"}.=substr($row[$j], 2, 1);
		}

		foreach (@indexY){
			$GTY{$_."1"}.=substr($row[$_], 0, 1);
			$GTY{$_."2"}.=substr($row[$_], 2, 1);
		}


	}
	close IN;
}

sub cal_Dxy{
	my ($X, $Y) = @_;
	my $count=0;
	my $empty=0;

	for(my $base=0; $base < length($X); $base++){
		my $baseX=substr($X, $base, 1);
		my $baseY=substr($Y, $base, 1);

		if($baseX eq "\." || $baseY eq "\."){
			$empty++;
			next;
		}else{
			$count += abs($baseX - $baseY);
		}
	}

	return $count;
}