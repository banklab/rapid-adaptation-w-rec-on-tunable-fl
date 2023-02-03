#!/usr/bin/perl

use strict;
use warnings;

my $fldir=$ARGV[0];
my $fixdir=$ARGV[1];
my $prefix=$ARGV[2]; ### parameter description
my $numFL=$ARGV[3];
my $flprefix=$prefix; $flprefix=~s/N[0-9]+_//g;

open FIX, "$fixdir/rmf_$prefix\_.dat" || die "$!";

#print("#landscape_id\trecombination_rate\tconf\tepistasisStep\n");

## number of loci;
my $L=$prefix;
$L=~s/L([0-9]+)_//g;
$L=$1;

my $tmp=<FIX>;

for(my $flid=0; $flid<$numFL; $flid++){
	open FL, "$fldir/rmf_landscape_$flprefix\_$flid.dat";
	### record the fitness landscape as hash (?). This might not the best way. 
	my %fl;
	while(<FL>){
		chomp;
		my @a=split;
		$fl{join("", @a[0..$L-1])}=$a[$L];
	}
	close FL;

	### we have five recombination rates now, 0, 0.001, 0.01. 0.1, 0.5. Each 100 populations at one fitness landscape. In the future, those code should be change for convinience. 
	for(my $conf=0; $conf<500; $conf++){
		my $tmp=<FIX>;
		my @a=split(/\s+/, $tmp);

		### get the final genotype
		my $final;
		for(my $i=4; $i<@a; $i+=2){
			$final.=$a[$i];
		}
		
		### fixation steps
		my @fixTime;

		for(my $i=3; $i<@a; $i+=2){
			push(@fixTime, $a[$i]);
		}
		my @steps=&fixStep(@fixTime);
		

		my @oneLoc;

		for(my $s=0; $s<$L; $s++){
			my $c1=substr($final, $steps[$s]-1, 1);
			my $a_c1=1-$c1;

			my $geno1=$final; substr($geno1, $steps[$s]-1, 1)=$a_c1;
			
			$oneLoc[$s]=log($fl{$final})-log($fl{$geno1});
			
		}
		print join("\t",(@a[0..2], @oneLoc)), "\n";
	}
}

close FIX;

sub fixStep(){
	my @fixTime=@_;

	my $L=@fixTime;
	my @order=1..@fixTime;
	for(my $i=0; $i<$L-1; $i++){
		for(my $j=$i+1; $j<$L; $j++){
			if($fixTime[$i]>$fixTime[$j]){
				my $tmp=$fixTime[$i];
				$fixTime[$i]=$fixTime[$j];
				$fixTime[$j]=$tmp;

				my $tmp1=$order[$i];
				$order[$i]=$order[$j];
				$order[$j]=$tmp1;
			}
		}
	}
	return(@order);
}











