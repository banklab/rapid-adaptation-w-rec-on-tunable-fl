#!/usr/bin/perl

use strict;
use warnings;

### check whether the fixed genotype is in the initial population. How much of the proportion?

open INTPOP, $ARGV[0] || die "$!"; ### initial population
#landscape_id   recombination_rate      conf    genotype        count
open FINAL, $ARGV[1] || die "$!"; ### final genotype
#landscape_id   recombination_rate   conf   fixation_time allele

my $popSize=$ARGV[2];

<INTPOP>;
<FINAL>;

while(<INTPOP>){
	chomp;
	my @int=split(/\s+/, $_);
	
	my $tmp=<FINAL>;
	chomp($tmp);
	my @a=split(/\s+/, $tmp);

	my $finalGenotype;
	for(my $i=4; $i<@a; $i+=2){
		$finalGenotype.=$a[$i];
	}

	my $finalGenoFrq=0;

	for(my $g=3; $g<@int; $g+=2){
		if($finalGenotype==$int[$g]){
			$finalGenoFrq=$int[$g+1];
		}
	}

	print $int[0], "\t", $int[1], "\t", $int[2], "\t", $finalGenotype, "\t", $finalGenoFrq, "\n";

}

close INTPOP;
close FINAL;

