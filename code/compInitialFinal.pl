#!/usr/bin/perl

use strict;
use warnings;

open INTPOP, $ARGV[0] || die "$!"; ### initial population
#landscape_id   recombination_rate      conf    genotype        count
open FINAL, $ARGV[1] || die "$!"; ### final genotype
#landscape_id   recombination_rate   conf   fixation_time allele
my $popSize=$ARGV[2];

<INTPOP>;
<FINAL>;

print("#landscape_id\trecombination_rate\tconf\tdistance\n");

while(<INTPOP>){
	chomp;
	my @a=split;
	my $tmp=<FINAL>;
	my @b=split(/\s+/, $tmp);
	


	my $final="";

	for(my $i=4; $i<@b; $i+=2){
		$final=$final.$b[$i];
	}
	

	my $dist=0;
	for(my $i=3; $i<@a; $i+=2){
		my $d=&compGenotype($final, $a[$i]);
		$dist+=$d*$a[$i+1];
	}

	$dist=$dist/$popSize;

	print join("\t", (@a[0..2], $dist)), "\n";
}

close INTPOP; close FINAL;
sub compGenotype(){
	my $s1=shift;
	my $s2=shift;
	my $l= length $s1;
	my $diff=0;

	for(my $i=0; $i<$l; $i++){
		my $c1 = substr $s1, $i, 1;
		my $c2 = substr $s2, $i, 1;
		if($c1 ne $c2) {
			$diff++;
		}
	}
	return $diff/$l;
}

