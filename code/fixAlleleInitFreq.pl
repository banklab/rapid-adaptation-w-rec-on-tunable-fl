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

print("#landscape_id\trecombination_rate\tconf\talleletype\tfixStep\tfixation_time\tinitialFreq\n");

while(<INTPOP>){
	chomp;
	my @frq=&initFreq($_);

	my $tmp=<FINAL>;
	my @a=split(/\s+/, $tmp);


#### extract fixation time
	my @fixTime;
	for(my $i=3; $i<@a; $i+=2){
		push(@fixTime, $a[$i]);
	}

#### order fixation time
	my @steps=&fixStep(@fixTime);

#### output 
	my @output=@a[0..2];
	for(my $i=4; $i<@a; $i+=2){
		my $j=($i-4)/2;
		# alleletype fixStep
		push(@output, ($a[$i], $steps[$j], $a[$i-1]));
		# initial frequency
		if($a[$i]==0){
			push(@output, 1-$frq[$j]/$popSize);
		}
		else{
			push(@output, $frq[$j]/$popSize);
		}
	}
	print join("\t", @output), "\n";

}

sub initFreq(){
	my $tmp=shift;
	my @a=split(/\s+/, $tmp);

	my @frq;
	for(my $i=3; $i<@a; $i+=2){
		my @b=split("", $a[$i]);
		my $j=$i+1;
		for(my $p=0; $p<@b; $p++){
			$frq[$p]+=$b[$p]*$a[$j];
		}
	}
	
	return(@frq);
}

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


close INTPOP; close FINAL;
