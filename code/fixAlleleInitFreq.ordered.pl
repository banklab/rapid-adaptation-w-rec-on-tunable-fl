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

#print("#landscape_id\trecombination_rate\tconf\talleletype\tfixStep\tfixation_time\tinitialFreq\n");

while(<INTPOP>){
	chomp;
	my @frq=&initFreq($_);

	my $tmp=<FINAL>;
	my @a=split(/\s+/, $tmp);


#### extract fixation time
	my @fixTime;
	my @frqStep;
	my @fixTstep;
	my @alleleStep;
	for(my $i=3; $i<@a; $i+=2){
		push(@fixTime, $a[$i]);
		push(@frqStep, 0);
		push(@fixTstep, 0);
		push(@alleleStep, 0);
	}

#### order fixation time;
	my @steps=&fixStep(@fixTime);

#### output allele type
	for(my $i=4; $i<@a; $i+=2){
		my $j=($i-4)/2;
		# initial frequency
		my $tmpFrq=0;
		
		if($a[$i]==0){
			$tmpFrq=1-$frq[$j]/$popSize;
		}
		else{
			$tmpFrq=$frq[$j]/$popSize;
		}
		$frqStep[$j]=$tmpFrq;
		# fixation time
		$fixTstep[$j]=$a[$i-1];
		# alleletype fixStep
		$alleleStep[$j]=$a[$i];
	}
	
		print $a[0], "\t", $a[1], "\t", $a[2];
	for(my $s=0; $s<@steps; $s++){
		my $orderVar=$steps[$s]-1;
		print "\t", $fixTstep[$orderVar], "\t", $alleleStep[$orderVar], "\t", $frqStep[$orderVar];
	}
	print "\n";
}

sub initFreq(){
	my $tmp=shift;
	my @genoInfo=split(/\s+/, $tmp);

	my @frq;
	for(my $i=3; $i<@genoInfo; $i+=2){
		my @b=split("", $genoInfo[$i]);
		my $j=$i+1;
		for(my $p=0; $p<@b; $p++){
			$frq[$p]+=$b[$p]*$genoInfo[$j];
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
				### swap elements
				my $tmp=$fixTime[$i];
				$fixTime[$i]=$fixTime[$j];
				$fixTime[$j]=$tmp;
				### swap orders
				my $tmp1=$order[$i];
				$order[$i]=$order[$j];
				$order[$j]=$tmp1;
			}
		}
	}
	return(@order);
}


close INTPOP; close FINAL;
