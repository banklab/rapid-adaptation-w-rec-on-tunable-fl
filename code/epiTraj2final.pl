#!/usr/bin/perl

use strict;
use warnings;

my $fldir=$ARGV[0];
my $fixdir=$ARGV[1];
my $prefix=$ARGV[2]; ### parameter description
my $flprefix=$prefix; $flprefix=~s/N[0-9]+_//g;

open FIX, "$fixdir/rmf_$prefix\_.dat" || die "$!";

#print("#landscape_id\trecombination_rate\tconf\tepistasisStep\n");

## number of loci;
my $L=$prefix;
$L=~s/L([0-9]+)_//g;
$L=$1;

my $tmp=<FIX>;

for(my $flid=0; $flid<100; $flid++){
	open FL, "$fldir/rmf_landscape_$flprefix\_$flid.dat";
	### record the fitness landscape as hash (?). This might not the best way. 
	my %fl;
	while(<FL>){
		chomp;
		my @a=split;
		$fl{join("", @a[0..$L-1])}=$a[$L];
	}
	close FL;

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
		

		my @epi;

		for(my $s=0; $s<$L-1; $s++){
			my $c1=substr($final, $steps[$s]-1, 1);
			my $a_c1=1-$c1;
			my $c2=substr($final, $steps[$s+1]-1, 1);
			my $a_c2=1-$c2;

			my $geno1=$final; substr($geno1, $steps[$s]-1, 1)=$a_c1;
			my $geno2=$final; substr($geno2, $steps[$s+1]-1, 1)= $a_c2;
			my $geno3=$geno2; substr($geno3, $steps[$s]-1, 1)= $a_c1;
			
			$epi[$s]=log($fl{$final})-(log($fl{$geno1})+log($fl{$geno2})-log($fl{$geno3}));
			
		}
		print join("\t",(@a[0..2], @epi)), "\n";
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

