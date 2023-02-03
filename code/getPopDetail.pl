#!/usr/bin/perl
use strict;
use warnings;

open IN,$ARGV[0] || die "$!";

###landscape_id  recombination_rate  conf  t  mean_fitness        variance_fitness         max_fitness         population_entropy  n_fixations  n_genotypes
###0              0.001               0

my $tmp=<IN>;
my $flID=0;
my $simID=0;
my $recomb=0.001;
my $t=1;
my $PopFitA=0;
my $PopFitZ=0;
my @b;

print "#landscape_id\trecombination_rate\tconf\tmean_fitness_initial\tmean_fitness_final\ttfixation\n";

while(<IN>){
	my @a=split;
	if($a[3]==1 && $t==1){
		$flID=$a[0];
		$recomb=$a[1];
		$simID=$a[2];
		$PopFitA=$a[4];
		$t=1;
	}
	elsif($a[3]==1 && $t!=1){
		print join("\t", ($flID, $recomb, $simID, $PopFitA, $b[4], $t)), "\n";
		$flID=$a[0];
		$recomb=$a[1];
		$simID=$a[2];
		$PopFitA=$a[4];
		$t=1;
	}
	$t=$a[3];
	@b=@a;
}

print join("\t", ($flID, $recomb, $simID, $PopFitA, $b[4], $t)), "\n";


close IN;
