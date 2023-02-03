#!/usr/bin/perl

use strict;
use warnings;

open FIX, $ARGV[0] || die "$!"; ### get fixed file
#landscape_id   recombination_rate      conf    genotype        count
my $popSize=$ARGV[1];

my $fldir=$ARGV[2];
my $prefix=$ARGV[3]; ### parameter description
my $numFL=$ARGV[4];
my $flprefix=$prefix; $flprefix=~s/N[0-9]+_//g;

my $L=$prefix;
$L=~s/L([0-9]+)_//g;
$L=$1;

my %globalGeno;

for(my $flid=0; $flid<$numFL; $flid++){
	my $fitness=0;
	open FL, "$fldir/rmf_landscape_$flprefix\_$flid.dat" || die "$!";
	while(<FL>){
		chomp;
		my @a=split;
		if($fitness < $a[$L]){
			$fitness=$a[$L];
			$globalGeno{$flid}=join("", @a[0..$L-1]);
		}
	}
	close FL;
}

print("#landscape_id\trecombination_rate\tconf\tdistFinalGlobal\n");

### Remove title
<FIX>;
while(<FIX>){
	chomp;
	my @a=split(/\s+/, $_);
	
	my $final="";

	for(my $i=4; $i<@a; $i+=2){
		$final=$final.$a[$i];
	}

	my $d=&compGenotype($globalGeno{$a[0]}, $final);

	print join("\t", (@a[0..2], $d)), "\n";
}

close FIX;

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

