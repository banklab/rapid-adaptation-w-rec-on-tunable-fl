#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

open FIX, $ARGV[0] ||  "$!";
#landscape_id   recombination_rate      conf    fixation_time allele


my @recomb=(0, 0.001, 0.01, 0.1, 0.5);

while(<FIX>){
	### fitness landscape ID
	for(my $i=0; $i<100;$i++){
		#### recombination ID
		for(my $r=0; $r<@recomb; $r++){
			my %genofrq;
			for(my $j=0; $j<100; $j++){
				my $tmp=<FIX>;
				chomp($tmp);
				my @a=split(/\s+/, $tmp);
				my $geno;
				for(my $allele=4; $allele<@a; $allele+=2){
					$geno.=$a[$allele];
				}
				if(defined $genofrq{$geno}){
					$genofrq{$geno}++;
				}
				else{
					$genofrq{$geno}=1;
				}
			}
			my $e=&entropy(values %genofrq);
			print "$i\t$recomb[$r]\t$e\n";
		}
	}
}



close FIX;

sub entropy(){
	my @freq=@_;

	my $total=sum(@freq);
	my $entropy=0;

	foreach my $f (@freq){
		$f=$f/$total;
		$entropy+=-$f*log($f)/(log(10));
	}

	return $entropy;
}
