
wdir=`pwd`
### Andre's github about FL in haploid population with recombination.
#hrM=/Users/lijuan/Desktop/work/fitnessLandscape/haploid/recomb/code/gitLab/recombination-main
hrM=$wdir/code/stunV0/recombination-main
data=$wdir/data

#:<<B
echo $data

for dt in 3
do
	[ -d $data/sfs$dt ] && rm -r $data/sfs$dt/
	mkdir $data/sfs$dt
	[ -d $wdir/landscapes/sfs$dt ] && rm -r $wdir/landscapes/sfs$dt/ 
	mkdir $wdir/landscapes/sfs$dt
	rm $wdir/landscapes/rmf*
	[ -d $data/data ] && rm -r $data/data/ $data/fixations/ $data/initial_population/ $data/population_details/
	mkdir $data/data $data/fixations  $data/initial_population $data/population_details

#:<<B
#### L: number of loci; add, additive effect on RMF fitness landscapes; epi: SD of epistasis on RMF landscapes;
	for L in 5 15
	do
		for add in 0.01
		do
			for epi in 0 0.001 0.01 0.05 0.1
			do
				##### generate fitness landscapes
				$hrM/generate_landscapes --landscapes 100 -L $L -m $add -s $epi
				for popSize in 500 5000 100
				do
					#### simulate adaptive dynamics.
					echo "$L\t$add\t$epi\t$popSize\t$dt"
					$hrM/haploid_recombination --landscapes 100 -L $L -m $add -s $epi -N $popSize -r 0,0.001,0.01,0.1,0.5 -c 100 --neutralsfs $dt &
				done
			done
		done
	done
	wait
	mv $data/data $data/fixations  $data/initial_population  $data/population_details $data/sfs$dt
	mv $wdir/landscapes/rmf_landscape* $wdir/landscapes/sfs$dt
#B
done
#B

exit 0
