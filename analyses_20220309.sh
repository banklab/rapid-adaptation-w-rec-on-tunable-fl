
wdir=`pwd`
code=$wdir/code
plot=$wdir/plot
data=$wdir/data
stat=$wdir/stat


### all the stats
#:<<A
for sim in sfs3
do
	for file in $(ls $data/$sim/population_details/* | xargs basename)
	do

		name=${file%_.dat}
		# get population size from the file name
		popSize=`perl -se 'my @a=split("_", $tmp); $a[2]=~s/N//g;print $a[2];' --  -tmp=$name`
		echo $popSize
#:<<B
#### fitness change between final genotype and initial populaiton
		perl $code/getPopDetail.pl $data/$sim/population_details/$file > $wdir/tmp/$sim.$name.popDetail
### genetic distance between final genotype and initial populaiton
		perl $code/compInitialFinal.pl $data/$sim/initial_population/$file $data/$sim/fixations/$file $popSize | cut -f4  >$wdir/tmp/$sim.$name.comp
		paste $wdir/tmp/$sim.$name.popDetail $wdir/tmp/$sim.$name.comp >$stat/$sim/$name.popDetail
#B
#:<<B
### Genetic distance between intial population and the global peak; The last parameter is number of fitness landscapes. Genetic distance between the final genotype and global peak
		prefix=${name#rmf_}
		perl $code/compInitialGlobal.pl $data/$sim/initial_population/$file $popSize $wdir/landscapes/$sim $prefix 100  >$wdir/tmp/$sim.$name.comp.ig
		perl $code/compFinalGlobal.pl $data/$sim/fixations/$file $popSize $wdir/landscapes/$sim $prefix 100 | cut -f4  >$wdir/tmp/$sim.$name.comp.fg
		paste $wdir/tmp/$sim.$name.comp.ig $wdir/tmp/$sim.$name.comp.fg >$stat/$sim/$name.dist2global
#B

#:<<B
		### fixation genotype diversity #### recomb
		perl $code/finalGenoDiversity.pl $data/$sim/fixations/$file >$stat/$sim/$name.finalGenoDiversity
		perl $code/fixAlleleInitFreq.pl $data/$sim/initial_population/$file $data/$sim/fixations/$file $popSize >$stat/$sim/$name.fixAlleleInitFreq
		perl $code/fixAlleleInitFreq.ordered.pl $data/$sim/initial_population/$file $data/$sim/fixations/$file $popSize >$stat/$sim/$name.fixAlleleInitFreq.ordered
		perl $code/fixGenoInit.pl $data/$sim/initial_population/$file $data/$sim/fixations/$file $popSize >$stat/$sim/$name.finalGenoIntFrq
B
#### Adaptive path
		prefix=${name#rmf_}
		echo $prefix
#:<<C
		perl $code/epiTraj2final.pl $wdir/landscapes/$sim $data/$sim/fixations $prefix 100 >$stat/$sim/$name.epiTraj2final
		perl $code/epiTraj2major.pl $wdir/landscapes/$sim $data/$sim/fixations $prefix 100 >$stat/$sim/$name.epiTraj2major
		perl $code/oneLocusEffect2major.pl $wdir/landscapes/$sim $data/$sim/fixations $prefix 100 >$stat/$sim/$name.oneLocEffect2major
		perl $code/oneLocusEffect2final.pl $wdir/landscapes/$sim $data/$sim/fixations $prefix 100 >$stat/$sim/$name.oneLocEffect2final
#C
	done
done
#A

### All the plots are here
#:<<B
for sim in sfs3
do
#### Scripts with "mt" generate figures in the main text; Other scripts generate figures in the supplementary.
#### Fixation time distribution
	Rscript $code/fixationTime.R $plot/$sim/popFL $data/$sim
	Rscript $code/fixationTime.mt.R $plot/$sim/popFL $data/$sim

#### Fitness change and fitness between fixed genotype and global peak
	Rscript $code/fitnessChange.R $plot/$sim/popFL $data/$sim  $stat/$sim
	Rscript $code/fitnessChange.mt.R $plot/$sim/popFL $data/$sim  $stat/$sim

##### Genetic distance
	Rscript $code/distance2initial.R $plot/$sim/popFL $data/$sim  $stat/$sim
	Rscript $code/distance2initial.mt.R $plot/$sim/popFL $data/$sim  $stat/$sim
	Rscript $code/distance2global.R $plot/$sim/popFL $data/$sim  $stat/$sim

##### The initial frequency of final genotype
	Rscript $code/fixGenoInit.R $plot/$sim/popFL $stat/$sim

#### Fixed alleles
	Rscript $code/fixMinorAllelePro.R $plot/$sim/popFL $stat/$sim 
	Rscript $code/initMinorAlleleFrq.R $plot/$sim/popFL $stat/$sim
	Rscript $code/initAlleleFrq.R $plot/$sim/popFL $stat/$sim
	Rscript $code/initAlleleFrq.mt.R $plot/$sim/popFL $stat/$sim
	Rscript $code/fixAlelleFrq.mt.R $plot/$sim/popFL $stat/$sim
	Rscript $code/fixAlelleFrq.mt.1.R $plot/$sim/popFL $stat/$sim

##### Adaptive path in the maintext.
	Rscript $code/epiTraj.mt.R $plot/$sim/popFL $stat/$sim
	Rscript $code/oneLocusEffect.mt.R $plot/$sim/popFL $stat/$sim
#### Adaptive path Figure S
	Rscript $code/oneLocusEffect2final.R $plot/$sim/popFL $stat/$sim
	Rscript $code/epiTraj2final.R $plot/$sim/popFL $stat/$sim
#### Adaptive path absolute values
	Rscript $code/oneLocusEffect2final.abs.R $plot/$sim/popFL $stat/$sim
	Rscript $code/epiTraj2final.abs.R $plot/$sim/popFL $stat/$sim
### Ratio on adaptive path
	Rscript $code/adaptivePath.ratio.R $plot/$sim/popFL $stat/$sim/
	Rscript $code/adaptivePath.ratio.mt.R $plot/$sim/popFL $stat/$sim/

done

### generate legend in figures
Rscript $code/legend.R 

### Estimate the ratio of one-locus effect and pairwise epsitasis to indicate the local ruggedness. 
python RMF_ruggedness_twoLocExpectation.py


#B
