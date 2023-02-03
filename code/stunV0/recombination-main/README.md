# Hi!
Sorry for creating a new repository but it was a mess and it was easier to me like this.

## Compiling
1. Install rust if it is not installed (https://www.rust-lang.org/tools/install)
2. run ./compile or cargo build --release (if you do the second option the files will be at target/release/)

(you may need to run chmod +x compile to be able to run that compile script)

Let me know if you have trouble or don't feel like compiling the programs. I'll send them to you already compiled.

This creates three programs: generate_landscapes, haploid_recombination and calculate_haplotypes

generate_landscapes: generates the RMF landscapes

haploid_recombination: runs the populations on those landscapes

calculate_haplotypes: calculate the diversity in starting populations

## Running
you need to run generate_landscapes first and then haploid_recombination

if you run ./generate_landscapes --help it will show you a help menu (same for haploid_recombination)

for example, to create 100 landscapes with 5 loci each and mu_a = 0.01, sigma_a = 0, sigma_b = 0.05 you would run
./generate_landscapes --landscapes 100 -L 5 -m 0.01 -S 0 -s 0.05 --neutralsfs

and then to run populations (with size 200 and 50 replicates of each, recombination rates 0.001, 0.01, 0.1) on them you would run
./haploid_recombination --landscapes 100 -L 5 -m 0.01 -S 0 -s 0.05 -N 200 -r 0.001,0.01,0.1 -c 50 --neutralsfs

(notice that recombination rates are only separated by commas, no spaces)

to use the neutral SFS as initial condition you should use the command --neutralsfs. If you want to use the other initial condition you should add -d 0.1 (where 0.1 is replaced with the probability of the minor allele)

the results will be stored in the folder data/

To calculate the diversity in the initial population you can use the program (where L receives a list of L's)
./calculate_haplotypes --neutralsfs -c 100 -N 50 -L 2,3,4,5,6,7,8,9,10,11,12,15,18,21,24


## SSWM
There's a python code now to calculate the expected times to fixation and fixation probabilities under SSWM assumptions. This is based on Claudia's PNAS (2016).

It can be run with
python3 SSWM.py filename.dat

where filename.dat is the file with the landscape (probability this will look something like landscapes/rmf_landscape_L4_mua0.01_sa0_sb0.05_0.dat, for example).

The output maybe it's a bit confusing for now, then we should think how to make it more readable. Ideas are welcome!

Note: the program requires numpy library. If not installed 'pip3 install numpy' should solve it
