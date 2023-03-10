# Rapid adaptation of recombining populations on tunable fitness landscapes
Accepted on Molecular Ecology in Feb.2023

Two main steps to generate our results presented in our paper:
## 1. Simulating dynamics of adaptation on Rough Mount Fuji (RMF) fitness landscapes. 

Simulations are summarized in bash_shell scripts 'fl_DT_20220309.sh'.

The main parameters are:

a. dt = 3 ### genetic drift threshold is set to 3, meaning the minimal allele frequency is 3.

b. L = 5 or 15 ### numbers of loci on fitness landscapes. L=5 represent small landscapes; L=15 represent large landscapes

c. add = 0.01 ### A constant value is assigned to each allele 1 as the additive term in RMF landscapes.

d. epi = 0 0.001 0.01 0.05 0.1 ### The variance of the normal distribution from which we drew random values and assigned values to epistatic terms of all genotypes in RMF landscapes. In the paper, we presented epi = 0.001, 0.01 and 0.1 to represent three levels of ruggedness, smooth, intermediate and high ruggedness.

e. popSize=100 500 5000 #### population size in the simulations. 

Note: All simulations were run by STUN (version 0) which contains two main functions, "generate_landscapes" and simulating "haploid_recombination". Now we have integrated those two functions into one efficient and user-friendly software, STUN (GitHub). See STUN details (user manual) and cite STUN (biorxiv). 

## 2. Analyzing features of the initial population and fitness landscapes, summaries of statistics, and generating figures and tables
Those analyses were coded in 'analyses_20220309.sh'.

Features of initial population and fitness landscapes include:

a. fitness change between final genotype and initial population;

b. genetic distance between final genotype and initial population;

c. genetic distance between the intial population and the global peak;

d. fixation genotype diversity

e. additive effect and epistasis on the adaptive path. 

These features were summarized by custom Perl scripts which can be found in the folder "code".

The figures and tables are generated by custom R scripts which can be found in the folder "code". Scripts for the main figures are with "mt" in the script file names. And all the rest of the R scripts are used to generate supplementary figures and table S1. 

In addition, we estimated the ratio of one-locus effect and pairwise epistasis to indicate the local ruggedness by a custom python script "RMF_ruggedness_twoLocExpectation.py" which can also be found in the folder "code". 
