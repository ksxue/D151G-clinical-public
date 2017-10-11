#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# This script sets up some of the main analyses 
# by preparing indexed reference genomes, compiled scripts, and so on.


##################################################################
# Index reference genomes and compile scripts.
##################################################################

module load bowtie2/2.2.3
bowtie2-build reference/H3N2-Victoria-2011.fasta reference/H3N2-Victoria-2011

module load mpc/0.8.2
module load mpfr/3.1.0
module load gmp/5.0.2
module load gcc/4.9.1
g++ -O scripts/SummarizeBAM.cpp -o bin/SummarizeBAM-1.21
g++ -O scripts/AnnotateVariants.cpp -o bin/AnnotateVariants-1.1
