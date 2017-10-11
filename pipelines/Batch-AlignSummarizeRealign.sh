#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# Script for batch submission of SRADownloadFASTQ jobs.

# Input parameters.
pipeline="pipelines/AlignSummarizeRealign.sh"
samplesheet="$1" # Give the path to a samplesheet listing rawread1, rawread2, trimmed1, trimmed2, samplename.
dir="$2" # Give the desired working directory for the trimmed FASTQ files.
reference="$3" # Give reference path WITHOUT .fasta or .fa extension.
refname="$4" # Give a short nickname for the reference sequence.
clean="$5" # If 0, then the script will run in its entirety, overwriting existing output.

# Extract the appropriate line of the list of run IDs to submit to the script.
fastq1=`cut -f3 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"` # Extract path of trimmed read 1.
fastq2=`cut -f4 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"` # Extract path of trimmed read 2.
sample=`cut -f5 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"` # Extract sample name.
${pipeline} ${fastq1} ${fastq2} ${sample} ${dir} ${reference} ${refname} ${clean}