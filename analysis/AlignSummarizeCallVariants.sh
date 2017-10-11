#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# This script trims Nextera adapters from sequencing reads,
# filters out reads that map to the human genome,
# aligns sequences to the designated reference genome,
# summarizes the alignment of the reads against each reference genome,
# creates a new FASTQ file for each sample representing the consensus sequences,
# re-aligns the sequencing reads against the sample-specific reference genome,
# summarizes the sample alignment,
# annotates each position in the genome,
# and calls variants in the resulting read summary file.

# File paths.
samplesheet="data/Samples.txt"
batchtrimadapters="pipelines/Batch-TrimAdaptersNextera.sh"
batchalignsummarize="pipelines/Batch-AlignSummarizeRealign.sh"
trimdir="nobackup/trimmed"
outdir="nobackup"
reference="reference/H3N2-Victoria-2011"
clean=1 # If set to 0, all analyses will rerun from the start, overwriting existing intermediates.

# Trim Nextera adapters from all reads.
# Calculate the number of biological samples to be analyzed.
# Note that this is different from the number of sequencing runs above.
numsamples="$(wc -l ${samplesheet} | cut -f1 -d' ')"
echo ${numsamples}

# Submit batch job to trim Nextera adapters from raw sequencing reads.
qsub -cwd -N TrimAdapters -l h_rt=48:00:00 -t 1-${numsamples} -tc 50 \
  -o nobackup/sge/ -e nobackup/sge/ \
  ${batchtrimadapters} ${samplesheet} ${trimdir} ${clean}
  
# Submit batch jobs to align reads to the Victoria-2011 reference.
qsub -cwd -N AlignVic11 -l h_rt=48:00:00 -l m_mem_free=4G \
  -t 1-${numsamples} -tc 50 -hold_jid TrimAdapters \
  -o nobackup/sge/ -e nobackup/sge/ \
  ${batchalignsummarize} ${samplesheet} ${outdir} ${reference} Vic11 \
  reference/H3N2-Victoria-2011.bed ${clean}
  
# Analyze replicability and call variants in the sequenced samples.
Rscript analysis/AnalyzeReplicability.R
Rscript analysis/CallVariants.R
