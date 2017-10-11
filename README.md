# D151G clinical sequencing

This repository analyzes deep-sequencing data from clinical H3N2 influenza samples. Clinical samples were obtained from the Washington State Public Health Laboratory (WSPHL).  The project aims in particular to determine whether the samples contain detectable frequencies of the D151G neuraminidase mutation. We have previously shown that the D151 and G151 viral variants are able to cooperate in cell culture and wish to determine whether this cooperation may also occur in clinical sequences. All scripts are meant to be run from the top level of this repository.

# repository organization

**reference** - This folder contains the A/Victoria/361/2011 (H3N2) reference sequence in FASTA format and the associated BED annotation file.

**scripts** - This folder contains custom C++ scripts written to process mapped reads. These scripts tally base frequencies and annotate variants within a gene.

**pipelines** - This folder contains shell scripts that process raw sequencing data by filtering out reads that map to the human genome, trimming adapters, mapping reads to a reference sequence, summarizing base frequencies at each position, extracting a consensus sequence for each sample, re-aligning reads to the consensus sequence for that sample, identifying common variants, and annotating them.

**analysis** - This folder contains custom R and shell scripts to perform relevant analyses. The "Setup.sh" script prepares genomes for mapping and compiles the C++ scripts used for custom analyses. The "AlignSumamrizeCallVariants.sh" script submits jobs to align all reads to the reference genomes, tally bases at each position, analyze sample replicates, and call variants in sequencing data.

**data** - This folder contains small data files, including a sample sheet listing the paths to all raw sequence reads. It also includes the GISAID acknowledgement table for the samples analyzed in order to find the samples sequenced in depth in this study, as well as the script, "ParseD151GFrequency.py", that was used to perform this analysis.