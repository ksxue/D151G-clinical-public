require(tidyverse)
require(cowplot)
require(foreach)
require(gridExtra)
require(RColorBrewer)

source("analysis/PlotThemes.R")

minfreq=0.03
mincoverage=100

# Find all annotated summary files in the given directory.
datapath <- "nobackup"
Files <- dir(datapath, pattern="*-annotated.summary")

# Read in all files, and combine them into a single data frame.
Data <- Files %>% 
  map_dfr(~ read_table2(file.path(datapath, .), col_names=FALSE))

# Assign column names to the data frame.
# Parse the column names for useful information.
colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                    "Sample")
Data <- Data %>% separate(Sample, into=c("Biosample","Replicate"),sep="-", remove=FALSE)

# Calculate base frequencies at each position.
Data <- Data %>% group_by(Sample, Chr, Pos) %>%
  mutate(Consensus=Base[which.max(Count)], Freq=Count/sum(Count),
         Coverage=sum(Count), Subtype="H3")

# Add HA codon numbering as well as a text label indicating gene and codon number.
Data$Codon <- as.numeric(Data$Codon)
Data <- Data %>% ungroup() %>%
  mutate(Codon=ifelse(Chr=="4-HA",
                        ifelse(Subtype=="H3", Codon-16, Codon-13),Codon))

# Analyze only genomic positions that are in the coding portions of genes.
Data <- Data %>% filter(Gene!="none") %>% filter(Codon > 0)

# Remove extraneous columns.
Data <- Data %>% dplyr::select(-Syn, -FourfoldSyn, -AvgQ, -AvgReadPos, -Subtype)

# Merge data for the two replicates.
DataMerged <- merge(Data %>% filter(Replicate==1), Data %>% filter(Replicate==2),
                    by=c("Biosample","GenomePos","Chr","Gene","Pos","Base","Codon"))

# Test the effect of different variant-calling thresholds
# on the proportion of variants called in one replicate that are
# also called in the other.
freqs <- rep(seq(0.01,0.1,0.01), each=3)
coverages <- rep(seq(100,200,50), 10)
FreqThresholds <- 
  foreach(freq=freqs, coverage=coverages, .combine='rbind') %do%
          {
           table((DataMerged %>% filter(Base!=Consensus.x) %>%
                    mutate(Variant.x=ifelse(Freq.x>freq & 
                                            Coverage.x>coverage,1,0),
                           Variant.y=ifelse(Freq.y>freq &
                                            Coverage.y>coverage,1,0),
                           VariantCount=Variant.x+2*Variant.y)
                  )$VariantCount)   
          }
FreqThresholds <- cbind(freqs, coverages, FreqThresholds)
colnames(FreqThresholds) <- c("Freq","Coverage","Var0","Var1","Var2","Var3")

# Calculate percentages of variants in replicates 1 and 2
# that are also called in the other replicate.
FreqThresholds <- as.data.frame(FreqThresholds) %>%
  mutate(PctDouble1=(Var3/(Var1+Var3)),
         PctDouble2=(Var3/(Var2+Var3))) %>%
  mutate(PctDoubleAvg=(PctDouble1+PctDouble2)/2)

# Plot the correspondance between replicates
# for different thresholds of minimum frequency and sequencing coverage.
p_thresholds <- ggplot(FreqThresholds) +
  geom_point(aes(x=Freq, y=PctDoubleAvg, shape=factor(Coverage))) +
  geom_line(aes(x=Freq, y=PctDoubleAvg, 
                group=factor(Coverage), linetype=factor(Coverage))) +
  ylim(0,1) + ylab("Percent correspondance between variant calls\nacross replicates") +
  scale_x_continuous(name="Minimum frequency",
                       limits=c(0,0.1), breaks=seq(0,0.1,0.01)) +
  scale_linetype_discrete(name="Minimum coverage") +
  scale_shape_discrete(name="Minimum coverage") +
  THEME_ALL
save_plot(paste(DIR,"thresholds.pdf", sep=""), p_thresholds)

# Using the frequency and coverage thresholds specified above,
# call variants in the sequencing data.
# Output a table of variants and their average frequencies between replicates.
Variants <- DataMerged %>% filter(Base!=Consensus.x) %>%
  filter(Freq.x>minfreq, Freq.y>minfreq, 
         Coverage.x>mincoverage, Coverage.y>mincoverage)

# Discard any variants in which the reference base differs between replicates,
# which would indicate some kind of analysis error.
Variants <- Variants %>% filter(RefAA.x==RefAA.y)

# Format the list of variants for printing.
Variants <- Variants %>% mutate(Frequency=(Freq.x+Freq.y)/2) %>%
  separate(Gene, into=c("GeneID","GeneShort"),sep="-") %>%
  mutate(Variant=paste0(GeneShort,'-',RefAA.x,Codon,AltAA.x)) %>%
  dplyr::select(Biosample, Variant, Frequency)
colnames(Variants) <- c("Sample","Variant","Frequency")

# Output table of variants.
save_plot(paste0(DIR,"variants.pdf"),
          grid.arrange(
            tableGrob(Variants, rows=NULL,
                      theme=ttheme_minimal(
                        base_size=10, padding=unit(c(2,2),"mm"),
                        core=list(fg_params=list(cex=0.5)),
                        colhead=list(fg_params=list(cex=0.5)),
                        rowhead=list(fg_params=list(cex=0.5, fontface=1)))
            ),
            nrow=1),
          base_width=2, base_height=4)

# Output table of variants in tab-delimited format.
write.table(format(Variants,digits=2), paste0(DIR,"variants.txt"),
            row.names=FALSE, quote=FALSE)

# Plot the number of variants called per clinical sample.
p_numvariants <- Variants %>% group_by(Sample) %>%
  summarize(NumVariants=n()) %>%
  ggplot() + geom_bar(aes(x=Sample, y=NumVariants), stat="identity") +
  THEME_ALL + theme(axis.text.x=element_text(angle=90)) +
  xlab("Sample") + ylab("Number of variants")
p_varfreq <- Variants %>% group_by(Sample) %>%
  separate(Variant,into=c("Gene","Site"), remove=FALSE) %>%
  ggplot() + 
  geom_point(aes(x=factor(Sample), y=Frequency, 
                 color=factor(Gene))) +
  scale_color_brewer(palette="Set2", name="Gene") +
  THEME_ALL + theme(axis.text.x=element_text(angle=90)) +
  xlab("Sample") + ylab("Variant frequency") +
  geom_text(aes(x=factor(Sample), y=(Frequency), label=Variant), 
            nudge_x=-0.15, nudge_y=0.01, size=2)
p_variants <- plot_grid(p_numvariants, p_varfreq, ncol=2)
p_variants
save_plot(paste0(DIR,"variantstats.pdf"), p_variants, ncol=2)

# Plot distribution of non-consensus amino acids at each codon site.
DataDist <- DataMerged %>% 
  mutate(Count=Count.x+Count.y, Coverage=Coverage.x+Coverage.y) %>%
  group_by(Biosample, Gene, Codon) %>% filter(Base!=Consensus.x) %>%
  summarize(Freq=sum(Count)/mean(Coverage))

# Plot the distribution of non-consensus amino acid counts
# at sites with sub-threshold variant frequencies.
p_distribution <- ggplot(DataDist) + geom_histogram(aes(x=Freq), binwidth=0.0005) +
  geom_vline(data=DataDist %>% filter(Gene=="6-NA", Codon==151),
             aes(xintercept=Freq), color="firebrick3", size=0.7) +
  xlim(0,minfreq) + ylab("Number of codon sites") + ylim(0,1000) +
  xlab("Frequency of minority variants") +
  facet_wrap(~Biosample, ncol=3, scales="free") +
  THEME_ALL
save_plot(paste0(DIR,"distribution.pdf"), p_distribution,
          base_height=5, base_width=5)
