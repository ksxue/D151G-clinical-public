require(tidyverse)
require(cowplot)

source("analysis/PlotThemes.R")

minfreq=0.01
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
  mutate(CodonHA=ifelse(Chr=="4-HA",
                        ifelse(Subtype=="H3", Codon-16, Codon-13),Codon))

# Analyze only genomic positions that are in the coding portions of genes.
Data <- Data %>% filter(Gene!="none")

# Ensure that positions in multiple genes are not counted multiple times.
Data <- Data %>% 
  dplyr::select(Sample, Biosample, Replicate,
                GenomePos, Chr, Pos, Base, Consensus, Count, Coverage, Freq) %>%
  distinct()

# Calculate the length and midpoint of each chromosome for plotting purposes.
ChrAnnotations <- Data %>% 
  dplyr::select(Chr,GenomePos) %>% distinct() %>%
  group_by(Chr) %>% summarize(Midpoint=as.integer(mean(GenomePos)),
                              GenomeMin=min(GenomePos),
                              GenomeMax=max(GenomePos))

# Plot coverage in 50bp bins across the genome.
p_coverage <- ggplot(Data %>% mutate(BinnedPos=floor(GenomePos/50)*50) %>%
         group_by(Biosample, Replicate, BinnedPos) %>%
         summarize(AvgCoverage=mean(Coverage))) +
  geom_rect(xmin=ChrAnnotations$GenomeMin[2], 
            xmax=ChrAnnotations$GenomeMax[2], 
            ymin=0, ymax=10e5, fill="gray80") +
  geom_rect(xmin=ChrAnnotations$GenomeMin[4], 
            xmax=ChrAnnotations$GenomeMax[4], 
            ymin=0, ymax=10e5, fill="gray80") +
  geom_rect(xmin=ChrAnnotations$GenomeMin[6], 
            xmax=ChrAnnotations$GenomeMax[6], 
            ymin=0, ymax=10e5, fill="gray80") +
  geom_rect(xmin=ChrAnnotations$GenomeMin[8], 
            xmax=ChrAnnotations$GenomeMax[8], 
            ymin=0, ymax=10e5, fill="gray80") +
  geom_line(aes(x=BinnedPos, y=AvgCoverage, 
                group=Replicate, linetype=factor(Replicate))) +
  facet_wrap(~Biosample, ncol=3, scales="free") +
  scale_y_log10(limits=c(10,10^5)) +
  scale_x_continuous(limits=c(0,13500), expand=c(0,0),
                     breaks=ChrAnnotations$Midpoint,
                     labels=c("PB2","PB1","PA","HA","NP","NA","M","NS")) +
  scale_linetype_discrete(name="Replicate") +
  xlab("Genome position (bp)") + ylab("Average coverage") +
  theme(strip.background=element_blank()) +
  THEME_ALL
save_plot(paste(DIR,"coverage.pdf", sep=""), p_coverage,
          base_height=5, base_width=7)

# Call variants that reach a certain frequency and coverage threshold.
Data <- Data %>%
  mutate(Variant=ifelse(Base!=Consensus & Freq>minfreq & Coverage>mincoverage,1,0))

# Split and then merge information for the two replicates.
DataMerged <- merge(Data %>% filter(Replicate==1), Data %>% filter(Replicate==2),
                    by=c("Biosample","GenomePos","Chr","Pos","Base"))

# Calculate the total number of variant calls at each site.
# 0 = site is called as a variant in no replicates
# 1 = site is called as a variant in replicate 1 only
# 2 = site is called as a variant in replicate 2 only
# 3 = site is called as a variant in both replicates
DataMerged <- DataMerged %>% mutate(VariantCalls=Variant.x + 2*Variant.y)
table(DataMerged$VariantCalls)

# Plot the frequencies of variants that are called in both samples.
p_replicability <- ggplot(DataMerged %>% filter(Variant.x==1, Variant.y==1)) +
  geom_point(aes(x=Freq.x, y=Freq.y)) +
  facet_wrap(~Biosample, ncol=3, scales="free") +
  xlim(0,1) + ylim(0,1) +
  theme(strip.background=element_blank()) +
  xlab("Mutation frequency, replicate 1") +
  ylab("Mutation frequency, replicate 2") +
  THEME_ALL
save_plot(paste(DIR,"replicability.pdf", sep=""), p_replicability,
          base_height=5, base_width=5)
