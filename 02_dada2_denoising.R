### NOTES ----

# Markers used:

# 16S
# Chiar16SF–Chiar16SR
# Forward: TARTYCAACATCGRGGTC
# Reverse: CYGTRCDAAGGTAGCATA
# Marquina et al. 2019

# CO1 (Bombus)
# BombusF-BombusR (~125 bp)
# Forward: AGWCAYCCTGGAATATGAA
# Reverse: GTGGRAAAGCTATATCAGG
# Milam et al. 2020

# CO1 (Jusino)
# ANML: LCO1490--CO1‐CFMRa (~180 bp)
# Forward: GGTCAACAAATCATAAAGATATTGG
# Reverse: GGWACTAATCAATTTCCAAATCC
# Folmer et al. 1994 and Jusino et al. 2019

# load libraries
library(dada2)
library(decontam)
library(tidyverse)

# Modified from: https://benjjneb.github.io/dada2/tutorial.html

# Import reads ----

# Set path to trimmed fastq files
path <- "path/to/trimmed"

# Read in list of full path names. Assumes forward and reverse fastq filenames
# have format: SAMPLE-NAME_R1.trimmed.fq.gz and SAMPLE-NAME_R2.trimmed.fq.gz
fnFs <- sort(list.files(path, pattern="_R1.trimmed.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.trimmed.fq.gz", full.names = TRUE))

# return sample names without _R1/_R2 or fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# Filter and trim ----

# Create subdirectory named "filtered" and place files in that directory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Run filterAndTrim step (settings dependent on marker)
# Jusino
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=c(2,2), maxN=0, truncLen=150, minQ=30, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
# bombus
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=c(2,2), maxN=0, truncLen=100, minQ=30, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
# 16S
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=c(2,2), maxN=0, truncLen=180, minQ=30, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Remove files with zero reads (otherwise they will cause an error during the
# inference step)
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

# Learn error rates ----

# Machine learning of error rates in the data set (this may take a little while...)
errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE, nbases=1e9, verbose=TRUE)

# Plot errors. Estimated error rates (black line) should roughly match observed
# rates (grey dots)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Infer number of true variants using pooled samples (see
# https://benjjneb.github.io/dada2/pseudo.html for more information)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = TRUE)

# Merge paired reads ----
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, trimOverhang = TRUE, verbose=TRUE)

# Make a sequence table from merged reads
seqtab <- makeSequenceTable(mergers)


# Remove contaminants ----
# See package `decontam` package here:
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Import metadata (must include a column with TRUE/FALSE indicating whether
# sample is a blank or is not)
meta <- read.table("path/to/meta_16S_bilsoda22.tsv", header = TRUE, sep="\t")

# Remove samples from metadata that were removed from analyses
samples <- as.data.frame(names(mergers))
colnames(samples) <- "sample_ID"
meta <- merge(samples, meta, by="sample_ID", all=FALSE)

# Extract TRUE/FALSE for blanks
vector_for_decontam <- meta$blank

# Identify contaminants based on blanks
contam_df <- isContaminant(seqtab, neg=vector_for_decontam, threshold=0.1)

# Pull out ASV sequences that are potential contaminants
contam_asvs <- as.matrix(row.names(contam_df[contam_df$contaminant == TRUE, ]))

# Export contaminants to fasta
contam_seq <- Biostrings::DNAStringSet(contam_asvs)
names(contam_seq) <- paste0("contam", 1:nrow(contam_asvs))
Biostrings::writeXStringSet(contam_seq, "ASV_contaminants.fasta", append=FALSE, compress=FALSE, format="fasta")

# Remove contaminants from data
seqtab.decontam <- seqtab[,!colnames(seqtab) %in% contam_asvs]

# Remove chimeric sequences ----
seqtab.nochim <- removeBimeraDenovo(seqtab.decontam, method="consensus", multithread=TRUE, verbose=TRUE)

# Remove ASVs outside the expected length range
# For CO1, this is +/- 1 codon from 181 bp (jusino) or 125 bp (bombus)
# For 16S, mean recorded length (from Marquina et al) is 348 bp, so did +/- 10% bp
# Note: for 16S, min before filtering = 163 and max = 348
# jusino
#seqtab.filtered <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 178:184]
# bombus
#seqtab.filtered <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 121:128]
# 16S
seqtab.filtered <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 313:383]

# Make stats table ----
getN <- function(x) sum(getUniques(x))

# Prepare data frames with statistics for merging
filtered <- as.data.frame(out)
denoisedF <- as.data.frame(sapply(dadaFs, getN))
denoisedR <- as.data.frame(sapply(dadaRs, getN))
merged <- as.data.frame(sapply(mergers, getN))
decontam <- as.data.frame(rowSums(seqtab.decontam))
nonchim <- as.data.frame(rowSums(seqtab.nochim))
sized <- as.data.frame(rowSums(seqtab.filtered))

# Make row names into a column
filtered$rn <- rownames(filtered)
denoisedF$rn <- rownames(denoisedF)
denoisedR$rn <- rownames(denoisedR)
merged$rn <- rownames(merged)
decontam$rn <- rownames(decontam)
nonchim$rn <- rownames(nonchim)
sized$rn <- rownames(sized)

# Clean up row names in the `filtered` data frame so they match the sample names
# in the rest of the data frames
filtered <- base::transform(filtered, rn = sub("_R1.trimmed.fq.gz", "", rn))

# Merge the data frames and clean up
track <- join_all(list(filtered,denoisedF,denoisedR,merged,decontam,nonchim,sized), by = 'rn', type = 'full')
rownames(track) <- track$rn
track <- subset(track, select = -rn)
colnames(track) <- c("Input", "After filtering", "After denoising (F)", "After denoising (R)", "After merging", "After decontamination", "After chimera removal", "After size filtering")
write.csv(track, "dada2_denoising_stats.csv", row.names = TRUE)

# Export data ----
# See also: https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# Extract and rename sequences
ASVs <- Biostrings::DNAStringSet(colnames(seqtab.filtered))
names(ASVs) = paste0("ASV", 1:ncol(seqtab.filtered))

# Export ASVs to fasta
Biostrings::writeXStringSet(ASVs, "ASVs.fasta", append=FALSE, compress=FALSE, format="fasta")

# Create a table of counts for each ASV sequence and export
colnames(seqtab.filtered) <- names(ASVs)
asv.tab <- t(seqtab.filtered)
write.table(asv.tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# GO TO 03_parse_BLAST_results.R
