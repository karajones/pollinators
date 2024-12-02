# COMMAND LINE ----

# produce a blastdatabase with the OTUs
makeblastdb -in ASVs.fasta -parse_seqids -dbtype nucl

#  blast the OTUs against the database
blastn -db ASVs.fasta -outfmt '6 qseqid sseqid pident' -out match_list_84.txt -qcov_hsp_perc 80 -perc_identity 84 -query ASVs.fasta
blastn -db ASVs.fasta -outfmt '6 qseqid sseqid pident' -out match_list_90.txt -qcov_hsp_perc 80 -perc_identity 90 -query ASVs.fasta
blastn -db ASVs.fasta -outfmt '6 qseqid sseqid pident' -out match_list_94.txt -qcov_hsp_perc 80 -perc_identity 94 -query ASVs.fasta
blastn -db ASVs.fasta -outfmt '6 qseqid sseqid pident' -out match_list_97.txt -qcov_hsp_perc 80 -perc_identity 97 -query ASVs.fasta

# END COMMAND LINE ----


# Curation with LULU ----
library(lulu) # https://github.com/tobiasgf/lulu

# curate using 84% min match (default)
matchlist <- read.table("match_list_84.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
lulu_result <- lulu(as.data.frame(asv.tab), matchlist, minimum_match = 84)
lulu_result_84_table <- lulu_result[["curated_table"]]

# curate using 90% match (suggested in Brandt et al. 2021)
matchlist_90 <- read.table("match_list_90.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
lulu_result_90 <- lulu(as.data.frame(asv.tab), matchlist_90, minimum_match = 90)
lulu_result_90_table <- lulu_result_90[["curated_table"]]

# curate using 94% match
matchlist_94 <- read.table("match_list_94.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
lulu_result_94 <- lulu(as.data.frame(asv.tab), matchlist_94, minimum_match = 94)
lulu_result_94_table <- lulu_result_94[["curated_table"]]

# curate using 97% match
matchlist_97 <- read.table("match_list_97.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
lulu_result_97 <- lulu(as.data.frame(asv.tab), matchlist_97, minimum_match = 97)
lulu_result_97_table <- lulu_result_97[["curated_table"]]

# Export count tables
write.table(lulu_result_84_table, "lulu/ASVs_counts_lulu_84.tsv", sep="\t", quote=F, col.names=NA)
write.table(lulu_result_90_table, "lulu/ASVs_counts_lulu_90.tsv", sep="\t", quote=F, col.names=NA)
write.table(lulu_result_94_table, "lulu/ASVs_counts_lulu_94.tsv", sep="\t", quote=F, col.names=NA)
write.table(lulu_result_97_table, "lulu/ASVs_counts_lulu_97.tsv", sep="\t", quote=F, col.names=NA)

# Export OTU maps
write.table(lulu_result[["otu_map"]], "lulu/otu_map_lulu_84.tsv", sep="\t", quote=F, col.names=NA)
write.table(lulu_result_90[["otu_map"]], "lulu/otu_map_lulu_90.tsv", sep="\t", quote=F, col.names=NA)
write.table(lulu_result_94[["otu_map"]], "lulu/otu_map_lulu_94.tsv", sep="\t", quote=F, col.names=NA)
write.table(lulu_result_97[["otu_map"]], "lulu/otu_map_lulu_97.tsv", sep="\t", quote=F, col.names=NA)

# Subset ASVs using vector of ASV names from lulu results
ASVs_lulu_84 <- ASVs[rownames(lulu_result_84_table)]
ASVs_lulu_90 <- ASVs[rownames(lulu_result_90_table)]
ASVs_lulu_94 <- ASVs[rownames(lulu_result_94_table)]
ASVs_lulu_97 <- ASVs[rownames(lulu_result_97_table)]

# Export ASVs to fasta
Biostrings::writeXStringSet(ASVs_lulu_84, "lulu/ASVs_lulu_84.fasta", append=FALSE, compress=FALSE, format="fasta")
Biostrings::writeXStringSet(ASVs_lulu_90, "lulu/ASVs_lulu_90.fasta", append=FALSE, compress=FALSE, format="fasta")
Biostrings::writeXStringSet(ASVs_lulu_94, "lulu/ASVs_lulu_94.fasta", append=FALSE, compress=FALSE, format="fasta")
Biostrings::writeXStringSet(ASVs_lulu_97, "lulu/ASVs_lulu_97.fasta", append=FALSE, compress=FALSE, format="fasta")


# Import LULU 97 results ----
# Note: do this part after running 03_parse_BLAST_results.R

# Import clustered LCA
taxa <- read.table("lulu/lulu_97_LCA.tsv", header = TRUE, sep="\t")
# Import clustered ASVs
ASVs <- Biostrings::readDNAStringSet("lulu/ASVs_lulu_97.fasta")
# Import ASV count table
asv.tab <- as.matrix(read.table("lulu/ASVs_counts_lulu_97.tsv", header = TRUE, sep="\t", row.names=1, check.names=FALSE))

# Grab ASV names from asv.tab to merge with taxa,
# otherwise phyloseq deletes ASVs without taxonomic assignments!
taxa_all_ASVs <- as_tibble(rownames(asv.tab)) %>% 
  # rename columm to ASV so it matches lca
  dplyr::rename("ASV" = value) %>%
  # join to LCA
  left_join(taxa, by = "ASV") %>%
  # add rownames
  column_to_rownames(var = "ASV") %>%
  as.matrix()

# Add sample names to metadata
rownames(meta) <- meta[,1]
meta$sample <- as.factor(meta$sample)

# make phyloseq object
ps <- phyloseq::phyloseq(otu_table(t(asv.tab), taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa_all_ASVs), refseq(ASVs))

