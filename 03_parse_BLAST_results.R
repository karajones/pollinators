library(tidyverse)

# BLAST to LCA ----

# This code assumes that BLAST is run in the command line and the results are
# available for import. To ensure that the blast results have the same columns
# as below, include the following in the blastn code: -outfmt="6 qseqid sseqid
# pident length qcovs mismatch gapopen evalue bitscore"


## Import and clean BLAST results ----

# column names for blast results
col_names <- c("ASV", "taxonomy", "Identity", "Length", "Coverage", "Mismatch", "Gap", "E-value", "Bitscore")

# read in BLAST results
blast <- read_table("blast_results.tsv", col_names = col_names, na = "NA")

# separate taxonomy out into separate columns using delimiter (semi-colon)
# automatically names the column by adding a number to end of "taxonomy"
blast <- separate_wider_delim(blast, taxonomy, delim = ";", names_sep = "")

# names for taxonomy columns
col_names_tax <- c(taxonomy1 = "Phylum", taxonomy2 = "Class", taxonomy3 = "Order", taxonomy4 = "Family", taxonomy5 = "Genus", taxonomy6 = "Species")

# clean up taxonomy
blast_clean <- blast %>% 
  # rename taxonomy columns using the vector of names stored in col_names_tax
  plyr::rename(col_names_tax) %>% 
  # remove any underscores and replace with spaces
  mutate(Phylum = str_replace_all(Phylum, "p__", "")) %>% 
  mutate(Class = str_replace_all(Class, "c__", "")) %>%
  mutate(Order = str_replace_all(Order, "o__", "")) %>%
  mutate(Family = str_replace_all(Family, "f__", "")) %>%
  mutate(Genus = str_replace_all(Genus, "g__", "")) %>%
  mutate(Species = str_replace_all(Species, "s__", "")) %>%
  mutate(Species = str_replace_all(Species, "_", " ")) %>%
  # remove genera from species column 
  # this replaces any values in the species column that do NOT have a space with NA
  mutate(Species = if_else(str_detect(Species, " "), Species, NA)) %>% 
  # change the columns that are numbers to numeric
  # otherwise the filtering below will not work correctly
  mutate_at(c("Identity", "Length", "Coverage","Mismatch","Gap","E-value","Bitscore"), as.numeric)


## Modified LCA ----
# subset species
# return only hits that are identified to species level with at least 99% identity and 90% coverage
species <- blast_clean %>%
  filter(!is.na(Species) & Identity >= 99 & Coverage >= 90) %>%
  dplyr::group_by(ASV) %>%
  dplyr::filter(Bitscore > max(Bitscore)-(max(Bitscore)*0.05)) %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species) %>%
  ungroup()

# Filter by bitscore
# Filter only bitscores within 5% of the maximum bitscore for each ASV
not_species <- blast_clean %>%
  filter(!ASV %in% unique(species$ASV)) %>%
  filter(!is.na(Species) & Identity >= 90 & Coverage >= 90) %>%
  dplyr::group_by(ASV) %>%
  dplyr::filter(Bitscore > max(Bitscore)-(max(Bitscore)*0.05)) %>%
  select(ASV, Phylum, Class, Order, Family, Genus) %>%
  ungroup()

# Combine species and not_species for LCA
final <- bind_rows(species, not_species)

# Run LCA function
final_lca <- aggregate(final[, 2:ncol(final)], by = list(final$ASV), FUN = function(taxonomy_levels) {
  lca <- Reduce(intersect, taxonomy_levels)
  return(ifelse(length(lca) == 0, NA, paste(lca, collapse = ""))) 
}) 

# Rename the first column back to ASV
final_lca <- dplyr::rename(final_lca, ASV = Group.1)

# Output taxa list
final_taxa_list <- final_lca %>%
  select(Phylum, Class, Order, Family, Genus, Species) %>%
  distinct(.keep_all = FALSE)

# Return list of species with multiple good hits
problematic_species <- blast_clean %>%
  # select top hits for species
  filter(!is.na(Species) & Identity >= 99 & Coverage >= 95) %>%
  group_by(ASV) %>%
  # return only hits that have two different species per ASV
  distinct(Species, .keep_all = TRUE) %>%
  # return ONLY ASVs that have multiple "best" blast hits
  filter(n()>1)

