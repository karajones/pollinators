# This script is modified from: https://pr2database.github.io/pr2database/articles/pr2_04_decipher.html

#### ----

library(DECIPHER)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(rio)
library(DT)

# Read data -------------------------------------------------------------

seqs <- readDNAStringSet("CO1_pollinator_final_database.fa")
seqs <- readDNAStringSet("16S_pollinator_final_database.dedup.fa")


# Taxonomic groups -------------------------------------------------------------

maxGroupSize <- Inf # max sequences per label (>= 1)
allowGroupRemoval <- FALSE

# 10 for 16s, 12 for CO1
maxIterations <- 10 # must be >= 1

# obtain the taxonomic assignments
groups <- names(seqs) # sequence names

# All taxos must start with "Root;"
groups <- str_c("Root;",groups)
names(seqs)<-groups

groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups
cat("Number of groups: ", length(u_groups), '\n') # number of group

taxid <- NULL

# Pruning group size ---

remove <- logical(length(seqs))

for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index), maxGroupSize)
  remove[index[-keep]] <- TRUE
}

cat("Number of sequences eliminated: ", sum(remove), "\n")

# Iteratively train classifier -----------------------------------------------

probSeqsPrev <- integer() # suspected problem sequences from prior iteration
df.problems <- list()

cat("Number of iterations:", maxIterations, "\n", sep=" ")

for (i in seq_len(maxIterations)) {
  
  cat("Training iteration: ", i, "\n", sep="")
  
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove], names(seqs)[!remove],taxid)
  
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  
  cat("Number of problem sequences: ", length(probSeqs), "\n", sep="")
  
  
  # Exit if no more problem sequences or same problems as previous or reach max Iter
  
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) && all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  
  if (i==maxIterations)
    break
  
  
  # remove any problem sequences
  
  probSeqsPrev <- probSeqs
  
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  
  df.problems[[i]] <- data.frame(index, trainingSet$problemSequences)
  
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}

cat("Total number of sequences eliminated: ", sum(remove), "\n") 
cat("Number of remaining problem sequences: ", length(probSeqs), "\n")

df.problems <- purrr::reduce(df.problems, bind_rows) %>% 
  select(-Index)

write.table(df.problems, "16S_pollinator_database_problem_seqs.tsv", sep="\t", quote=F, row.names = FALSE, col.names=TRUE)

save.image("16S_pollinator_database_idtaxa.RData")

### END DECIPHER ----
