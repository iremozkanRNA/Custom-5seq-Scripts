rm(list=ls())
lapply(c("dplyr", "tibble"), library, character.only = T)

fh1 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/KsgData_fRanked.csv")
fh1 <- fh1 %>% add_column(top.antisense = "", comp.antisense = "",.after = "termUTR")

genes <- read.delim("Stuff/AMeyer_Lab/Experiment_results/3seq/RNAFOLD/RNAfold_DRIVE/Spn/Chl/Nofilter/NC_003028.v3.17.ncrna.genes", header = T, sep = "")

colnames(genes) <- c("genome", "from", "to", "strand", "name", "old.name", "length", "protein.Name", "symbol")

genes.top <- genes %>% filter(strand == "+")
genes.comp <- genes %>% filter(strand == "-")

top <- fh1 %>% filter(SuperStrand == "+")
comp <- fh1 %>% filter(SuperStrand == "-")
cols_to_numeric <- c("from", "to")
# Loop through the selected columns and convert them to numeric
genes.comp[cols_to_numeric] <- lapply(genes.comp[cols_to_numeric], as.numeric)
genes.top[cols_to_numeric] <- lapply(genes.top[cols_to_numeric], as.numeric)

# Look for antisense
for (i in 1:nrow(top)) {
  for (k in 1:nrow(genes.comp)) {
    if (!is.na(genes.comp$to[k]) && !is.na(genes.comp$from[k + 1])) {
      #first annotate 5'UTR
      if (
        top$SuperPos[i] < genes.comp$from[k+1] + 50 &&
        top$SuperPos[i] > genes.comp$to[k] - 10 &&
        top$SuperPos[i] <= genes.comp$to[k] + 550
      ) {
        top$comp.antisense[i] <- "5'UTR"
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      } else {
        #Annotate everything within genes as genic if the distance is <500 bp
        if (
          top$SuperPos[i] > genes.comp$from[k] &&
          top$SuperPos[i] < genes.comp$to[k]
        ) {
          top$comp.antisense[i] <- "genic"
          break  # Exit the inner loop once a 5' UTR annotation is found
        }
      }
    }
}


# Look for antisense
for (i in 1:nrow(comp)) {
  for (k in 1:nrow(genes.top)) {
    if (!is.na(genes.top$to[k]) && !is.na(genes.top$from[k + 1])) {
      #first annotate 5'UTR
      if (
        comp$SuperPos[i] < genes.top$from[k+1] + 10 &&
        comp$SuperPos[i] > genes.top$to[k] - 50 &&
        comp$SuperPos[i] <= genes.top$from[k] - 550
      ) {
        comp$top.antisense[i] <- "5'UTR"
        break  # Exit the inner loop once a 5' UTR annotation is found
      }
    } else {
      #Annotate everything within genes as genic if the distance is <500 bp
      if (
        comp$SuperPos[i] > genes.top$from[k] &&
        comp$SuperPos[i] < genes.top$to[k]
      ) {
        comp$top.antisense[i] <- "genic"
        break  # Exit the inner loop once genic annotation is found
      }
    }
  }
}

fh1 <- rbind(top,comp) %>% arrange(SuperPos)

# Fixed antisense checks (handles both NA and empty strings)
fh1$has_antisense <- as.integer(
  (nzchar(fh1$top.antisense, keepNA = TRUE) | 
     nzchar(fh1$comp.antisense, keepNA = TRUE)) 
)

# Fixed combined column (ignores empty strings)
fh1$combined_antisense <- ifelse(
  nzchar(fh1$top.antisense, keepNA = TRUE),
  fh1$top.antisense,
  fh1$comp.antisense
)
write.csv(fh1, file = "~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/KsgData_fRanked_antisenseKsg.csv")

