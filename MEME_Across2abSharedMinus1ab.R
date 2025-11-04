rm(list = ls())
lapply(c("dplyr", "tibble"), library, character.only = TRUE)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput")

Data <- read.csv("ChlData_fRanked.csv", header = TRUE)

# 1. Define genome names
ChlNDC <- "NDCt60"
Chl1q <- "Chl1q_t60"
Chl3q <- "Chl3q_t60"

# 2. Find SuperPos shared between ChlNDC_t60and Chl1q_t60
shared_superpos <- Data %>%
  filter(Genome %in% c(ChlNDC, Chl1q)) %>%
  group_by(SuperPos) %>%
  filter(n_distinct(Genome) == 2) %>%
  pull(SuperPos) %>%
  unique()

# 3. Identify SuperPos present in Chl3q
superpos_3q <- Data %>%
  filter(Genome == Chl3q) %>%
  pull(SuperPos) %>%
  unique()

# 4. Filter out those found in ChlNDC
final_superpos <- setdiff(shared_superpos, superpos_3q)

# 5. Select primary TSS entries matching filtered SuperPos (from ChlNDC or Chl1q_t60)
primary <- Data %>%
  filter(rank == "primary", TSSUTR != "", Locus == "", Genome %in% c(ChlNDC, Chl1q), SuperPos %in% final_superpos) %>%
  distinct(SuperPos, .keep_all = TRUE)

# 6. Select columns for FASTA output (adjust column numbers as needed)
fastaformat <- primary[, c(1, 2, 3, 4, 5, 8, 15)]

# 7. Add sequence number before SuperPos
fastaformat <- fastaformat %>%
  add_column(number = 1:nrow(fastaformat), .before = "SuperPos")

# 8. Duplicate rows for FASTA spacing
fastaformat2 <- as.data.frame(fastaformat[, 1])
fastaformat2[, 2:8] <- ""
colnames(fastaformat2) <- colnames(fastaformat)

# 9. Combine and arrange by number
fastaformat <- rbind(fastaformat, fastaformat2) %>%
  arrange(number)

# 10. Format headers with '>ID' prefix, remove spaces
fastaformat$number <- paste(">ID", fastaformat$number)
fastaformat$number <- gsub(" ", "", fastaformat$number)

# 11. Reset rownames
rownames(fastaformat) <- 1:nrow(fastaformat)

# 12. Clear headers where SuperPos is not empty
for (i in 1:nrow(fastaformat)) {
  if (fastaformat$SuperPos[i] != "") {
    fastaformat$number[i] <- ""
  }
}

# 13. Remove double quotes from headers and sequences
sequence_column <- "Sequence..50.nt.upstream...TSS..51nt."
fastaformat$number <- gsub("\"", "", fastaformat$number)
fastaformat[[sequence_column]] <- gsub("\"", "", fastaformat[[sequence_column]])

# 14. Merge headers and sequences into FASTA lines
fasta_lines <- character(nrow(fastaformat) * 2)
fasta_lines[c(TRUE, FALSE)] <- fastaformat$number
fasta_lines[c(FALSE, TRUE)] <- fastaformat[[sequence_column]]

# 15. Remove empty FASTA lines
fasta_lines <- fasta_lines[fasta_lines != ""]

# 16. Write FASTA lines to file
writeLines(fasta_lines, con = "~/Desktop/ChlNDC1Q_shared.fasta")
