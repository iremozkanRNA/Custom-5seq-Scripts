rm(list = ls())
lapply(c("dplyr", "tibble"), library, character.only = TRUE)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput")

Data <- read.csv("ChlData_fRanked.csv", header = TRUE)

# Extract SuperPos values from NDCt60 and Ksg1q_t60 genomes
exclude_superpos <- Data %>%
  filter(Genome %in% c("Chl1q_t60", "NDCt60")) %>%
  pull(SuperPos) %>%
  unique()

# Filter primary TSS for Ksg1q_t60, exclude SuperPos in other genomes
primary <- Data %>%
  filter(rank == "primary", TSSUTR != "", Locus == "", Genome == "Chl3q_t60") %>%
  filter(!SuperPos %in% exclude_superpos)%>%
  distinct(SuperPos, .keep_all = TRUE) 

# Select columns for FASTA output
fastaformat <- primary[, c(1, 2, 3, 4, 5, 8, 15)]

# Add sequence number before SuperPos column
fastaformat <- fastaformat %>%
  add_column(number = 1:nrow(fastaformat), .before = "SuperPos")

# Duplicate rows for FASTA format spacing
fastaformat2 <- as.data.frame(fastaformat[, 1])
fastaformat2[, 2:8] <- ""
colnames(fastaformat2) <- colnames(fastaformat)

# Combine and arrange by number
fastaformat <- rbind(fastaformat, fastaformat2) %>%
  arrange(number)

# Format headers with ">ID" prefix and remove spaces
fastaformat$number <- paste(">ID", fastaformat$number)
fastaformat$number <- gsub(" ", "", fastaformat$number)

# Reset row names
rownames(fastaformat) <- 1:nrow(fastaformat)

# Clear headers for rows where SuperPos is not empty
for (i in 1:nrow(fastaformat)) {
  if (fastaformat$SuperPos[i] != "") {
    fastaformat$number[i] <- ""
  }
}

# Remove double quotes from headers and sequences
sequence_column <- "Sequence..50.nt.upstream...TSS..51nt."
fastaformat$number <- gsub("\"", "", fastaformat$number)
fastaformat[[sequence_column]] <- gsub("\"", "", fastaformat[[sequence_column]])

# Merge headers and sequences into FASTA lines
fasta_lines <- character(nrow(fastaformat) * 2)
fasta_lines[c(TRUE, FALSE)] <- fastaformat$number
fasta_lines[c(FALSE, TRUE)] <- fastaformat[[sequence_column]]

# Remove empty FASTA lines
fasta_lines <- fasta_lines[fasta_lines != ""]

# Write FASTA lines to file
writeLines(fasta_lines, con = "~/Desktop/Chl3q.fasta")

