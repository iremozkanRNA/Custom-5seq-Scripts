# Load required libraries
rm(list = ls())
lapply(c("Biostrings", "dplyr", "tidyr"), library, character.only = TRUE)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")

# Function to reverse complement a sequence
strReverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
}

# Read and process FASTA file
temp <- readDNAStringSet(filepath = "TIGR4.fasta")
fasta <- data.frame(
  seqname = names(temp),
  sequence = paste(temp)
)

# Function to process data and extract sequences
process_data <- function(data_file, genome_prefix) {
  # Read data
  data <- read.csv(file = data_file)
  data <- data %>% add_column(Sequence_200nt = "", TSS_1nt = "")
  
  # Separate strands
  datatop <- data %>% filter(SuperStrand == "+")
  datacomp <- data %>% filter(SuperStrand == "-")
  
  # Extract sequences for + strand
  for (i in seq_along(datatop$SuperPos)) {
    datatop$Sequence_200nt[i] <- substr(fasta$sequence, datatop$SuperPos[i] - 100, datatop$SuperPos[i] + 100)
    datatop$TSS_1nt[i] <- substr(fasta$sequence, datatop$SuperPos[i] + 1, datatop$SuperPos[i] + 1)
  }
  
  # Extract sequences for - strand
  for (j in seq_along(datacomp$SuperPos)) {
    datacomp$Sequence_200nt[j] <- substr(fasta$sequence, datacomp$SuperPos[j] + 100, datacomp$SuperPos[j] - 100)
    datacomp$TSS_1nt[j] <- substr(fasta$sequence, datacomp$SuperPos[j], datacomp$SuperPos[j])
  }
  
  # Reverse complement the - strand sequences
  datacomp$Sequence_200nt <- chartr("ATCG", "TAGC", datacomp$Sequence_200nt)
  datacomp$Sequence_200nt <- strReverse(datacomp$Sequence_200nt)
  
  # Combine and adjust genome names
  data <- rbind(datatop, datacomp) %>% arrange(SuperPos)
  
  for (i in seq_along(data$Genome)) {
    if (data$Genome[i] == "NDCt60") {
      data$Genome[i] <- paste0(genome_prefix, "NDCt60")
    }
  }
  
  return(data)
}

# Process Chl and Ksg datasets
Chl <- process_data("ChlData_fRanked_metadata_UTRlength_RNAseqDEseqadded.csv", "Chl")
Ksg <- process_data("KsgData_fRanked_metadata_UTRlength_RNAseqDEseqadded.csv", "Ksg")

# Combine datasets
data <- rbind(Chl, Ksg) %>% arrange(SuperPos)

# Filter rows where TSSUTR is not empty
filtered_data <- data %>% filter(TSSUTR != "")

# Function to calculate nucleotide distribution
calculate_nt_distribution <- function(data) {
  # Split sequences into individual nucleotides
  nt_matrix <- do.call(rbind, strsplit(data$Sequence_200nt, ""))
  
  # Create a dataframe for positions and nucleotides
  nt_distribution <- data.frame(
    Sample = rep(1:nrow(nt_matrix), each = ncol(nt_matrix)),
    Position = rep(seq(-100, 100), times = nrow(nt_matrix)),
    Nucleotide = as.vector(t(nt_matrix))
  )
  
  # Calculate nucleotide counts at each position
  nt_summary <- nt_distribution %>%
    group_by(Position, Nucleotide) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Nucleotide, values_from = Count, values_fill = list(Count = 0))
  
  return(nt_summary)
}

# Generate nucleotide distribution for filtered data
nt_distribution_filtered <- calculate_nt_distribution(filtered_data)

# Write outputs to files
write.csv(filtered_data, file = "Filtered_Data_TSSUTR_NotEmpty_with_200nt.csv", row.names = FALSE)
write.csv(nt_distribution_filtered, file = "Nucleotide_Distribution_TSSUTR_NotEmpty_-100_to_+100.csv", row.names = FALSE)

