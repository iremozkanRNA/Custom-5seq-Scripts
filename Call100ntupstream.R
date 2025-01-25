rm(list=ls())
lapply(c("Biostrings", "dplyr"), library, character.only=T)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")

# Arrange your fasta
temp <- readDNAStringSet(filepath = "TIGR4.fasta")
fasta <- data.frame(
  seqname = names(temp),
  sequence = paste(temp)
)

data <- read.csv(file = "ChlData_fRanked_metadata_UTRlength_RNAseqDEseqadded.csv")
data <- data %>% add_column(Sequence...100nt..upstream = "", TSS..1nt.. = "")
datatop <- data %>% filter(SuperStrand == "+")
datacomp <- data %>% filter(SuperStrand == "-")

# Call 100 nt upstream
for(i in seq_along(datatop$SuperPos)){
  for(j in seq_along(datacomp$SuperPos)){
    datatop$Sequence...100nt..upstream[i] <- substr(fasta$sequence, datatop$SuperPos[i] - 100, datatop$SuperPos[i])
    datatop$TSS..1nt..[i] <- substr(fasta$sequence, datatop$SuperPos[i] +1, datatop$SuperPos[i]+1)
    datacomp$Sequence...100nt..upstream[j] <- substr(fasta$sequence, datacomp$SuperPos[j], datacomp$SuperPos[j] + 100)
    datacomp$TSS..1nt..[j] <- substr(fasta$sequence, datacomp$SuperPos[j]-1, datacomp$SuperPos[j]-1)
  }
}
#reverse complement the - strand
datacomp$Sequence...100nt..upstream <- chartr("ATCG", "TAGC", datacomp$Sequence...100nt..upstream)
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL),rev),paste, collapse="")

datacomp$Sequence...100nt..upstream <- strReverse(datacomp$Sequence...100nt..upstream)
datacomp$TSS..1nt.. <- strReverse(datacomp$TSS..1nt..)
data <- rbind(datatop, datacomp) %>% arrange(SuperPos)

for(i in seq_along(data$Genome)){
  if(data$Genome[i] == "NDCt60"){
    data$Genome[i] <- "ChlNDCt60"
  }
}

Chl <- data

data <- read.csv(file = "KsgData_fRanked_metadata_UTRlength_RNAseqDEseqadded.csv")
data <- data %>% add_column(Sequence...100nt..upstream = "", TSS..1nt.. = "")
datatop <- data %>% filter(SuperStrand == "+")
datacomp <- data %>% filter(SuperStrand == "-")

# Call 100 nt upstream
for(i in seq_along(datatop$SuperPos)){
  for(j in seq_along(datacomp$SuperPos)){
    datatop$Sequence...100nt..upstream[i] <- substr(fasta$sequence, datatop$SuperPos[i] - 100, datatop$SuperPos[i])
    datatop$TSS..1nt..[i] <- substr(fasta$sequence, datatop$SuperPos[i] + 1, datatop$SuperPos[i] +1)
    datacomp$Sequence...100nt..upstream[j] <- substr(fasta$sequence, datacomp$SuperPos[j], datacomp$SuperPos[j] + 100)
    datacomp$TSS..1nt..[j] <- substr(fasta$sequence, datacomp$SuperPos[j]-1, datacomp$SuperPos[j]-1)
  }
}
#reverse complement the - strand
datacomp$Sequence...100nt..upstream <- chartr("ATCG", "TAGC", datacomp$Sequence...100nt..upstream)
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL),rev),paste, collapse="")

datacomp$Sequence...100nt..upstream <- strReverse(datacomp$Sequence...100nt..upstream)
datacomp$TSS..1nt.. <- strReverse(datacomp$TSS..1nt..)
data <- rbind(datatop, datacomp) %>% arrange(SuperPos)

for(i in seq_along(data$Genome)){
  if(data$Genome[i] == "NDCt60"){
    data$Genome[i] <- "KsgNDCt60"
  }
}

data <- rbind(Chl, data) %>% arrange(SuperPos)

write.csv(data, file = "ChlKsg_Data_fRanked_metadata_UTRlength_RNAseqDEseqadded_100ntupstream.csv", row.names = F, col.names = T)
