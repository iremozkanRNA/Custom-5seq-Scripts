# Clear environment and load required packages
rm(list = ls())
lapply(c("dplyr", "tibble", "ggVennDiagram"), library, character.only = TRUE)

# Set working directory
setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")

# Read data
Data <- read.csv("ChlData_fRanked.csv", header = TRUE)
Data2 <- read.csv("KsgData_fRanked.csv", header = TRUE)

#######################################################
################ TSS MOTIF Analysis ###################
#######################################################

# Filter for primary, genic TSSs
primary <- Data %>% filter(rank == "primary", TSSUTR != "" & Locus == "", Genome == "NDCt60")
primary2 <- Data2 %>% filter(rank == "primary", TSSUTR != "" & Locus == "", Genome == "NDCt60")

# Standardize Genome names
primary$Genome[primary$Genome == "NDCt60"] <- "ChlNDCt60"
primary2$Genome[primary2$Genome == "NDCt60"] <- "KsgNDC60"

# Find SuperPos unique to each dataset
unique_SuperPos1 <- setdiff(primary$SuperPos, primary2$SuperPos) # Unique to Data
unique_SuperPos2 <- setdiff(primary2$SuperPos, primary$SuperPos) # Unique to Data2

# --- For ChlData (Data) unique SuperPos ---
primary_unique <- primary %>% filter(SuperPos %in% unique_SuperPos1)
fastaformat1 <- primary_unique[,c(1,2,3,4,5,8,15)]
fastaformat1 <- fastaformat1 %>% 
  add_column(number = 1:nrow(fastaformat1), .before = "SuperPos")
fastaformat2_1 <- as.data.frame(fastaformat1[,c(1)]) 
fastaformat2_1[,c(2:8)] <- ""
colnames(fastaformat2_1) <- colnames(fastaformat1)
fastaformat1 <- rbind(fastaformat1, fastaformat2_1)
fastaformat1 <- fastaformat1 %>% arrange(number)
fastaformat1$number <- paste(">ID", fastaformat1$number)
fastaformat1$number <- sub(" ", "", fastaformat1$number)
rownames(fastaformat1) <- 1:nrow(fastaformat1)
for(i in 1:nrow(fastaformat1)){
  if(!fastaformat1$SuperPos[i] == ""){
    fastaformat1$number[i] <- ""
  }
}
fastaformat1$Sequence..50.nt.upstream...TSS..51nt. <- paste(fastaformat1$number, fastaformat1$Sequence..50.nt.upstream...TSS..51nt.)
fastaformat1 <- fastaformat1[,-c(1)]
rownames(fastaformat1) <- 1:nrow(fastaformat1)
fastaformatPOS1 <- fastaformat1
fastaformatPOS1 <- fastaformat1[-c(nrow(fastaformatPOS1)),]
write.csv(fastaformatPOS1, file = "~/Desktop/Manuscript/Tables/MEME_Chl_primary_TSSUTR_SuperPos_UNIQUE.csv", row.names = FALSE, quote = FALSE)

# --- For KsgData (Data2) unique SuperPos ---
primary2_unique <- primary2 %>% filter(SuperPos %in% unique_SuperPos2)
fastaformat2 <- primary2_unique[,c(1,2,3,4,5,8,15)]
fastaformat2 <- fastaformat2 %>% 
  add_column(number = 1:nrow(fastaformat2), .before = "SuperPos")
fastaformat2_2 <- as.data.frame(fastaformat2[,c(1)]) 
fastaformat2_2[,c(2:8)] <- ""
colnames(fastaformat2_2) <- colnames(fastaformat2)
fastaformat2 <- rbind(fastaformat2, fastaformat2_2)
fastaformat2 <- fastaformat2 %>% arrange(number)
fastaformat2$number <- paste(">ID", fastaformat2$number)
fastaformat2$number <- sub(" ", "", fastaformat2$number)
rownames(fastaformat2) <- 1:nrow(fastaformat2)
for(i in 1:nrow(fastaformat2)){
  if(!fastaformat2$SuperPos[i] == ""){
    fastaformat2$number[i] <- ""
  }
}
fastaformat2$Sequence..50.nt.upstream...TSS..51nt. <- paste(fastaformat2$number, fastaformat2$Sequence..50.nt.upstream...TSS..51nt.)
fastaformat2 <- fastaformat2[,-c(1)]
rownames(fastaformat2) <- 1:nrow(fastaformat2)
fastaformatPOS2 <- fastaformat2
fastaformatPOS2 <- fastaformat2[-c(nrow(fastaformatPOS2)),]
write.csv(fastaformatPOS2, file = "~/Desktop/Manuscript/Tables/MEME_Ksg_primary_TSSUTR_SuperPos_UNIQUE.csv", row.names = FALSE, quote = FALSE)
