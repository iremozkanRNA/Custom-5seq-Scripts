rm(list = ls())
lapply(c("dplyr", "tibble"), library, character.only=T)
####################################
#Read your TSSpreditor file output.
####################################
setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")
fRank <- read.csv("KsgData_fRanked_metadata.csv") %>% filter(rank == "primary")
top <- fRank %>% filter((SuperStrand == "+"))# & (Genome != "NDCt60"))
comp <- fRank %>% filter(SuperStrand == "-" ) #& (Genome != "NDCt60"))

####################################
#Read your annotation file
#################################### 
# Read your annotation file
genes <- read.delim("NC_003028.v3.17.ncrna.genes", header = F)
genes <- genes[complete.cases(genes),]
genes.top <- genes %>% filter(V4 == "+")
colnames(genes.top) <- c("genome", "from", "to", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )
genes.comp <- genes %>% filter(V4 == "-")
colnames(genes.comp) <- c("genome", "to", "from", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )

# List of columns to convert to numeric
cols_to_numeric <- c("from", "to")
# Loop through the selected columns and convert them to numeric
genes.comp[cols_to_numeric] <- lapply(genes.comp[cols_to_numeric], as.numeric)
genes.top[cols_to_numeric] <- lapply(genes.top[cols_to_numeric], as.numeric)
comp$SuperPos <- as.numeric(comp$SuperPos)
top$SuperPos <- as.numeric(top$SuperPos)
############################################
#   Calculate relative ORF positions - TOP strand
############################################ 

top <- top %>% add_column(relORFpos = "", .after = "SuperPos")
genes.top <- genes.top %>% add_column(length = genes.top$to - genes.top$from, .after = "to")
for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if(top$SuperPos[i] > genes.top$from[j] &
       top$SuperPos[i] <= genes.top$to[j]){
      top$relORFpos[i] <- (top$SuperPos[i] - genes.top$from[j])*100/genes.top$length[j]
    }
  }
}
top$relORFpos <- as.numeric(top$relORFpos)
top <- top %>% filter(
  (!is.na(top$relORFpos) & top$TSSUTR == "")) 

############################################
#   Calculate relative ORF positions - COMP strand
############################################ 

comp <- comp %>% add_column(relORFpos = "", .after = "SuperPos")
genes.comp <- genes.comp %>% add_column(length = genes.comp$from - genes.comp$to, .after = "to")

for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if(comp$SuperPos[i] > genes.comp$to[j] &
       comp$SuperPos[i] <= genes.comp$from[j]){
      comp$relORFpos[i] <- (genes.comp$from[j] - comp$SuperPos[i])*100/genes.comp$length[j]
    }
  }
}
comp$relORFpos <- as.numeric(comp$relORFpos)
comp <- comp %>% filter(!is.na(comp$relORFpos) & comp$TSSUTR == "")


# Combine top comp strands
 
fh <- rbind(top, comp) %>% distinct(Locus, .keep_all = T) 

fh <- fh %>% filter(relORFpos >= 50, Genome != "NDCt60")

write.csv(fh, file="Ksgfiltered_relORF.csv", row.names = F, col.names = T)
