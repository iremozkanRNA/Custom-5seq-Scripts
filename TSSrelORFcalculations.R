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
#   Retrieve relORF positions - TOP Strand
############################################ 

graph <- as.data.frame(matrix(ncol = 4, nrow = length(seq(0,100, by=2))))
colnames(graph) <- c("ORFpos", "TSSnumberNDC", "TSSnumber1Q", "TSSnumber3Q")
graph$ORFpos <- seq(from = 0, to=100, by = 2)

# Initialize columns in graph
graph[c(1:51),c(2:4)] <- numeric(nrow(graph))

# Filter top for the specific genome
top_NDCt60 <- top[top$Genome == "NDCt60", ]

for (i in 1:(nrow(graph) - 1)) {
  n <- 0
  for (j in 1:nrow(top_NDCt60)) {
    k <- ifelse(top_NDCt60$relORFpos[j] > graph$ORFpos[i] &
                  top_NDCt60$relORFpos[j] <= graph$ORFpos[i + 1], 1, 0)
    n <- n + k
  }
  graph$TSSnumberNDC[i] <- n
}

top1Q <- top[top$Genome == "Ksg1q_t60",]

for (i in 1:(nrow(graph) - 1)) {
  n <- 0
  for (j in 1:nrow(top1Q)) {
    k <- ifelse(top1Q$relORFpos[j] > graph$ORFpos[i] &
                  top1Q$relORFpos[j] <= graph$ORFpos[i + 1], 1, 0)
    n <- n + k
  }
  graph$TSSnumber1Q[i] <- n
}

top3Q <- top[top$Genome == "Ksg3q_t60",]

for (i in 1:(nrow(graph) - 1)) {
  n <- 0
  for (j in 1:nrow(top1Q)) {
    k <- ifelse(top1Q$relORFpos[j] > graph$ORFpos[i] &
                  top1Q$relORFpos[j] <= graph$ORFpos[i + 1], 1, 0)
    n <- n + k
  }
  graph$TSSnumber3Q[i] <- n
}

graphtop <- graph

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

############################################
#   Retrieve relORF positions - COMP Strand
############################################ 

graph[c(1:51), c(2:4)] <- numeric(nrow(graph))

compNDC <- comp[comp$Genome == "NDCt60",]

for(i in 1:nrow(graph) - 1){
  n <- 0
  for(j in 1:nrow(compNDC)){
    k <- ifelse(compNDC$relORFpos[j] > graph$ORFpos[i] &
                  compNDC$relORFpos[j] <= graph$ORFpos[i+1],1,0)
    n <- n + k
  }
  graph$TSSnumberNDC[i] <- n
}

comp1Q <- comp[comp$Genome == "Ksg1q_t60",]

for(i in 1:nrow(graph) - 1){
  n <- 0
  for(j in 1:nrow(comp1Q)){
    k <- ifelse(comp1Q$relORFpos[j] > graph$ORFpos[i] &
                  comp1Q$relORFpos[j] <= graph$ORFpos[i+1],1,0)
    n <- n + k
  }
  graph$TSSnumber1Q[i] <- n
}

comp3Q <- comp[comp$Genome == "Ksg3q_t60",]

for(i in 1:nrow(graph) - 1){
  n <- 0
  for(j in 1:nrow(comp3Q)){
    k <- ifelse(comp3Q$relORFpos[j] > graph$ORFpos[i] &
                  comp3Q$relORFpos[j] <= graph$ORFpos[i+1],1,0)
    n <- n + k
  }
  graph$TSSnumber3Q[i] <- n
}

graphcomp <- graph

############################################
#  Combine Strands ORF positions
############################################ 

peaks <- rbind(top, comp) %>% arrange(SuperPos)
write.csv(peaks, file = "Chl_genic_relORF.csv", row.names = F)

graph <- as.data.frame(matrix(ncol = 10, nrow = length(seq(0,100, by = 2))))
graph[c(1:51), c(1)] <- seq(0,100, by=2)
graphcomp <- graphcomp[,-c(1)]
comb <- cbind(graphtop,graphcomp)
graph[c(1:51),c(1:7)] <- comb
colnames(graph) <- c("ORFpos", "TSS_NDCtop", "TSS_1Qtop", "TSS_3Qtop",
                     "TSS_NDCcomp", "TSS_1Qcomp", "TSS_3Qcomp", 
                     "TSS_NDC", "TSS_1Q", "TSS_3Q")

graph$TSS_NDC <- graph$TSS_NDCtop + graph$TSS_NDCcomp
graph$TSS_1Q <- graph$TSS_1Qtop + graph$TSS_1Qcomp
graph$TSS_3Q <- graph$TSS_3Qtop + graph$TSS_3Qcomp

write.csv(graph, file = "Ksg_relORF_BinOf2_primary.csv", col.names = T, row.names = F)
