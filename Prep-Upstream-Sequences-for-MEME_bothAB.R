rm(list = ls())
lapply(c("dplyr", "tibble", "ggVennDiagram"), library ,character.only=T)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")

Data <- read.csv("ChlData_fRanked.csv", header = T)
Data2 <- read.csv("KsgData_fRanked.csv", header = T)
#######################################################

################ TSS MOTIF Analysis #################

#######################################################

primary <- Data %>% filter(rank == "primary",  TSSUTR == "" & Locus != "" , Genome =="NDCt60")
primary$Genome[primary$Genome == "NDCt60"] <- "ChlNDCt60"
primary2 <- Data2 %>% filter(rank == "primary",  TSSUTR == "" & Locus != "" ,Genome =="NDCt60" )
primary2$Genome[primary2$Genome == "NDCt60"] <- "KsgNDC60"
primary <- rbind(primary, primary2)
primary <- primary[!duplicated(primary$SuperPos),]
fastaformat <- primary[,c(1,2,3,4,5,8,15)]
fastaformat <- fastaformat %>% 
  add_column(number = 1:nrow(fastaformat), .before = "SuperPos")
fastaformat2 <- as.data.frame(fastaformat[,c(1)]) 
fastaformat2[,c(2:8)] <- ""
colnames(fastaformat2) <- colnames(fastaformat)
fastaformat <- rbind(fastaformat, fastaformat2)
fastaformat <- fastaformat %>% arrange(number)
fastaformat$number <- paste( ">ID", fastaformat$number)
fastaformat$number <- sub(" ", "", fastaformat$number)
rownames(fastaformat) <- 1:nrow(fastaformat)

for(i in 1:nrow(fastaformat)){
  if(!fastaformat$SuperPos[i] == ""){
    fastaformat$number[i] <- ""
  }
}

fastaformat$Sequence..50.nt.upstream...TSS..51nt. <- paste(fastaformat$number, fastaformat$Sequence..50.nt.upstream...TSS..51nt.)
fastaformat <- fastaformat[,-c(1)]
rownames(fastaformat) <- 1:nrow(fastaformat)
fastaformatPOS <- fastaformat
fastaformatPOS <- fastaformat[-c(nrow(fastaformatPOS)),]

write.csv(fastaformatPOS, file = "~/Desktop/Manuscript/Tables/Table_S4_MEME_ChlKsg_primaryTSS_Locus_NDC.csv", sep = "", row.names = F, col.names = F)
