rm(list = ls())
lapply(c("dplyr", "tibble", "ggVennDiagram"), library ,character.only=T)

setwd("~/path/to/your/file/")

Data <- read.csv("ChlData_fRanked.csv", header = T)

#######################################################

################ TSS MOTIF Analysis #################

#######################################################

primary <- Data %>% filter(rank == "secondary",  TSSUTR != "" & Locus == "" & Genome== "NDCt60")

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

write.csv(fastaformatPOS, file = "yourfilename.csv", sep = "", row.names = F, col.names = F)
