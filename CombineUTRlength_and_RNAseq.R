rm(list = ls())

library(dplyr)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")

data <- read.csv("KsgData_fRanked_metadata_UTRlength.csv")
rnaseq1Q <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/RNA-Seq Feb 24,2022/ShinyomicsDerivedSignidicanceFiles/T4_Ksg_RNAseq_1Q.csv")
rnaseq3Q <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/RNA-Seq Feb 24,2022/ShinyomicsDerivedSignidicanceFiles/T4_Ksg_RNAseq_3Q.csv")


rnaseq3Q_subset <- rnaseq3Q %>% select(Gene, Value)
data_subset <- data %>% subset(Genome == "Ksg3q_t60")
data_subset <-  data_subset %>% left_join(rnaseq3Q_subset, by = c("TSSUTR" = "Gene"))

rnaseq1Q_subset <- rnaseq1Q %>% select(Gene, Value)
data_subset1 <-  data %>% subset(Genome == "Ksg1q_t60")
data_subset1 <-  data_subset1 %>% left_join(rnaseq1Q_subset, by = c("TSSUTR" = "Gene"))

data_subset2 <-  data %>% subset(Genome == "NDCt60")

data <- bind_rows(data_subset, data_subset1, data_subset2) %>% arrange(SuperPos)
write.csv(data,"KsgData_fRanked_metadata_UTRlength_RNAseqDEseqadded.csv", row.names = F, col.names = T)
