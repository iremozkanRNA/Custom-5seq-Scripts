rm(list=ls())
lapply(c("dplyr", "UpSetR"), library, character.only=T)

fh1 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/ChlData_fRanked.csv") %>% filter(rank == "primary",Genome == "NDCt60", TSSUTR != "")
fh2 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/KsgData_fRanked.csv") %>% filter(rank == "primary", Genome == "NDCt60", TSSUTR != "")

tss1 <- unique(fh1$SuperPos)
tss2 <- unique((fh2$SuperPos))

tss_list <- list(Chl = tss1, Ksg = tss2)

upset( fromList(tss_list),
       sets = c("Chl", "Ksg"),
       order.by="degree",
       sets.bar.color = c("#108080", "#7F0440"))

