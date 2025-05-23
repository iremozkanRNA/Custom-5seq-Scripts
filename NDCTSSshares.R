rm(list=ls())
lapply(c("dplyr", "UpSetR"), library, character.only=T)

fh1 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/ChlData_fRanked.csv") %>% filter(rank == "primary",Genome == "NDCt60")
fh2 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/KsgData_fRanked.csv") %>% filter(rank == "primary", Genome == "NDCt60")

tss1 <- unique(fh1$Locus)
tss2 <- unique((fh2$Locus))

tss_list <- list(Chl = tss1, Ksg = tss2)

upset( fromList(tss_list),
       sets = c("Chl", "Ksg"),
       order.by="freq",
       sets.bar.color = c("lightgreen", "pink"))
