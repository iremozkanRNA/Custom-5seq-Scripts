rm(list = ls())
lapply(c("dplyr", "ggplot2", "purrr", "ggVennDiagram"), library, character.only = T )

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")

Data <- read.csv("ChlData_fRanked.csv", header = T) %>% 
  filter(rank == "primary")

list2env(set_names(lapply(c("NDCt60", "Chl1q_t60", "Chl3q_t60"), function(df) 
  Data %>% filter(Genome == df)), 
  c("NDC", "Chl1Q", "Chl3Q")), envir = .GlobalEnv)

#Function to create Venn Diagram
create_venn <- function(NDC, Chl1Q, Chl3Q) {
  
  #Create list for Venn Diagram
  venn_list <- list(NDC = NDC$SuperPos, Chl1Q = Chl1Q$SuperPos, Chl3Q = Chl3Q$SuperPos) 

  #Create Venn diagram
  venn_plot <- ggVennDiagram(venn_list) +
    labs(title = paste("TSS Primary Overlap Between Antibiotic Conditions by Position"))
  
  return(venn_plot)
}
venn <- create_venn(NDC, Chl1Q, Chl3Q)
print(venn)

