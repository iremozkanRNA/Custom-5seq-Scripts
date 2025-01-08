# Clear the environment
rm(list = ls())

# Load necessary libraries
lapply(c("dplyr","grid","fuzzyjoin"), library, character.only = TRUE)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/")
# Read in the datasets
set1 <- read.csv("ChlData_fRanked.csv") %>%
  filter(rank == "secondary", Genome == "NDCt60")

set2 <- read.csv("KsgData_fRanked.csv") %>%
  filter(rank == "secondary", Genome == "NDCt60")

# Separate into + and - strands for both datasets
Chlpos <- set1 %>% filter(SuperStrand == "+")
Chlneg <- set1 %>% filter(SuperStrand == "-")
Ksgpos <- set2 %>% filter(SuperStrand == "+")
Ksgneg <- set2 %>% filter(SuperStrand == "-")

# Perform an intersection with a tolerance of ±10
common <- difference_inner_join(
  Chlpos, Ksgpos,
  by = c("SuperPos" = "SuperPos"),
  max_dist = 20
)

# Convert the result to a data frame if needed
common <- as.data.frame(common)

# View the result
print(common)

common2 <- difference_inner_join(
  Chlneg, Ksgneg,
  by = c("SuperPos" = "SuperPos"),
  max_dist = 20
)

common2 <- as.data.frame(common2)
commoncomb <- rbind(common,common2) %>% arrange("SuperPos")

Chlunique <- nrow(set1) - nrow(commoncomb)
Ksgunique <- nrow(set2) - nrow(commoncomb)
ChluniquePerc <-  Chlunique/(Chlunique + nrow(commoncomb) + Ksgunique)*100
KsguniquePerc <- Ksgunique/(Chlunique + nrow(commoncomb) + Ksgunique)*100
commoncombPerc <- nrow(commoncomb)/(Chlunique + nrow(commoncomb) + Ksgunique)*100
