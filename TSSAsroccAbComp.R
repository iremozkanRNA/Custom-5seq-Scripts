# Clean environment and load libraries
rm(list=ls())
lapply(c("dplyr", "UpSetR"), library, character.only = TRUE)

# Read in data
fh1 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/ChlData_fRanked.csv")
fh2 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/KsgData_fRanked.csv")

# Filter for each genome and extract unique SuperPos
ndc1 <- fh1 %>% filter(rank == "primary", Genome == "NDCt60", TSSUTR != "")
ndc2 <- fh2 %>% filter(rank == "primary", Genome == "NDCt60", TSSUTR != "")
chl1 <- fh1 %>% filter(rank == "primary", Genome == "Chl3q_t60", TSSUTR != "")
ksg1 <- fh2 %>% filter(rank == "primary", Genome == "Ksg3q_t60", TSSUTR != "")

superpos_ndc1 <- unique(ndc1$SuperPos)
superpos_ndc2 <- unique(ndc2$SuperPos)
superpos_chl1 <- unique(chl1$SuperPos)
superpos_ksg1 <- unique(ksg1$SuperPos)

# Combine NDC SuperPos as NDCall
NDCall <- union(superpos_ndc1, superpos_ndc2)

# Get all unique SuperPos across all sets
all_superpos <- union(union(NDCall, superpos_chl1), superpos_ksg1)

# Build the binary matrix for UpSetR
ndc_matrix <- data.frame(
  SuperPos = all_superpos,
  NDCall = as.integer(all_superpos %in% NDCall),
  Chl3q_t60 = as.integer(all_superpos %in% superpos_chl1),
  Ksg3q_t60 = as.integer(all_superpos %in% superpos_ksg1)
)
rownames(ndc_matrix) <- ndc_matrix$SuperPos
ndc_matrix$SuperPos <- NULL
# Set colors for the sets in the same order as in the upset() call
set_colors <- c("gray", "#7F0440","#108080")
# Draw the UpSet plot
upset(ndc_matrix,
      sets = c("NDCall", "Chl3q_t60", "Ksg3q_t60"),
      order.by = "freq",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "SuperPos per Group",
      sets.bar.color = set_colors)
