rm(list=ls())
lapply(c("dplyr", "UpSetR"), library, character.only = T)

fh1 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/ChlData_fRanked.csv") %>% filter(rank == "primary")

# Get unique Genome types
genomes <- unique(fh1$Genome)

# Create binary matrix
all_superpos <- unique(fh1$SuperPos)
binary_matrix <- sapply(genomes, function(g) all_superpos %in% fh1$SuperPos[fh1$Genome == g])
binary_matrix <- as.data.frame(binary_matrix)
colnames(binary_matrix) <- genomes
rownames(binary_matrix) <- all_superpos

# Assign colors (customize as needed)
set_colors <- c("gray", "lightgreen", "darkgreen")[seq_along(genomes)]
binary_matrix[] <- lapply(binary_matrix, as.integer)

# Plot
upset(
  binary_matrix,
  nsets = length(genomes),
  sets = genomes,
  sets.bar.color = set_colors,
  keep.order = TRUE
)
print(genomes)
length(genomes)
str(binary_matrix)
colnames(binary_matrix)
print(set_colors)
length(set_colors)
