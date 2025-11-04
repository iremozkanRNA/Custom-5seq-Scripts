# Clear environment
rm(list = ls())

# Load dplyr
library(dplyr)

# Read BED files without headers
fh1 <- read.delim("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/PIPETS/5seq/Processed/316-5PTappKsg3quarterMICt60.bed", header = FALSE)
fh2 <- read.delim("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/PIPETS/5seq/Unprocessed/316-5NTapuKsg3quarterMICt60.bed", header = FALSE)

# Read CSV with HighestPeak filter
fh3 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/PIPETS/5seq/Processed/Ksg/316-5PTappKsg3quarterMICt60.bed_TopStrandResults.csv", header = TRUE) %>%
  filter(HighestPeak > 776000 & HighestPeak < 785000)

# Remove unwanted columns (4th and 5th)
fh1 <- fh1[, -c(4, 5)]
fh2 <- fh2[, -c(4, 5)]

# Separate based on strand (+ or -)
fh1Minus <- fh1 %>% filter(V6 == "-")
fh1Plus  <- fh1 %>% filter(V6 == "+")
fh2Minus <- fh2 %>% filter(V6 == "-")
fh2Plus  <- fh2 %>% filter(V6 == "+")

# Create count tables: frequency of coordinates
fh1MinusCounts <- as.data.frame(table(fh1Minus$V3))
fh1PlusCounts  <- as.data.frame(table(fh1Plus$V2))
fh2MinusCounts <- as.data.frame(table(fh2Minus$V3))
fh2PlusCounts  <- as.data.frame(table(fh2Plus$V2))

# Convert coordinate factors to numeric
fh1MinusCounts$Var1 <- as.numeric(as.character(fh1MinusCounts$Var1))
fh2MinusCounts$Var1 <- as.numeric(as.character(fh2MinusCounts$Var1))

# Define coordinate range from fh3 HighestPeak
min_peak <- min(fh3$HighestPeak)
max_peak <- max(fh3$HighestPeak)

# Filter counts within region and counts greater than 10
fh1_filtered <- fh1MinusCounts %>%
  filter(Var1 >= min_peak & Var1 <= max_peak & Freq > 10)

fh2_filtered <- fh2MinusCounts %>%
  filter(Var1 >= min_peak & Var1 <= max_peak & Freq > 10)

# Merge data frames by coordinate (Var1)
merged_counts <- merge(fh1_filtered, fh2_filtered, by = "Var1", suffixes = c("_fh1", "_fh2"))

# Calculate ratio per coordinate
merged_counts <- merged_counts %>%
  mutate(ratio = Freq_fh1 / Freq_fh2)

# View results
print(merged_counts)

write.csv(merged_counts,file = "~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/PIPETS/5seq/5seq_processedoverunprocessedratio_for_purineregion_Ksg.csv", row.names = F, col.names = T)
