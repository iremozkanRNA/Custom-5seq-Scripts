# Load libraries
library(ggplot2)
library(plotly)

# Read data
fh1 <- read.csv("~/Desktop/AllRNAseq/t60/316-Chl1quarterMICt60-316-Chlndct60-DGE-results.csv")
fh2 <- read.csv("~/Desktop/AllRNAseq/t60/316-Kas1quarterMICt60-316-Kasndct60-DGE-results.csv")
fh3 <- read.csv("~/Desktop/AllRNAseq/t60/316-Tet1quarterMICt60-316-Tetndct60-DGE-results.csv")

# Merge all three by 'Gene'
merged12 <- merge(fh1, fh2, by = "Gene", suffixes = c(".fh1", ".fh2"))
merged <- merge(merged12, fh3, by = "Gene")

# Ensure proper column names for log2FC and pvalue from fh3
# If needed, adjust these to match your fh3 column names
colnames(merged)[which(colnames(merged) == "Value")] <- "Value"
colnames(merged)[which(colnames(merged) == "pvalue")] <- "pvalue"

# Remove non-finite values
merged <- merged[
  is.finite(merged$Value.fh1) &
    is.finite(merged$Value.fh2) &
    is.finite(merged$Value), 
]

# Logical vectors for significance in each antibiotic
sig_chl <- abs(merged$Value.fh1) > 1 & merged$pvalue.fh1 < 0.05
sig_ksg <- abs(merged$Value.fh2) > 1 & merged$pvalue.fh2 < 0.05
sig_tet <- abs(merged$Value) > 1 & merged$pvalue < 0.05

# Count number of antibiotics each gene is significant in
sig_count <- sig_chl + sig_ksg + sig_tet

# Assign new significance categories
merged$SigCategory <- NA
merged$SigCategory[sig_count == 1 & sig_chl] <- "Chl only"
merged$SigCategory[sig_count == 1 & sig_ksg] <- "Ksg only"
merged$SigCategory[sig_count == 1 & sig_tet] <- "Tet only"
merged$SigCategory[sig_count == 2] <- "Significant in 2 Ab"
merged$SigCategory[sig_count == 3] <- "Significant in all Ab"

# Filter to plot only significant genes
merged_sig <- subset(merged, !is.na(SigCategory))

# Define custom colors
custom_colors <- c(
  "Chl only" = "lightgreen",
  "Ksg only" = "pink",
  "Tet only" = "royalblue",
  "Significant in 2 Ab" = "goldenrod",
  "Significant in all Ab" = "cyan"
)

# 3D Plotly scatter plot
plot_ly(
  merged_sig,
  x = ~Value.fh1,
  y = ~Value.fh2,
  z = ~Value,
  color = ~SigCategory,
  colors = custom_colors,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5, opacity = 0.75)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "log2FC (1/4 xMIC Chl)"),
      yaxis = list(title = "log2FC (1/4 xMIC Ksg)"),
      zaxis = list(title = "log2FC (1/4 xMIC Tet)")
    ),
    title = "Gene-wise 3D Value Comparison at T60"
  )
