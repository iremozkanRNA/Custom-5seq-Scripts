rm(list=ls())
library(dplyr)
library(purrr)

# Read and filter the data
fh1 <- read.csv("~/Desktop/Stuff/Desktop/antisenseKsg.csv") %>%
  filter(has_antisense == 1) %>%
  rename(SuperStrand = strand)

fh2 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/KsgData_fRanked.csv")

genes <- read.delim("Stuff/AMeyer_Lab/Experiment_results/3seq/RNAFOLD/RNAfold_DRIVE/Spn/Chl/Nofilter/NC_003028.v3.17.ncrna.genes", header = TRUE, sep = "")
colnames(genes) <- c("genome", "from", "to", "strand", "name", "old.name", "length", "protein.Name", "symbol")

# Split genes by strand
genes.top <- genes %>% filter(strand == "+")
genes.comp <- genes %>% filter(strand == "-")

cols_to_numeric <- c("from", "to")
genes.comp[cols_to_numeric] <- lapply(genes.comp[cols_to_numeric], as.numeric)
genes.top[cols_to_numeric] <- lapply(genes.top[cols_to_numeric], as.numeric)

strand_levels <- c("+", "-")
genome_levels <- unique(fh2$Genome)

# Function to check overlap and classify asRNA
classify_asRNA <- function(row, genes_opposite) {
  if (nrow(genes_opposite) == 0) return(FALSE)
  start <- row$from
  end <- row$to
  length <- abs(end - start) + 1
  overlaps <- genes_opposite %>%
    filter(!(to < start | from > end)) %>%
    mutate(
      overlap_start = pmax(from, start),
      overlap_end = pmin(to, end),
      overlap_len = pmax(0, overlap_end - overlap_start + 1)
    ) %>%
    filter(overlap_len > 0)
  if (nrow(overlaps) == 0) return(FALSE)
  max_overlap <- max(overlaps$overlap_len)
  if (max_overlap >= 100) return(TRUE)
  if (length < 200 && max_overlap >= 0.5 * length) return(TRUE)
  return(FALSE)
}

results_list <- map(genome_levels, function(genome) {
  fh2_genome <- fh2 %>% filter(Genome == genome)
  fh2_split <- lapply(strand_levels, function(s) fh2_genome %>% filter(SuperStrand == s))
  names(fh2_split) <- strand_levels
  
  fh1_split <- lapply(strand_levels, function(s) fh1 %>% filter(SuperStrand == s))
  names(fh1_split) <- strand_levels
  
  genes_split <- list(
    "+" = genes.comp, # opposite of +
    "-" = genes.top   # opposite of -
  )
  
  map2_df(
    fh1_split,
    fh2_split,
    ~ {
      if(nrow(.x) == 0) return(NULL)
      strand <- unique(.x$SuperStrand)[1]
      genes_opposite <- genes_split[[strand]]
      
      .x <- .x %>%
        rowwise() %>%
        mutate(
          from = PeakCoord,  # If you have real start/end, use them here
          to = PeakCoord,
          asRNA = classify_asRNA(cur_data(), genes_opposite)
        ) %>%
        ungroup()
      
      if(nrow(.y) == 0) {
        .x %>%
          mutate(closest_TSS = NA_real_, TSS_distance = NA_real_, Genome_in_fh2 = genome)
      } else {
        if(strand == "+") {
          .x %>%
            mutate(
              closest_TSS = map_dbl(
                PeakCoord,
                function(pos) {
                  candidates <- .y$SuperPos[.y$SuperPos < pos]
                  if(length(candidates) == 0) return(NA_real_)
                  candidates[which.max(candidates)]
                }
              ),
              TSS_distance = PeakCoord - closest_TSS,
              Genome_in_fh2 = genome
            )
        } else {
          .x %>%
            mutate(
              closest_TSS = map_dbl(
                PeakCoord,
                function(pos) {
                  candidates <- .y$SuperPos[.y$SuperPos > pos]
                  if(length(candidates) == 0) return(NA_real_)
                  candidates[which.min(candidates)]
                }
              ),
              TSS_distance = closest_TSS - PeakCoord,
              Genome_in_fh2 = genome
            )
        }
      }
    },
    .id = "strand"
  )
})

fh1_processed_all_genomes <- bind_rows(results_list)
write.csv(fh1_processed_all_genomes, file = "~/Desktop/Ksg_antisense.csv", row.names = F, col.names = T)
