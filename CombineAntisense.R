rm(list=ls())
lapply(c("dplyr", "purrr"), library, character.only = TRUE)

# Read and filter the data
fh1 <- read.csv("~/Desktop/Stuff/Desktop/antisenseChl.csv") %>%
  filter(Sig_Control == 1, has_antisense == 1)

fh2 <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/ChlData_fRanked.csv") %>%
  filter(Genome == "NDCt60") %>%
  distinct(SuperPos, .keep_all = TRUE)

# Ensure consistent strand column naming
fh1 <- fh1 %>% rename(SuperStrand = strand)

# Define strand levels
strand_levels <- c("+", "-")

# Split both datasets by strand, ensuring all levels are present
fh1_split <- lapply(strand_levels, function(s) fh1 %>% filter(SuperStrand == s))
names(fh1_split) <- strand_levels

fh2_split <- lapply(strand_levels, function(s) fh2 %>% filter(SuperStrand == s))
names(fh2_split) <- strand_levels

# For each strand, match according to biological convention
fh1_processed <- map2_df(
  fh1_split,
  fh2_split,
  ~ {
    if(nrow(.x) == 0) return(NULL)
    if(nrow(.y) == 0) {
      .x %>%
        mutate(closest_antisense = NA_real_, length = NA_real_)
    } else {
      if(unique(.x$SuperStrand)[1] == "+") {
        # For + strand: closest SuperPos < PeakCoord (upstream)
        .x %>%
          mutate(
            closest_antisense = map_dbl(
              PeakCoord,
              function(pos) {
                candidates <- .y$SuperPos[.y$SuperPos < pos]
                if(length(candidates) == 0) return(NA_real_)
                candidates[which.max(candidates)]
              }
            ),
            length = PeakCoord - closest_antisense
          )
      } else {
        # For - strand: closest SuperPos > PeakCoord (downstream)
        .x %>%
          mutate(
            closest_antisense = map_dbl(
              PeakCoord,
              function(pos) {
                candidates <- .y$SuperPos[.y$SuperPos > pos]
                if(length(candidates) == 0) return(NA_real_)
                candidates[which.min(candidates)]
              }
            ),
            length = closest_antisense - PeakCoord
          )
      }
    }
  },
  .id = "strand"
)
