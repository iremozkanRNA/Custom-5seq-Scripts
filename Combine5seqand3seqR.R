rm(list = ls())
library(dplyr)

fh1 <- read.csv("Stuff/AMeyer_Lab/Experiment_results/3seq/RNAFOLD/QuinnLast/Spn/Chl/Nofilter/outputresults.csv") 
fh2 <- read.csv("5seq/Unprocessed/Chl/316-5NTapuChlNDCt60_BothStrands.csv")
fh3 <- read.csv("5seq/Unprocessed/Chl/316-5NTapuChl1quarterMICt60_BothStrands.csv")
fh4 <- read.csv("5seq/Unprocessed/Chl/316-5NTapuChl3quarterMICt60_BothStrands.csv")
fh1$unprocessed5seq <- ""

##### NDC
for(i in 1:nrow(fh1)) {
  # Check your locus and Sig_Control conditions once per row i
  if(fh1$Locus[i] != "" & fh1$Sig_Control[i] == 1) {
    for(j in 1:nrow(fh2)) {
      # Check if PeakCoord[i] is within 5 nt of HighestPeak[j]
      if(abs(fh1$PeakCoord[i] - fh2$HighestPeak[j]) <= 55) {
        fh1$unprocessed5seq[i] <- "unprocessed 5'seq in NDC"
        break  # no need to continue looping j for this i once condition is met
      }
    }
  }
}

##### 1Q
for(i in 1:nrow(fh1)) {
  # Check your locus and Sig_Control conditions once per row i
  if(fh1$Locus[i] != "" & fh1$Sig_Cond1[i] == 1) {
    for(j in 1:nrow(fh3)) {
      # Check if PeakCoord[i] is within 5 nt of HighestPeak[j]
      if(abs(fh1$PeakCoord[i] - fh3$HighestPeak[j]) <= 55) {
        fh1$unprocessed5seq[i] <- "unprocessed 5'seq in 1Q"
        break  # no need to continue looping j for this i once condition is met
      }
    }
  }
}
##### 3Q
for(i in 1:nrow(fh1)) {
  # Check your locus and Sig_Control conditions once per row i
  if(fh1$Locus[i] != "" & fh1$Sig_Cond2[i] == 1) {
    for(j in 1:nrow(fh4)) {
      # Check if PeakCoord[i] is within 5 nt of HighestPeak[j]
      if(abs(fh1$PeakCoord[i] - fh4$HighestPeak[j]) <= 55) {
        fh1$unprocessed5seq[i] <- "unprocessed 5'seq in 3Q"
        break  # no need to continue looping j for this i once condition is met
      }
    }
  }
}

sum(fh1$unprocessed5seq != "")

write.csv(fh1, file = "Stuff/Desktop/Chl_RNAfold_Unprocessed5seqInfoAdded.csv", row.names = F, col.names = T)


