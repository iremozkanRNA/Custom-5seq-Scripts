rm(list=ls())

data <- read.csv("path/to/your/file/that/includes/UTRlenghts") %>% filter(Genome == "NDCt60" & TSSUTR != "")

total <- nrow(data)
count <- sum(!is.na(data$UTRlength != "") & data$UTRlength > -50 & data$UTRlength < 0)

frequency <- count/total * 100
frequency 
