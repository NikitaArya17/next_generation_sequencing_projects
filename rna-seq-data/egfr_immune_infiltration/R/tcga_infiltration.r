#Load the required packages
library(readr)
library(tidyverse)
library(stringr)

#Read in the data
leuk_est <- read_tsv("TCGA_all_leuk_estimate.masked.20170107.tsv")
full_ids <- read_tsv("TCGA.Kallisto.fullIDs.cibersort.relative.tsv")

#Explore the data
head(full_ids)
glimpse(full_ids)

#Demonstrate that the second file contains relative immune cell abundance
df_cell_types <- full_ids 
                    %>% select(contains(c("cells", "cytes", "phages", "phils")))
glimpse(df_cell_types)
