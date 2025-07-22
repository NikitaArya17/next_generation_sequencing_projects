#Load the required packages
library(readr)
library(tidyverse)
library(ggplot2)

#Read in the data
leuk_est <- read_tsv("TCGA_all_leuk_estimate.masked.20170107.tsv", col_names = FALSE)
immune_est <- read_tsv("TCGA.Kallisto.fullIDs.cibersort.relative.tsv")

#Explore the data
head(immune_est)
glimpse(immune_est)

#Identify and remove duplicated rows
sum(duplicated(immune_est))
immune_est[duplicated(immune_est), ]
immune_est <- distinct(immune_est)
sum(duplicated(immune_est))

#There are multiple aliquots from the same sample.
sum(duplicated(immune_est$SampleID))

#Demonstrate that the second file contains relative immune cell abundance
df_cell_types <- select(immune_est, contains(c("cells", "cytes", "phages", "phils")))
glimpse(df_cell_types)

df_cell_types$cell_type_fractions <- rowSums(df_cell_types)
colnames(df_cell_types)

#As all immune cell type fractions sum to 1,
#we have proved that this file contains relative immune cell abundance.
unique(df_cell_types$cell_type_fractions)

#Cleaning the leukaemia fraction dataset
#to create a common column on which to join the two files.
glimpse(leuk_est)
colnames(leuk_est) <- c("CancerType", "SampleID", "Leukocytes")

colnames(leuk_est) %in% colnames(immune_est)

head(immune_est$SampleID)

colnames(leuk_est)
head(leuk_est)

