# Install and load the necessary packages
#BiocManager::install("DESeq2")
#install.packages("pheatmap")

library(DESeq2)
library(tidyverse)
library(pheatmap)

# Read in the metadata file and the gene expression matrix
meta_filepath <- "/Users/nikit/sequencing data analysis projects in R/metadata.csv"

meta_data <- read.csv(meta_filepath)
gene_matrix <- read.csv("https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv")

# Exploring the metadata file
head(meta_data)
str(meta_data)
summary(meta_data)
colnames(meta_data)

meta_data <- meta_data %>% select(!ends_with("_color"))
