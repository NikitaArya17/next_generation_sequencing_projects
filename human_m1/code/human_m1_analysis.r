# Install and load the necessary packages
BiocManager::install("DESeq2")
install.packages("pheatmap")

library(DESeq2)
library(tidyverse)
library(pheatmap)

# Read in the metadata file and the gene expression matrix
meta_data <- read.csv("https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv")
exprs_mat <- read.csv("https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv")

# Explore the data files
head(meta_data)
head(exprs_mat)
