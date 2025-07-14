# Install and load the necessary packages
#BiocManager::install(c("DESeq2", "biomaRt"))
#install.packages(c("pheatmap", "data.table"))

library(DESeq2)
library(biomaRt)
library(tidyverse)
library(pheatmap)

# Read in the metadata file and the gene expression matrix
meta_filepath <- "/Users/nikit/sequencing data analysis projects in R/metadata.csv"
matrix_link <- "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv"

meta_data <- read.csv(meta_filepath)
gene_matrix <- data.table::fread(matrix_link)

# Exploring the metadata file
head(meta_data)
str(meta_data)
summary(meta_data)
colnames(meta_data)

meta_data <- meta_data %>% select(!ends_with("_color") & !ends_with("_order"))

# Exploring the gene matrix file
head(gene_matrix)
str(gene_matrix)
