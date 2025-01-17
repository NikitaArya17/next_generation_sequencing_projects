#install.packages("pheatmap")
#BiocManager::install("DESeq2")

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Reading in the raw count data
smoc2_rawcounts <- read.csv("fibrosis_smoc2_rawcounts_unordered.csv")
head(smoc2_rawcounts)

# Creating the metadata dataframe
genotype <- c("smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe")
condition <- c("fibrosis", "fibrosis", "fibrosis", "fibrosis", "normal", "normal", "normal")

smoc2_metadata <- data.frame(genotype, condition)
rownames(smoc2_metadata) <- c("smoc2_fibrosis1", "smoc2_fibrosis2", "smoc2_fibrosis3", "smoc2_fibrosis4", "smoc2_normal1", "smoc2_normal3", "smoc2_normal4")
smoc2_metadata

# Organising the data for DESeq2
# Use the match() function to reorder the columns of the raw counts 
#so it is in the same order as the rows of the metadata
reorder_idx <- match(rownames(smoc2_metadata), colnames(smoc2_rawcounts))

# Reorder the columns of the count data
reordered_smoc2_rawcounts <- smoc2_rawcounts[, reorder_idx]

# Create a DESeq2 object in a format that can be used for analysis
dds_smoc2 <- DESeqDataSetFromMatrix(countData =  reordered_smoc2_rawcounts,
                              colData = smoc2_metadata,
                              design = ~ condition)

#Quality control on samples
#To generate normalised counts 
#to accurately compare gene expression between samples

# Determine the size factors to use for normalization
dds_smoc2 <- estimateSizeFactors(dds_smoc2)

# Extract the normalized counts
smoc2_normalized_counts <- counts(dds_smoc2, normalized = TRUE)

#Plotting a hierarchical heatmap allows us to
#visualise the variance between the samples

# Log-transform the normalized counts
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Extract the matrix of transformed counts
vsd_mat_smoc2 <- assay(vsd_smoc2)

# Compute the correlation values between samples
vsd_cor_smoc2 <- cor(vsd_mat_smoc2)

# Plot the heatmap
pheatmap(vsd_cor_smoc2, annotation = select(smoc2_metadata, condition))
     
#PCA - Principal Component Analysis
#Plots your samples and plots lines between them known as Principal Components.
#These basically tell you how much of the variation in your data can be explained 
#by the fact that there are different samples
#PCA is a method of dimension reduction.
#It reduces a multi-dimensional dataset to two dimensions.

# Transform the normalized counts
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Plot the PCA of PC1 and PC2
plotPCA(vsd_smoc2, intgroup="condition")

#We now convert the samples matrix into a DESeq2 object as before, so it can be analysed
#We then finally run the analysis on it

# Create DESeq2 object
dds_smoc2 <- DESeqDataSetFromMatrix(countData = reordered_smoc2_rawcounts,
                 colData = smoc2_metadata,
                 design = ~ condition) #The design argument specifies the variables you want to control for

# Run the DESeq2 analysis
dds_smoc2_ <- DESeq(dds_smoc2)

 
#In addition to running an analysis as above, 
#we can plot the data using a negative binomial model
#And we can assess the fit of the data to the model

# Plot dispersions
plotDispEsts(dds_smoc2_)


#We have explored dispersions and fit our data to the DESeq2 model
#Now we need to extract the results

# Extract the results of the differential expression analysis
smoc2_res <- results(dds_smoc2_,
                contrast = c("condition", "fibrosis", "normal"),  
                alpha = 0.05)#To extract the results for fibrosis relative to normal


#We can improve the fold change estimates for the data by shrinking the log2 fold changes
#Log2 fold change = log2(fibrosis samples/normal samples)

# Shrink the log2 fold change estimates to be more accurate
smoc2_res <- lfcShrink(dds_smoc2_ ,
                    type = "normal",
                    contrast =  c("condition", "fibrosis", "normal"),
                    res = smoc2_res)

#We filter out genes not likely to be biologically relevant due to their being outliers
#We also want to reduce the chances of throwing out genes that are biologically relevant in the process
#We due this by changing the log2 fold change threshold

# Explore the results() function
#?results

# Extract results
smoc2_res <- results(dds_smoc2_,
                contrast = c("condition", "fibrosis", "normal"),
                alpha = 0.05,
                lfcThreshold = 0.32) #log changes threhold of 1.25 fold = log2 0.32

# Shrink the log2 fold changes
smoc2_res <- lfcShrink(dds_smoc2_,
                    type = "normal",
                    contrast = c("condition", "fibrosis", "normal"),
                    res = smoc2_res)

#We can now explore our extracted results
#and locate the differentially expressed genes according to oue alpha

# Get an overview of the results
summary(smoc2_res)
     

#We now extract significant results

# Save results as a data frame
smoc2_res_all <- data.frame(smoc2_res)

# Subset the results to only return the significant genes
#with p-adjusted values less than 0.05
smoc2_res_sig <- subset(smoc2_res_all, padj < 0.05)

#The results can be adequately explored through visualisations

# Create MA plot
plotMA(smoc2_res)

# Generate logical column
smoc2_res_all <- data.frame(smoc2_res) %>% mutate(threshold = padj < 0.05)

# Create the volcano plot
ggplot(smoc2_res_all) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
        xlab("log2 fold change") +
        ylab("-log10 adjusted p-value") +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))

#Continued
#Convert normalized counts to a dataframe
smoc2_normalized_counts <- data.frame(smoc2_normalized_counts)

# Subset normalized counts to significant genes
sig_norm_counts_smoc2 <- smoc2_normalized_counts[rownames(smoc2_res_sig), ]

# Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

# Plot heatmap
pheatmap(sig_norm_counts_smoc2,
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = select(smoc2_metadata, condition),
         scale = "row")

smoc2_res_sig <- smoc2_res_sig %>% 
                  data.frame() %>% 
                  rownames_to_column(var = "geneID") %>% 
                  arrange(padj) %>% 
                  select(geneID, padj) %>% 
                  head()

print(smoc2_res_sig)