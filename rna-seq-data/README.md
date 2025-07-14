The R script contains an analysis conducted as part of the 'RNA-Seq with Bioconductor in R' course offered by Mary Piper on [DataCamp](https://www.datacamp.com/courses/rna-seq-with-bioconductor-in-r). The R script contains the basic workflow for the analysis of RNA-sequencing data, using data from a study of fibrosis in mouse kidney cells as an example. The code for this basic workflow was built off of the code that was taught in the course.

The script also contains code to convert Entrez Gene IDs to their corresponding gene symbols, at the end of the analysis. This code was built off of the documentation for the `biomaRt` library, available on [the Bioconductor webpage](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html).


# RNA-Seq Differential Expression Analysis Workflow Summary

Check that all of the samples are in the same order in the metadata and count data

`all(rownames(all_metadata) %in% colnames(all_rawcounts))`

## DESeq object to test for the effect of fibrosis regardless of genotype

`dds_all <- DESeqDataSetFromMatrix(countData = all_rawcounts, colData = all_metadata, design = ~ condition + genotype)`

## DESeq object to test for the effect of genotype on the effect of fibrosis

`dds_complex <- DESeqDataSetFromMatrix(countData = all_rawcounts, colData = all_metadata, design = ~ genotype + condition + genotype:condition)`

## Log transform counts for QC (Quality Control)

`vsd_all <- vst(dds_all, blind = TRUE)`

## Create heatmap of sample correlation values

`vsd_all %>% assay() %>% cor() %>% pheatmap(annotation = select(all_metadata, c("genotype", "condition")))`

## Create the PCA plot for PC1 and PC2 and color by condition

`plotPCA(vsd_all, intgroup = "condition")`

## Create the PCA plot for PC1 and PC2 and color by genotype

`plotPCA(vsd_all, intgroup = "genotype")`

## Select significant genes with padj < 0.05

`smoc2_sig <- subset(res_all, padj < 0.05) %>% data.frame() %>% rownames_to_column(var = "geneID")`

## Extract the top 6 genes with padj values

`smoc2_sig %>% arrange(padj) %>% select(geneID, padj) %>% head()`
