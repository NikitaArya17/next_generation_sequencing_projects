# Single Cell ATAC Sequencing

This is a sequencing method that allows us to find regions of chromatin that are not condensed or wound around histone proteins. In other words, these regions are available to bind with DNA binding proteins.

The data obtained from a single cell ATAC seq experiment can be analysed with the `Signac` and `Seurat` packages in R. A demonstration has been provided in this notebook, using a dataset of human peripheral blood mono-nuclear cells provided by 10X Genomics, available from the [Stuart Lab](https://stuartlab.org/signac/articles/pbmc_vignette#non-linear-dimension-reduction-and-clustering).

## Preprocessing

The first step in the workflow is to read in and preprocess the data. This is done by creating a `Seurat` object and adding gene annotations to it.


## Quality Control

There are five metrics that are considered for quality control of scATAC-seq data:

1. Nucleosome signal and nucleosome banding pattern
2. Transcription Start Site (TSS) Enrichment
3. Total number of fragments in peaks
4. Fraction of fragments in peaks
5. Ratio reads in genomic blacklist regions

These metrics will be examined in detail along with the code used to perform quality control.


## Normalisation and Linear Dimension Reduction

Normalisation is essential to correct for differences in cellular sequencing depth across cells and peaks. In other words, it corrects for technical differences between cells.

Dimension reduction is the process of selecting and retaining only relevant features of the dataset, those that explain the highest amount of variability in the dataset.


## Clustering

The cells are now embedded in a low-dimensional space. Methods for graph-based clustering and non-linear dimension reduction for visualisation can be applied.


## Downstream Analysis

The data is now preprocessed and ready for downstream analysis. These include procedures such as cell identity annotation.
