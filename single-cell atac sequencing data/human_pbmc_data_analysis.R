if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install(c("EnsDb.Hsapiens.v75", "biovizBase"))
install.packages(c("Seurat", "Matrix", "irlba", "hdf5r", "remotes", "tidyverse"), type = "source")
remotes::install_github("stuart-lab/signac", ref = "develop")

#To update outdated packages
BiocManager::install(c(
  "bit", "cpp11", "GenomeInfoDb", "Hmisc", "parallelly", "RSQLite", "Signac"
), update = TRUE, ask = FALSE, force = TRUE)

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(biovizBase)
library(GenomicRanges)
library(tidyverse)
library(hdf5r)
library(Matrix)
library(irlba)

# QUALITY CONTROL AND DATA PREPROCESSING

# The fragment file contains all unique fragments from each single cell.
#frag.file <- read.delim("10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz", header = F, nrows = 10)
#head(frag.file)

# The rows of the counts matrix contain the regions on the chromosome that are unwound.
#counts <- Read10X_h5("10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")
#counts[1:10, 1:10]

# Both of the files mentioned previously are used to create a chromatin assay.
# This allows the use of specialised functions to analyse single-cell genomic data.

chrom_assay <- CreateChromatinAssay(
  counts = Read10X_h5("10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5"),
  sep = c(":", "-"),
  fragments = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)


str(chrom_assay)

# The metadata contains essential information about the data
#metadata <- read.csv("10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv", header = T, row.names = 1)

View(metadata)

# The metadata and the chromatin assay are used to create the Seurat object.
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = read.csv("10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv", header = T, row.names = 1),
  assay = "ATAC"
)

str(pbmc)


#There are currently no added annotations, so this track will be empty.
pbmc$ATAC$annotation

# We can add gene annotations for the human genome.
# This will allow downstream functions to obtain the gene annotation information directly from the Seurat object.
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
annotations

seqlevels(annotations) <- paste("chr", seqlevels(annotations))

Annotation(pbmc) <- annotations

#The annotations from the Ensembl Database have been added.
pbmc$ATAC$annotation

# 1. Nucleosome Signal
# The primary goal of scATAC-seq is to detect regions of unwound DNA, where DNA binding proteins are free to bind.
# The signal helps us to detect regions of DNA that are condensed into nucleosomes and thus cannot bind to DNA binding proteins.
pbmc <- NucleosomeSignal(pbmc)


#2. TSS Enrichment
# DNA that is unwound is available to be transcribed.
# Therefore, we are interested in locating peaks found on the transcription start sites of the DNA.
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

#3. Blacklist Ratio
# This metric helps us to locate and exclude certain regions of DNA that pose problems for high-throughput sequencing analyses.
# These regions include unstructured or repetitive regions.
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

#4. Fraction of Reads in Peaks
# Peaks with low fractions of reads are considered of low quality and should be filtered out.
# Generally, peaks with < 15-20% reads are filtered out.
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100


View(pbmc@meta.data)

# Visualising the Quality Control Metrics helps us decide our cut-off thresholds.
colnames(pbmc@meta.data)
a1 <- DensityScatter(pbmc, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)

a2 <- DensityScatter(pbmc, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

a1 | a2

VlnPlot(object = pbmc,
        features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", "blacklist_ratio", "fraction_reads_in_peaks"),
        pt.size = 0.1,
        ncol = 6)

# Filtering out low quality cells according to our selected thresholds.
pbmc <- subset(x = pbmc,
               subset = nCount_ATAC > 3000 &
                 fraction_reads_in_peaks > 15 &
                 blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 3)


#Normalisation for differences in cellular sequencing depth across cells and peaks

pbmc <- RunTFIDF(pbmc) # Normalisation using Term Frequency - Inverse Document Frequency
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') # selecting top features
pbmc <- RunSVD(pbmc) # Dimensionality reduction by Singular Value Decomposition


DepthCor(pbmc) #correlates sequencing depth with reduced dimensions

# Clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30) #The first component is excluded as it captures the technical variation
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

#Visualise the data
DimPlot(object = pbmc, label = TRUE) + NoLegend()

# DOWNSTREAM ANALYSIS

# Create a gene activity matrix to integrate the scATAC data with scRNA-seq data.

gene.activities <- GeneActivity(pbmc)
gene.activities[1:10, 1:10]

