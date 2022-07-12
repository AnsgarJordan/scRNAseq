
library(dplyr)
library(Seurat)
library(patchwork)

cnc.data <- Read10X(data.dir = "filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
cnc <- CreateSeuratObject(counts = cnc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
cnc
cnc[["percent.mt"]] <- PercentageFeatureSet(cnc, pattern = "^MT-")

# prune out the unnecessary cells (Quality Control (QC))
# cnc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cnc <- scaleData(cnc)

cnc <- JackStraw(cnc, num.replicate = 100)
cnc <- ScoreJackStraw(cnc, dims = 1:15)

# JackStrawPlot(cnc, dims = 1:15)
ElbowPlot(cnc)

## use the same number of PC's as prescribed by the elbowplot

cnc <- FindNeighbors(cnc, dims = 1:12)
cnc <- FindClusters(cnc, resolution = 0.5)

cnc <- RunUMAP(cnc, dims = 1:12)
DimPlot(cnc, reduction = "umap")