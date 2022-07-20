
library(dplyr)
library(Seurat)
library(patchwork)

cnc.data <- Read10X(data.dir = "filtered_feature_bc_matrix")

cnc <- CreateSeuratObject(counts = cnc.data, project = "cnc3k", min.cells = 3, min.features = 200)


# identify the mitochondrial counts for each cell
cnc[["percent.mt"]] <- PercentageFeatureSet(cnc, pattern = "^MT-")

# choose a subset of the data set depending on parameters
cnc <- subset(cnc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# normalization
cnc <- NormalizeData(cnc)


# Feature Selection 
cnc <- FindVariableFeatures(cnc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cnc), 10)

# scale with all genes
all.genes <- rownames(cnc)
cnc <- ScaleData(cnc, features = all.genes)

# scale with just the 2000 most variable genes
cnc <- ScaleData(cnc)

# PCA
cnc <- RunPCA(cnc, features = VariableFeatures(object = cnc))

VizDimLoadings(cnc, dims = 1:2, reduction = "pca")

# visualize PCA
# DimPlot(cnc, reduction = "pca")

DimHeatmap(cnc, dims = 1:15, cells = 500, balanced = TRUE)


# cnc <- JackStraw(cnc, num.replicate = 100)
# cnc <- ScoreJackStraw(cnc, dims = 1:20)
# JackStrawPlot(cnc, dims = 1:15)

ElbowPlot(cnc)


cnc <- FindNeighbors(cnc, dims = 1:12)
cnc <- FindClusters(cnc, resolution = 0.5)
# using UMAP
cnc <- RunUMAP(cnc, dims = 1:12)
DimPlot(cnc, reduction = "umap")




## create CSV file to transfer from Loupe -> Seurat and back
umap = cbind("Barcode" = rownames(Embeddings(object = cnc, reduction = "umap")), Embeddings(object = cnc, reduction = "umap"))
write.table(umap, file="umap_rna.csv", sep = ",", quote = F, row.names = F, col.names = T)



