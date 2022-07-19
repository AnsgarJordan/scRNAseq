
library(dplyr)
library(Seurat)
library(patchwork)

cnc.data <- Read10X(data.dir = "filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
cnc <- CreateSeuratObject(counts = cnc.data, project = "cnc3k", min.cells = 3, min.features = 200)
# cnc
cnc[["percent.mt"]] <- PercentageFeatureSet(cnc, pattern = "^MT-")

# prune out the unnecessary cells (Quality Control (QC))
cnc <- subset(cnc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Feature Selection
cnc <- FindVariableFeatures(cnc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cnc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cnc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scaling - always do this before doing PCA

cnc <- ScaleData(cnc)

# PCA
cnc <- RunPCA(cnc, features = VariableFeatures(object = cnc))

# visualize PCA
DimPlot(cnc, reduction = "pca")


cnc <- FindNeighbors(cnc, dims = 1:10)
cnc <- FindClusters(cnc, resolution = 0.5)
# using UMAP
cnc <- RunUMAP(cnc, dims = 1:10)
DimPlot(cnc, reduction = "umap")




## 
umap = cbind("Barcode" = rownames(Embeddings(object = cnc, reduction = "umap")), Embeddings(object = cnc, reduction = "umap"))
write.table(umap, file="umap_rna.csv", sep = ",", quote = F, row.names = F, col.names = T)



