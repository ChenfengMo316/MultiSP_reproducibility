library(reticulate)
library(dplyr)
library(Seurat)
library(patchwork)
library(Signac)
library(anndata)

setwd("./data/P22 mouse brain/")

#RNA
adata <- read_h5ad("P22_mouse_brain_adata_RNA.h5ad")
rna <- CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)
print(rna)

rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^mt-")

rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))

#ATAC
adata_atac <- read_h5ad("P22_mouse_brain_adata_ATAC.h5ad")
atac <- CreateSeuratObject(counts = t(adata_atac$X), meta.data = adata_atac$obs,assay='peaks')
print(atac)

atac <- RunTFIDF(atac, assay = "peaks")  #normalization
atac <- FindTopFeatures(atac, min.cutoff = 'q0', assay = "peaks")
atac <- RunSVD(atac, assay = "peaks")

find_optimal_clusters <- function(combined, graph.name = "wsnn", target_clusters, key_added = "optimal_clusters", 
                                  range_min = 0, range_max = 3, max_steps = 30, tolerance = 0) {
  this_step <- 0
  this_min <- range_min
  this_max <- range_max
  
  while (this_step < max_steps) {
    this_resolution <- this_min + ((this_max - this_min) / 2)
    
    # Run FindClusters with the current resolution
    combined <- FindClusters(combined, graph.name = graph.name, verbose = FALSE, resolution = this_resolution)
    this_clusters <- length(unique(Idents(combined)))
    
    if (this_clusters > target_clusters + tolerance) {
      this_max <- this_resolution
    } else if (this_clusters < target_clusters - tolerance) {
      this_min <- this_resolution
    } else {
      print(sprintf("Succeed to find %d clusters at resolution %.3f", target_clusters, this_resolution))
      combined[[key_added]] <- Idents(combined)
      return(combined)
    }
    
    this_step <- this_step + 1
  }
  
  print('Cannot find the number of clusters')
  combined[[key_added]] <- Idents(combined)
  return(combined)
}

#Combined
combined <- rna
combined@reductions$lsi=atac@reductions$lsi

combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "lsi"), dims.list = list(1:10, 2:10))
combined <-find_optimal_clusters(combined,target_clusters=11)
combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")








