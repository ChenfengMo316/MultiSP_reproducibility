library(sceasy)
library(Seurat)
library(dplyr)
library(SeuratDisk)
library(anndata)

setwd("./data/human lymph node/")

adata <- read_h5ad("human lymph node adata_RNA.h5ad")
rna <- CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)
print(rna)


adata_adt <- read_h5ad("human lymph node adata_Protein.h5ad")
adt <- CreateSeuratObject(counts = t(adata_adt$X), meta.data = adata_adt$obs)
print(adt)


#RNA
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^mt-")
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))

#ADT
adt <- NormalizeData(adt, normalization.method = 'CLR', margin = 2)
adt <- FindVariableFeatures(adt, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(adt)
adt <- ScaleData(adt, features = all.genes)
adt <- RunPCA(adt,features = VariableFeatures(object =adt),reduction.name = 'apca' )

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

#####################
#Combined
combined <- rna
combined@reductions$apca=adt@reductions$apca

combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "apca"), dims.list = list(1:10, 2:10))
combined <-find_optimal_clusters(combined,target_clusters=11)
combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
DimPlot(combined , reduction = "wnn.umap", label = T) 


