args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

library(Seurat)
library(tidyverse)
library(qs)
library(corral)
library(glmGamPoi)

# Define corral DR function
RunCorral <- function(obj, counts = T, alt = NULL, n = 30, residual = 'standardized', reduction.name = 'corral'){
  if (counts) {
    dimred <- corral(seur_obj@assays$RNA@counts, ncomp = n, rtype = residual)
  }  else {
    dimred <- corral(alt, ncomp = n, rtype = residual)
  }
  seur_obj[[reduction.name]] <- CreateDimReducObject(embeddings = dimred$v, loadings = dimred$u, key = reduction.name)
  return(seur_obj)
}

# Load data and create seurat object
data <- read_csv(input_file) %>% 
  replace(is.na(.), 0) %>% 
  rowid_to_column('index') %>%
  column_to_rownames('index')

metadata <- as.data.frame((data[1:9]))
genexp <- as.data.frame(data[10:ncol(data)])

seur_obj <- CreateSeuratObject(counts = t(genexp),
    project = 'spatialtx', meta.data = metadata)

cell_coords <- seur_obj@meta.data %>% mutate(coords_1 = x, coords_2 = y) %>% select(c(coords_1, coords_2))

seur_obj[['spatial']] <- CreateDimReducObject(as.matrix(cell_coords), key = 'spatial_')


seur_obj <- RunCorral(obj = seur_obj)

seur_obj <- NormalizeData(seur_obj)
seur_obj <- ScaleData(seur_obj)
seur_obj <- RunPCA(seur_obj, features = rownames(seur_obj), reduction.name = 'lognorm_pca')


seur_obj <- SCTransform(seur_obj)
seur_obj <- RunPCA(seur_obj, features = rownames(seur_obj), reduction.name = 'sct_poisson_pca', assay = 'SCT')

seur_obj <- SCTransform(seur_obj, method = 'glmGamPoi')
seur_obj <- RunPCA(seur_obj, features = rownames(seur_obj), reduction.name = 'sct_glmGamPoi_pca', assay = 'SCT')

qsave(seur_obj, output_file)
