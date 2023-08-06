args <- commandArgs(TRUE)
st_input <- args[1]
output_file <- args[2]

library(Seurat)
library(tidyverse)
library(qs)
library(corral)
library(cluster)
library(Metrics)
library(glmGamPoi)
library(fields)

RunCorral <- function(seur_obj, counts = T, alt = NULL, n = 30, residual = 'standardized', reduction.name = 'corral'){
  if (counts) {
    dimred <- corral(seur_obj@assays$RNA@counts, ncomp = n, rtype = residual)
  }  else {
    dimred <- corral(alt, ncomp = n, rtype = residual)
  }
  seur_obj[[reduction.name]] <- CreateDimReducObject(embeddings = dimred$v, loadings = dimred$u, key = reduction.name)
  return(seur_obj)
}

DimReduc <- function(data_obj, dimreduc = 'lognorm_pca', knn = 20, ndims = 30){
  DefaultAssay(data_obj) <- 'RNA'
  data_obj <- data_obj[,unname(which(colSums(GetAssayData(data_obj))!=0))]
  if (dimreduc %in% c('raw_pca', 'lognorm_pca', 'sct_pca', 'qp_sct_pca', 'glm_sct_pca')){
    if (dimreduc == 'lognorm_pca'){
      data_obj <- NormalizeData(data_obj)
    }
    if (dimreduc != 'sct_pca'){
      data_obj <- FindVariableFeatures(data_obj, nfeatures = 100)
      data_obj <- ScaleData(data_obj)
    }
    if (dimreduc == 'sct_pca'){
      data_obj <- SCTransform(data_obj)
    }
    if (dimreduc == 'qp_sct_pca'){
      data_obj <- SCTransform(data_obj, method = 'qpoisson')
    }
    if (dimreduc == 'glm_sct_pca'){
      data_obj <- SCTransform(data_obj, method = 'glmGamPoi')
    }
    data_obj <- RunPCA(data_obj, reduction.name = 'reduction')
  }
  if (dimreduc == 'lognorm_corral'){
    data_obj <- NormalizeData(data_obj)
    data_obj <- RunCorral(data_obj, counts = F, alt = data_obj@assays$RNA@data, reduction.name = 'reduction')
  }
  if (dimreduc == 'raw_corral'){
    data_obj <- RunCorral(data_obj, reduction.name = 'reduction')
  }
  if (dimreduc == 'sct_corral'){
    data_obj <- SCTransform(data_obj)
    data_obj <- RunCorral(data_obj, counts = F, alt = data_obj@assays$SCT@data, reduction.name = 'reduction')
  }
  data_obj <- FindNeighbors(data_obj, reduction = 'reduction', k.param = knn+1,
                                  return.neighbor = T, graph.name = 'nn_obj', dims = 1:ndims)
  return(data_obj)
}

FindNNIndices <- function(data_obj){
  nn_graph <- data_obj@neighbors$nn_obj@nn.idx[, -1]
  return(nn_graph)
}

geneRMSE <- function(data_obj){
  expr_data <- as.matrix(data_obj@assays$RNA@counts)
  nn_graph <- FindNNIndices(data_obj)
  RMSE <- c()
  for (gene in 1:nrow(expr_data)){
    actual_exp <- expr_data[gene, ]
    predicted_exp <- c()
    for (cell in 1:ncol(expr_data)){
      predicted_exp[cell] <- mean(expr_data[gene, nn_graph[cell,]]) 
    }
    RMSE[gene] <- rmse(actual_exp, predicted_exp)
  }
  return(mean(RMSE))
}

cellRMSE <- function(data_obj){
  expr_data <- as.matrix(data_obj@assays$RNA@counts)
  nn_graph <- FindNNIndices(data_obj)
  RMSE <- c()
  for (cell in 1:nrow(nn_graph)){
    actual_exp <- expr_data[, cell]
    predicted_exp <- rowMeans(expr_data[, nn_graph[cell, ]])
    RMSE[cell] <- rmse(actual_exp, predicted_exp)
  }
  return(mean(RMSE))
}

ConsistencyTest <- function(data_obj, method, knn, ndims, rand_seed = T) {
  if (rand_seed){
    rm(.Random.seed, envir=globalenv())
  }
  indices <- sample(1:nrow(data_obj), size = nrow(data_obj)/2, replace = FALSE)
  subset1 <- data_obj[indices, ]
  subset2 <- data_obj[-indices, ]
  subset1 <- DimReduc(subset1, method, knn, ndims)
  subset1 <- DimReduc(subset2, method, knn, ndims)
  graph1 <- FindNNIndices(subset1)
  graph2 <- FindNNIndices(subset2)
  nn_list = c()
  for (i in 1:nrow(graph1)){
    nn_list[i] <- length(intersect(graph1[i,], graph2[i,]))
  }
  return(mean(nn_list)/knn)
}

SilhouetteScore <- function(data_obj, knn, ndims){
  data_obj <- RunUMAP(data_obj, reduction = 'reduction', dims = 1:ndims)
  data_obj <- FindNeighbors(data_obj, reduction = 'umap', k.param = knn+1, graph.name = 'nn_graph', dims = 1:2)
  data_obj <- FindClusters(data_obj, graph.name = 'nn_graph', resolution = 0.5)
  
  dist.matrix <- dist(x = Embeddings(object = data_obj[['umap']])[, 1:2])
  clusters <- data_obj$seurat_clusters
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  sl_df <- as.data.frame(sil[, c(1,3)]) %>% group_by(cluster) %>% summarise(mean_coeff = mean(sil_width), cluster_size = n())
  return(mean(sl_df$mean_coeff))
}

SpatialNeighbourhood <- function(data_obj){
  nn_graph <- data_obj@neighbors$nn_obj@nn.idx[, -1]
  distmat <- rdist(data.frame(x = data_obj@meta.data$x, y = data_obj@meta.data$y))
  mean_dists <- c()
  for (cell in 1:nrow(nn_graph)){
    mean_dists[cell] <- mean(distmat[cell, nn_graph[cell, ]])
  }
  return(mean(mean_dists))
}

# Load spatial data and create Seurat object
st_data <- read_csv(st_input) %>%
  replace(is.na(.), 0) %>%
  filter(area > quantile(area, 0.1),
  area < quantile(area, 0.9),
  n_transcripts > quantile(n_transcripts, 0.1)) %>% 
  rowid_to_column('index') %>%
  column_to_rownames('index') %>% replace(is.na(.), 0)
st_metadata <- as.data.frame((st_data[1:6]))
st_genexp <- as.data.frame(st_data[7:ncol(st_data)])

st_seur <- CreateSeuratObject(counts = t(st_genexp),
    project = 'spatial', meta.data = st_metadata)

methods <- c('raw_pca','lognorm_pca', 'sct_pca', 'raw_corral', 'lognorm_corral', 'sct_corral', 'qp_sct_pca', 'glm_sct_pca')
dr_dims <- c(10,20,30)
k_neighs <- c(10,20,30,50)

all_rows <- list()
i <- 1
for (method in methods){
  for (ndim in dr_dims){
    for (k in k_neighs){
      dr_data <- DimReduc(st_seur, dimreduc = method, ndims = ndim, knn = k)
      cell_RMSE <- cellRMSE(dr_data)
      gene_RMSE <- geneRMSE(dr_data)
      consistency_score <- ConsistencyTest(dr_data, method = method, knn = k, ndims = ndim)
      silhouette_coefficient <- SilhouetteScore(dr_data, knn = k, ndims = ndim)
      neighbourhood_score <- SpatialNeighbourhood(dr_data)
      all_rows[[i]] <- c(method, ndim, k, cell_RMSE, gene_RMSE, consistency_score, silhouette_coefficient, neighbourhood_score)
      i <- i+1
      i
    }
  }
}

results_df <- as.data.frame(do.call(rbind, all_rows))
colnames(results_df) <- c('method', 'dims', 'k', 'cell_RMSE', 'gene_RMSE', 'consistency_score', 'silhouette_coefficient', 'neighbourhood_score')

qsave(results_df, output_file)
