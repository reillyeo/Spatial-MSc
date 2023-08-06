args <- commandArgs(TRUE)
st_input <- args[1]
sc_input <- args[2]
output_file <- args[3]

library(Seurat)
library(SeuratWrappers)
library(qs)
library(corral)
library(tidyverse)
library(HGNChelper)

st.data <- qread(st_input)
sc.data <- qread(sc_input)

set_prep <- function(db, cell_type){
    cell_markers = db
    cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    # correct gene symbols from the given DB (up-genes)
    cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
      
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
        paste0(markers_all, collapse=",")
      } else {
        ""
      }
    })
    
    # correct gene symbols from the given DB (down-genes)
    cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
      
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
        paste0(markers_all, collapse=",")
      } else {
        ""
      }
    })
    
    cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
     
    gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
    gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
    
    list(gs_positive = gs, gs_negative = gs2)
}


sc.data <- RunALRA(sc.data)

Idents(object = sc.data) <- sc.data@meta.data$reflabels

sc_markers <- FindAllMarkers(sc.data,assay = 'alra')
top10 <- sc_markers %>% group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)
bot10 <- sc_markers %>%  group_by(cluster) %>% 
  top_n(n = -5, wt = avg_log2FC)

marker_lists <- list()
for (label in unique(top10$cluster)){
  top <- top10 %>% filter(cluster == label)
  bottom <- bot10 %>% filter(cluster == label)
  marker_lists[[label]] <- c("Arc", label, (toString(top$gene)), (toString(bottom$gene)))
}


marker_db <- t(data.frame(marker_lists))
colnames(marker_db) <- c('tissueType', 'cellName', 'geneSymbolmore1', "geneSymbolmore2")
rownames(marker_db) <- NULL

genelist <- as.data.frame(marker_db)

# load libraries and functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from database
gs_list = set_prep(genelist) # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

st.data <- RunUMAP(st.data, reduction = 'corral', dims = 1:30, reduction.name = 'corral_umap')
st.data <- FindNeighbors(st.data, reduction = 'corral', k.param = 15, dims = 1:30)
st.data <- FindClusters(st.data, resolution = 1.5, graph.name = 'RNA_nn')

# assign cell types
st.data <- RunALRA(st.data)
data_matrix <- as.matrix(st.data@assays$alra@data)
es.max = sctype_score(scRNAseqData = data_matrix, scaled = F, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(st.data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(st.data@meta.data[st.data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(st.data@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

st.data@meta.data$label1 = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  st.data@meta.data$label1[st.data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  st.data@meta.data$score1[st.data@meta.data$seurat_clusters == j] = as.numeric(cl_type$scores[1])
}

st.data <- RunUMAP(st.data, reduction = 'sct_glmGamPoi_pca', dims = 1:30, reduction.name = 'sct_umap')
st.data <- FindNeighbors(st.data, reduction = 'sct_glmGamPoi_pca', k.param = 20, dims = 1:30, graph.name = 'SCT_snn')
st.data <- FindClusters(st.data, resolution = 1.5, graph.name = 'SCT_snn')

# assign cell types
es.max = sctype_score(scRNAseqData = data_matrix, scaled = F, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(st.data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(st.data@meta.data[st.data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(st.data@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

st.data@meta.data$label2 = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  st.data@meta.data$label2[st.data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  st.data@meta.data$score2[st.data@meta.data$seurat_clusters == j] = as.numeric(cl_type$scores[1])
}

st.data@meta.data <- st.data@meta.data %>% mutate(final_label = ifelse(score1 > score2, label1, label2))
            

qsave(st.data, output_file)