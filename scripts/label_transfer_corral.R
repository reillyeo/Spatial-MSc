args <- commandArgs(TRUE)
st_input <- args[1]
sc_input <- args[2]
output_file <- args[3]

library(Seurat)
library(qs)
library(corral)
library(tidyverse)

RunCorral <- function(seur_obj, counts = T, alt = NULL, n = 30, residual = 'standardized', reduction.name = 'corral'){
  if (counts) {
    dimred <- corral(seur_obj@assays$RNA@counts, ncomp = n, rtype = residual)
  }  else {
    dimred <- corral(alt, ncomp = n, rtype = residual)
  }
  seur_obj[[reduction.name]] <- CreateDimReducObject(embeddings = dimred$v, loadings = dimred$u, key = reduction.name)
  return(seur_obj)
}

st.query <- qread(st_input)
sc.reference <- qread(sc_input)

# find common genes between reference and query datasets and subset objects to only these genes
common_genes <- intersect(rownames(sc.reference), rownames(st.query))
sc.reference <- subset(sc.reference, features = common_genes)
st.query <- subset(st.query, features = common_genes)

# remove cells with zero counts
st.query=st.query[,unname(which(colSums(GetAssayData(st.query))!=0))]
sc.reference=sc.reference[,unname(which(colSums(GetAssayData(sc.reference))!=0))]

# run corral preprocessing on both objects
st_preproc <- corral_preproc(st.query@assays$RNA@counts)
sc_preproc <- corral_preproc(sc.reference@assays$RNA@counts)
sc.reference[['corral_pp']] <- CreateAssayObject(counts = sc_preproc)

st.query[['corral_pp']] <- CreateAssayObject(counts = st_preproc)

# add pre-processed data to scale.data slot and make default assay
DefaultAssay(sc.reference) <- 'corral_pp'
DefaultAssay(st.query) <- 'corral_pp'

# run dimension reduction
st.query <- RunCorral(st.query)

anchorset <- FindTransferAnchors(reference = sc.reference, query = st.query, scale = T, reference.assay = DefaultAssay(sc.reference),
                                 query.assay = DefaultAssay(st.query), reduction = 'cca', dims = 1:20, features = common_genes)

# create celltype prediction for spatial data based on transfer anchors and add as metadata
celltype.predictions <- TransferData(anchorset, refdata = sc.reference$celltype_labels,
                                   weight.reduction = st.query[['corral']], dims = 1:20)

celltype.predictions <- celltype.predictions %>%
  mutate(predictions_corral = predicted.id, prediction_scores_corral = prediction.score.max) %>%
  select(c(predictions_corral, prediction_scores_corral))

st.query <- AddMetaData(st.query, metadata = celltype.predictions)



qsave(st.query, output_file)