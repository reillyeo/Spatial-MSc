args <- commandArgs(TRUE)
st_input <- args[1]
sc_input <- args[2]
output_file <- args[3]

library(Seurat)
library(qs)
library(corral)
library(tidyverse)

st.query <- qread(st_input)
sc.reference <- qread(sc_input)

# find common genes between reference and query datasets and subset objects to only these genes
common_genes <- intersect(rownames(sc.reference), rownames(st.query))
sc.reference <- subset(sc.reference, features = common_genes)
st.query <- subset(st.query, features = common_genes)

# remove cells with zero counts
st.query=st.query[,unname(which(colSums(GetAssayData(st.query))!=0))]
sc.reference=sc.reference[,unname(which(colSums(GetAssayData(sc.reference))!=0))]

# run SCTransform on both objects
st.query <- SCTransform(st.query, method = 'glmGamPoi')
sc.reference <- SCTransform(sc.reference, method = 'glmGamPoi')

# run dimension reduction
st.query <- RunPCA(st.query, features = rownames(st.query))

anchorset <- FindTransferAnchors(reference = sc.reference, query = st.query, scale = T, reference.assay = 'SCT', recompute.residuals = T,query.assay = 'SCT', reduction = 'cca', dims = 1:20, features = common_genes, normalization.method = 'SCT')

# create celltype prediction for spatial data based on transfer anchors and add as metadata
celltype.predictions <- TransferData(anchorset, refdata = sc.reference$celltype_labels,
                                   weight.reduction = st.query[['pca']], dims = 1:30)

celltype.predictions <- celltype.predictions %>%
  mutate(predictions_sct = predicted.id, prediction_scores_sct = prediction.score.max) %>%
  select(c(predictions_sct, prediction_scores_sct))

st.query <- AddMetaData(st.query, metadata = celltype.predictions)

qsave(st.query, output_file)