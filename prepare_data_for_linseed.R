library(Seurat)
library(dplyr)
library(Matrix)

load('data/SRA779509_SRS3805254.sparse.RData')

bone_m <- CreateSeuratObject(counts = sm)
bone_m <- NormalizeData(object = bone_m)
bone_m <- FindVariableFeatures(object = bone_m, nfeatures = 5000)

bone_m_raw <- CreateSeuratObject(counts = sm)
bone_m_raw <- FindVariableFeatures(object = bone_m_raw, nfeatures = 5000)

var_genes <- rownames(GetAssayData(object = bone_m, slot = "scale.data"))

matrix_linseed_raw <- bone_m[rownames(GetAssayData(object = bone_m, slot = "counts")) %in% var_genes]
matrix_linseed_raw <- GetAssayData(object = matrix_linseed_raw, slot = "counts")


matrix_linseed_norm <- bone_m_raw[rownames(GetAssayData(object = bone_m_raw, slot = "counts")) %in% var_genes]
matrix_linseed_norm <- GetAssayData(object = matrix_linseed_norm, slot = "counts")

saveRDS(matrix_linseed_norm, 'results/normalized_bone_marrow.rds')
saveRDS(matrix_linseed_raw, 'results/raw_bone_marrow.rds')
