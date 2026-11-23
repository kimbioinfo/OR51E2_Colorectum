setwd('/mnt/nas-nepi/nepi_SW/plc/single-cell/')
library(dplyr)
library(stringr)
library(Seurat)
library(SingleCellExperiment)
library(qs)
library(scran)
library(scater)
library(harmony)
library(scDblFinder)
library(BiocParallel)
options(scipen=10)

### Loading seurat objects
gse183916 <- qread('lenos/gse183916_seurat_object_filtered_unprocessed.qs')
gse132465 <- qread('gse132465/gse132465_all_samples_seurat_object_filtered_unprocessed.qs')
gse161277 <- qread('gse161277/gse161277_all_samples_seurat_object_filtered_unprocessed.qs')

### Merging datasets - GSE132465, GSE161277, GSE183916
## Modifying metadata
# GSE183916
gse183916.new <- CreateSeuratObject(counts=CreateAssayObject(gse183916@assays$RNA$counts), project='GSE183916', meta.data=NULL)
gse183916.new$condition <- gse183916$condition
gse183916.new$study <- 'GSE183916'
gse183916.new$sample.id <- paste0(gse183916.new$study,':',gse183916$sample.id)
gse183916.new$percent.mt <- PercentageFeatureSet(gse183916.new, pattern='MT-')
gse183916.new@assays$RNA@misc <- list()

gse183916.new$orig.ident <- as.factor(gse183916.new$orig.ident)
gse183916.new$condition <- as.factor(gse183916.new$condition)
gse183916.new$study <- as.factor(gse183916.new$study)
gse183916.new$sample.id <- as.factor(gse183916.new$sample.id)

# GSE132465
gse132465.new <- CreateSeuratObject(counts=CreateAssayObject(gse132465@assays$RNA$counts), project='GSE132465', meta.data=NULL)
gse132465.new@assays$RNA@key <- 'rna_'
gse132465.new@assays$RNA@misc <- list()
gse132465.new$orig.ident <- as.factor('GSE132465')

gse132465$orig.ident <- as.character(gse132465$orig.ident)
gse132465.new$condition <- '.'
gse132465.new$condition[which(endsWith(gse132465$orig.ident, 'T'))] <- 'Tumor'
gse132465.new$condition[which(endsWith(gse132465$orig.ident, 'N'))] <- 'Normal'
gse132465.new$study <- 'GSE132465'
gse132465.new$sample.id <- paste0(gse132465.new$study,':',gse132465$orig.ident)

gse132465.new$orig.ident <- as.factor(gse132465.new$orig.ident)
gse132465.new$condition <- as.factor(gse132465.new$condition)
gse132465.new$study <- as.factor(gse132465.new$study)
gse132465.new$sample.id <- as.factor(gse132465.new$sample.id)

# GSE161277
gse161277.new <- CreateSeuratObject(counts=CreateAssayObject(gse161277@assays$RNA$counts), project='GSE161277', meta.data=NULL)
gse161277.new$condition <- gse161277$condition
gse161277.new$study <- 'GSE161277'
gse161277.new$sample.id <- paste0(gse161277.new$study,':',gse161277$sample.id)
gse161277.new$percent.mt <- PercentageFeatureSet(gse161277.new, pattern='MT-')
gse161277.new@assays$RNA@misc <- list()

gse161277.new$orig.ident <- as.factor(gse161277.new$orig.ident)
gse161277.new$condition <- as.factor(gse161277.new$condition)
gse161277.new$study <- as.factor(gse161277.new$study)
gse161277.new$sample.id <- as.factor(gse161277.new$sample.id)

## Merging
npt.sr <- merge(gse132465.new, list(gse161277.new, gse183916.new))

# Deleting unused objects
rm(gse132465, gse161277, gse183916, gse132465.new, gse161277.new, gse183916.new)
gc()

### Filtering cells
npt.sr[["percent.mt"]] <- PercentageFeatureSet(npt.sr, pattern = "^MT-")
FeatureScatter(npt.sr, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position='none')
FeatureScatter(npt.sr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position='none')
VlnPlot(npt.sr, c('nCount_RNA','nFeature_RNA','percent.mt'), pt.size=0)

qsave(npt.sr, 'merged/Processing/oldPublic_merged_seuratObject_filtered_unprocessed.qs')
