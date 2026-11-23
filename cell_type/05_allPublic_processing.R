setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(qs)
library(dplyr)
library(stringr)
library(BiocParallel)


### Loading a seurat object
all.sr <- qread('allPublic_merged_seuratObject_filtered_unprocessed.qs')


### Filtering
all.sr <- subset(all.sr, subset=nCount_RNA>1000 & nFeature_RNA>500 & percent.mt<20)
all.sr <- subset(all.sr, subset=sample.id%in%names(which(table(all.sr$sample.id)>=100)))
all.sr$study.condition <- paste0(all.sr$study,'_',all.sr$condition)
gc()

### Processing
## Splitting into layers
all.sr[['RNA']] <- split(all.sr[['RNA']], f=all.sr[['study.condition']][,'study.condition'])
gc()

## Dimension reduction
all.sr <- NormalizeData(all.sr)
all.sr <- FindVariableFeatures(all.sr, nfeatures=3000)
all.sr <- ScaleData(all.sr)
all.sr <- RunPCA(all.sr)
gc()

ElbowPlot(all.sr, reduction='pca', ndims=30)

## Clustering without integration
all.sr <- FindNeighbors(all.sr, reduction='pca', dims=1:23)
all.sr <- FindClusters(all.sr, resolution=1)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='pca', dims=1:23, reduction.name='umap.unintegrated')
gc()

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.unintegrated')
DimPlot(all.sr, reduction='umap.unintegrated', group.by='study')
DimPlot(all.sr, reduction='umap.unintegrated', group.by='condition')
gc()


### Integration
## RPCA
all.sr <- IntegrateLayers(all.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
all.sr[['RNA']] <- JoinLayers(all.sr[['RNA']]) 
gc()

## Clustering with integration
all.sr <- FindNeighbors(all.sr, reduction='integrated.rpca', dims=1:30)
all.sr <- FindClusters(all.sr, resolution=1)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_afterRPCA_beforeSCDBL.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', group.by='study')
DimPlot(all.sr, reduction='umap.rpca', group.by='condition')
gc()

qsave(all.sr, '../Processing/allPublic_all_sr_temp.qs')
gc()


### Removing doublets
## To an SCE object
all.sce <- as.SingleCellExperiment(all.sr)
all.sce <- scDblFinder(all.sce, clusters='seurat_clusters', samples='sample.id')
all.sr$scDblFinder.score <- all.sce$scDblFinder.score
all.sr$scDblFinder.class <- all.sce$scDblFinder.class
gc()
qsave(all.sr, '../Processing/allPublic_all_afterRPCA_beforeSCDBL.qs')
gc()

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')

FeaturePlot(all.sr, reduction='umap.rpca', 'scDblFinder.score', raster=F) + ggtitle(NULL)
DimPlot(all.sr, reduction='umap.rpca', group.by='scDblFinder.class', raster=F) + ggtitle(NULL)  + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byDoublet_afterRPCA_beforeSCDBL.png', width=700, height=500, unit='px', bg='transparent', scale=3)
VlnPlot(all.sr, c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), group.by='scDblFinder.class', pt.size=0, raster=F)
ggsave('allPublic_all_vlnplot_doublet_features_beforeSCDBL.png', width=800, height=400, unit='px', bg='transparent', scale=3)

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
all.sr <- subset(all.sr, cells=colnames(all.sr)[all.sr$scDblFinder.class=='singlet'])
qsave(all.sr, 'allPublic_all_afterRPCA_afterSCDBL.qs')
gc()


#--------------- 2nd analysis ----------------


setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
all.sr <- qread('allPublic_all_afterRPCA_afterSCDBL.qs')

### Processing
## Splitting into layers
all.sr[['RNA']] <- split(all.sr[['RNA']], f=all.sr[['study.condition']][,'study.condition'])

## Dimension reduction
all.sr <- NormalizeData(all.sr)
all.sr <- FindVariableFeatures(all.sr, nfeatures=3000)
all.sr <- ScaleData(all.sr)
all.sr <- RunPCA(all.sr)
gc()

ElbowPlot(all.sr, reduction='pca', ndims=30)

## Clustering without integration
all.sr <- FindNeighbors(all.sr, reduction='pca', dims=1:22)
all.sr <- FindClusters(all.sr, resolution=0.8)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='pca', dims=1:22, reduction.name='umap.unintegrated')
gc()

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_beforeRPCA_afterSCDBL.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byStudy_beforeRPCA_afterSCDBL.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCondition_beforeRPCA_afterSCDBL.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
all.sr <- IntegrateLayers(all.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
all.sr[['RNA']] <- JoinLayers(all.sr[['RNA']]) 
gc()

## Clustering with integration
all.sr <- FindNeighbors(all.sr, reduction='integrated.rpca', dims=1:30)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_afterRPCA_afterSCDBL.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', group.by='study')
ggsave('allPublic_all_umap_byStudy_afterRPCA_afterSCDBL.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.rpca', group.by='condition')
ggsave('allPublic_all_umap_byCondition_afterRPCA_afterSCDBL.png', width=450, height=350, unit='px', bg='transparent', scale=4)


### Removing doublets
## To an SCE object
all.sce <- as.SingleCellExperiment(all.sr)
gc()
all.sce <- scDblFinder(all.sce, clusters='seurat_clusters', samples='sample.id')
all.sr$scDblFinder.score <- all.sce$scDblFinder.score
all.sr$scDblFinder.class <- all.sce$scDblFinder.class
gc()
#qsave(all.sr, 'allPublic_all_afterRPCA_beforeSCDBLx2.qs')

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')

#FeaturePlot(all.sr, reduction='umap.rpca', 'scDblFinder.score', raster=F) + ggtitle(NULL)
DimPlot(all.sr, reduction='umap.rpca', group.by='scDblFinder.class', raster=F) + ggtitle(NULL)  + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byDoublet_afterRPCA_beforeSCDBLx2.png', width=700, height=500, unit='px', bg='transparent', scale=3)
VlnPlot(all.sr, c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), group.by='scDblFinder.class', pt.size=0, raster=F)
ggsave('allPublic_all_vlnplot_doublet_features_beforeSCDBLx2.png', width=800, height=400, unit='px', bg='transparent', scale=3)

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
all.sr <- subset(all.sr, cells=colnames(all.sr)[all.sr$scDblFinder.class=='singlet'])
qsave(all.sr, 'allPublic_all_afterRPCA_afterSCDBLx2.qs')



#--------------- 3rd analysis ----------------


setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
all.sr <- qread('allPublic_all_afterRPCA_afterSCDBLx2.qs')
gc()

### Processing
## Splitting into layers
all.sr[['RNA']] <- split(all.sr[['RNA']], f=all.sr[['study.condition']][,'study.condition'])
gc()

## Dimension reduction
all.sr <- NormalizeData(all.sr)
all.sr <- FindVariableFeatures(all.sr, nfeatures=3000)
all.sr <- ScaleData(all.sr)
all.sr <- RunPCA(all.sr)
gc()

ElbowPlot(all.sr, reduction='pca', ndims=30)

## Clustering without integration
all.sr <- FindNeighbors(all.sr, reduction='pca', dims=1:18)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='pca', dims=1:18, reduction.name='umap.unintegrated')
gc()

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_beforeRPCA_afterSCDBLx2.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byStudy_beforeRPCA_afterSCDBLx2.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCondition_beforeRPCA_afterSCDBLx2.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
all.sr <- IntegrateLayers(all.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
all.sr[['RNA']] <- JoinLayers(all.sr[['RNA']]) 
gc()

## Clustering with integration
all.sr <- FindNeighbors(all.sr, reduction='integrated.rpca', dims=1:30)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
qsave(all.sr, 'allPublic_all_afterRPCA_afterSCDBLx2_readyToUse.qs')

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_afterRPCA_afterSCDBLx2.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', raster=F, label=T) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_withLabel_afterRPCA_afterSCDBLx2.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', group.by='study')
ggsave('allPublic_all_umap_byStudy_afterRPCA_afterSCDBLx2.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.rpca', group.by='condition')
ggsave('allPublic_all_umap_byCondition_afterRPCA_afterSCDBLx2.png', width=450, height=350, unit='px', bg='transparent', scale=4)


### Removing doublets
## To an SCE object
all.sce <- as.SingleCellExperiment(all.sr)
gc()
all.sce <- scDblFinder(all.sce, clusters='seurat_clusters', samples='sample.id')
all.sr$scDblFinder.score <- all.sce$scDblFinder.score
all.sr$scDblFinder.class <- all.sce$scDblFinder.class
gc()

qsave(all.sr, 'allPublic_all_afterRPCA_beforeSCDBLx3.qs')

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')

FeaturePlot(all.sr, reduction='umap.rpca', 'scDblFinder.score', raster=F) + ggtitle(NULL)
DimPlot(all.sr, reduction='umap.rpca', group.by='scDblFinder.class', raster=F) + ggtitle(NULL)  + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byDoublet_afterRPCA_beforeSCDBLx3.png', width=700, height=500, unit='px', bg='transparent', scale=3)
VlnPlot(all.sr, c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), group.by='scDblFinder.class', pt.size=0, raster=F)
ggsave('allPublic_all_vlnplot_doublet_features_beforeSCDBLx3.png', width=800, height=400, unit='px', bg='transparent', scale=3)

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
all.sr <- subset(all.sr, cells=colnames(all.sr)[all.sr$scDblFinder.class=='singlet'])
qsave(all.sr, 'allPublic_all_afterRPCA_afterSCDBLx3.qs')



#--------------- 4th analysis ----------------


setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
gc()

### Processing
## Splitting into layers
all.sr[['RNA']] <- split(all.sr[['RNA']], f=all.sr[['study.condition']][,'study.condition'])
gc()

## Dimension reduction
all.sr <- NormalizeData(all.sr)
all.sr <- FindVariableFeatures(all.sr, nfeatures=3000)
all.sr <- ScaleData(all.sr)
all.sr <- RunPCA(all.sr)
gc()

ElbowPlot(all.sr, reduction='pca', ndims=30)

## Clustering without integration
all.sr <- FindNeighbors(all.sr, reduction='pca', dims=1:21)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='pca', dims=1:21, reduction.name='umap.unintegrated')
gc()

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_beforeRPCA_afterSCDBLx3.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byStudy_beforeRPCA_afterSCDBLx3.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCondition_beforeRPCA_afterSCDBLx3.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
all.sr <- IntegrateLayers(all.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
all.sr[['RNA']] <- JoinLayers(all.sr[['RNA']]) 
gc()

## Clustering with integration
all.sr <- FindNeighbors(all.sr, reduction='integrated.rpca', dims=1:30)
all.sr <- FindClusters(all.sr, resolution=1.4)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()


all.sr <- subset(all.sr, subset=ident!=48)

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
gc()

### Processing
## Splitting into layers
all.sr[['RNA']] <- split(all.sr[['RNA']], f=all.sr[['study.condition']][,'study.condition'])
gc()

## Dimension reduction
all.sr <- NormalizeData(all.sr)
all.sr <- FindVariableFeatures(all.sr, nfeatures=3000)
all.sr <- ScaleData(all.sr)
all.sr <- RunPCA(all.sr)
gc()

ElbowPlot(all.sr, reduction='pca', ndims=30)

## Clustering without integration
all.sr <- FindNeighbors(all.sr, reduction='pca', dims=1:23)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='pca', dims=1:23, reduction.name='umap.unintegrated')
gc()

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_beforeRPCA_afterSCDBLx3.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byStudy_beforeRPCA_afterSCDBLx3.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCondition_beforeRPCA_afterSCDBLx3.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
all.sr <- IntegrateLayers(all.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
all.sr[['RNA']] <- JoinLayers(all.sr[['RNA']]) 
gc()

## Clustering with integration
all.sr <- FindNeighbors(all.sr, reduction='integrated.rpca', dims=1:30)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()


setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
qsave(all.sr, 'allPublic_all_afterRPCA_afterSCDBLx3_readyToUse.qs')
#all.sr <- qread('allPublic_all_afterRPCA_afterSCDBLx3_readyToUse.qs')

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_afterRPCA_afterSCDBLx3.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', raster=F, label=T) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_withLabel_afterRPCA_afterSCDBLx3.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', group.by='study')
ggsave('allPublic_all_umap_byStudy_afterRPCA_afterSCDBLx3.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.rpca', group.by='condition')
ggsave('allPublic_all_umap_byCondition_afterRPCA_afterSCDBLx3.png', width=450, height=350, unit='px', bg='transparent', scale=4)


all.sr <- subset(all.sr, subset=ident!=23)

setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
gc()

### Processing
## Splitting into layers
all.sr[['RNA']] <- split(all.sr[['RNA']], f=all.sr[['study.condition']][,'study.condition'])
gc()

## Dimension reduction
all.sr <- NormalizeData(all.sr)
all.sr <- FindVariableFeatures(all.sr, nfeatures=3000)
all.sr <- ScaleData(all.sr)
all.sr <- RunPCA(all.sr)
gc()

ElbowPlot(all.sr, reduction='pca', ndims=30)

## Clustering without integration
all.sr <- FindNeighbors(all.sr, reduction='pca', dims=1:23)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='pca', dims=1:23, reduction.name='umap.unintegrated')
gc()

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_beforeRPCA_afterSCDBLx3.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byStudy_beforeRPCA_afterSCDBLx3.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCondition_beforeRPCA_afterSCDBLx3.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
all.sr <- IntegrateLayers(all.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
all.sr[['RNA']] <- JoinLayers(all.sr[['RNA']]) 
gc()

## Clustering with integration
all.sr <- FindNeighbors(all.sr, reduction='integrated.rpca', dims=1:30)
all.sr <- FindClusters(all.sr, resolution=0.6)
all.sr$ident <- all.sr$seurat_clusters
all.sr <- RunUMAP(all.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()


setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
qsave(all.sr, 'allPublic_all_afterRPCA_afterSCDBLx3_readyToUse.qs')
#all.sr <- qread('allPublic_all_afterRPCA_afterSCDBLx3_readyToUse.qs')

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(all.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_afterRPCA_afterSCDBLx3.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', raster=F, label=T) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_all_umap_byCluster_withLabel_afterRPCA_afterSCDBLx3.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(all.sr, reduction='umap.rpca', group.by='study')
ggsave('allPublic_all_umap_byStudy_afterRPCA_afterSCDBLx3.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(all.sr, reduction='umap.rpca', group.by='condition')
ggsave('allPublic_all_umap_byCondition_afterRPCA_afterSCDBLx3.png', width=450, height=350, unit='px', bg='transparent', scale=4)



FeaturePlot(all.sr, reduction='umap.rpca', features='OR51E2', raster=F, order=T)
FeaturePlot(all.sr, reduction='umap.rpca', features='RSPO2', raster=F, order=T)
FeaturePlot(all.sr, reduction='umap.rpca', features='PLN', raster=F, order=T)
FeaturePlot(all.sr, reduction='umap.rpca', features='AKT3', raster=F, order=T)
FeaturePlot(all.sr, reduction='umap.rpca', features='AKTIP', raster=F, order=T)


### Annotation - simple
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output/')
all.sr@reductions$umap.rpca@key <- 'UMAP_'
colnames(all.sr@reductions$umap.rpca@cell.embeddings) <- c('UMAP_1','UMAP_2')

## OR51E2
FeaturePlot(all.sr, reduction='umap.rpca', features='OR51E2', raster=F, order=T)
ggsave('allPublic_all_featplot_or51e2.png', width=700, height=600, unit='px', bg='transparent', scale=3)
VlnPlot(all.sr, features='OR51E2') + theme(legend.position='none')

## Red blood cell
FeaturePlot(all.sr, reduction='umap.rpca', features=c('HBB','HBA1'), raster=T, order=T)
ggsave('allPublic_all_featplot_rbc.png', width=1300, height=600, unit='px', bg='transparent', scale=2)
VlnPlot(all.sr, features=c('HBB','HBA1'), pt.size=0) + theme(legend.position='none')
ggsave('allPublic_all_vlnplot_rbc.png', width=800, height=400, unit='px', bg='transparent', scale=3)
# No cluster

## Epithelial cells
FeaturePlot(all.sr, reduction='umap.rpca', features=c('EPCAM','CDH1'), raster=T, order=T)
ggsave('allPublic_all_featplot_epithelial.png', width=1300, height=600, unit='px', bg='transparent', scale=2)
VlnPlot(all.sr, features=c('EPCAM','CDH1'), pt.size=0) + theme(legend.position='none')
ggsave('allPublic_all_vlnplot_epithelial.png', width=1200, height=400, unit='px', bg='transparent', scale=4)
# Cluster 0,2,4,5,7,8,10,11,14,18,21,22,23,24,26,28

## Endothelial cells
FeaturePlot(all.sr, reduction='umap.rpca', features=c('PECAM1','CDH5'), raster=T, order=T)
ggsave('allPublic_all_featplot_endothelial.png', width=1300, height=600, unit='px', bg='transparent', scale=2)
VlnPlot(all.sr, features=c('PECAM1','CDH5'), pt.size=0) + theme(legend.position='none')
ggsave('allPublic_all_vlnplot_endothelial.png', width=1200, height=400, unit='px', bg='transparent', scale=4)
# Cluster 19

## Immune cells
# General markers
FeaturePlot(all.sr, reduction='umap.rpca', features=c('PTPRC'), raster=F, order=T)
# T cell markers
FeaturePlot(all.sr, reduction='umap.rpca', features=c('CD3D','CD3E'), raster=F, order=T)
# B cell markers
FeaturePlot(all.sr, reduction='umap.rpca', features=c('MS4A1','CD79A'), raster=F, order=T)
# Plasma cell markers
FeaturePlot(all.sr, reduction='umap.rpca', features=c('SDC1','MZB1'), raster=F, order=T)
# Mast cell markers
FeaturePlot(all.sr, reduction='umap.rpca', features=c('GATA2','CPA3','HPGDS'), raster=F, order=T)
# Macrophage markers
FeaturePlot(all.sr, reduction='umap.rpca', features=c('CD163','C1QB','C1QC'), raster=F, order=T)
# All
FeaturePlot(all.sr, reduction='umap.rpca', features=c('PTPRC','CD3D','CD3E','MS4A1','CD79A','SDC1','MZB1','CPA3','HPGDS','CD163','C1QB','C1QC'), raster=T, order=T)
ggsave('allPublic_all_featplot_immune.png', width=2600, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(all.sr, features=c('PTPRC','CD3D','CD3E','MS4A1','CD79A','SDC1','MZB1','CPA3','HPGDS','CD163','C1QB','C1QC'), pt.size=0) + theme(legend.position='none')
ggsave('allPublic_all_vlnplot_immune.png', width=2400, height=1200, unit='px', bg='transparent', scale=4)
# Cluster 1,3,6,9,13,16,17,29,30,31,32,33

## Neural cells
FeaturePlot(all.sr, reduction='umap.rpca', features=c('SCN7A','GPR155','SPP1','LGI4'), raster=F, order=T)
ggsave('allPublic_all_featplot_neural.png', width=1300, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(all.sr, features=c('SCN7A','GPR155','SPP1','LGI4'), pt.size=0) + theme(legend.position='none')
ggsave('allPublic_all_vlnplot_neural.png', width=1800, height=800, unit='px', bg='transparent', scale=4)
# Cluster 27

## Mesenchymal cells
FeaturePlot(all.sr, reduction='umap.rpca', features=c('VIM','PDGFRA','PDGFRB','ACTA2'), raster=T, order=T)
ggsave('allPublic_all_featplot_mesenchymal.png', width=1300, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(all.sr, features=c('VIM','PDGFRA','PDGFRB','ACTA2'), pt.size=0) + theme(legend.position='none')
ggsave('allPublic_all_vlnplot_mesenchymal.png', width=1800, height=800, unit='px', bg='transparent', scale=4)
# Cluster 12,15,20,25

### Annotation
all.sr$celltype.simple <- '.'
all.sr$celltype.simple[which(all.sr$ident%in%c(0,2,4,5,7,8,10,11,14,18,21,22,23,24,26,28))] <- 'Epithelial'
all.sr$celltype.simple[which(all.sr$ident%in%c(19))] <- 'Endothelial'
all.sr$celltype.simple[which(all.sr$ident%in%c(1,3,6,9,13,16,17,29,30,31,32,33))] <- 'Immune'
all.sr$celltype.simple[which(all.sr$ident%in%c(27))] <- 'Neural'
all.sr$celltype.simple[which(all.sr$ident%in%c(12,15,20,25))] <- 'Mesenchymal'


## Plotting
DimPlot(all.sr, reduction='umap.rpca', group.by='celltype.simple', raster=F) + ggtitle(NULL)
ggsave('allPublic_all_umap_byCelltype_simple.png', width=700, height=500,unit='px', bg='transparent', scale=3)

## Saving an object
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
qsave(all.sr, 'allPublic_all_afterRPCA_afterRemovingDoublet_readyToUse.qs')


### Subsetting
epi.sr <- subset(all.sr, subset=celltype.simple=='Epithelial')
mes.sr <- subset(all.sr, subset=celltype.simple=='Mesenchymal')
immune.sr <- subset(all.sr, subset=celltype.simple=='Immune')

## Saving an object
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
qsave(all.sr, 'allPublic_all_afterRPCA_afterRemovingDoublet_simply_annotated.qs')
qsave(epi.sr, 'allPublic_epithelial.qs')
qsave(mes.sr, 'allPublic_mesenchymal.qs')
qsave(immune.sr, 'allPublic_immune.qs')

