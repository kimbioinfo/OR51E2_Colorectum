setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(qs)
library(dplyr)
library(stringr)
library(BiocParallel)


### Loading an object
epi.sr <- qread('allPublic_epithelial.qs')
sort(table(epi.sr$study.condition))
head(sort(table(epi.sr$sample.id)))


### Processing
## Splitting into layers
epi.sr[['RNA']] <- split(epi.sr[['RNA']], f=epi.sr[['study.condition']][,'study.condition'])

## Dimension reduction
epi.sr <- NormalizeData(epi.sr)
epi.sr <- FindVariableFeatures(epi.sr, nfeatures=3000)
epi.sr <- ScaleData(epi.sr)
epi.sr <- RunPCA(epi.sr)

ElbowPlot(epi.sr, reduction='pca', ndims=30)
gc()

## Clustering without integration
epi.sr <- FindNeighbors(epi.sr, reduction='pca', dims=1:21)
epi.sr <- FindClusters(epi.sr, resolution=0.8)
epi.sr$ident <- epi.sr$seurat_clusters
epi.sr <- RunUMAP(epi.sr, reduction='pca', dims=1:21, reduction.name='umap.unintegrated')

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(epi.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_epi_umap_byCluster_beforeRPCA.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(epi.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_epi_umap_byStudy_beforeRPCA.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(epi.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_epi_umap_byCondition_beforeRPCA.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
epi.sr <- IntegrateLayers(epi.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
epi.sr[['RNA']] <- JoinLayers(epi.sr[['RNA']]) 
gc()

## Clustering with integration
epi.sr <- FindNeighbors(epi.sr, reduction='integrated.rpca', dims=1:20)
epi.sr <- FindClusters(epi.sr, resolution=3)
epi.sr$ident <- epi.sr$seurat_clusters
epi.sr <- RunUMAP(epi.sr, reduction='integrated.rpca', dims=1:20, reduction.name='umap.rpca')
gc()

epi.sr@reductions$umap.rpca@key <- 'UMAP_'
colnames(epi.sr@reductions$umap.rpca@cell.embeddings) <- c('UMAP_1','UMAP_2')

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(epi.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_epi_umap_byCluster_afterRPCA.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(epi.sr, reduction='umap.rpca', raster=F, label=T) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_epi_umap_byCluster_withLabel_afterRPCA.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(epi.sr, reduction='umap.rpca', group.by='study')
ggsave('allPublic_epi_umap_byStudy_afterRPCA.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(epi.sr, reduction='umap.rpca', group.by='condition')
ggsave('allPublic_epi_umap_byCondition_afterRPCA.png', width=450, height=350, unit='px', bg='transparent', scale=4)


### OR51E2 distribution
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output/')
FeaturePlot(epi.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=F, order=T)
ggsave('allPublic_epi_featplot_OR51E2_noLabel.png', width=700, height=600, unit='px', bg='transparent', scale=2)
FeaturePlot(epi.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=T, order=T)
ggsave('allPublic_epi_featplot_OR51E2_withLabel.png', width=700, height=600, unit='px', bg='transparent', scale=2)
FeaturePlot(epi.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=T, order=T, split.by='condition')
ggsave('allPublic_epi_featplot_OR51E2_splitbyCondition_withLabel.png', width=2000, height=600, unit='px', bg='transparent', scale=2)
VlnPlot(epi.sr, features='OR51E2') + theme(legend.position='none')
ggsave('allPublic_epi_vlnplot_OR51E2.png', width=700, height=400, unit='px', bg='transparent', scale=2)
VlnPlot(epi.sr, features='OR51E2', split.by='condition') + theme(legend.position='none')
ggsave('allPublic_epi_vlnplot_OR51E2_splitByCondition.png', width=700, height=400, unit='px', bg='transparent', scale=2)


### Cell type annotation
## Enteroendocrine
eec.genes <- c('CHGA','CHGB','NEUROD1')
FeaturePlot(epi.sr, reduction='umap.rpca', features=eec.genes, raster=F, label=F, order=T)
ggsave('allPublic_epi_featplot_endocrine.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(epi.sr, features=eec.genes, pt.size=0)
ggsave('allPublic_epi_vlnplot_endocrine.png', width=1800, height=400, unit='px', bg='transparent', scale=3)
# cluster 53


### Annotation
epi.sr$celltype.detailed <- '.'
epi.sr$celltype.detailed[which(epi.sr$ident%in%c(53))] <- 'Enteroendocrine'
epi.sr$celltype.detailed[which(!(epi.sr$ident%in%c(53)))] <- 'Others'

## Plotting
DimPlot(epi.sr, reduction='umap.rpca', group.by='celltype.detailed', raster=F) + ggtitle(NULL)
ggsave('allPublic_epi_umap_byCelltype_detailed.png', width=700, height=500,unit='px', bg='transparent', scale=3)

## Saving an object
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing/')
qsave(epi.sr, 'allPublic_epi_afterRPCA_fully_annotated.qs')
