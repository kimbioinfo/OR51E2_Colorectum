setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(qs)
library(dplyr)
library(stringr)
library(BiocParallel)
library(RColorBrewer)


### Loading an object
mes.sr <- qread('allPublic_mesenchymal.qs')
sort(table(mes.sr$study.condition))
sort(table(mes.sr$sample.id))
above10.samples <- names(table(mes.sr$sample.id)[table(mes.sr$sample.id)>10])
mes.sr <- subset(mes.sr, subset=sample.id%in%above5.samples)
#mes.sr <- subset(mes.sr, subset=study.condition!='GSE161277_Polyp')
#mes.sr <- subset(mes.sr, subset=study.condition!='GSE161277_Tumor')


### Processing
## Splitting into layers
mes.sr[['RNA']] <- split(mes.sr[['RNA']], f=mes.sr[['study.condition']][,'study.condition'])

## Dimension reduction
mes.sr <- NormalizeData(mes.sr)
mes.sr <- FindVariableFeatures(mes.sr, nfeatures=3000)
mes.sr <- ScaleData(mes.sr)
mes.sr <- RunPCA(mes.sr)

ElbowPlot(mes.sr, reduction='pca', ndims=30)
gc()

## Clustering without integration
mes.sr <- FindNeighbors(mes.sr, reduction='pca', dims=1:18)
mes.sr <- FindClusters(mes.sr, resolution=0.8)
mes.sr$ident <- mes.sr$seurat_clusters
mes.sr <- RunUMAP(mes.sr, reduction='pca', dims=1:18, reduction.name='umap.unintegrated')

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(mes.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCluster_beforeRPCA.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(mes.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byStudy_beforeRPCA.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(mes.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCondition_beforeRPCA.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
mes.sr <- IntegrateLayers(mes.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
mes.sr[['RNA']] <- JoinLayers(mes.sr[['RNA']]) 
gc()

## Clustering with integration
mes.sr <- FindNeighbors(mes.sr, reduction='integrated.rpca', dims=1:30)
mes.sr <- FindClusters(mes.sr, resolution=3)
mes.sr$ident <- mes.sr$seurat_clusters
mes.sr <- RunUMAP(mes.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()

mes.sr@reductions$umap.rpca@key <- 'UMAP_'
colnames(mes.sr@reductions$umap.rpca@cell.embeddings) <- c('UMAP_1','UMAP_2')

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(mes.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCluster_afterRPCA.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(mes.sr, reduction='umap.rpca', raster=F, label=T) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCluster_withLabel_afterRPCA.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(mes.sr, reduction='umap.rpca', group.by='study')
ggsave('allPublic_mes_umap_byStudy_afterRPCA.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(mes.sr, reduction='umap.rpca', group.by='condition')
ggsave('allPublic_mes_umap_byCondition_afterRPCA.png', width=450, height=350, unit='px', bg='transparent', scale=4)



### Distribution of OR51E2 gene
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output/')
## OR51E2
FeaturePlot(mes.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=F, order=T)
ggsave('allPublic_mes_featplot_OR51E2_noLabel.png', width=700, height=600, unit='px', bg='transparent', scale=2)
FeaturePlot(mes.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=T, order=T)
ggsave('allPublic_mes_featplot_OR51E2_withLabel.png', width=700, height=600, unit='px', bg='transparent', scale=2)
FeaturePlot(mes.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=T, order=T, split.by='condition')
ggsave('allPublic_mes_featplot_OR51E2_splitbyCondition_withLabel.png', width=2000, height=600, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features='OR51E2') + theme(legend.position='none')
ggsave('allPublic_mes_vlnplot_OR51E2.png', width=700, height=400, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features='OR51E2', split.by='condition') + theme(legend.position='none')
ggsave('allPublic_mes_vlnplot_OR51E2_splitByCondition.png', width=700, height=400, unit='px', bg='transparent', scale=2)


### Annotation
## Stromal cells
# Stromal 1
stromal1.genes <- c('ADAMDEC1','ADAM28','CCL8','CCL13','SFTA1P')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal1.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal1.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal1.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal1.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 1,2,26,36

# Stromal 2
stromal2.genes <- c('CDH13','PLK2','F3','PTPRE')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal2.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal2.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal2.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal2.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 12,20

# Stromal 3
stromal3.genes <- c('C7','KCNN3','JAM2')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal3.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal3.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal3.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal3.png', width=1800, height=400, unit='px', bg='transparent', scale=3)
# cluster 10

# Transitional Stromal 3 (C3+)
stromal3.trans.genes <- c('C3','MFAP5','FBN1','SFRP2')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal3.trans.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal3_trans.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal3.trans.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal3_trans.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 3,5,7,8,9,11,15,19,28,31,32,33,34,37

## Myofibroblast
# Myofibroblast (general)
myof.genes <- c('ACTA2','TAGLN','DCN','NPNT','HHIP','SOSTDC1')
FeaturePlot(mes.sr, reduction='umap.rpca', features=myof.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_myof.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=myof.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_myof.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 25

# Myofibroblast (RSPO2+)
myof.genes3 <- c('ACTA2','RSPO2','HHIP','NPNT')
FeaturePlot(mes.sr, reduction='umap.rpca', features=myof.genes3, raster=F, order=T)
ggsave('allPublic_mes_featplot_myof_rspo2.png', width=1500, height=1200, unit='px', bg='transparent', scale=3)
VlnPlot(mes.sr, features=myof.genes3)
ggsave('allPublic_mes_vlnplot_myof_rspo2.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 24

## Pericyte
# Pericyte (general)
peri.genes <- c('NOTCH3','PDGFRB','CD36','ARHGDIB','STEAP4')
FeaturePlot(mes.sr, reduction='umap.rpca', features=peri.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_pericyte.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=peri.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_pericyte.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 13

# Contractile pericyte
contractile.genes <- c('NOTCH3','PDGFRB','PLN','NET1','RERGL','KCNA5')
FeaturePlot(mes.sr, reduction='umap.rpca', features=contractile.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_contractile_pericyte.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=contractile.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_contractile_pericyte.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 14,30

## SMC
smc.genes <- c('ACTA2','MYH11','TAGLN','DES','ACTG2','CHRM2','PART1','LBR','ENAH','HACD1','MTCL1','NCAM1')
FeaturePlot(mes.sr, reduction='umap.rpca', features=smc.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_smc.png', width=2900, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=smc.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_smc.png', width=2400, height=1200, unit='px', bg='transparent', scale=3)
# cluster 0,4,6,6,16,17,18,21,22,27

## ICC
icc.genes <- c('ETV1','FHL2','KIT','SNCA')
FeaturePlot(mes.sr, reduction='umap.rpca', features=icc.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_icc.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=icc.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_icc.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 35

## T reticular & mLTo
tret.mlto.genes <- c('GREM1','PTGDS','CCL21','CCL19','TNFSF11','SLC22A3')
FeaturePlot(mes.sr, reduction='umap.rpca', features=tret.mlto.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_tret_mlto.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=tret.mlto.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_tret_mlto.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 23

## FDC
fdc.genes <- c('CR2','SRGN','FDCSP','CR1')
FeaturePlot(mes.sr, reduction='umap.rpca', features=fdc.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_fdc.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=fdc.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_fdc.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 29


### Annotation
mes.sr$celltype.detailed <- '.'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(1,2,26,36))] <- 'Stromal 1'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(12,20))] <- 'Stromal 2'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(10))] <- 'Stromal 3'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(3,5,7,8,9,11,15,19,28,32,33))] <- 'Transitional Stromal 3'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(25))] <- 'Myofibroblast'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(21))] <- 'Myofibroblast (RSPO2+)'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(13))] <- 'Pericyte'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(15,34))] <- 'Contractile pericyte'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(0,2,3,13,14,18,22,24,28))] <- 'Smooth muscle cell'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(35))] <- 'ICC'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(23))] <- 'T reticular & mLTo'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(29))] <- 'FDC'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(31,34,37))] <- 'Doublet'

## Plotting
DimPlot(mes.sr, reduction='umap.rpca', group.by='celltype.detailed', raster=F) + ggtitle(NULL)
ggsave('allPublic_mes_umap_byCelltype_detailed.png', width=700, height=500,unit='px', bg='transparent', scale=3)

### Removing doublets
## Manually
library(tibble)
# cluster 34
markers.c34 <- FindMarkers(mes.sr, ident.1=34, only.pos=T, logfc.threshold=5) %>% 
  rownames_to_column('geneName')
head(markers.c34, 20)
markers.c34.pct <- FindMarkers(mes.sr, ident.1=34, only.pos=T, logfc.threshold=5) %>% 
  rownames_to_column('geneName') %>% 
  arrange(desc(pct.1))
FeaturePlot(mes.sr, reduction='umap.rpca', features=markers.c34.pct$geneName[1:12], raster=F, order=T)

# cluster 41
markers.c41 <- FindMarkers(mes.sr, ident.1=41, only.pos=T, logfc.threshold=5) %>% 
  rownames_to_column('geneName')
head(markers.c41, 20)
markers.c41.pct <- FindMarkers(mes.sr, ident.1=41, only.pos=T, logfc.threshold=5) %>% 
  rownames_to_column('geneName') %>% 
  arrange(desc(pct.1))
FeaturePlot(mes.sr, reduction='umap.rpca', features=markers.c41.pct$geneName[1:12], raster=F, order=T)

# cluster 37
markers.c37 <- FindMarkers(mes.sr, ident.1=37, only.pos=T, logfc.threshold=5) %>% 
  rownames_to_column('geneName')
head(markers.c37, 20)
markers.c37.pct <- FindMarkers(mes.sr, ident.1=37, only.pos=T, logfc.threshold=5) %>% 
  rownames_to_column('geneName') %>% 
  arrange(desc(pct.1))
FeaturePlot(mes.sr, reduction='umap.rpca', features=markers.c37.pct$geneName[1:12], raster=F, order=T)

mes.sr <- subset(mes.sr, subset=ident!=34)
mes.sr <- subset(mes.sr, subset=ident!=35)
mes.sr <- subset(mes.sr, subset=ident!=40)
mes.sr <- subset(mes.sr, subset=ident!=41)
mes.sr <- subset(mes.sr, subset=ident!=37)


#----------------- 2nd analysis -----------------


sort(table(mes.sr$study.condition))


### Processing
## Splitting into layers
mes.sr[['RNA']] <- split(mes.sr[['RNA']], f=mes.sr[['study.condition']][,'study.condition'])

## Dimension reduction
mes.sr <- NormalizeData(mes.sr)
mes.sr <- FindVariableFeatures(mes.sr, nfeatures=3000)
mes.sr <- ScaleData(mes.sr)
mes.sr <- RunPCA(mes.sr)

ElbowPlot(mes.sr, reduction='pca', ndims=30)
gc()

## Clustering without integration
mes.sr <- FindNeighbors(mes.sr, reduction='pca', dims=1:23)
mes.sr <- FindClusters(mes.sr, resolution=0.8)
mes.sr$ident <- mes.sr$seurat_clusters
mes.sr <- RunUMAP(mes.sr, reduction='pca', dims=1:23, reduction.name='umap.unintegrated')

## Checking batch effects
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(mes.sr, reduction='umap.unintegrated', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCluster_beforeRPCA_afterRemovingDoublet.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(mes.sr, reduction='umap.unintegrated', group.by='study') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byStudy_beforeRPCA_afterRemovingDoublet.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(mes.sr, reduction='umap.unintegrated', group.by='condition') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCondition_beforeRPCA_afterRemovingDoublet.png', width=450, height=350, unit='px', bg='transparent', scale=4)
gc()


### Integration
## RPCA
mes.sr <- IntegrateLayers(mes.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
mes.sr[['RNA']] <- JoinLayers(mes.sr[['RNA']]) 
gc()

## Clustering with integration
mes.sr <- FindNeighbors(mes.sr, reduction='integrated.rpca', dims=1:30)
mes.sr <- FindClusters(mes.sr, resolution=3)
mes.sr$ident <- mes.sr$seurat_clusters
mes.sr <- RunUMAP(mes.sr, reduction='integrated.rpca', dims=1:30, reduction.name='umap.rpca')
gc()

mes.sr@reductions$umap.rpca@key <- 'UMAP_'
colnames(mes.sr@reductions$umap.rpca@cell.embeddings) <- c('UMAP_1','UMAP_2')

## Checking results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output')
DimPlot(mes.sr, reduction='umap.rpca', raster=F) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCluster_afterRPCA_afterRemovingDoublet.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(mes.sr, reduction='umap.rpca', raster=F, label=T) + theme(legend.position='none') + xlab('UMAP_1') + ylab('UMAP_2')
ggsave('allPublic_mes_umap_byCluster_withLabel_afterRPCA_afterRemovingDoublet.png', width=700, height=600, unit='px', bg='transparent', scale=3)
DimPlot(mes.sr, reduction='umap.rpca', group.by='study')
ggsave('allPublic_mes_umap_byStudy_afterRPCA_afterRemovingDoublet.png', width=500, height=350, unit='px', bg='transparent', scale=4)
DimPlot(mes.sr, reduction='umap.rpca', group.by='condition')
ggsave('allPublic_mes_umap_byCondition_afterRPCA_afterRemovingDoublet.png', width=450, height=350, unit='px', bg='transparent', scale=4)


### OR51E2 distribution
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Output/')
FeaturePlot(mes.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=F, order=T)
ggsave('allPublic_mes_featplot_OR51E2_noLabel.png', width=700, height=600, unit='px', bg='transparent', scale=2)
FeaturePlot(mes.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=T, order=T)
ggsave('allPublic_mes_featplot_OR51E2_withLabel.png', width=700, height=600, unit='px', bg='transparent', scale=2)
FeaturePlot(mes.sr, reduction='umap.rpca', features='OR51E2', raster=F, label=T, order=T, split.by='condition')
ggsave('allPublic_mes_featplot_OR51E2_splitbyCondition_withLabel.png', width=2000, height=600, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features='OR51E2') + theme(legend.position='none')
ggsave('allPublic_mes_vlnplot_OR51E2.png', width=700, height=400, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features='OR51E2', split.by='condition') + theme(legend.position='none')
ggsave('allPublic_mes_vlnplot_OR51E2_splitByCondition.png', width=700, height=400, unit='px', bg='transparent', scale=2)


### Cell type annotation
## Stromal cells
# Stromal 1
stromal1.genes <- c('ADAMDEC1','ADAM28','CCL8','CCL13','SFTA1P')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal1.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal1.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal1.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal1.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 0,13,24,28,31,39

# Stromal 2
stromal2.genes <- c('CDH13','PLK2','F3','PTPRE')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal2.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal2.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal2.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal2.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 15,23

# Stromal 3
stromal3.genes <- c('C7','KCNN3','JAM2')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal3.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal3.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal3.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal3.png', width=1800, height=400, unit='px', bg='transparent', scale=3)
# cluster 8

# Transitional Stromal 3 (C3+)
stromal3.trans.genes <- c('C3','MFAP5','FBN1','SFRP2')
FeaturePlot(mes.sr, reduction='umap.rpca', features=stromal3.trans.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_stromal3_trans.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=stromal3.trans.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_stromal3_trans.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 3,4,5,6,9,12,14,20,29,30,33,34,37

## Myofibroblast
# Myofibroblast (general)
myof.genes <- c('ACTA2','TAGLN','DCN','NPNT','HHIP','SOSTDC1')
FeaturePlot(mes.sr, reduction='umap.rpca', features=myof.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_myof.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=myof.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_myof.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 27

# Myofibroblast (RSPO2+)
myof.genes3 <- c('ACTA2','RSPO2','HHIP','NPNT')
FeaturePlot(mes.sr, reduction='umap.rpca', features=myof.genes3, raster=F, order=T)
ggsave('allPublic_mes_featplot_myof_rspo2.png', width=1500, height=1200, unit='px', bg='transparent', scale=3)
VlnPlot(mes.sr, features=myof.genes3)
ggsave('allPublic_mes_vlnplot_myof_rspo2.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 21

## Pericyte
# Pericyte (general)
peri.genes <- c('NOTCH3','PDGFRB','CD36','ARHGDIB','STEAP4')
FeaturePlot(mes.sr, reduction='umap.rpca', features=peri.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_pericyte.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=peri.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_pericyte.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 11

# Contractile pericyte
contractile.genes <- c('NOTCH3','PDGFRB','PLN','NET1','RERGL','KCNA5')
FeaturePlot(mes.sr, reduction='umap.rpca', features=contractile.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_contractile_pericyte.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=contractile.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_contractile_pericyte.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 16

## SMC
smc.genes <- c('ACTA2','MYH11','TAGLN','DES','ACTG2','CHRM2','PART1','LBR','ENAH','HACD1','MTCL1','NCAM1')
FeaturePlot(mes.sr, reduction='umap.rpca', features=smc.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_smc.png', width=2900, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=smc.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_smc.png', width=2400, height=1200, unit='px', bg='transparent', scale=3)
# cluster 1,2,7,10,17,18,19,22,25,32,36

## ICC
icc.genes <- c('ETV1','FHL2','KIT','SNCA')
FeaturePlot(mes.sr, reduction='umap.rpca', features=icc.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_icc.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=icc.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_icc.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 38

## T reticular & mLTo
tret.mlto.genes <- c('GREM1','PTGDS','CCL21','CCL19','TNFSF11','SLC22A3')
FeaturePlot(mes.sr, reduction='umap.rpca', features=tret.mlto.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_tret_mlto.png', width=1500, height=1800, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=tret.mlto.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_tret_mlto.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 26

## FDC
fdc.genes <- c('CR2','SRGN','FDCSP','CR1')
FeaturePlot(mes.sr, reduction='umap.rpca', features=fdc.genes, raster=F, order=T)
ggsave('allPublic_mes_featplot_fdc.png', width=1500, height=1200, unit='px', bg='transparent', scale=2)
VlnPlot(mes.sr, features=fdc.genes, pt.size=0)
ggsave('allPublic_mes_vlnplot_fdc.png', width=1800, height=800, unit='px', bg='transparent', scale=3)
# cluster 35


### Annotation
mes.sr$celltype.detailed <- '.'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(0,13,24,28,31,39))] <- 'Stromal 1'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(15,23))] <- 'Stromal 2'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(8))] <- 'Stromal 3'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(3,4,5,6,9,12,14,20,29,30,33,34,37))] <- 'Transitional Stromal 3'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(27))] <- 'Myofibroblast'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(21))] <- 'Myofibroblast (RSPO2+)'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(11))] <- 'Pericyte'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(16))] <- 'Contractile pericyte'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(1,2,7,10,17,18,19,22,25,32,36))] <- 'Smooth muscle cell'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(38))] <- 'ICC'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(26))] <- 'T reticular & mLTo'
mes.sr$celltype.detailed[which(mes.sr$ident%in%c(35))] <- 'FDC'

## Plotting
DimPlot(mes.sr, reduction='umap.rpca', group.by='celltype.detailed', raster=F) + ggtitle(NULL)
ggsave('allPublic_mes_umap_byCelltype_detailed.png', width=700, height=500,unit='px', bg='transparent', scale=3)

## Saving an object
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing/')
qsave(mes.sr, 'allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')



### OR51E2 cells
three.celltypes <- c('Smooth muscle cell','Contractile pericyte','Myofibroblast (RSPO2+)')
or51e2 <- WhichCells(mes.sr, expression=OR51E2>0)
total.np.clusters <- table(mes.sr@meta.data[["condition"]], mes.sr$celltype.detailed)[,three.celltypes]
or51e2.np.clusters <- table(mes.sr@meta.data[or51e2,"condition"], mes.sr@meta.data[or51e2,"celltype.detailed"])[,three.celltypes]
round(or51e2.np.clusters/total.np.clusters*100, 1)












### Cell type proportions
## Loading an object
mes.sr <- qread('allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')

## Total cell numbers
prop.total.condition <- as.data.frame(table(mes.sr$condition))
prop.total.study <- as.data.frame(table(mes.sr$study))

ggplot(prop.total.condition, aes(x='Total', y=Freq, fill=Var1)) +
  theme_classic(base_size=15) +
  geom_col(position='fill', width=0.5) +
  xlab(NULL) + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(size=0))

ggplot(prop.total.study, aes(x='Total', y=Freq, fill=Var1)) +
  theme_classic(base_size=15) +
  geom_col(position='fill', width=0.5) +
  xlab(NULL) + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(size=0))

#####
## Per sample
prop.samples <- as.data.frame(table(mes.sr$celltype.detailed, mes.sr$sample.id))
ggplot(prop.samples, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=9) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

# Considering conditions
samples.normal <- unique(mes.sr@meta.data[mes.sr@meta.data$condition=='Normal', 6])
samples.polyp <- unique(mes.sr@meta.data[mes.sr@meta.data$condition=='Polyp', 6])
samples.tumor <- unique(mes.sr@meta.data[mes.sr@meta.data$condition=='Tumor', 6])

prop.samples.normal <- prop.samples[prop.samples$Var2%in%samples.normal,]
prop.samples.polyp <- prop.samples[prop.samples$Var2%in%samples.polyp,]
prop.samples.tumor <- prop.samples[prop.samples$Var2%in%samples.tumor,]

emf('', emfPlus = FALSE)
ggplot(prop.samples.normal, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=8) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

ggplot(prop.samples.polyp, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=8) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

ggplot(prop.samples.tumor, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=8) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

# GSE132465 Matched
gse132465.normal <- unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE132465_Normal', 6])[-9]
gse132465.tumor <- unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE132465_Tumor', 6])[1:9]

prop.gse132465.normal <- prop.samples[prop.samples$Var2%in%gse132465.normal,]
prop.gse132465.tumor <- prop.samples[prop.samples$Var2%in%gse132465.tumor,]

ggplot(prop.gse132465.normal, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  scale_x_discrete(labels=as.data.frame(str_split(gse132465.normal,':',simplify=T))[,2]) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

ggplot(prop.gse132465.tumor, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  scale_x_discrete(labels=as.data.frame(str_split(gse132465.tumor,':',simplify=T))[,2]) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

# GSE201348 Matched
gse201348.normal <- sort(unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE201348_Normal', 6])[1:15])
gse201348.polyp <- sort(unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE201348_Polyp', 6])[-c(20:23)])

prop.gse201348.normal <- prop.samples[prop.samples$Var2%in%gse201348.normal,] %>% arrange(Var2)
prop.gse201348.polyp <- prop.samples[prop.samples$Var2%in%gse201348.polyp,] %>% arrange(Var2)

ggplot(prop.gse201348.normal, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  scale_x_discrete(labels=as.data.frame(str_split(gse201348.normal,':',simplify=T))[,2]) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

ggplot(prop.gse201348.polyp, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  scale_x_discrete(labels=as.data.frame(str_split(gse201348.polyp,':',simplify=T))[,2]) +
  geom_col(position='fill', width=0.5) +
  xlab('Sample') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

## Per study
prop.studies <- as.data.frame(table(mes.sr$celltype.detailed, mes.sr$study))
ggplot(prop.studies, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  geom_col(position='fill', width=0.5) +
  xlab('Study') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

# Considering conditions
prop.studies.temp <- as.data.frame(table(mes.sr$celltype.detailed, mes.sr$study.condition))
prop.studies.normal <- prop.studies.temp[grep('Normal',prop.studies.temp$Var2),]
prop.studies.polyp <- prop.studies.temp[grep('Polyp',prop.studies.temp$Var2),]
prop.studies.tumor <- prop.studies.temp[grep('Tumor',prop.studies.temp$Var2),]

ggplot(prop.studies.normal, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  geom_col(position='fill', width=0.5) +
  xlab('Study') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

ggplot(prop.studies.polyp, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  geom_col(position='fill', width=0.5) +
  xlab('Study') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

ggplot(prop.studies.tumor, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  geom_col(position='fill', width=0.5) +
  xlab('Study') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))

## Per condition
prop.condition <- as.data.frame(table(mes.sr$celltype.detailed, mes.sr$condition))
ggplot(prop.condition, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=15) +
  geom_col(position='fill', width=0.5) +
  xlab('Condition') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))


sr <- qread('allPublic_all_afterRPCA_afterRemovingDoublet_readyToUse.qs')
prop.studies <- as.data.frame(table(sr$condition, sr$study))
ggplot(prop.studies, aes(x=Var2, y=Freq, fill=Var1)) +
  theme_classic(base_size=12) +
  geom_col(position='fill', width=0.5) +
  xlab('Study') + ylab('Proportion') +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  theme(legend.title=element_blank(), axis.text.x=element_text(angle=-90, vjust=0.5))


for (i in 1:length(unique(prop.gse201348.normal$Var2))) {
  print(sum(prop.gse201348.normal[prop.gse201348.normal$Var2==unique(prop.gse201348.normal$Var2)[i],3]))
}

for (i in 1:length(unique(prop.gse201348.polyp$Var2))) {
  print(sum(prop.gse201348.polyp[prop.gse201348.polyp$Var2==unique(prop.gse201348.polyp$Var2)[i],3]))
}
