setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/gutcellatlas/Processing/')
library(Seurat)
library(SeuratDisk)
library(zellkonverter)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)
library(qs)

### Loading materials
## H5AD to SeuratObject
Convert('gutcellatlas_spacetime_all_raw_counts.h5ad', dest = "h5seurat", overwrite = TRUE)
gca.sr <- LoadH5Seurat('gutcellatlas_spacetime_all_raw_counts.h5seurat')
colnames(gca.sr@meta.data)[c(1,4,5)] <- c('donor','sample.id','region_code')

### Subsetting subjects
gca.sr <- subset(gca.sr, subset=Diagnosis=='Healthy adult')
gca.sr <- subset(gca.sr, subset=Region%in%c('LargeInt','REC'))
gca.sr <- subset(gca.sr, subset=sample.id%in%names(which(table(gca.sr$sample.id)>=100)))


gca.sr$sample.id <- as.factor(as.character(gca.sr$sample.id))


gca.sr[['RNA']] <- split(gca.sr[['RNA']], f=gca.sr[['sample.id']][,'sample.id'])

gca.sr <- NormalizeData(gca.sr)
gca.sr <- FindVariableFeatures(gca.sr, nfeatures=3000)
gca.sr <- ScaleData(gca.sr)
gca.sr <- RunPCA(gca.sr)

ElbowPlot(gca.sr, reduction='pca', ndims=30)

gca.sr <- FindNeighbors(gca.sr, reduction='pca', dims=1:20)
gca.sr <- FindClusters(gca.sr, resolution=0.8)
gca.sr$ident <- gca.sr$seurat_clusters
gca.sr <- RunUMAP(gca.sr, reduction='pca', dims=1:20)

UMAPPlot(gca.sr)
UMAPPlot(gca.sr, group.by='donor')
UMAPPlot(gca.sr, group.by='sample.id') + theme(legend.position='none')

gca.sr <- IntegrateLayers(gca.sr, method=RPCAIntegration, orig.reduction='pca', new.reduction='integrated.rpca')
gca.sr[['RNA']] <- JoinLayers(gca.sr[['RNA']])

gca.sr <- FindNeighbors(gca.sr, reduction='integrated.rpca', dims=1:30)
gca.sr <- FindClusters(gca.sr, resolution=0.8)
gca.sr$ident <- gca.sr$seurat_clusters
gca.sr <- RunUMAP(gca.sr, reduction='integrated.rpca', dims=1:30)

DimPlot(gca.sr, reduction='umap', raster=F, label=T)
DimPlot(gca.sr, reduction='umap', group.by='sample.id') + theme(legend.position='none')


gca.sce <- as.SingleCellExperiment(gca.sr)
gca.sce <- scDblFinder(gca.sce, clusters='seurat_clusters', samples='sample.id', BPPARAM=MulticoreParam(4))
gca.sr$scDblFinder.score <- gca.sce$scDblFinder.score
gca.sr$scDblFinder.class <- gca.sce$scDblFinder.class

FeaturePlot(gca.sr, 'scDblFinder.score', raster=F) + ggtitle(NULL)
UMAPPlot(gca.sr, group.by='scDblFinder.class', raster=F) + ggtitle(NULL)
VlnPlot(gca.sr, c('nCount_RNA', 'nFeature_RNA', 'pct_counts_mt'),
        group.by='scDblFinder.class', pt.size=0, raster=F)

gca.sr <- subset(gca.sr, cells=colnames(gca.sr)[gca.sr$scDblFinder.class=='singlet'])

qsave(gca.sr, 'gutcellatlas_subset_seuratObject_filtered_unprocessed.qs')
