setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
library(Seurat)
library(ggplot2)
library(qs)


### Loading merged seurat objects
gse201348 <- readRDS('gse201348_merged_seuratObject_filtered_unprocessed.rds')
gca <- qread('gutcellatlas_subset_seuratObject_filtered_unprocessed.qs')
oldpub <- qread('oldPublic_merged_seuratObject_filtered_unprocessed.qs')


### Modifying metadata
## GSE201348
colnames(gse201348@meta.data)[6] <- 'sample.id'
gse201348$sample.id <- paste0(gse201348$study,':',gse201348$sample.id)
gse201348[['percent.mt']] <- PercentageFeatureSet(gse201348, pattern='MT-')
gse201348@assays$RNA@misc <- list()

## GutCellAtlas
gca.new <- CreateSeuratObject(counts=CreateAssayObject(gca@assays$RNA$counts), project='GutCellAtlas', meta.data=NULL)
gca.new$condition <- 'Normal'
gca.new$study <- 'GutCellAtlas'
gca.new$sample.id <- paste0(gca.new$study,':',gca$sample.id)
gca.new$percent.mt <- PercentageFeatureSet(gca.new, pattern='MT-')
gca.new@assays$RNA@misc <- list()

## oldPublic
oldpub@assays$RNA@misc <- list()

### Merging
all <- merge(gse201348, list(gca.new, oldpub))

all@assays$RNA@misc <- list()
all@assays$RNA@key <- 'rna_'
unique(all$condition)
unique(all$study)
unique(all$sample.id)

### Saving
qsave(all, 'allPublic_merged_seuratObject_filtered_unprocessed.qs')
