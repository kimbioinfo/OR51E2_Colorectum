setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/gse201348/Processing')
library(Seurat)
library(SeuratObject)
library(qs)

### Loading previously processed seurat objects
gse201348.n.depri <- qread('seurat_n_simply_annotated.qs')
gse201348.p1.depri <- qread('seurat_p1_simply_annotated.qs')
gse201348.p2.depri <- qread('seurat_p2_simply_annotated.qs')
gse201348.c.depri <- qread('seurat_c_simply_processed.qs')


### Creating new seurat objects without layers
gse201348.n.sr <- CreateSeuratObject(counts=CreateAssayObject(counts=gse201348.n.depri@assays$originalexp$counts), project='GSE201348', meta.data=NULL)
gse201348.p1.sr <- CreateSeuratObject(counts=CreateAssayObject(counts=gse201348.p1.depri@assays$originalexp$counts), project='GSE201348', meta.data=NULL)
gse201348.p2.sr <- CreateSeuratObject(counts=CreateAssayObject(counts=gse201348.p2.depri@assays$originalexp$counts), project='GSE201348', meta.data=NULL)
gse201348.c.sr <- CreateSeuratObject(counts=CreateAssayObject(counts=gse201348.c.depri@assays$RNA$counts), project='GSE201348', meta.data=NULL)


### Adding metadata
gse201348.n.sr@assays$RNA@misc <- list()
gse201348.n.sr$condition <- 'Normal'
gse201348.n.sr$study <- 'GSE201348'
gse201348.n.sr$sample.name <- gse201348.n.depri$id

gse201348.p1.sr@assays$RNA@misc <- list()
gse201348.p1.sr$condition <- 'Polyp'
gse201348.p1.sr$study <- 'GSE201348'
gse201348.p1.sr$sample.name <- gse201348.p1.depri$id

gse201348.p2.sr@assays$RNA@misc <- list()
gse201348.p2.sr$condition <- 'Polyp'
gse201348.p2.sr$study <- 'GSE201348'
gse201348.p2.sr$sample.name <- gse201348.p2.depri$id

gse201348.c.sr@assays$RNA@misc <- list()
gse201348.c.sr$condition <- 'Tumor'
gse201348.c.sr$study <- 'GSE201348'
gse201348.c.sr$sample.name <- gse201348.c.depri$sample.id

rm(gse201348.n.depri, gse201348.p1.depri, gse201348.p2.depri, gse201348.c.depri)
gc()


### Merging
gse201348 <- merge(gse201348.n.sr, list(gse201348.p1.sr, gse201348.p2.sr, gse201348.c.sr))

rm(gse201348.n.sr, gse201348.p1.sr, gse201348.p2.sr, gse201348.c.sr)
gc()


### Saving
qsave(gse201348, 'gse201348_merged_seuratObject_filtered_unprocessed.qs')
