setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing')
library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(stringr)
library(tibble)
library(qs)
library(ggplot2)
library(patchwork)
library(glue)

set.seed(214)
enableWGCNAThreads(nThreads=12)

### Loading a seurat object
mes <- qread('allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')
mes$celltype.simple <- as.factor(mes$celltype.simple)
mes$celltype.detailed <- as.factor(mes$celltype.detailed)
mes$sample.id <- as.factor(mes$sample.id)
mes$study.condition <- as.factor(mes$study.condition)
mes$sampleID.celltype <- paste0(mes$sample.id,'_',mes$celltype.detailed)

## Modifying celltypes
mes$celltype.or51e2 <- '.'
mes$celltype.or51e2[which(mes$celltype.detailed%in%c('Contractile pericyte','Smooth muscle cell','Myofibroblast (RSPO2+)'))] <- 'Positive'
mes$celltype.or51e2[which(!(mes$celltype.detailed%in%c('Contractile pericyte','Smooth muscle cell','Myofibroblast (RSPO2+)')))] <- 'Negative'

## Subsetting
mes <- subset(mes, subset=celltype.or51e2=='Positive')
sort(table(mes$sample.id))

mes.or <- subset(mes, cells=WhichCells(mes, expression=OR51E2>0))
or.sample <- as.data.frame(table(mes.or$sample.id, mes.or$celltype.detailed)) %>% arrange(desc(Freq))
wgcna.sample <- paste0(or.sample[which(or.sample$Freq>0),1],'_',or.sample[which(or.sample$Freq>0),2])

mes.wgcna <- subset(mes, subset=sampleID.celltype%in%wgcna.sample)   # 4394 cells
mes.wgcna <- subset(mes.wgcna, subset=sampleID.celltype%in%names(which(table(as.character(mes.wgcna$sampleID.celltype))>5)))   # 4376 cells

DimPlot(mes.wgcna, reduction='umap.rpca', group.by='celltype.detailed')
DimPlot(mes.wgcna, reduction='umap.rpca', group.by='condition')
DimPlot(mes.wgcna, reduction='umap.rpca', group.by='study')
FeaturePlot(mes.wgcna, features='OR51E2', reduction='umap.rpca', raster=F, order=T)
FeaturePlot(mes.wgcna, features='OR51E2', reduction='umap.rpca', raster=F, order=T)
VlnPlot(mes.wgcna, features='OR51E2', group.by='study')
VlnPlot(mes.wgcna, features='OR51E2', group.by='condition')


#### WGCNA with PSEUDOBULK ####
## Processing
mes.wgcna <- SetupForWGCNA(mes.wgcna,
                        gene_select='fraction',
                        fraction=0.03,
                        wgcna_name='wgcna')

pseudo <- ConstructPseudobulk(mes.wgcna,
                              group.by='celltype.or51e2',
                              replicate_col='sample.id',
                              assay='RNA',
                              slot='counts',
                              wgcna_name='wgcna',
                              min_reps=10,
                              label_col='celltype.detailed')

cpm <- t(apply(pseudo, 1, function(x){
  y <- (x) / sum(x) * 1000000
  log2(y + 1)
}))

mes.pseudo <- SetDatExpr(mes.wgcna, mat=cpm)

mes.pseudo <- TestSoftPowers(mes.pseudo, networkType='unsigned')
plot_list <- PlotSoftPowers(mes.pseudo)
wrap_plots(plot_list, ncol=2)   # sftP=5

mes.pseudo.unsigned <- ConstructNetwork(mes.pseudo,
                               tom_name='pseudobulk',
                               overwrite_tom=T,
                               soft_power=5,
                               minModuleSize=30,
                               mergeCutHeight=0.4,
                               networkType='unsigned',
                               randomSeed=214)

jpeg('~/plc/Figures/Suppl_Fig5/suppl_fig5a_hdWGCNA_dendrogram.png', width=180, height=120, units='mm', res=500)
PlotDendrogram(mes.pseudo.unsigned)
dev.off()

mod.pseudo.unsigned <- GetModules(mes.pseudo.unsigned)
sort(table(mod.pseudo.unsigned$module), decreasing=T)
mod.pseudo.unsigned[c('OR51E2','YY1','ATP5F1A', 'ATP5F1B', 'ATP5MC3', 'ATP5MD', 'ATP5PD', 'ATP5PF'),2]
write.table(mod.pseudo.unsigned[,c(1,2)], '~/plc/Figures/Suppl_Table/WGCNA_result_all_modules_geneList.txt', col.names=F, row.names=F, sep='\t', quote=F)

or.module <- mod.pseudo.unsigned[mod.pseudo.unsigned$module=='red',]


### Saving results
setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/network/')
write.table(or.module, 'or.module_red.txt', col.names=T, row.names=F, sep='\t', quote=F)
qsave(or.module, 'or51e2_module_red.qs')
qsave(mes.pseudo.unsigned, 'allPublic_or51e2-positive_hdWGCNAed_seuratObject_red.qs')

exp.pseudo <- mes.pseudo.unsigned@misc[["wgcna"]][["datExpr"]] %>% as.data.frame()
or.posCor <- NULL
or.negCor <- NULL
for (gene in setdiff(or.module$gene_name, 'OR51E2')) {
  cor.temp <- cor(exp.pseudo[,'OR51E2'], exp.pseudo[,gene])
  if (cor.temp > 0) {or.posCor <- append(or.posCor, gene)}
  if (cor.temp < 0) {or.negCor <- append(or.negCor, gene)}
}   # All genes correlate positively with OR51E2

writeLines(or.posCor, 'or.module_red_posCor.txt')
writeLines(or.negCor, 'or.module_red_negCor.txt')

## Other modules for validation
other.modules <- mod.pseudo.unsigned[mod.pseudo.unsigned$module!='red',1:2]
rownames(other.modules) <- NULL

write.module <- function(module) {
  temp <- data.frame(gene=unique(other.modules[other.modules$module==module,'gene_name']))
  write.table(temp, glue('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/network/modules_fromRed/{module}_genes.txt'),
              col.names=F, row.names=F, sep='\t', quote=F)
}

for (module in names(sort(table(mod.pseudo.unsigned$module), decreasing=T))) {
  write.module(module)
}


