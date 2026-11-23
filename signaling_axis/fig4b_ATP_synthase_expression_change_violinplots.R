library(qs)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devEMF)


### Loading
setwd('~/plc/single-cell/merged/Processing')
mes <- qread('allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')


### Changes in OR51E2 expression
## Processing seurat objects
mes@meta.data[mes@meta.data=='Tumor'] <- 'Cancer'
mes$isNormal <- '.'
mes@meta.data[mes@meta.data$condition=='Normal', 'isNormal'] <- 'Normal'
mes@meta.data[mes@meta.data$condition%in%c('Polyp','Cancer'), 'isNormal'] <- 'Polyp/Cancer'
Idents(mes) <- 'isNormal'

mes.or <- subset(mes, subset=celltype.detailed%in%c('Contractile pericyte','Smooth muscle cell','Myofibroblast (RSPO2+)'))
Idents(mes.or) <- 'isNormal'


### ATP synthase genes
## Loading 14 ATP synthase genes
setwd('~/plc/single-cell/network')
or.module <- read.table('or.module_red_posCor.txt', header=F)[,'V1']
atp.genes <- sort(or.module[grep('^ATP5', or.module)])
atp.genes <- setdiff(atp.genes, 'ATP5IF1')

## Calculating module scores
mes <- AddModuleScore(mes, features=list(atp.genes), name='ATP_synthase',ctrl = 5)
mes.or <- AddModuleScore(mes.or, features=list(atp.genes), name='ATP_synthase',ctrl = 5)

FeaturePlot(mes, reduction='umap.rpca', features='ATP_synthase1', raster=F, order=T)
FeaturePlot(mes.or, reduction='umap.rpca', features='ATP_synthase1', raster=F, order=T)

FeaturePlot(mes, reduction='umap.rpca', features='ATP_synthase1', split.by='isNormal', raster=F, order=T)
FeaturePlot(mes.or, reduction='umap.rpca', features='ATP_synthase1', split.by='isNormal', raster=F, order=T)

## Violin plots for ATP synthases
# Across all cell types
VlnPlot(mes, features='ATP_synthase1', group.by='isNormal', pt.size=0) + geom_boxplot(width=0.1, fill="white")
VlnPlot(mes, features='ATP_synthase1', group.by='celltype.detailed', split.by='isNormal', pt.size=0)

# Across COIs
atp.vln.coi <- VlnPlot(mes.or, features='ATP_synthase1', group.by='isNormal', pt.size=0) + geom_boxplot(width=0.1, fill="white") + ggtitle('COIs') + xlab(NULL)
atp.vln.cp <- VlnPlot(subset(mes.or, subset=celltype.detailed=='Contractile pericyte'), features='ATP_synthase1', group.by='isNormal', pt.size=0) + geom_boxplot(width=0.1, fill="white") + ggtitle('CP') + xlab(NULL)
atp.vln.smc <- VlnPlot(subset(mes.or, subset=celltype.detailed=='Smooth muscle cell'), features='ATP_synthase1', group.by='isNormal', pt.size=0) + geom_boxplot(width=0.1, fill="white") + ggtitle('SMC') + xlab(NULL)
atp.vln.myof <- VlnPlot(subset(mes.or, subset=celltype.detailed=='Myofibroblast (RSPO2+)'), features='ATP_synthase1', group.by='isNormal', pt.size=0) + geom_boxplot(width=0.1, fill="white") + ggtitle('MF (RSPO2+)') + xlab(NULL)

normal.coi.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Normal','ATP_synthase1']
tumor.coi.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Polyp/Cancer','ATP_synthase1']
wilcox.test(x=tumor.coi.atp, y=normal.coi.atp, alternative='two.sided')   # Wilcoxon rank sum test (= Mann-Whitney test)

normal.cp.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Normal' & mes.or@meta.data$celltype.detailed=='Contractile pericyte', 'ATP_synthase1']
normal.myof.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Normal' & mes.or@meta.data$celltype.detailed=='Myofibroblast (RSPO2+)', 'ATP_synthase1']
normal.smc.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Normal' & mes.or@meta.data$celltype.detailed=='Smooth muscle cell', 'ATP_synthase1']
tumor.cp.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Polyp/Cancer' & mes.or@meta.data$celltype.detailed=='Contractile pericyte', 'ATP_synthase1']
tumor.myof.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Polyp/Cancer' & mes.or@meta.data$celltype.detailed=='Myofibroblast (RSPO2+)', 'ATP_synthase1']
tumor.smc.atp <- mes.or@meta.data[mes.or@meta.data$isNormal=='Polyp/Cancer' & mes.or@meta.data$celltype.detailed=='Smooth muscle cell', 'ATP_synthase1']
wilcox.test(x=tumor.cp.atp, y=normal.cp.atp, alternative='two.sided')   # Wilcoxon rank sum test (= Mann-Whitney test)
wilcox.test(x=tumor.myof.atp, y=normal.myof.atp, alternative='two.sided')   # Wilcoxon rank sum test (= Mann-Whitney test)
wilcox.test(x=tumor.smc.atp, y=normal.smc.atp, alternative='two.sided')   # Wilcoxon rank sum test (= Mann-Whitney test)


### Saving plots
emf('~/plc/Figures/Fig4/fig4b_ATP_synthase_expression_change_violinplots.emf', width=350, height=90, units='mm', coordDPI=500, emfPlus=FALSE)
atp.vln.coi + atp.vln.cp + atp.vln.myof + atp.vln.smc + plot_layout(ncol=4, guides="collect")
dev.off()
