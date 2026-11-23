library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(devEMF)


### Loading a seurat object
setwd('~/plc/single-cell/merged/Processing')
mes <- qread('allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')


### Subsetting and modifying
mes@meta.data[mes@meta.data=='Tumor'] <- 'Cancer'
mes$isNormal <- '.'
mes@meta.data[mes@meta.data$condition=='Normal', 'isNormal'] <- 'Normal'
mes@meta.data[mes@meta.data$condition%in%c('Polyp','Cancer'), 'isNormal'] <- 'Polyp/Cancer'
mes.coi <- subset(mes, subset=celltype.detailed%in%c('Contractile pericyte','Smooth muscle cell','Myofibroblast (RSPO2+)'))

mes.coi$expressOR <- '.'
mes.coi@meta.data[WhichCells(mes.coi, expression=OR51E2>0), 'expressOR'] <- 'Yes'
mes.coi@meta.data[WhichCells(mes.coi, expression=OR51E2==0), 'expressOR'] <- 'No'
mes.coi.meta <- mes.coi@meta.data


### Calculating percentages of OR51E2-expressing cells
normal.all.percYes <- length(which(mes.coi.meta$isNormal=='Normal' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Normal'))*100
tumor.all.percYes <- length(which(mes.coi.meta$isNormal=='Polyp/Cancer' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Polyp/Cancer'))*100

normal.cp.percYes <- length(which(mes.coi.meta$isNormal=='Normal' & mes.coi.meta$celltype.detailed=='Contractile pericyte' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Normal' & mes.coi.meta$celltype.detailed=='Contractile pericyte'))*100
tumor.cp.percYes <- length(which(mes.coi.meta$isNormal=='Polyp/Cancer' & mes.coi.meta$celltype.detailed=='Contractile pericyte' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Polyp/Cancer' & mes.coi.meta$celltype.detailed=='Contractile pericyte'))*100
normal.myof.percYes <- length(which(mes.coi.meta$isNormal=='Normal' & mes.coi.meta$celltype.detailed=='Myofibroblast (RSPO2+)' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Normal' & mes.coi.meta$celltype.detailed=='Myofibroblast (RSPO2+)'))*100
tumor.myof.percYes <- length(which(mes.coi.meta$isNormal=='Polyp/Cancer' & mes.coi.meta$celltype.detailed=='Myofibroblast (RSPO2+)' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Polyp/Cancer' & mes.coi.meta$celltype.detailed=='Myofibroblast (RSPO2+)'))*100
normal.smc.percYes <- length(which(mes.coi.meta$isNormal=='Normal' & mes.coi.meta$celltype.detailed=='Smooth muscle cell' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Normal' & mes.coi.meta$celltype.detailed=='Smooth muscle cell'))*100
tumor.smc.percYes <- length(which(mes.coi.meta$isNormal=='Polyp/Cancer' & mes.coi.meta$celltype.detailed=='Smooth muscle cell' & mes.coi.meta$expressOR=='Yes'))/length(which(mes.coi.meta$isNormal=='Polyp/Cancer' & mes.coi.meta$celltype.detailed=='Smooth muscle cell'))*100

coi.percYes.df <- data.frame(celltype=c('CP','CP','MFB (RSPO2+)','MFB (RSPO2+)','SMC','SMC'),
                             status=rep(c('Normal','Polyp/Cancer'),3),
                             percentage=c(normal.cp.percYes, tumor.cp.percYes, normal.myof.percYes,
                                          tumor.myof.percYes, normal.smc.percYes, tumor.smc.percYes))

### Plotting
emf('~/plc/Figures/Fig3/fig3f_coi_or51e2_expressing_cell_percentage.emf', width=120, height=80, units='mm', coordDPI=500, emfPlus=FALSE)
ggplot(coi.percYes.df) +
  geom_bar(aes(x=celltype, y=percentage, fill=status), stat='identity', position='dodge') +
  coord_cartesian(ylim=c(0,7.75)) +
  theme_classic(base_size=14) +
  labs(x=NULL)
dev.off()
