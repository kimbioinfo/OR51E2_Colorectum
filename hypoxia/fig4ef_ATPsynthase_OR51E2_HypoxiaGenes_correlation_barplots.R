setwd('~/plc/single-cell/merged/Processing')
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(qs)
library(glue)
library(Seurat)
library(devEMF)
library(edgeR)

### Loading seurat object
## Seurat object with all kinds of cells
all.sr <- qread('allPublic_all_afterRPCA_afterRemovingDoublet_simply_annotated.qs')
mes.sr <- qread('allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')
table(all.sr$celltype.simple)
table(all.sr$sample.id)

## Transfer of mesenchymal detail cell types
all(rownames(mes.sr@meta.data)%in%rownames(all.sr@meta.data))   # TRUE
all.sr$celltype.detailed <- all.sr$celltype.simple
all.sr$celltype.detailed <- as.character(all.sr$celltype.detailed)
mes.sr$celltype.detailed <- as.character(mes.sr$celltype.detailed)
all.sr@meta.data[colnames(mes.sr), 'celltype.detailed'] <- mes.sr@meta.data$celltype.detailed
all.sr$celltype.coi <- all.sr$celltype.detailed
all.sr@meta.data[all.sr@meta.data$celltype.detailed%in%c('Smooth muscle cell', 'Contractile pericyte', 'Myofibroblast (RSPO2+)'),'celltype.coi'] <- 'COI'
all.sr@meta.data[all.sr@meta.data$celltype.simple=='Mesenchymal' & !(all.sr@meta.data$celltype.detailed%in%c('Smooth muscle cell', 'Contractile pericyte', 'Myofibroblast (RSPO2+)')),'celltype.coi'] <- 'Non-COI'

mes.sr$celltype.detailed <- as.factor(mes.sr$celltype.detailed)
all.sr$celltype.detailed <- as.factor(all.sr$celltype.detailed)
all.sr$celltype.coi <- as.factor(all.sr$celltype.coi)
table(all.sr$celltype.detailed)
table(all.sr$celltype.coi)

## Subset with OR51E2 expression in mesenchymal cells
coi.sr <- subset(all.sr, subset=celltype.coi=='COI')
coi.sample <- unique(coi.sr$sample.id)
all.coi <- subset(all.sr, subset=sample.id%in%coi.sample)
all.coi <- subset(all.coi, subset=condition!='Normal')

# Manually setting input tissue types (Healthy normal tisseus OR Adjacent normal tissues)
#all.coi <- subset(all.coi, subset=condition=='Normal')
#healthy <- c(unique(all.coi$sample.id[grep('^GSE201348:B',all.coi$sample.id)]), unique(all.coi$sample.id[grep('^GutCellAtlas',all.coi$sample.id)]))
#adjacent <- c(paste0('GSE132465:SMC0',seq(0,9),'-N'),'GSE201348:A001-C-124','GSE201348:A002-C-024','GSE201348:A002-C-121',
#              'GSE201348:A002-C-121-R0','GSE201348:A002-C-212','GSE201348:A014-C-054','GSE201348:A014-C-111','GSE201348:A014-C-114',
#              'GSE201348:A015-C-010')
#all.coi <- subset(all.coi, subset=sample.id%in%healthy)
#all.coi <- subset(all.coi, subset=sample.id%in%adjacent)

coi.or <- subset(coi.sr, subset=OR51E2>0)
or.sample <- unique(coi.or$sample.id)
all.or <- subset(all.sr, subset=sample.id%in%or.sample)
all.or <- subset(all.or, subset=condition!='Normal')

#all.or <- subset(all.or, subset=condition=='Normal')
#all.or <- subset(all.or, subset=sample.id%in%healthy)
#all.or <- subset(all.or, subset=sample.id%in%adjacent)

rm(all.sr, mes.sr)
gc()


### Pseudobulk samples
## Gene lists
setwd('~/plc/single-cell/network')
or.module <- read.table('or.module_red.txt', header=T)[,'gene_name']

atp.genes <- sort(or.module[grep('^ATP5', or.module)])
atp.genes <- setdiff(atp.genes, 'ATP5IF1')

hif.Immune.genes <- c('ITGB2','IL22','CXCL2','CCL20','TREM1','TLR4','IL10','CCR9','ALDH1A2','CXCL8','TNF','FOXP3','ALDOA','PGK1','PKM','IFNG','CXCR2','CXCR3','CCR5','KCNK5')

## Pseudobulk object
pseudo.or <- AggregateExpression(all.or, assays='RNA', return.seurat=F, group.by=c('sample.id','celltype.coi')) %>% as.data.frame()
colnames(pseudo.or) <- gsub('RNA.','',colnames(pseudo.or))
pseudo.coi <- AggregateExpression(all.coi, assays='RNA', return.seurat=F, group.by=c('sample.id','celltype.coi')) %>% as.data.frame()
colnames(pseudo.coi) <- gsub('RNA.','',colnames(pseudo.coi))

pseudo.or <- cpm(pseudo.or) %>% as.data.frame()
pseudo.coi <- cpm(pseudo.coi) %>% as.data.frame()

rm(all.or, all.coi, coi.or, coi.sr)
gc()

# Manual iterative execution with different 'ref' variables
ref <- 'or51e2'
ref <- 'atp'
pseudo.sr <- get(ifelse(ref=='or51e2', 'pseudo.or', 'pseudo.coi'))

pseudo.goi <- pseudo.sr[c('OR51E2',hif.Immune.genes),]
pseudo.atp <- apply(pseudo.sr[c(atp.genes),],2,mean) %>% t() %>% as.data.frame()
rownames(pseudo.atp) <- 'ATP_synthase'
pseudo.goi <- rbind(pseudo.goi, pseudo.atp)

## Processing dataframe
sample.ct <- as.data.frame(str_split(colnames(pseudo.goi), '_', simplify=T))
sample.ct <- data.frame(sample=paste0(sample.ct$V1,'_',sample.ct$V2), celltype=sample.ct$V3)

eachGene <- function(gene, ct) {
  samples <- sample.ct[sample.ct$celltype==ct,'sample']
  if (ref=='or51e2') {
    temp <- data.frame(or51e2=as.numeric(pseudo.goi['OR51E2',paste0(samples,'_COI')]),
                       goi=as.numeric(pseudo.goi[gene,paste0(samples,'_',ct)]), row.names=samples)
  }
  if (ref=='atp') {
    temp <- data.frame(atp=as.numeric(pseudo.goi['ATP_synthase',paste0(samples,'_COI')]),
                       goi=as.numeric(pseudo.goi[gene,paste0(samples,'_',ct)]), row.names=samples)
  }
  assign(glue('{gene}.{ct}'), temp, globalenv())
}

eachGene(hif.Immune.genes, 'Immune')

## Testing correlation
getCoef <- function(gene, ct, method) {
  temp <- get(glue('{gene}.{ct}'))
  if (ref=='or51e2') {
    result <- cor.test(temp$or51e2, temp$goi, method=method)
  }
  if (ref=='atp') {
    result <- cor.test(temp$atp, temp$goi, method=method)
  }
  return(result$estimate)
}

getPval <- function(gene, ct, method) {
  temp <- get(glue('{gene}.{ct}'))
  if (ref=='or51e2') {
    result <- cor.test(temp$or51e2, temp$goi, method=method)
  }
  if (ref=='atp') {
    result <- cor.test(temp$atp, temp$goi, method=method)
  }
  return(result$p.value)
}

getCor <- function(ct, method) {
  genes <- get(glue('hif.{ct}.genes'))
  temp <- data.frame(coef=rep('.',length(genes)), P=rep('.',length(genes)), row.names=genes)
  for (gene in genes) {
    temp[gene,'coef'] <- getCoef(gene, ct, method)
    temp[gene,'P'] <- getPval(gene, ct, method)
  }
  temp$coef <- as.numeric(temp$coef)
  temp$P <- as.numeric(temp$P)
  assign(glue('cor.{ct}'), temp, globalenv())
}

getCor('Immune', 'spearman')

# Visualization
isSignif <- function(ct) {
  cor <- get(glue('cor.{ct}'))
  cut <- 0.05
  cor$isSignif <- ifelse(cor$P<cut, glue('P < {cut}'), glue('P >= {cut}'))
  cor$isSignif <- factor(cor$isSignif, levels=c(glue('P < {cut}'), glue('P >= {cut}')))
  cor <- cor %>% rownames_to_column('gene')
  assign(glue('cor.{ct}'), cor, globalenv())
}

getPlot <- function(ct) {
  cor <- get(glue('cor.{ct}'))
  if (length(which(is.na(cor)))!=0) {
    cor <- cor[-which(is.na(cor$coef)),]
  }
  cut <- 0.05
  if (all(unique(as.character(cor$isSignif))==glue('P >= {cut}'))) {
    plot <- ggplot(cor, aes(x=coef, y=reorder(gene, coef), fill=isSignif)) +
      theme_bw() +
      geom_bar(stat='identity') +
      scale_fill_manual(values=c('grey40')) +
      theme(axis.title.y=element_blank()) +
      xlab('Correlation coefficient')
  } else {
    plot <- ggplot(cor, aes(x=coef, y=reorder(gene, coef), fill=isSignif)) +
      theme_bw() +
      geom_bar(stat='identity') +
      scale_fill_manual(values=c('#fb736d','grey40')) +
      theme(axis.title.y=element_blank()) +
      xlab('Correlation coefficient')
  }
  return(plot)
}

for (ct in celltypes) {
  isSignif(ct)
}


library(devEMF)
## When ref=='atp'
cordf.atp.immune <- cor.Immune
corplot.atp.immune <- ggplot(cordf.atp.immune, aes(x=coef, y=reorder(gene, coef), fill=isSignif)) +
  theme_bw() +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('#fb736d','grey40')) +
  theme(axis.title.y=element_blank()) +
  xlab('Correlation coefficient') +
  theme(text=element_text(size=9))

qsave(cordf.atp.immune, 'spearman_correlation_ATPsynthase_hypoxiaGenes_immune_tumor.qs')   # Tumor
#qsave(cordf.atp.immune, 'spearman_correlation_ATPsynthase_hypoxiaGenes_immune_adjacentNormal.qs')   # Adjacent normal
#qsave(cordf.atp.immune, 'spearman_correlation_ATPsynthase_hypoxiaGenes_immune_healthyNormal.qs')   # Healthy normal

## When ref=='or51e2'
# Immune cells
cordf.or.immune <- cor.Immune
corplot.or.immune <- ggplot(cordf.or.immune, aes(x=coef, y=reorder(gene, cordf.atp.immune$coef), fill=isSignif)) +
  theme_bw() +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('#fb736d','grey40')) +
  theme(axis.title.y=element_blank()) +
  xlab('Correlation coefficient') +
  theme(text=element_text(size=9))

qsave(cordf.or.immune, 'spearman_correlation_OR51E2_hypoxiaGenes_immune_tumor.qs')   # Tumor
#qsave(cordf.or.immune, 'spearman_correlation_OR51E2_hypoxiaGenes_immune_adjacentNormal.qs')   # Adjacent normal
#qsave(cordf.or.immune, 'spearman_correlation_OR51E2_hypoxiaGenes_immune_healthyNormal.qs')   # Healthy normal


##### --------------------------------------------------------------------

### Merged plot
## Merging correlation results
atp.tumor <- qread('spearman_correlation_ATPsynthase_hypoxiaGenes_immune_tumor.qs') %>% mutate(tissueType='Tumor') %>% arrange(coef, gene)
atp.adjacent <- qread('spearman_correlation_ATPsynthase_hypoxiaGenes_immune_adjacentNormal.qs') %>% mutate(tissueType='AdjacentNormal')
atp.healthy <- qread('spearman_correlation_ATPsynthase_hypoxiaGenes_immune_healthyNormal.qs') %>% mutate(tissueType='HealthyNormal')
atp.healthy[atp.healthy=='CCL2'] <- 'CCL20'

or.tumor <- qread('spearman_correlation_OR51E2_hypoxiaGenes_immune_tumor.qs') %>% mutate(tissueType='Tumor') %>% arrange(coef, gene)
or.adjacent <- qread('spearman_correlation_OR51E2_hypoxiaGenes_immune_adjacentNormal.qs') %>% mutate(tissueType='AdjacentNormal')
or.healthy <- qread('spearman_correlation_OR51E2_hypoxiaGenes_immune_healthyNormal.qs') %>% mutate(tissueType='HealthyNormal')
or.healthy$isSignif <- ifelse(or.healthy$P<0.05, 'P < 0.05', 'P >= 0.05')
or.healthy[or.healthy=='CCL2'] <- 'CCL20'

atp.merged <- rbind(atp.tumor, atp.adjacent, atp.healthy)
or.merged <- rbind(or.tumor, or.adjacent, or.healthy)

atp.merged$gene <- factor(atp.merged$gene, levels=atp.tumor$gene)
or.merged$gene <- factor(or.merged$gene, levels=atp.tumor$gene)

atp.merged$isSignif <- factor(atp.merged$isSignif, levels=c('P < 0.05','P >= 0.05'))
or.merged$isSignif <- factor(or.merged$isSignif, levels=c('P < 0.05','P >= 0.05'))

atp.merged$tissueType <- factor(atp.merged$tissueType, levels=c('Tumor','AdjacentNormal','HealthyNormal'))
or.merged$tissueType <- factor(or.merged$tissueType, levels=c('Tumor','AdjacentNormal','HealthyNormal'))

## Heatmap
# ATP synthase genes
heat.atp.merged <- ggplot(atp.merged) +
  geom_tile(aes(x=gene, y=rev(tissueType), fill=coef)) +
  geom_point(aes(x=gene, y=rev(tissueType), alpha=isSignif), shape=42, size=4, color='white') +
  scale_alpha_manual(values=c(1,0)) +
  scale_fill_gradientn(limits=c(-1,1), breaks=c(-1,-0.5,0,0.5,1), colors=c('#012d5e','#4F96E9','white','#EC736E','#9a130f')) +
  theme_void() +
  theme(axis.text.x=element_text(angle=-90, v=0.5, hjust=0))

heat.or.merged <- ggplot(or.merged) +
  geom_tile(aes(x=gene, y=rev(tissueType), fill=coef)) +
  geom_point(aes(x=gene, y=rev(tissueType), alpha=isSignif), shape=42, size=4, color='white') +
  scale_alpha_manual(values=c(1,0)) +
  scale_fill_gradientn(limits=c(-1,1), breaks=c(-1,-0.5,0,0.5,1), colors=c('#012d5e','#4F96E9','white','#EC736E','#9a130f')) +
  theme_void() +
  theme(axis.text.x=element_text(angle=-90, v=0.5, hjust=0))

## Saving
emf('~/plc/Figures/Fig5/fig5e_ATPsynthase_Hypoxia_genes_Immune_heatmap.emf', width=240, height=50, units='mm', coordDPI=500, emfPlus=FALSE)
heat.atp.merged
dev.off()

emf('~/plc/Figures/Fig5/fig5e_OR51E2_Hypoxia_genes_Immune_heatmap.emf', width=240, height=50, units='mm', coordDPI=500, emfPlus=FALSE, )
heat.or.merged
dev.off()

qsave(heat.atp.merged, 'heatmap_ATPsynthase_hypoxiaGenes_immune.qs')
qsave(heat.or.merged, 'heatmap_OR51E2_hypoxiaGenes_immune.qs')
