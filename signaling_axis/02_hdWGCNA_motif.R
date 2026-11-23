setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/network/')
library(Seurat)
library(Signac)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(stringr)
library(tibble)
library(qs)
library(ggplot2)
library(glue)
library(patchwork)
library(TFBSTools)
library(dbplyr)
library(JASPAR2020)
library(EnsDb.Hsapiens.v86)
library(GeneOverlap)


### Loading materials
## hdWGCNAed Seurat object
mes.or <- qread('allPublic_or51e2-positive_hdWGCNAed_seuratObject_red.qs')

## TF motif dataset based on ChIP-seq evidence
pfm_core <- getMatrixSet(x=JASPAR2020,
                         opts = list(collection='CORE', type='ChIP-seq', tax_group='vertebrates', species=9606, all_versions=F))
pfm_core.list <- as.list(pfm_core)
tf.list <- purrr::map(1:length(pfm_core.list), function(i){pfm_core.list[[i]]@name}) %>% unlist()
dups <- which(duplicated(tf.list))   # No dups


### TF motif enrichment
## Putting TF motif information into the seurat object
mes.or <- MotifScan(mes.or,
                    species_genome='hg38',
                    pfm=pfm_core,
                    EnsDb=EnsDb.Hsapiens.v86)
dim(GetMotifMatrix(mes.or))
motif.matrix <- GetMotifMatrix(mes.or) %>% as.data.frame()

## Updating TF target gene names
# Updating whole gene set
target_genes <- GetMotifTargets(mes.or)
motif.names <- mes.or@misc$motifs$motif_df[,1:2] %>% unique()
rownames(motif.names) <- NULL
motif.names <- motif.names %>% column_to_rownames('motif_ID')
names(target_genes) <- motif.names[names(target_genes),'motif_name']
names(target_genes)[2] <- 'EWSR1-FLI1'

target.uniq <- unique(purrr::map(1:length(target_genes), function(i){target_genes[[i]]}) %>% unlist())
target.uniq.update <- UpdateSymbolList(target.uniq)

# Finding and Fixing trouble points
isSame <- target.uniq==target.uniq.update
isSame[4301:4500]   # 'False' starting from 4402
target.uniq[4402]   # LTBP2
target.uniq[4403]   # FCF1
target.uniq.update[4402]   # FCF1
'LTBP2'%in%rownames(mes.or)   # TRUE
target.uniq.update <- c(target.uniq.update[1:4401],'LTBP2',target.uniq.update[4402:length(target.uniq.update)])

isSame <- target.uniq==target.uniq.update
isSame[12801:12900]
target.uniq[12837]   # GABPB1
target.uniq[12838]  # USP50
target.uniq.update[12837]   # USP50
'GABPB1'%in%rownames(mes.or)   # TRUE
target.uniq.update <- c(target.uniq.update[1:12836],'GABPB1',target.uniq.update[12837:length(target.uniq.update)])

isSame <- target.uniq==target.uniq.update
isSame[13901:14000]
target.uniq[13926]   # SLC23A2
target.uniq[13927]   # BMP2
target.uniq.update[13926]   # BMP2
'SLC23A2'%in%rownames(mes.or)
target.uniq.update <- c(target.uniq.update[1:13925],'SLC23A2',target.uniq.update[13926:length(target.uniq.update)])

isSame <- target.uniq==target.uniq.update
isSame[19501:19604]   # All OK

# Updating TF motif object
target_genes.update <- target_genes
for (i in 1:length(target_genes.update)) {
  for (j in 1:length(target_genes.update[[i]])) {
    temp <- target_genes.update[[i]][j]
    target_genes.update[[i]][j] <- target.uniq.update[which(target.uniq==temp)]
  }
  print(glue('{i} Done'))
}

mes.or@misc$motifs$motif_targets <- target_genes.update

mes.or <- OverlapModulesMotifs(mes.or)

motif.allModules <- as.data.frame(mes.or@misc$wgcna$motif_module_overlaps)
motif.allModules <- motif.allModules[motif.allModules$Significance=='***',]

motif.red <- motif.allModules[motif.allModules$module=='red',] %>% arrange(fdr)
rownames(motif.red) <- NULL

yy1.target <- mes.or@misc[["motifs"]][["motif_targets"]]$YY1
red.genes <- mes.or@misc$wgcna$wgcna_modules[mes.or@misc$wgcna$wgcna_modules$module=='red', 'gene_name']
yy1.target.red <- intersect(yy1.target, red.genes)

sort(yy1.target)
sort(yy1.target.red)   # 9 of 14 ATP synthase genes

### Recording target genes for each TF
motif.red$target_genes <- '.'
for (tf in motif.red$tf) {
  temp.target <- mes.or@misc[["motifs"]][["motif_targets"]][[tf]]
  temp.target.red <- sort(intersect(temp.target, red.genes))
  motif.red[motif.red$tf==tf,'target_genes'] <- paste(temp.target.red, collapse=' | ')
}

qsave(mes.or, 'allPublic_or51e2-positive_hdWGCNAed_Motifed_seuratObject_red.qs')
write.table(motif.red, 'motifs_for_or51e2_module_red.txt', col.names=T, row.names=F, sep='\t', quote=F)

