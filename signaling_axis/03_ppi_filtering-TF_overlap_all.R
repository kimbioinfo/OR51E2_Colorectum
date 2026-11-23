setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/ppi')
library(Seurat)
library(qs)
library(dplyr)
library(stringr)
library(tibble)


### Loading an seurat object
mes.or <- qread('../network/allPublic_or51e2-positive_hdWGCNAed_Motifed_seuratObject_red.qs')
DimPlot(mes.or, reduction='umap.rpca', group.by='celltype.detailed')
#ggsave('Output/OR51E2-positive_umap.png', width=950, height=600, unit='px', bg='transparent', scale=2)

### TFs for each cell type
## Loading motif analysis results
motif.red <- read.table('../network/motifs_for_or51e2_module_red.txt', header=T, sep='\t') %>% arrange(fdr)
motif.red <- motif.red[motif.red$Significance=='***',]

## Refining TF names
df.temp1 <- as.data.frame(str_split(motif.red$tf, pattern='\\(', simplify=T))
tf.temp1 <- df.temp1$V1
df.temp2 <- as.data.frame(str_split(tf.temp1, pattern='::', simplify=T))
tf.temp2 <- c(df.temp2$V1, df.temp2$V2, df.temp2$V3)
unique(tf.temp2)
tf.temp2 <- tf.temp2[nchar(tf.temp2)>0]
tf.red <- unique(toupper(tf.temp2))
tf.red <- intersect(tf.red, rownames(mes.or@assays$RNA$counts))   # 144 of 145 TFs

rm(df.temp1, df.temp2, tf.temp1, tf.temp2)

## Cell proportion
prop.tf <- DotPlot(mes.or, features=tf.red, group.by='celltype.detailed')$data[,c(3,4,2)]
prop.tf.cp <- prop.tf[prop.tf$id=='Contractile pericyte',]
prop.tf.myo <- prop.tf[prop.tf$id=='Myofibroblast (RSPO2+)',]
prop.tf.smc <- prop.tf[prop.tf$id=='Smooth muscle cell',]

## Specifying TFs for cell types
tf.cp <- as.character(prop.tf.cp[prop.tf.cp$pct.exp>=5,1])
tf.myo <- as.character(prop.tf.myo[prop.tf.myo$pct.exp>=5,1])
tf.smc <- as.character(prop.tf.smc[prop.tf.smc$pct.exp>=5,1])
tf.all <- intersect(intersect(tf.cp, tf.myo), tf.smc)


#-------------------------------------------

### 1st filtering
## Candidates before filtering
before1 <- read.csv('All/ppi_1st_all_beforeFilter.csv', header=T)[,3]
before1 <- before1[before1!='OR51E2']

## Checking expression
# Cell proportion
prop.before1 <- DotPlot(mes.or, features=before1, group.by='celltype.detailed')$data[,c(3,4,2)]
prop.before1.cp <- prop.before1[prop.before1$id=='Contractile pericyte',]
prop.before1.myo <- prop.before1[prop.before1$id=='Myofibroblast (RSPO2+)',]
prop.before1.smc <- prop.before1[prop.before1$id=='Smooth muscle cell',]

## Classfying
# Contractile pericyte
after1.pos.cp <- as.character(prop.before1.cp[prop.before1.cp$pct.exp>=5,1])

# Myofibroblast (RSPO2+)
after1.pos.myo <- as.character(prop.before1.myo[prop.before1.myo$pct.exp>=5,1])

# Smooth muscle cell
after1.pos.smc <- as.character(prop.before1.smc[prop.before1.smc$pct.exp>=5,1])

# For all cell types
after1.pos.all <- intersect(intersect(after1.pos.cp, after1.pos.myo), after1.pos.smc)

## Filtering with coexpression proportion
first <- data.frame(genes=after1.pos.all)
first$coex.prop <- '.'

expr <- FetchData(mes.or, vars='OR51E2')
or.num <- length(WhichCells(mes.or[,expr>0]))
for (i in 1:nrow(first)) {
  expr2 <- FetchData(mes.or, vars=first[i,1])
  co.num <- length(WhichCells(mes.or[,expr>0 & expr2>0]))
  first[i,2] <- co.num/or.num*100
}
first$coex.prop <- as.numeric(first$coex.prop)
after1.pos.all <- first[first$coex.prop>=5,1]

## Results
# Visualization
#FeaturePlot(mes.or, reduction='umap.rpca', features=after1.pos.all, order=T)
#ggsave('Output/OR51E2-positive_featplot_1st_CP-positive.png', width=2200, height=1800, unit='px', bg='transparent', scale=2)

# Saving
first <- data.frame(genes=after1.pos.all)

write.table(first, 'All/ppi_1st_afterFilter.txt', col.names=T, row.names=F, sep='\t', quote=F)


#-------------------------------------------

### 2nd filtering
## Candidates before filtering
before2 <- read.csv('All/ppi_2nd_all_beforeFilter.csv', header=T)[,3]
before2 <- setdiff(before2, c(after1.pos.all, 'OR51E2'))
before2[before2=='MARCHF7'] <- 'MARCH7'


## Checking expression
# Cell proportion
prop.before2 <- DotPlot(mes.or, features=before2, group.by='celltype.detailed')$data[,c(3,4,2)]
prop.before2.cp <- prop.before2[prop.before2$id=='Contractile pericyte',]
prop.before2.myo <- prop.before2[prop.before2$id=='Myofibroblast (RSPO2+)',]
prop.before2.smc <- prop.before2[prop.before2$id=='Smooth muscle cell',]

## Classfying
# Contractile pericyte
after2.pos.cp <- as.character(prop.before2.cp[prop.before2.cp$pct.exp>=5,1])

# Myofibroblast (RSPO2+)
after2.pos.myo <- as.character(prop.before2.myo[prop.before2.myo$pct.exp>=5,1])

# Smooth muscle cell
after2.pos.smc <- as.character(prop.before2.smc[prop.before2.smc$pct.exp>=5,1])

# For all cell types
after2.pos.all <- intersect(intersect(after2.pos.cp, after2.pos.myo), after2.pos.smc)

## Filtering with coexpression proportion
second <- data.frame(genes=after2.pos.all)
second$coex.prop <- '.'

expr <- FetchData(mes.or, vars='OR51E2')
or.num <- length(WhichCells(mes.or[,expr>0]))
for (i in 1:nrow(second)) {
  expr2 <- FetchData(mes.or, vars=second[i,1])
  co.num <- length(WhichCells(mes.or[,expr>0 & expr2>0]))
  second[i,2] <- co.num/or.num*100
}
second$coex.prop <- as.numeric(second$coex.prop)
after2.pos.all <- second[second$coex.prop>=5,1]

## Results
# Visualization
FeaturePlot(mes.or, reduction='umap.rpca', features=after2.pos.all, order=T)
ggsave('Output/OR51E2-positive_featplot_2nd_CP-positive.png', width=2200, height=1800, unit='px', bg='transparent', scale=2)

# Saving
second <- data.frame(genes=after2.pos.all)

write.table(second, 'All/ppi_2nd_afterFilter.txt', col.names=T, row.names=F, sep='\t', quote=F)

#-------------------------------------------

### 3rd filtering
## Candidates before filtering
before3 <- read.csv('All/ppi_3rd_all_beforeFilter.csv', header=T)[,3]
before3 <- setdiff(before3, c(after1.pos.all, after2.pos.all, 'OR51E2'))
before3[before3=='MARCHF7'] <- 'MARCH7'

## Checking expression
# Cell proportion
prop.before3 <- DotPlot(mes.or, features=before3, group.by='celltype.detailed')$data[,c(3,4,2)]
prop.before3.cp <- prop.before3[prop.before3$id=='Contractile pericyte',]
prop.before3.myo <- prop.before3[prop.before3$id=='Myofibroblast (RSPO2+)',]
prop.before3.smc <- prop.before3[prop.before3$id=='Smooth muscle cell',]

## Classfying
# Contractile pericyte
after3.pos.cp <- as.character(prop.before3.cp[prop.before3.cp$pct.exp>=5,1])

# Myofibroblast (RSPO2+)
after3.pos.myo <- as.character(prop.before3.myo[prop.before3.myo$pct.exp>=5,1])

# Smooth muscle cell
after3.pos.smc <- as.character(prop.before3.smc[prop.before3.smc$pct.exp>=5,1])

# For all cell types
after3.pos.all <- intersect(intersect(after3.pos.cp, after3.pos.myo), after3.pos.smc)

## Filtering with coexpression proportion
third <- data.frame(genes=after3.pos.all)
third$coex.prop <- '.'

expr <- FetchData(mes.or, vars='OR51E2')
or.num <- length(WhichCells(mes.or[,expr>0]))
for (i in 1:nrow(third)) {
  expr2 <- FetchData(mes.or, vars=third[i,1])
  co.num <- length(WhichCells(mes.or[,expr>0 & expr2>0]))
  third[i,2] <- co.num/or.num*100
}

third$coex.prop <- as.numeric(third$coex.prop)
after3.pos.all <- third[third$coex.prop>=5,1]

## Results
# Visualization
FeaturePlot(mes.or, reduction='umap.rpca', features=after3.pos.all, order=T)
ggsave('Output/OR51E2-positive_featplot_3rd_CP-positive.png', width=2200, height=1800, unit='px', bg='transparent', scale=2)

# Saving
third <- data.frame(genes=after3.pos.all)

write.table(third, 'All/ppi_3rd_afterFilter.txt', col.names=T, row.names=F, sep='\t', quote=F)


#-------------------------------------------

### 4th filtering
## Candidates before filtering
before4 <- read.csv('All/ppi_4th_all_beforeFilter.csv', header=T)[,3]
before4 <- setdiff(before4, c(after1.pos.all, after2.pos.all, after3.pos.all, 'OR51E2'))

## Checking expression
# Cell proportion
prop.before4 <- DotPlot(mes.or, features=before4, group.by='celltype.detailed')$data[,c(3,4,2)]
prop.before4.cp <- prop.before4[prop.before4$id=='Contractile pericyte',]
prop.before4.myo <- prop.before4[prop.before4$id=='Myofibroblast (RSPO2+)',]
prop.before4.smc <- prop.before4[prop.before4$id=='Smooth muscle cell',]

## Classfying
# Contractile pericyte
after4.pos.cp <- as.character(prop.before4.cp[prop.before4.cp$pct.exp>=5,1])

# Myofibroblast (RSPO2+)
after4.pos.myo <- as.character(prop.before4.myo[prop.before4.myo$pct.exp>=5,1])

# Smooth muscle cell
after4.pos.smc <- as.character(prop.before4.smc[prop.before4.smc$pct.exp>=5,1])

# For all cell types
after4.pos.all <- intersect(intersect(after4.pos.cp, after4.pos.myo), after4.pos.smc)

## Filtering with coexpression proportion
fourth <- data.frame(genes=after4.pos.all)
fourth$coex.prop <- '.'

expr <- FetchData(mes.or, vars='OR51E2')
or.num <- length(WhichCells(mes.or[,expr>0]))
for (i in 1:nrow(fourth)) {
  expr2 <- FetchData(mes.or, vars=fourth[i,1])
  co.num <- length(WhichCells(mes.or[,expr>0 & expr2>0]))
  fourth[i,2] <- co.num/or.num*100
}
fourth$coex.prop <- as.numeric(fourth$coex.prop)

after4.pos.all <- fourth[fourth$coex.prop>=5,1]

## Results
# Visualization
FeaturePlot(mes.or, reduction='umap.rpca', features=after4.pos.all, order=T)
ggsave('Output/OR51E2-positive_featplot_4th_CP-positive.png', width=2200, height=1800, unit='px', bg='transparent', scale=2)

# Saving
fourth <- data.frame(genes=after4.pos.all)

write.table(fourth, 'All/ppi_4th_afterFilter.txt', col.names=T, row.names=F, sep='\t', quote=F)


#-------------------------------------------

### 5th filtering
## Candidates before filtering
before5 <- read.csv('All/ppi_5th_all_beforeFilter.csv', header=T)[,3]
before5 <- setdiff(before5, c(after1.pos.all, after2.pos.all, after3.pos.all, after4.pos.all, 'OR51E2'))

## Checking expression
# Cell proportion
prop.before5 <- DotPlot(mes.or, features=before5, group.by='celltype.detailed')$data[,c(3,4,2)]
prop.before5.cp <- prop.before5[prop.before5$id=='Contractile pericyte',]
prop.before5.myo <- prop.before5[prop.before5$id=='Myofibroblast (RSPO2+)',]
prop.before5.smc <- prop.before5[prop.before5$id=='Smooth muscle cell',]

## Classfying
# Contractile pericyte
after5.pos.cp <- as.character(prop.before5.cp[prop.before5.cp$pct.exp>=5,1])

# Myofibroblast (RSPO2+)
after5.pos.myo <- as.character(prop.before5.myo[prop.before5.myo$pct.exp>=5,1])

# Smooth muscle cell
after5.pos.smc <- as.character(prop.before5.smc[prop.before5.smc$pct.exp>=5,1])

# For all cell types
after5.pos.all <- intersect(intersect(after5.pos.cp, after5.pos.myo), after5.pos.smc)

## Filtering with coexpression proportion
fifth <- data.frame(genes=after5.pos.all)
fifth$coex.prop <- '.'

expr <- FetchData(mes.or, vars='OR51E2')
or.num <- length(WhichCells(mes.or[,expr>0]))
for (i in 1:nrow(fifth)) {
  expr2 <- FetchData(mes.or, vars=fifth[i,1])
  co.num <- length(WhichCells(mes.or[,expr>0 & expr2>0]))
  fifth[i,2] <- co.num/or.num*100
}
fifth$coex.prop <- as.numeric(fifth$coex.prop)

after5.pos.all <- fifth[fifth$coex.prop>=5,1]

## Results
# Visualization
FeaturePlot(mes.or, reduction='umap.rpca', features=after5.pos.all, order=T)
ggsave('Output/OR51E2-positive_featplot_5th_CP-positive.png', width=2200, height=1800, unit='px', bg='transparent', scale=2)

# Saving
fifth <- data.frame(genes=after5.pos.all)

write.table(fifth, 'All/ppi_5th_afterFilter.txt', col.names=T, row.names=F, sep='\t', quote=F)


#-------------------------------------------

### 6th filtering
## Candidates before filtering
before6 <- read.csv('All/ppi_6th_all_beforeFilter.csv', header=T)[,3]
before6 <- setdiff(before6, c(after1.pos.all, after2.pos.all, after3.pos.all, after4.pos.all, after5.pos.all, 'OR51E2'))

## Checking expression
# Cell proportion
prop.before6 <- DotPlot(mes.or, features=before6, group.by='celltype.detailed')$data[,c(3,4,2)]
prop.before6.cp <- prop.before6[prop.before6$id=='Contractile pericyte',]
prop.before6.myo <- prop.before6[prop.before6$id=='Myofibroblast (RSPO2+)',]
prop.before6.smc <- prop.before6[prop.before6$id=='Smooth muscle cell',]

## Classfying
# Contractile pericyte
after6.pos.cp <- as.character(prop.before6.cp[prop.before6.cp$pct.exp>=5,1])

# Myofibroblast (RSPO2+)
after6.pos.myo <- as.character(prop.before6.myo[prop.before6.myo$pct.exp>=5,1])

# Smooth muscle cell
after6.pos.smc <- as.character(prop.before6.smc[prop.before6.smc$pct.exp>=5,1])

# For all cell types
after6.pos.all <- intersect(intersect(after6.pos.cp, after6.pos.myo), after6.pos.smc)

## Filtering with coexpression proportion
sixth <- data.frame(genes=after6.pos.all)
sixth$coex.prop <- '.'

expr <- FetchData(mes.or, vars='OR51E2')
or.num <- length(WhichCells(mes.or[,expr>0]))
for (i in 1:nrow(sixth)) {
  expr2 <- FetchData(mes.or, vars=sixth[i,1])
  co.num <- length(WhichCells(mes.or[,expr>0 & expr2>0]))
  sixth[i,2] <- co.num/or.num*100
}
sixth$coex.prop <- as.numeric(sixth$coex.prop)

after6.pos.all <- sixth[sixth$coex.prop>=5,1]

## Results
# Visualization
FeaturePlot(mes.or, reduction='umap.rpca', features=after6.pos.all, order=T)
ggsave('Output/OR51E2-positive_featplot_6th_CP-positive.png', width=2200, height=1800, unit='px', bg='transparent', scale=2)

# Saving
sixth <- data.frame(genes=after6.pos.all)

write.table(sixth, 'All/ppi_6th_afterFilter.txt', col.names=T, row.names=F, sep='\t', quote=F)


#-------------------------------------------

### Integration
## Merging all PPI results
ppi.all1 <- c('OR51E2', after1.pos.all) %>% unique()
ppi.all2 <- c('OR51E2', after1.pos.all, after2.pos.all) %>% unique()
ppi.all3 <- c('OR51E2', after1.pos.all, after2.pos.all, after3.pos.all) %>% unique()
ppi.all4 <- c('OR51E2', after1.pos.all, after2.pos.all, after3.pos.all, after4.pos.all)
ppi.all5 <- c('OR51E2', after1.pos.all, after2.pos.all, after3.pos.all, after4.pos.all, after5.pos.all)
ppi.all6 <- c('OR51E2', after1.pos.all, after2.pos.all, after3.pos.all, after4.pos.all, after5.pos.all, after6.pos.all)

## Specifying putative effector TFs
tf.ppi.all2 <- intersect(ppi.all2, tf.all)
tf.ppi.all3 <- intersect(ppi.all3, tf.all)


## Saving
write.table(data.frame(genes=ppi.all3), 'All/ppi_allStep3_all_afterFilter.txt', col.names=F, row.names=F, sep='\t', quote=F)
write.table(data.frame(genes=ppi.all4), 'All/ppi_allStep4_all_afterFilter.txt', col.names=F, row.names=F, sep='\t', quote=F)
write.table(data.frame(genes=ppi.all5), 'All/ppi_allStep5_all_afterFilter.txt', col.names=F, row.names=F, sep='\t', quote=F)
write.table(data.frame(genes=ppi.all6), 'All/ppi_allStep6_all_afterFilter.txt', col.names=F, row.names=F, sep='\t', quote=F)

