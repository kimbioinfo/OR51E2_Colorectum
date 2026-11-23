setwd('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/network/')
library(Seurat)
library(dplyr)
library(stringr)
library(tibble)
library(qs)
library(ggplot2)
library(glue)
library(patchwork)
library(devEMF)


### Loading materials
mes.or <- qread('allPublic_or51e2-positive_hdWGCNAed_Motifed_seuratObject_red.qs')
motif.red <- read.table('motifs_for_or51e2_module_red.txt', header=T, sep='\t') %>% arrange(fdr)


### Filtering with expression
## Refining TF names
df.temp1 <- as.data.frame(str_split(motif.red$tf, pattern='\\(', simplify=T))
tf.temp1 <- df.temp1$V1
motif.red$tf <- tf.temp1

df.temp2 <- as.data.frame(str_split(tf.temp1, pattern='::', simplify=T))
tf.temp2 <- c(df.temp2$V1, df.temp2$V2, df.temp2$V3)
unique(tf.temp2)
tf.temp2 <- tf.temp2[nchar(tf.temp2)>0]
tf.red <- unique(toupper(tf.temp2))
tf.red <- intersect(tf.red, rownames(mes.or@assays$RNA$counts))   # 167 of 168 TFs

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

## Selecting only TFs with expression
which.include <- NULL

for (i in 1:nrow(motif.red)) {
  temp.tf <- motif.red$tf[i]
  if (str_detect(temp.tf, '::', negate=T)) {
    if (temp.tf%in%tf.all) {
      which.include <- append(which.include, i)
    }
  }
  if (str_detect(temp.tf, '::', negate=F)) {
    temp.tfs <- str_split(temp.tf,'::')[[1]]
    if (all(temp.tfs%in%tf.all)) {
      which.include <- append(which.include, i)
    }
  }
}

motif.red.filt <- motif.red[which.include,]
rownames(motif.red.filt) <- NULL


### Saving whole TF statistics
motif.red.stat <- motif.red
motif.red.stat$COI_expression <- 'No'
motif.red.stat[which.include, 'COI_expression'] <- 'Yes'
motif.red.stat <- motif.red.stat[,c(2,4:9,11,10)]
write.table(motif.red.stat, 'motifs_for_or51e2_module_red_withCOIexpression.txt', col.names=T, row.names=F, sep='\t', quote=F)

### Visualization
emf('~/plc/Figures/Fig4/fig4c_TF_motifs_red_genes_barplot_lineplot.emf', width=180, height=100, units='mm', coordDPI=500, emfPlus=FALSE)
ggplot(motif.red.filt[1:30,]) +
  geom_bar(aes(x=reorder(tf, fdr), y=-log10(fdr)), fill='#fb736d',stat='identity', position='dodge') +
  coord_cartesian(ylim=c(0,33)) +
  xlab(NULL) + ylab('-log10(FDR)') +
  geom_line(aes(x=reorder(tf, fdr), y=size_intersection/6, group=module), color='grey25') +
  #geom_bar(aes(x=reorder(tf, fdr), y=size_intersection/6), fill='grey25', stat='identity', position='dodge', width=0.07) +
  geom_point(aes(x=reorder(tf, fdr), y=size_intersection/6, group=module), shape=19, cex=1.2, color='grey25') +
  scale_y_continuous(sec.axis = sec_axis(trans=~.*6, name='Number of target genes')) +
  theme_bw(base_size=14) +
  theme(axis.text.x=element_text(angle=60, hjust=1.1))
dev.off()
