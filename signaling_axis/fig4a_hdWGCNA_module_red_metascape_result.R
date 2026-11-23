setwd('~/plc/single-cell/network/')
library(dplyr)
library(tibble)
library(ggplot2)
library(devEMF)

pos <- read.csv('metascape_final_results_red_posCor.csv') %>%
  select(GO,Description,LogP,X.GeneInGO,X.GeneInGOAndHitList,Log.q.value.,GROUP_ID) %>% 
  dplyr::rename(InList=X.GeneInGO, InTerm=X.GeneInGOAndHitList) %>% 
  mutate(negLogP=-LogP, negLogQ=-Log.q.value., GeneRatio=InTerm/InList)
pos.top.qval <- pos %>% group_by(GROUP_ID) %>% summarise(negLogQ=max(negLogQ))
pos.top <- right_join(pos, pos.top.qval)
pos.top.filtered <- pos.top %>% filter(negLogQ>-log10(0.05)) %>% arrange(desc(negLogQ))
pos.top.filtered$direction <- 'Positive'

neg <- read.csv('metascape_final_results_red_negCor.csv') %>%
  select(GO,Description,LogP,X.GeneInGO,X.GeneInGOAndHitList,Log.q.value.,GROUP_ID) %>% 
  dplyr::rename(InList=X.GeneInGO, InTerm=X.GeneInGOAndHitList) %>% 
  mutate(negLogP=-LogP, negLogQ=-Log.q.value., GeneRatio=InTerm/InList)
neg.top.qval <- neg %>% group_by(GROUP_ID) %>% summarise(negLogQ=max(negLogQ))
neg.top <- right_join(neg, neg.top.qval)
neg.top.filtered <- neg.top %>% filter(negLogQ>-log10(0.05)) %>% arrange(desc(negLogQ))
neg.top.filtered$direction <- 'Negative'

top.filtered <- rbind(pos.top.filtered, neg.top.filtered)
top.filtered <- top.filtered[-which(duplicated(paste0(top.filtered$direction,'_',top.filtered$GROUP_ID))),]

### Dot plot
setwd('/mnt/nepi3_22TB/nepi_SW/plc/Figures/Fig4')

emf('fig4a_hdWGCNA_module_red_metascape_result.emf', width=230, height=90, units='mm', coordDPI=500, emfPlus=FALSE)
ggplot(top.filtered, aes(x=negLogQ, y=reorder(Description, negLogQ), color=direction, size=GeneRatio)) +
  geom_point() +
  labs(x="-log10(q-value)", y=NULL, title="GO Biological Processes", color='Correlation', size='Gene ratio') +
  coord_cartesian(xlim=c(0,15)) +
  scale_size_continuous(limits=c(0,0.2), breaks=c(0.05,0.10,0.15,0.20)) +
  theme_bw(base_size=16)
dev.off()

