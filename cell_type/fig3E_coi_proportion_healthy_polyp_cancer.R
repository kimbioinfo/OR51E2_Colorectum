setwd('~/plc/single-cell/merged/Processing')
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(qs)
library(glue)
library(Seurat)
library(devEMF)


### Loading seurat object
## Seurat object with all kinds of cells
mes.sr <- qread('allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')


### Among mesenchymal cells
## Processing
mes.sr.coi <- mes.sr
mes.sr.coi@meta.data[mes.sr.coi@meta.data$celltype.detailed%in%c('Smooth muscle cell','Contractile pericyte','Myofibroblast (RSPO2+)'), 'celltype.detailed'] <- 'COI'
prop.samples.coi <- as.data.frame(table(mes.sr.coi$celltype.detailed, mes.sr.coi$sample.id))

GSE201348.normal <- sort(unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE201348_Normal',6]))
GSE201348.normal <- sort(c(GSE201348.normal[which(startsWith(GSE201348.normal, c('GSE201348:A')))]))
GSE201348.polyp <- sort(unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE201348_Polyp',6]))
GSE201348.cancer <- sort(unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE201348_Tumor',6]))
GSE201348.healthy <- sort(unique(mes.sr@meta.data[mes.sr@meta.data$study.condition=='GSE201348_Normal',6]))
GSE201348.healthy <- sort(setdiff(GSE201348.healthy, GSE201348.normal))

prop.GSE201348.normal <- prop.samples.coi[prop.samples.coi$Var2%in%GSE201348.normal,] %>% arrange(Var2)
prop.GSE201348.polyp <- prop.samples.coi[prop.samples.coi$Var2%in%GSE201348.polyp,] %>% arrange(Var2)
prop.GSE201348.cancer <- prop.samples.coi[prop.samples.coi$Var2%in%GSE201348.cancer,] %>% arrange(Var2)
prop.GSE201348.healthy <- prop.samples.coi[prop.samples.coi$Var2%in%GSE201348.healthy,] %>% arrange(Var2)


prop.GSE201348.normal.coi <- left_join(prop.GSE201348.normal %>% group_by(Var2) %>% summarise(Total=sum(Freq)),
                                       prop.GSE201348.normal[prop.GSE201348.normal$Var1=='COI',]) %>% mutate(Perc=Freq/Total*100, Type='Unaffected') %>% as.data.frame()
prop.GSE201348.polyp.coi <- left_join(prop.GSE201348.polyp %>% group_by(Var2) %>% summarise(Total=sum(Freq)),
                                      prop.GSE201348.polyp[prop.GSE201348.polyp$Var1=='COI',]) %>% mutate(Perc=Freq/Total*100, Type='Polyp') %>% as.data.frame()
prop.GSE201348.cancer.coi <- left_join(prop.GSE201348.cancer %>% group_by(Var2) %>% summarise(Total=sum(Freq)),
                                       prop.GSE201348.cancer[prop.GSE201348.cancer$Var1=='COI',]) %>% mutate(Perc=Freq/Total*100, Type='Cancer') %>% as.data.frame()
prop.GSE201348.healthy.coi <- left_join(prop.GSE201348.healthy %>% group_by(Var2) %>% summarise(Total=sum(Freq)),
                                        prop.GSE201348.healthy[prop.GSE201348.healthy$Var1=='COI',]) %>% mutate(Perc=Freq/Total*100, Type='Healthy') %>% as.data.frame()

## Plotting
# Per sample (Overlapping individuals across and within tissue types)
prop.GSE201348.all.coi <- rbind(prop.GSE201348.healthy.coi, prop.GSE201348.polyp.coi, prop.GSE201348.cancer.coi)# %>% dplyr::filter(Total>5)
prop.GSE201348.all.coi$Type <- factor(prop.GSE201348.all.coi$Type, levels=c('Healthy','Polyp','Cancer'))
prop.GSE201348.all.coi$Donor <- apply(str_split(prop.GSE201348.all.coi$Var2,'-',simplify=T)[,1:2], 1, paste, collapse='-')
table(prop.GSE201348.all.coi$Type)

# Per individual
prop.GSE201348.indiv.coi <- prop.GSE201348.all.coi[-which(prop.GSE201348.all.coi$Donor%in%c('GSE201348:A001-C','GSE201348:A002-C') &
                                                            prop.GSE201348.all.coi$Type=='Polyp'),] %>% group_by(Donor) %>% summarise(Type=unique(Type), Total=sum(Total), Freq=sum(Freq)) %>% mutate(Perc=Freq/Total*100)
table(prop.GSE201348.indiv.coi$Type)   # Healthy: 2 / Polyp: 4 / Cancer: 4

emf('~/plc/Figures/Fig3/fig3e_coi_proportion_healthy_polyp_cancer_GSE201348.emf', width=120, height=70, units='mm', coordDPI=500, emfPlus=FALSE)
ggplot(prop.GSE201348.indiv.coi) +
  geom_boxplot(aes(x=Type, y=Perc, fill=Type)) +
  theme_classic(base_size=17) +
  labs(x=NULL, y='Percentage (%)', fill='Condition')
dev.off()


### Statistical significance
## Preprocessing data
prop.GSE201348.indiv.coi[prop.GSE201348.indiv.coi=='GSE201348:F072B-'] <- 'GSE201348:F072B'
prop.GSE201348.indiv.coi$Donor <- as.factor(prop.GSE201348.indiv.coi$Donor)

## Fitting the model
fit <- lm(data=prop.GSE201348.indiv.coi, formula=Perc~Type)
summary(fit)

library(DHARMa)
test <- simulateResiduals(fittedModel = fit)
plot(test)
DHARMa::plotQQunif(test, testUniformity=F, testDispersion=F, testOutliers=F)

## ANOVA
bartlett.test(prop.GSE201348.indiv.coi$Perc~prop.GSE201348.indiv.coi$Type)   # Passed
anova <- aov(data=prop.GSE201348.indiv.coi, formula=Perc~Type)
summary(anova)

## Post-hoc tests
shapiro.test(prop.GSE201348.indiv.coi$Perc)
pairwise.t.test(prop.GSE201348.indiv.coi$Perc, prop.GSE201348.indiv.coi$Type, p.adjust.method='bonferroni')

pairwise.wilcox.test(prop.GSE201348.indiv.coi$Perc, prop.GSE201348.indiv.coi$Type, p.adjust.method='bonferroni')

TukeyHSD(anova)

