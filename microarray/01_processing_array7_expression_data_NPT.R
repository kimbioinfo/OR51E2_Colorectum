setwd('/mnt/nepi3_22TB/nepi_SW/plc/network/array_7/Rawdata')
library(dplyr)
library(stringr)
library(tibble)
library(affy)
library(tidyverse)
library(GEOquery)
library(hthgu133pluspm.db)
library(ggplot2)
library(ggpubr)
library(qs)


### Loading CEL files
cel.npt <- ReadAffy(celfile.path=c('Normal_matched/','Polyp_matched/','Tumor_matched/'))


### normalization
## RMA
exp.npt <- rma(cel.npt)


### Getting expression matrices
## From RMA set
exp.npt <- as.data.frame(exprs(exp.npt))


### Matching probe IDs and gene symbols
## Retrieving a feature set
features <- as.data.frame(hthgu133pluspmSYMBOL)
colnames(features)[1] <- 'ID'
rownames(features) <- NULL

## Matching with expression data
exp.npt <- exp.npt %>%
  rownames_to_column(var='ID') %>% 
  inner_join(., features, by='ID')
exp.npt <- exp.npt[,c(1,173,2:172)]

## Modifying column names
colnames(exp.npt)[-c(1,2)] <- str_split(colnames(exp.npt)[-c(1,2)], '_', simplify=T)[,1]

## Probe filtering
# Only probes with 'PM_at'
features.at <- rbind(features[endsWith(features$ID, 'PM_at'),], features[43009:43018,])
features.sat <- features[setdiff(rownames(features), rownames(features.at)),]

exp.npt.sat <- exp.npt %>% inner_join(., features.sat)
exp.npt <- exp.npt %>% inner_join(., features.at)


### Checking gene expression patterns
## Between groups
t.test(as.numeric(apply(exp.npt[exp.npt$symbol=='OR51E2', c(3:59)], 2, mean)),
       as.numeric(apply(exp.npt[exp.npt$symbol=='OR51E2', c(60:116)], 2, mean)))
t.test(as.numeric(apply(exp.npt[exp.npt$symbol=='OR51E2', c(60:116)], 2, mean)),
       as.numeric(apply(exp.npt[exp.npt$symbol=='OR51E2', c(117:173)], 2, mean)))

# One-way ANOVA
or51e2.match.melt <- data.frame(exp=c(or51e2.normal, or51e2.polyp, or51e2.tumor),
                                tissue=c(rep('Normal',57), rep('Polyp',57), rep('Tumor',57)))
aov.or <- aov(exp~tissue, or51e2.match.melt)
summary(aov.or)
TukeyHSD(aov.or)
plot(TukeyHSD(aov.or))
boxplot(or51e2.normal, or51e2.polyp, or51e2.tumor, main='OR51E2', xlab='Tissue type', ylab='Normalized level', names=c('Normal','Polyp','Tumor'))
summary(aov.or, split = list(tissue=list(linear=1,quadratic=2)))

# p for trend test (Linear mixed model)
library(lmerTest)
normal.patient.id <- as.character(as.data.frame(str_split(match$Name_normal, ' ', simplify=T))[,6])
polyp.patient.id <- as.character(as.data.frame(str_split(match$Name_polyp, ' ', simplify=T))[,5])
tumor.patient.id <- as.character(as.data.frame(str_split(match$Name_tumor, ' ', simplify=T))[,5])

or51e2.match.melt$patientID <- c(normal.patient.id, polyp.patient.id, tumor.patient.id)
or51e2.match.melt$patientID <- as.factor(or51e2.match.melt$patientID)

or51e2.match.melt$tissue <- as.character(or51e2.match.melt$tissue)
or51e2.match.melt[or51e2.match.melt=='Normal'] <- 1
or51e2.match.melt[or51e2.match.melt=='Polyp'] <- 2
or51e2.match.melt[or51e2.match.melt=='Tumor'] <- 3
or51e2.match.melt$tissue <- as.numeric(or51e2.match.melt$tissue)

lme.or <- lmer(data=or51e2.match.melt, formula=exp~tissue+(1|patientID))
summary(lme.or)

ggplot(or51e2.match.melt, aes(x=as.character(tissue), y=exp, fill=tissue)) + geom_bar(stat='summary', fun='mean') +
  stat_summary(geom='errorbar', width=0.15) +
  xlab('Type of tissue') + ylab('Expression level') +
  scale_x_discrete(labels=c('Normal','Polyp','Tumor')) +
  theme(text=element_text(size=18))


