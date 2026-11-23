setwd('~/plc/single-cell/merged/Processing')
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(Seurat)
library(devEMF)
library(qs)


### Fig 3AB (Whole cells)
## Loading a seurat object
all.sr <- qread('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing/allPublic_all_afterRPCA_afterRemovingDoublet_simply_annotated.qs')

## Plotting
# Cell types
cols <- c('#FFCF9D', '#219C90', '#FFB000', '#D83F31', '#004225')
all.sr$celltype.simple <- factor(all.sr$celltype.simple, levels = c('Endothelial', 'Epithelial', 'Immune', 'Mesenchymal', 'Neural'))

p <- DimPlot(all.sr, reduction='umap.rpca', group.by='celltype.simple', raster=F, cols = cols) +
  theme_void() + 
  ggtitle(NULL) + 
  NoLegend()
ggsave('~/plc/Figures/Fig3/fig3a_umap_all_cells_with_celltypes.jpg', p, width = 4.5, height = 4.5, scale = 0.8, dpi=300)

# OR51E2
p <- FeaturePlot(all.sr, reduction='umap.rpca', features='OR51E2', raster=F, order=T)+
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = 4),
        legend.text = element_text(size = 20))

p <- FeaturePlot(all.sr, reduction='umap.rpca', features='OR51E2', raster=F, order=T)+
  theme_void() +
  ggtitle(NULL) +
  theme(legend.text = element_text(size = 20))
ggsave('~/plc/Figures/Fig3/fig3b_umap_all_cells_with_or51e2.jpg', p, width = 5.3, height = 4.5, scale = 0.8, dpi=300)


### Fig 3CD (Mesenchymal cells)
## Loading a seurat object
mes.sr <- qread('/mnt/nepi3_22TB/nepi_SW/plc/single-cell/merged/Processing/allPublic_mes_afterRPCA_afterRemovingDoublet_fully_annotated.qs')

## Plotting
# Cell types
cols <- c('#AC233B', '#F67F98', '#F78E3B', '#F46052', '#F68F81', '#81D2F3', '#C3E7F9', '#70CBD3', '#1AA0AE', '#1FA4DD', '#9278B5', '#AAD54E')
mes.sr$celltype.detailed <- factor(mes.sr$celltype.detailed,
                                   levels = c('Contractile pericyte',
                                              'Pericyte',
                                              'Myofibroblast',
                                              'Myofibroblast (RSPO2+)',
                                              'Smooth muscle cell',
                                              'Stromal 1',
                                              'Stromal 2',
                                              'Stromal 3',
                                              'Transitional Stromal 3',
                                              'T reticular & mLTo',
                                              'FDC',
                                              'ICC'))

p <- DimPlot(mes.sr, reduction='umap.rpca', group.by='celltype.detailed', raster=F, cols = cols) +
  theme_void() + 
  ggtitle(NULL) + 
  NoLegend()
ggsave('~/plc/Figures/Fig3/fig3c_umap_mesenchymal_cells_with_celltypes.jpg', p, width = 4.5, height = 4.5, scale = 0.8, dpi=300)

# OR51E2
p <- FeaturePlot(mes.sr, reduction='umap.rpca', features='OR51E2', raster=F, order=T) +
  theme_void() +
  ggtitle(NULL) +
  theme(legend.text = element_text(size = 20))
ggsave('~/plc/Figures/Fig3/fig3d_umap_mesenchymal_cells_with_or51e2.jpg', p, width = 5.3, height = 4.5, scale = 0.8, dpi=300)

