library(tidyverse)
library(Seurat)
library(viridis)

seuratobj = readRDS("seuratObj_diet_BD007_2024.rds") # spleen

# Simplify annotation

Idents(seuratobj) = "annotated_clusters"
seuratobj = RenameIdents(object = seuratobj, #remove immature and proliferating cells, simplify annotation
                         # "Neutrophils" = "", 
                         "Activated Neutrophils" = "Neutrophils", 
                         # "Basophils" = "", 
                         # "Mast cells" = "", 
                         "Monocytes 1"  = "Ly6c2hi Monocytes", 
                         "Monocytes 2" = "Ly6c2low Monocytes", 
                         "Proliferating monocytes" = "Other",
                         "Macrophages 1" = "Macrophages", 
                         "Macrophages 2" = "Macrophages", 
                         # "cDC1s" = "", 
                         "Proliferating cDC1s" = "Other",
                         "pre-cDC2s" = "Other", 
                         # "cDC2s" = "",  
                         "Proliferating cDC2s" = "Other",
                         # "Migratory DCs" = "", 
                         "Plasmacytoid cells" = "pDCs", 
                         "Immature B cells" = "Other", 
                         "Mature B cells" = "B cells", 
                         # "Plasma cells" = "", 
                         "Naive T cells" = "T cells", 
                         "T helper cells or ILCs" = "T cells", 
                         # "Tregs" = "",
                         "Activated Th2 cells" = "T cells", 
                         "Memory CD8 T cells" = "T cells", 
                         "Gamma delta T cells" = "T cells", 
                         "T cells/NK cells" = "T cells", 
                         # "NK cells" = "", 
                         "Eythroid cells" = "Erythroid cells",
                         "Unknown cells" = "Eosinophils")
Idents(seuratobj) = factor(Idents(seuratobj), levels = c(
  "Neutrophils", "Eosinophils", "Basophils", "Mast cells", "Ly6c2hi Monocytes" , "Ly6c2low Monocytes", "Proliferating monocytes", 
  "Macrophages", "cDC1s", "Proliferating cDC1s", "pre-cDC2s", "cDC2s",  "Proliferating cDC2s", "Migratory DCs" , 
  "pDCs", "B cells", "Plasma cells", "T cells", "Tregs", "NK cells", "Erythroid cells", "Other"))
table(seuratobj@meta.data$annotation_pooled)
table(Idents(seuratobj))
seuratobj@meta.data$annotation_pooled = Idents(seuratobj)

Fc_receptor_atlas_colors = c(
  
  # cell populations
  "neutros" = "#7f0000", "neutrophils" = "#7f0000", "granulocyte" = "#7f0000", "Neutrophils" = "#7f0000", 
  "eos" = "#d40000", "Eosinophils" = "#d40000",
  "basos" = "#ff8585", "basophil" = "#ff8585", "basophils" = "#ff8585", "Basophils" = "#ff8585",
  "mast cells" = "#7f3f00", "Mast cells" = "#7f3f00",
  "Monos" = "#cd9966", "monocyte" = "#cd9966",  "monocytes" = "#cd9966",
  "class monos" = "#cd9966",
  "Ly6Chi monos" = "#cd9966", "Ly6c2hi Monocytes" = "#cd9966",
  "non-class monos" = "#e6a666", "Ly6c2low Monocytes" = "#e6a666",
  "patr monos" = "#e6a666",
  "macs" = "#ffc891", "macrophage" = "#ffc891", "macrophages" = "#ffc891", "Macrophages" = "#ffc891",# same as RPMs, AMs, KCs, LPMs
  "cDC1s" = "#828200", "cDCs" = "#828200",
  "cDC2s" = "#acac00",
  "Migratory DCs" = "#c2c202",
  "pDCs" = "#d4d400",
  "B cells" = "#3f8000", "B cell" = "#3f8000",
  "plasma cells" = "#52990e", "Plasma cells" = "#52990e",
  "CD8+ T cells" = "#008241",
  "CD4+ T cells" = "#00d66b",
  "T cells" = "#00d66b",   "T cell" = "#00d66b",
  "Tregs" = "#48ffff",
  "NK cells" = "#003f7f", "natural killer cell" = "#003f7f",
  "platelets" = "purple",
  "Erythroid cells" = "#653780",
  "endothelial" = "grey25",
  "epithelial" = "grey40",
  "stromal" = "grey55",
  "parenchymal" = "grey70",
  "FMO" = "#4c4c4c",
  "Other" = "grey80",
  
  # Fc receptors 
  "FcγRI" = "#e00480", # shared human and mouse
  "FcγRIIA" = "#22a395", # human only
  "FcγRIIB" = "#ff9a00", # shared human and mouse
  "FcγRIIC" = "#dac600", # human only
  "FcγRIII" = "#0045b9", # mouse only, likely related to human FcγRIIA
  "FcγRIIIA" = "#006131", # human only
  "FcγRIIIB" = "#009700", # human only
  "FcγRIV" = "#01bb00", # mouse only, ortholog of human FcγRIIIA/B
  "FcRn" = "#ff5c27", # shared human and mouse
  
  # groups
  "wild-type" = "#0080ff", # alpha = 0.6
  "hFcγR" = "#ff7f00", # alpha = 0.6
  "FcεR1γ-/-" = "#4c4c4c" # alpha = 0.6
)

DimPlot(seuratobj, label = TRUE, repel = TRUE, label.size = 3.7) +
  NoLegend() +
  NoAxes() + scale_color_manual(values = Fc_receptor_atlas_colors)

FeaturePlot(seuratobj, "Jchain")
saveRDS(seuratobj, "sorted_splenocytes_wt_females.RDS")


cluster_markers = FindAllMarkers(subset(seuratobj, downsample = 100), logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE)
# write.csv(cluster_markers, "mouse_enriched_splenocytes_findallmarkers.csv")

topn = cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DoHeatmap(subset(seuratobj, downsample=10), features = topn$gene, 
          # assay = "RNA", 
          size = 3,
          slot = "counts",
          label = TRUE, disp.min = -0.5, draw.lines = TRUE, lines.width = 1) + 
  NoLegend() +
  theme(text = element_text(size = 9)) +
  scale_fill_binned(type = "viridis", na.value = "white")

# Plot Fc gamma receptors

fcreceptors = c("Fcgr1", "Fcgr2b", "Fcgr3", "Fcgr4", "Fcer1g", "Fcgrt", "Trim21") # murine Fc receptors
fcreceptors = c("Chia1", "Chil1", "Chil3", "Chil4", "Chit1") # murine Fc receptors
fcreceptors = c("Fcgr1", "Fcgr2b", "Fcgr3", "Fcgr4", "Fcer1g", "Fcgrt", "Trim21") # murine Fc receptors

plotList <- lapply(
  fcreceptors,
  function(i) {
    plot = FeaturePlot(seuratobj, i, cols = c("grey", viridis(10)), 
                       order = T,
                       max.cutoff = "q95") + NoAxes() + theme(title = element_text(face = "plain"))
    print(plot)
  }
)
allplots <- ggpubr::ggarrange(plotlist=plotList, ncol = 3, nrow = 3, align = "hv", common.legend = TRUE, legend = "right")
allplots

DotPlot(seuratobj, features = rev(fcreceptors), 

        cols= "RdYlBu",
        scale.by = "radius",
        scale = T) + 
  # ggtitle("Lung cells of COVID-19 patients") +
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) + 
  theme(axis.text=element_text(size=9)) +
  theme(axis.title=element_blank()) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=10)) +
  scale_size(breaks = c(0, 40, 80)) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0)
  ) + 
  scale_y_discrete(position = "right") +
  coord_flip()

plotList <- lapply(
  fcreceptors,
  function(i) {
    df = FetchData(seuratobj, vars = c(i, "annotation_pooled"))
    colnames(df) <- c("parameter", "celltype")
    df$celltype = factor(df$celltype, levels = rev(c(
      "Neutrophils", "Eosinophils", "Basophils", "Mast cells", "Ly6c2hi Monocytes" , "Ly6c2low Monocytes", "Proliferating monocytes", 
      "Macrophages", "cDC1s", "Proliferating cDC1s", "pre-cDC2s", "cDC2s",  "Proliferating cDC2s", "Migratory DCs" , 
      "pDCs", "B cells", "Plasma cells", "T cells", "Tregs", "NK cells", "Erythroid cells", "Other")))
    plot = ggplot(df, aes(y=celltype, x=parameter, fill=celltype)) + 
      ggridges::geom_density_ridges(aes(height =..ndensity..), scale = 1.5, alpha = 0.5,
                                    color = "grey70", size = 0.3) + 
      theme(axis.title = element_blank()) +
      NoLegend() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), #axis.text.y.left = element_blank(),
            axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text = element_text(colour = "grey30"), axis.text.y = element_text(vjust = -1),
            plot.title = element_text(hjust = 0.5)) + 
      scale_fill_manual(values = Fc_receptor_atlas_colors) +
      ggtitle(i)
    print(plot)
  }
)
allplots <- ggpubr::ggarrange(plotlist=plotList, ncol = 7, nrow = 1)
allplots