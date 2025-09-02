library(Seurat)
library(viridis)
library(tidyverse)

##### Data downloading and conversion to seurat #####

# Data downloaded via https://figshare.com/projects/Tabula_Sapiens/100973
library(anndata)
reticulate::install_python(version = '<version>')
data <- read_h5ad("ts_stromal.h5ad")
ts_stromal <- CreateSeuratObject(counts = MatrixExtra::t_shallow(data$X), meta.data = data$obs)
saveRDS(ts_stromal, "ts_stromal.rds")
rm(data)

data <- read_h5ad("ts_epithelial.h5ad")
ts_epithelial <- CreateSeuratObject(counts = MatrixExtra::t_shallow(data$X), meta.data = data$obs)
saveRDS(ts_epithelial, "ts_epithelial.rds")
rm(data)

data <- read_h5ad("ts_endothelial.h5ad")
ts_endothelial <- CreateSeuratObject(counts = MatrixExtra::t_shallow(data$X), meta.data = data$obs)
saveRDS(ts_endothelial, "ts_endothelial.rds")
rm(data)

data <- read_h5ad("ts_immune.h5ad")
ts_immune <- CreateSeuratObject(counts = MatrixExtra::t_shallow(data$X), meta.data = data$obs)
saveRDS(ts_immune, "ts_immune.rds")
rm(data)

##### Dataset reannotation and merging ####
ts_stromal = DietSeurat(readRDS("ts_stromal.rds"))
ts_epithelial = DietSeurat(readRDS("ts_epithelial.rds"))
ts_endothelial = DietSeurat(readRDS("ts_endothelial.rds"))
ts_immune = DietSeurat(readRDS("ts_immune_preproc_pooled_annot.rds"))

# Immune cell reannotation
Idents(ts_immune) = "cell_ontology_class"
names(table(Idents(ts_immune)))
table(Idents(ts_immune))
ts_immune = RenameIdents(object = ts_immune,
                         "b cell" = "B cells",
                         "basophil" = "basophils",
                         "cd1c-positive myeloid dendritic cell" = "cDCs",
                         "cd4-positive alpha-beta t cell" = "CD4+ T cells",
                         "cd4-positive, alpha-beta memory t cell" = "CD4+ T cells",
                         "cd4-positive, alpha-beta t cell" = "CD4+ T cells",
                         "cd4-positive helper t cell" = "CD4+ T cells",
                         "cd8-positive alpha-beta t cell" = "CD8+ T cells",
                         "cd8-positive, alpha-beta cytokine secreting effector t cell" = "CD8+ T cells",
                         "cd8-positive, alpha-beta cytotoxic t cell" = "CD8+ T cells",
                         "cd8-positive, alpha-beta memory t cell" = "CD8+ T cells",
                         "cd8-positive, alpha-beta t cell" = "CD8+ T cells",
                         "cd8b-positive nk t cell" = "NK cells",
                         "cd24 neutrophil" = "neutrophils",
                         "cd141-positive myeloid dendritic cell" = "cDCs",
                         "classical monocyte" = "monocytes",
                         "dendritic cell" = "cDCs",
                         "dn1 thymic pro-t cell" = "other",
                         "dn3 thymocyte" = "other",
                         "dn4 thymocyte" = "other",
                         "double-positive, alpha-beta thymocyte" = "other",
                         "erythrocyte" = "other",
                         "erythroid lineage cell" = "other",
                         "erythroid progenitor" = "other",
                         "granulocyte" = "neutrophils",
                         "hematopoietic stem cell" = "other",
                         "immature natural killer cell" = "other",
                         "immune cell" = "other",
                         "innate lymphoid cell" = "other",
                         "intermediate monocyte" = "monocytes",
                         "langerhans cell" = "cDCs",
                         "leucocyte" = "other",
                         "liver dendritic cell" = "cDCs",
                         "macrophage" = "macrophages",
                         "mast cell" = "mast cells",
                         "mature conventional dendritic cell" = "cDCs",
                         "mature nk t cell" = "other",
                         "memory b cell" = "B cells",
                         "mesenchymal stem cell" = "other",
                         "microglial cell" = "macrophages",
                         "monocyte" = "monocytes",
                         "myeloid cell" = "macrophages",
                         "myeloid dendritic cell" = "cDCs",
                         "myeloid progenitor" = "other",
                         "naive b cell" = "B cells",
                         "naive regulatory t cell" = "CD4+ T cells",
                         "naive thymus-derived cd4-positive, alpha-beta t cell" = "CD4+ T cells",
                         "naive thymus-derived cd8-positive, alpha-beta t cell" = "CD8+ T cells",
                         "nampt neutrophil" = "neutrophils",
                         "neutrophil" = "neutrophils",
                         "nk cell" = "NK cells",
                         "nkt cell" = "other",
                         "non-classical monocyte" = "monocytes",
                         "plasma cell" = "plasma cells",
                         "plasmablast" = "plasma cells",
                         "plasmacytoid dendritic cell" = "pDCs",
                         "platelet" = "platelets",
                         "regulatory t cell" = "Tregs",
                         "t cell" = "CD4+ T cells",
                         "t follicular helper cell" = "CD4+ T cells",
                         "thymocyte" = "other",
                         "type i nk t cell" = "other"
)
names(table(Idents(ts_immune)))
ts_immune@meta.data$annotation_pooled = Idents(ts_immune)

# Merging
seuratobj = merge(ts_stromal, y = c(ts_epithelial, ts_endothelial, ts_immune))

# Add Metadata
endothelial = WhichCells(ts_endothelial)
epithelial = WhichCells(ts_epithelial)
stromal = WhichCells(ts_stromal)

annot = as_tibble(seuratobj@meta.data$annotation_pooled)

rownames(annot) = colnames(seuratobj)
annot$value[which(rownames(annot) %in% stromal)] = "stromal"
rownames(annot) = colnames(seuratobj)
annot$value[which(rownames(annot) %in% endothelial)] = "endothelial"
rownames(annot) = colnames(seuratobj)
annot$value[which(rownames(annot) %in% epithelial)] = "epithelial"

table(annot$value)

seuratobj@meta.data$annotation_pooled = annot$value
table(seuratobj@meta.data$annotation_pooled)

seuratobj@meta.data$annotation_pooled = factor(seuratobj@meta.data$annotation_pooled, levels = c("neutrophils", "basophils", "mast cells", "monocytes", "macrophages", "cDCs", "pDCs",
                                                                                                 "B cells", "plasma cells", "CD8+ T cells", "CD4+ T cells", "Tregs", "NK cells", "platelets", "other",
                                                                                                 "endothelial", "epithelial", "stromal"))
Idents(seuratobj) = "annotation_pooled"
saveRDS(seuratobj, "ts_pooled.RDS")

rm(annot)
rm(annot_df)
rm(endothelial)
rm(stromal)
rm(epithelial)
rm(ts_endothelial)
rm(ts_epithelial)
rm(ts_stromal)
rm(ts_immune)

##### Preprocessing ####
library(Seurat)
library(viridis)
library(tidyverse)
seuratobj = readRDS("ts_pooled.rds")


Idents(seuratobj) = "annotation_pooled"

seuratobj[["percent.mt"]] <- PercentageFeatureSet(seuratobj, pattern = "^MT")
# VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0 , ncol = 1) + NoLegend()
VlnPlot(seuratobj, features = "nFeature_RNA", pt.size = 0) + NoLegend()
VlnPlot(seuratobj, features = "nCount_RNA", pt.size = 0) + NoLegend()
VlnPlot(seuratobj, features = "percent.mt", pt.size = 0) + NoLegend()

seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst"#, nfeatures = 2000
)
head(VariableFeatures(seuratobj), 10)

fcreceptors = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", 
                "FCER1G", "FCGRT", "TRIM21", "FCRL3", "FCRL4", "FCRL5", 
                "FCER1A", "FCER2", "FCAR", "FCAMR", "FCMR", "PIGR")

fcreceptors[fcreceptors %in% VariableFeatures(seuratobj)]
seuratobj <- ScaleData(seuratobj, features = #rownames(seuratobj)
                         c(VariableFeatures(seuratobj), fcreceptors)
                      )
fcreceptors[fcreceptors %in% VariableFeatures(seuratobj)]
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
fcreceptors[fcreceptors %in% rownames(seuratobj@assays[["RNA"]]@scale.data)]
# length(fcreceptors[fcreceptors %in% rownames(seuratobj@assays[["RNA"]]@scale.data)])
DimPlot(seuratobj, reduction = "pca", label = TRUE) # + NoLegend()
DimHeatmap(seuratobj, dims = 1:15, cells = 200, balanced = TRUE)
DimHeatmap(seuratobj, dims = 16:30, cells = 200, balanced = TRUE) 
ElbowPlot(seuratobj, ndims = 50) # Elbow at +- 30
seuratobj <- RunUMAP(seuratobj, dims = 1:30, verbose = TRUE)

Fc_receptor_atlas_colors = c(
  # cell populations
  "neutros" = "#7f0000", "neutrophils" = "#7f0000", "granulocyte" = "#7f0000",
  "eos" = "#d40000",
  "basos" = "#ff8585", "basophil" = "#ff8585", "basophils" = "#ff8585",
  "mast cells" = "#7f3f00",
  "Monos" = "#cd9966", "monocyte" = "#cd9966",  "monocytes" = "#cd9966",
  "class monos" = "#cd9966",
  "Ly6Chi monos" = "#cd9966",
  "non-class monos" = "#e6a666",
  "patr monos" = "#e6a666",
  "macs" = "#ffc891", "macrophage" = "#ffc891", "macrophages" = "#ffc891", # same as RPMs, AMs, KCs, LPMs
  "cDC1s" = "#828200", "cDCs" = "#828200",
  "cDC2s" = "#acac00",
  "pDCs" = "#d4d400",
  "B cells" = "#3f8000", "B cell" = "#3f8000",
  "plasma cells" = "#52990e",
  "CD8+ T cells" = "#008241",
  "CD4+ T cells" = "#00d66b",
  "T cells" = "#00d66b",   "T cell" = "#00d66b",
  "Tregs" = "#48ffff",
  "NK cells" = "#003f7f", "natural killer cell" = "#003f7f",
  "platelets" = "purple",
  "endothelial" = "grey25",
  "epithelial" = "grey40",
  "stromal" = "grey55",
  "parenchymal" = "grey70",
  "FMO" = "#4c4c4c",
  
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

Idents(seuratobj) = "annotation_pooled"
p1 = DimPlot(seuratobj, reduction = "umap", label = TRUE, label.size = 4) + NoAxes() +
  scale_color_manual(values = Fc_receptor_atlas_colors) + theme(legend.position = "bottom")
p2 = DimPlot(seuratobj, reduction = "umap", label = TRUE, label.size = 4, group.by = "organ_tissue") + NoAxes() + theme(legend.position = "bottom")
p3 = DimPlot(seuratobj, reduction = "umap", label = TRUE, label.size = 1.5, group.by = "cell_ontology_class") + NoAxes() + NoLegend()
ggpubr::ggarrange(p1, p2, p3, ncol = 3, align = "hv")
rm(p1) + rm(p2) + rm(p3)

saveRDS(seuratobj, "ts_pooled.rds")

##### Visualisation ####
fcreceptors = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", 
                "FCER1G", "FCGRT", "TRIM21", "FCRL3", "FCRL4", "FCRL5", 
                "FCER1A", "FCER2", "FCAR", "FCAMR", "FCMR", "PIGR")


# FeaturePlots
FeaturePlot(seuratobj, "MKI67", cols = c("grey", viridis(10)), order = T, max.cutoff = "q95") + NoAxes()
# i = "FCGRT"
plotList <- lapply(
  fcreceptors,
  function(i) {
    plot = FeaturePlot(seuratobj, i, cols = c("grey", viridis(10)), order = T, max.cutoff = "q95") + NoAxes() + NoLegend()
    print(plot)
  }
)
allplots <- ggpubr::ggarrange(plotlist=plotList, ncol = 6, nrow = 3, legend = "bottom", common.legend = TRUE)
allplots
rm(allplots) + rm(plotList) + rm (plot) + rm(i)

# Heatmap
DoHeatmap(subset(seuratobj, downsample = 100),
          slot = "scale.data",
          features = fcreceptors, size = 3.5, #disp.min = -0.5, disp.max = 3.5, 
          draw.lines = TRUE, lines.width = 4) + 
  NoLegend() +
  theme(text = element_text(size = 12)) +
  scale_fill_binned(type = "viridis", na.value = "white") +
  ggtitle("Tabula Sapiens (scale.data)")

# DotPlot
seuratobj@meta.data$annotation_pooled = factor(seuratobj@meta.data$annotation_pooled, 
                                               levels = rev(c("neutrophils", "basophils", "mast cells", "monocytes", "macrophages", "cDCs", "pDCs",
                                                          "B cells", "plasma cells", "CD8+ T cells", "CD4+ T cells", "Tregs", "NK cells", "platelets", "other",
                                                          "endothelial", "epithelial", "stromal")))
Idents(seuratobj) = "annotation_pooled"

DotPlot(subset(seuratobj, downsample = 1000), 
        features = fcreceptors,
        # assay = "RNA",
        # scale = FALSE,
        group.by = "annotation_pooled", 
        # col.min = -1,
        # col.max = 2,
        # dot.min = 0,
        # dot.scale = 5,
        cols= "RdYlBu",
        scale.by = "radius",
        # idents = cells_more_than_250
        ) + 
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
  scale_y_discrete(position = "left") +
  scale_x_discrete(position = "top") +
  coord_flip()

# Ridgeplot
seuratobj@meta.data$annotation_pooled = factor(seuratobj@meta.data$annotation_pooled, 
                                               levels = c("neutrophils", "basophils", "mast cells", "monocytes", "macrophages", "cDCs", "pDCs",
                                                              "B cells", "plasma cells", "CD8+ T cells", "CD4+ T cells", "Tregs", "NK cells", "platelets", "other",
                                                              "endothelial", "epithelial", "stromal"))
Idents(seuratobj) = "annotation_pooled"
seuratobj_down = subset(seuratobj, downsample = 1000)
# i = "FCGR3A"

plotList <- lapply(
  fcreceptors,
  function(i) {
    df = FetchData(seuratobj_down, vars = c(i, "annotation_pooled"), slot = "scale.data")
    colnames(df) <- c("parameter", "celltype")
    df$celltype = factor(df$celltype, levels = rev(c("neutrophils", "basophils", "mast cells", "monocytes", "macrophages", "cDCs", "pDCs",
                                                     "B cells", "plasma cells", "CD8+ T cells", "CD4+ T cells", "Tregs", "NK cells", "platelets", "other",
                                                     "endothelial", "epithelial", "stromal")))
    plot = ggplot(df, aes(y=celltype, x=parameter, fill=celltype)) + 
      ggridges::geom_density_ridges(aes(height =..ndensity..), scale = 1.5, alpha = 0.5,
                                    color = "grey70", size = 0.3) + 
      theme(axis.title = element_blank()) +
      NoLegend() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.text.y.left = element_blank(),
            axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text = element_text(colour = "grey30"), axis.text.y = element_text(vjust = -0.5),
            plot.title = element_text(hjust = 0.5)) + 
      scale_fill_manual(values = Fc_receptor_atlas_colors) +
      ggtitle(i)
    print(plot)
  }
)
allplots <- ggpubr::ggarrange(plotlist=plotList, ncol = 6, nrow = 3)
allplots

# Specific ridgeplots

# Endothelial, epithelial and stromal
Idents(seuratobj) = "annotation_pooled"
table(Idents(seuratobj))
seuratobj_compartment = subset(seuratobj, idents = c("endothelial", "epithelial", "stromal"), downsample = 5000)

df = FetchData(seuratobj_compartment, vars = c("FCGRT", "annotation_pooled"), slot = "scale.data")
colnames(df) <- c("parameter", "celltype")
df$celltype = factor(df$celltype, levels = rev(c("endothelial", "stromal", "epithelial")))
ggplot(df, aes(y=celltype, x=parameter, color=celltype, fill = celltype)) + 
  ggridges::geom_density_ridges(aes(height =..ndensity..), scale = 1.5, alpha = 0.1, size = 0.3) + 
  theme(axis.title = element_blank()) +
  NoLegend() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        # axis.text.y.left = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text = element_text(colour = "grey30"), axis.text.y = element_text(vjust = -0.5),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("endothelial" = "grey5",
                                "epithelial" = "#300030",
                                "stromal" = "#400000")) +
  scale_fill_manual(values = c("endothelial" = "grey5",
                               "epithelial" = "#300030",
                               "stromal" = "#400000")) +
  xlim(NA, 3)

# Rigdeplot for stromal cells (grouped fibroblasts, smooth muscle, skeletal muscle)
Idents(seuratobj) = "annotation_pooled"
seuratobj_compartment = subset(seuratobj, idents = "stromal")
table(seuratobj_compartment@meta.data$cell_ontology_class, seuratobj_compartment@meta.data$organ_tissue)
seuratobj_stromal_ridgeplot = subset(seuratobj_compartment,
                                     idents = c(
                                       "alveolar fibroblast", "fibroblast", "fibroblast of breast", "fibroblast of cardiac tissue", "myofibroblast cell", "pancreatic stellate cell", 
                                       "bronchial smooth muscle cell", "smooth muscle cell", "vascular associated smooth muscle cell", 
                                       "fast muscle cell", "slow muscle cell"))
fibroblasts = WhichCells(seuratobj_stromal_ridgeplot, 
                         idents = c("alveolar fibroblast", "fibroblast", "fibroblast of breast", "fibroblast of cardiac tissue", "myofibroblast cell", "pancreatic stellate cell"))
smooth_muscle = WhichCells(seuratobj_stromal_ridgeplot, 
                           idents = c("bronchial smooth muscle cell", "smooth muscle cell", "vascular associated smooth muscle cell"))
skeletal_muscle = WhichCells(seuratobj_stromal_ridgeplot, 
                           idents = c("fast muscle cell", "slow muscle cell"))
annot = as_tibble(as.character(seuratobj_stromal_ridgeplot@meta.data$annotation_pooled))
rownames(annot) = colnames(seuratobj_stromal_ridgeplot)
annot$value[which(rownames(annot) %in% fibroblasts)] = "fibroblasts"
rownames(annot) = colnames(seuratobj_stromal_ridgeplot)
annot$value[which(rownames(annot) %in% smooth_muscle)] = "smooth muscle"
rownames(annot) = colnames(seuratobj_stromal_ridgeplot)
annot$value[which(rownames(annot) %in% skeletal_muscle)] = "skeletal muscle"
table(annot$value)
seuratobj_stromal_ridgeplot@meta.data$annotation_pooled = annot$value
table(seuratobj_stromal_ridgeplot@meta.data$annotation_pooled)

df = FetchData(seuratobj_stromal_ridgeplot, vars = c("FCGRT", "annotation_pooled"), slot = "scale.data")
colnames(df) <- c("parameter", "celltype")
df$celltype = factor(df$celltype, levels = rev(c("fibroblasts", "smooth muscle", "skeletal muscle")))
ggplot(df, aes(y=celltype, x=parameter, color=celltype, fill = celltype)) + 
  ggridges::geom_density_ridges(aes(height =..ndensity..), scale = 1.5, alpha = 0.1, size = 0.3) + 
  theme(axis.title = element_blank()) +
  NoLegend() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        # axis.text.y.left = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text = element_text(colour = "grey30"), axis.text.y = element_text(vjust = -0.5),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("fibroblasts" = "#400000",
                                "smooth muscle" = "#bf0000",
                                "skeletal muscle" = "grey40")) + 
  scale_fill_manual(values = c("fibroblasts" = "#400000",
                                "smooth muscle" = "#bf0000",
                                "skeletal muscle" = "#999999")) +
  xlim(NA, 3)

# Ridgeplots for epithelial (intestynal crypt stem cell, type ii pneumocyte and epithelial cell in skin)
Idents(seuratobj) = "annotation_pooled"
seuratobj_compartment = subset(seuratobj, idents = "epithelial")
Idents(seuratobj_compartment) = "organ_tissue"
seuratobj_epithelial_ridgeplot = subset(seuratobj_compartment, idents = c("Skin", "Lung", "Small_Intestine", "Large_Intestine"))
Idents(seuratobj_epithelial_ridgeplot) = "cell_ontology_class"
seuratobj_epithelial_ridgeplot = subset(seuratobj_epithelial_ridgeplot,
                                     idents = c("epithelial cell", "type ii pneumocyte", "intestinal crypt stem cell"))
table(seuratobj_epithelial_ridgeplot@meta.data$cell_ontology_class, seuratobj_epithelial_ridgeplot@meta.data$organ_tissue)

df = FetchData(seuratobj_epithelial_ridgeplot, vars = c("FCGRT", "cell_ontology_class"), slot = "scale.data")
colnames(df) <- c("parameter", "celltype")
df$celltype = factor(df$celltype, levels = rev(c("intestinal crypt stem cell", "type ii pneumocyte", "epithelial cell")))
ggplot(df, aes(y=celltype, x=parameter, color=celltype, fill = celltype)) + 
  ggridges::geom_density_ridges(aes(height =..ndensity..), scale = 1.5, alpha = 0.1, size = 0.3) + 
  theme(axis.title = element_blank()) +
  NoLegend() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        # axis.text.y.left = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text = element_text(colour = "grey30"), axis.text.y = element_text(vjust = -0.5),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("intestinal crypt stem cell" = "#300030",
                                "type ii pneumocyte" = "#680840",
                                "epithelial cell" = "grey40")) + 
  scale_fill_manual(values = c("intestinal crypt stem cell" = "#300030",
                               "type ii pneumocyte" = "#680840",
                               "epithelial cell" = "#999999")) +
  xlim(NA, 3)

# Ridgeplots for endothelial (LSECs) vs non-LSECs
Idents(seuratobj) = "annotation_pooled"
seuratobj_compartment = subset(seuratobj, idents = "endothelial")
Idents(seuratobj_compartment) = "cell_ontology_class"
LSECs = WhichCells(seuratobj_compartment, 
                         idents = c("endothelial cell of hepatic sinusoid"))
non_LSECs = WhichCells(seuratobj_compartment, 
                       idents = c("endothelial cell of hepatic sinusoid"), invert = TRUE)
annot = as_tibble(as.character(seuratobj_compartment@meta.data$cell_ontology_class))
rownames(annot) = colnames(seuratobj_compartment)
annot$value[which(rownames(annot) %in% LSECs)] = "LSECs"
rownames(annot) = colnames(seuratobj_compartment)
table(annot$value)
seuratobj_compartment@meta.data$annotation_pooled = annot$value

df = FetchData(seuratobj_compartment, vars = c("FCGR2B", "annotation_pooled"), slot = "scale.data")
colnames(df) <- c("parameter", "celltype")
# df$celltype = factor(df$celltype, levels = rev(c("intestinal crypt stem cell", "type ii pneumocyte", "epithelial cell")))
ggplot(df, aes(y=celltype, x=parameter, color=celltype, fill = celltype)) + 
  ggridges::geom_density_ridges(aes(height =..ndensity..), scale = 1.5, alpha = 0.1, size = 0.3) + 
  theme(axis.title = element_blank()) +
  NoLegend() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        # axis.text.y.left = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text = element_text(colour = "grey30"), axis.text.y = element_text(vjust = -0.5),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("non_LSECs" = "#999999",
                                "LSECs" = "grey20")) +
  scale_fill_manual(values = c("non_LSECs" = "#999999",
                               "LSECs" = "grey20"))