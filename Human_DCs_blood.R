##### Data integration and processing ####
library(Seurat)
library(tidyverse)

# Downloaded data from: https://singlecell.broadinstitute.org/single_cell/study/SCP43/atlas-of-human-blood-dendritic-cells-and-monocytes
expr_matrix_tpm = read.delim("expression_matrix_tpm.txt", header = TRUE)
rownames(expr_matrix_tpm) = expr_matrix_tpm[,1]
expr_matrix_tpm_named = expr_matrix_tpm[,-1]

seuratobj<-CreateSeuratObject(counts = expr_matrix_tpm_named, min.cells = 0, min.features = 0, project = "DC_Atlas")

metadata <- read.delim("metadata.txt")
clusters = metadata$CellType[-1]

Idents(seuratobj) = clusters
seuratobj@meta.data$celltype_novel = Idents(seuratobj)
table(Idents(seuratobj))

fcreceptors = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", 
                "FCER1G", "FCGRT", "TRIM21", "FCRL3", "FCRL4", "FCRL5", 
                "FCER1A", "FCER2", "FCAR", "FCAMR", "FCMR", "PIGR")

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

Idents(seuratobj) = "orig.ident"
seuratobj = RenameIdents(object = seuratobj,
             "CD141" = "cDC1s",
             "CD1C" = "cDC2s",
             "DoubleNeg" = "CD1C- CD141-",
             "Mono" = "Monos",
             "pDC" = "pDCs")
seuratobj@meta.data$celltype = Idents(seuratobj)
             
Idents(seuratobj) = "celltype_novel"
celltype = Idents(seuratobj)
celltype = factor(celltype, levels = c("DC1", "DC2", "DC3", "DC4", "DC5", "DC6", 
                                          "Mono1", "Mono2", "Mono3", "Mono4"))
seuratobj@meta.data$celltype_novel = celltype

seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst"#, nfeatures = 2000
)
head(VariableFeatures(seuratobj), 10)
fcreceptors[fcreceptors %in% VariableFeatures(seuratobj)]
seuratobj <- ScaleData(seuratobj, features = #rownames(seuratobj)
                         c(VariableFeatures(seuratobj), fcreceptors)
)
fcreceptors[fcreceptors %in% VariableFeatures(seuratobj)]

saveRDS(seuratobj, "dc_atlas.rds")
seuratobj = readRDS("dc_atlas.rds")

fcreceptors = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", 
                "FCER1G", "FCGRT", "TRIM21", "FCRL3", "FCRL4", "FCRL5", 
                "FCER1A", "FCER2", "FCAR", "FCAMR", #"FCMR", 
                "PIGR")

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

##### Data visualisation ####
library(Seurat)
library(tidyverse)
library(ggpubr)
seuratobj = readRDS("dc_atlas.rds")

Idents(seuratobj) = "celltype"
table(Idents(seuratobj))

Idents(seuratobj) = factor(Idents(seuratobj), levels = rev(c("cDC1s", "cDC2s", "CD1C- CD141-", "pDCs", "Monos")))
Idents(seuratobj) = factor(Idents(seuratobj), levels = c("cDC1s", "cDC2s", "CD1C- CD141-", "pDCs", "Monos"))

# i = "PIGR"
plotList <- lapply(
  fcreceptors,
  function(i) {
    df = FetchData(subset(seuratobj, idents = c("cDC1s", "cDC2s", "pDCs", "Monos")), vars = c(i, "celltype"), slot = "scale.data")
    colnames(df) <- c("parameter", "celltype")
    df$celltype = factor(df$celltype, levels = rev(c("cDC1s", "cDC2s", "pDCs", "Monos")))
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
allplots <- ggpubr::ggarrange(plotlist=plotList, ncol = 8, nrow = 3)
allplots

# Heatmap

DoHeatmap(subset(seuratobj, idents = rev(c("cDC1s", "cDC2s", "pDCs", "Monos")), downsample = 50), 
          slot = "counts", features = fcreceptors, size = 3.5, disp.min = -0.5, disp.max = 3.5, draw.lines = TRUE, lines.width = 2) + 
  NoLegend() +
  theme(text = element_text(size = 12)) +
  scale_fill_binned(type = "viridis", na.value = "white")

# Dotplot

metadata = as_tibble(seuratobj@meta.data[["celltype"]])
metadata$celltype_novel = seuratobj@meta.data[["celltype_novel"]]
combined = metadata %>% 
  mutate(combined = case_when(celltype_novel == "Mono1" ~ "Mono1",
                              celltype_novel == "Mono2" ~ "Mono2",
                              .default = value)) %>%
  pull(combined)

seuratobj@meta.data$celltype_combined = combined

Idents(seuratobj) = "celltype_combined"
Idents(seuratobj) = factor(Idents(seuratobj),
                           levels = c("cDC1s", "cDC2s", "CD1C- CD141-", "pDCs", "Mono1", "Mono2", "Monos"))

table(Idents(seuratobj))

DotPlot(seuratobj, 
        idents = c("cDC1s", "cDC2s", "CD1C- CD141-", "pDCs", "Mono1", "Mono2"),
        features = c(rev(fcreceptors)), 
        # group.by = "celltype", 
        col.min = -1,
        col.max = 2,
        dot.min = 0,
        dot.scale = 5,
        cols= "RdYlBu",
        scale.by = "radius",
        # idents = cells_more_than_250, 
        scale = T) + 
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


##### Clustering heatmap ####
Idents(seuratobj) = "celltype_combined"
table(Idents(seuratobj))

# Idents(seuratobj) = droplevels(Idents(seuratobj))
table(Idents(seuratobj))

cluster.markers <- lapply(levels(Idents(seuratobj)[[1]]),
                          function(x)FindMarkers(seuratobj,ident.1 = x,min.pct = 0.40)) 

sapply(0:(length(cluster.markers)-1),function(x) {
  cluster.markers[[x+1]]$gene <<- rownames(cluster.markers[[x+1]])
  cluster.markers[[x+1]]$cluster <<- x})
as_tibble(do.call(rbind,cluster.markers)) %>% arrange(p_val_adj) -> cluster.markers
cluster.markers
saveRDS(cluster.markers, "cluster_markers_dc_celltype_annot.rds")

top10 = cluster.markers %>%
  filter(cluster != 6) %>%
  group_by(cluster) %>%
  slice_max(n = 8, order_by = avg_log2FC)

# Idents(seuratobj) = factor(Idents(seuratobj),
                           # levels = c("cDC1s", "cDC2s", "CD1C- CD141-", "pDCs", "Mono1", "Mono2", "Monos"))

DoHeatmap(subset(seuratobj, downsample=30,
                 idents = c("cDC1s", "cDC2s", "CD1C- CD141-", "pDCs", "Mono1", "Mono2")
                 ), 
          features = top10$gene, assay = "RNA", size = 3,
          slot = "scale.data", label = TRUE, disp.min = -1, draw.lines = TRUE, lines.width = 1) + 
  NoLegend() +
  theme(text = element_text(size = 8)) +
  scale_fill_binned(type = "viridis", na.value = "white")

##### Ridgeplots for new classification of DCs #####
library(Seurat)
library(tidyverse)
library(ggpubr)
library(viridis)
seuratobj = readRDS("dc_atlas.rds")
source(file.path("D:/Experiments/20230906 Tabula Sapiens and Muris/colors", 'Cell_pop_colors.R'))


fcreceptors = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", 
                "FCER1G", "FCGRT" #, "TRIM21", "FCRL3", "FCRL4", "FCRL5", 
                #"FCER1A", "FCER2", "FCAR", "FCAMR", #"FCMR", 
                #"PIGR"
                )

Idents(seuratobj) = "celltype_novel"

DC_atlas_colors = c(
  "DC1" = "#cc6a4e",
  "DC2" = "#aebedf",
  "DC3" = "#d4dded",
  "DC4" = "#b384b8",
  "DC5" = "#efb77f",
  "DC6" = "#e07e73",
  "Mono1" = "#8ec78f",
  "Mono2" = "#9fd3f4",
  "Mono3" = "#fcb4c0",
  "Mono4" = "#aa8cc4"
)

plotList <- lapply(
  fcreceptors,
  function(i) {
    df = FetchData(seuratobj, vars = c(i, "celltype_novel"), slot = "scale.data")
    colnames(df) <- c("parameter", "celltype")
    df$celltype = factor(df$celltype, levels = rev(c("DC1", "DC2", "DC3", "DC4", "DC5", "DC6", "Mono1", "Mono2", "Mono3", "Mono4")))
    plot = ggplot(df, aes(y=celltype, x=parameter, fill=celltype)) + 
      ggridges::geom_density_ridges(aes(height =..ndensity..), scale = 1.5, alpha = 0.5,
                                    color = "grey70", size = 0.3) + 
      theme(axis.title = element_blank()) +
      NoLegend() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.text.y.left = element_blank(),
            axis.ticks = element_blank(), axis.text.x = element_text(size = 8), axis.text = element_text(colour = "grey30"), axis.text.y = element_text(vjust = -0.5),
            plot.title = element_text(hjust = 0.5)) + 
      scale_fill_manual(values = DC_atlas_colors) +
      ggtitle(i)
    print(plot)
  }
)
allplots <- ggpubr::ggarrange(plotlist=plotList, ncol = 8, nrow = 1)
allplots

FeaturePlot(seuratobj, "FCGR1A", cols = c("grey", viridis(10)), 
            order = T,
            max.cutoff = "q95"
            ) + NoAxes() + theme(title = element_text(face = "plain"))

VlnPlot(seuratobj, "FCGR3B") + NoLegend() +     
  scale_fill_manual(values = cell_colors)

