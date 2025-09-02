source(file.path("D:/Experiments/20230906 Tabula Sapiens and Muris/colors", 'Cell_pop_colors.R'))
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

seuratobj = readRDS("immune_dictionary_SeuratV5.RDS")
seuratobj = readRDS("immune_dictionary_SeuratV5_preproc.RDS")

Idents(seuratobj) = "celltype"
table(seuratobj@meta.data[["celltype"]])
cytokines = unique(seuratobj@meta.data[["sample"]])[!unique(seuratobj@meta.data[["sample"]]) == "PBS"]

seuratobj[["percent.mt"]] <- PercentageFeatureSet(seuratobj, pattern = "^mt-")
VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) & 
  scale_fill_manual(values = cell_colors)

seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
"Tnfaip3" %in% VariableFeatures(seuratobj)
seuratobj <- ScaleData(seuratobj, features = c(VariableFeatures(seuratobj), "TNFAIP3"))
seuratobj <- RunPCA(seuratobj, features = c(VariableFeatures(seuratobj), "TNFAIP3"))
DimPlot(seuratobj, reduction = "pca", label = TRUE) # + NoLegend()
DimHeatmap(seuratobj, dims = 1:15, cells = 200, balanced = TRUE)
DimHeatmap(seuratobj, dims = 16:30, cells = 200, balanced = TRUE) 
ElbowPlot(seuratobj, ndims = 50) # Elbow at +- 30
# seuratobj <- FindNeighbors(seuratobj, dims = 1:30)
# seuratobj <- FindClusters(seuratobj, resolution = 0.5)
seuratobj <- RunUMAP(seuratobj, dims = 1:30, verbose = TRUE)

DimPlot(seuratobj, label = T, label.size = 5, repel = T, reduction = "umap") + 
  NoLegend() + NoAxes() +
  scale_color_manual(values = cell_colors) + ggtitle("Populations in murine LN")

saveRDS(seuratobj, "immune_dictionary_SeuratV5_preproc.RDS")
seuratobj = readRDS("immune_dictionary_SeuratV5_preproc.RDS")

##### Fc gamma receptor analysis ----

fcreceptors = c("Fcgr1", "Fcgr2b", "Fcgr3", "Fcgr4", "Fcer1g", "Fcgrt")

seuratobj = readRDS("immune_dictionary_SeuratV5_preproc.RDS")
gene_data <- FetchData(seuratobj, vars = c(fcreceptors, "celltype", "sample"))

i = "Fcgr1"
data_Fcgr1 <- gene_data %>%
  select(celltype, sample, i) %>%
  group_by(celltype, sample) %>%
  filter(n() >= 3) %>%
  group_by(celltype, sample) %>%
  summarise(mean = mean(Fcgr1, na.rm = TRUE), 
            expression_list = list(Fcgr1),
            .groups = 'drop') %>%
  mutate(fcreceptor = i)
i = "Fcgr2b"
data_Fcgr2b <- gene_data %>%
  select(celltype, sample, i) %>%
  group_by(celltype, sample) %>%
  filter(n() >= 3) %>%
  group_by(celltype, sample) %>%
  summarise(mean = mean(Fcgr2b, na.rm = TRUE), 
            expression_list = list(Fcgr2b),
            .groups = 'drop') %>%
  mutate(fcreceptor = i)
i = "Fcgr3"
data_Fcgr3 <- gene_data %>%
  select(celltype, sample, i) %>%
  group_by(celltype, sample) %>%
  filter(n() >= 3) %>%
  group_by(celltype, sample) %>%
  summarise(mean = mean(Fcgr3, na.rm = TRUE), 
            expression_list = list(Fcgr3),
            .groups = 'drop') %>%
  mutate(fcreceptor = i)
i = "Fcgr4"
data_Fcgr4 <- gene_data %>%
  select(celltype, sample, i) %>%
  group_by(celltype, sample) %>%
  filter(n() >= 3) %>%
  group_by(celltype, sample) %>%
  summarise(mean = mean(Fcgr4, na.rm = TRUE), 
            expression_list = list(Fcgr4),
            .groups = 'drop') %>%
  mutate(fcreceptor = i)
i = "Fcgrt"
data_Fcgrt <- gene_data %>%
  select(celltype, sample, i) %>%
  group_by(celltype, sample) %>%
  filter(n() >= 3) %>%
  group_by(celltype, sample) %>%
  summarise(mean = mean(Fcgrt, na.rm = TRUE), 
            expression_list = list(Fcgrt),
            .groups = 'drop') %>%
  mutate(fcreceptor = i)
i = "Fcer1g"
data_Fcer1g <- gene_data %>%
  select(celltype, sample, i) %>%
  group_by(celltype, sample) %>%
  filter(n() >= 3) %>%
  group_by(celltype, sample) %>%
  summarise(mean = mean(Fcer1g, na.rm = TRUE), 
            expression_list = list(Fcer1g),
            .groups = 'drop') %>%
  mutate(fcreceptor = i)

data = rbind(data_Fcgr1, data_Fcgr2b, data_Fcgr3, data_Fcgr4, data_Fcgrt, data_Fcer1g)

results <- list()
for (fcr in unique(data$fcreceptor)) {

  data_fcr = data %>%
    filter(fcreceptor == fcr)
  
  for (ct in unique(data$celltype)) {
    sub_data <- data_fcr %>%
      filter(celltype == ct)
    
    control_expr <- sub_data$expression_list[sub_data$sample == "PBS"][[1]]
    
    for (i in seq_along(sub_data$sample)) {
      if (sub_data$sample[i] != "PBS") {
        treatment_expr <- sub_data$expression_list[[i]]
        p_value <- wilcox.test(treatment_expr, control_expr, alternative = "two.sided", exact = FALSE)$p.value
        log2FC <- log2(mean(treatment_expr, na.rm = TRUE) / mean(control_expr, na.rm = TRUE))
        
        results[[paste(ct, sub_data$sample[i], fcr, sep = "_")]] <- list(celltype = ct, cytokine = sub_data$sample[i], log2FC = log2FC, p_value = p_value, fcreceptor = fcr)
      }
    }
  }
}

results = bind_rows(results)  %>%
  mutate(direction = case_when(
    log2FC > 0 ~ "positive",
    log2FC < 0 ~ "negative")) %>%
  filter(!is.infinite(log2FC)) %>%
  filter(!is.na(log2FC))

write_csv(results, "fcreceptor_cytokine_influence_immune_dictionary.csv")

results_sign = results %>%
  filter(p_value < 0.05) %>%
  filter(!is.infinite(log2FC)) %>%
  mutate(celltype_receptor = paste(fcreceptor, celltype)) %>%
  filter(abs(log2FC) > 1)

sort(table(results_sign$celltype))


deg_spec_gene = results_sign %>% 
  filter(!celltype %in% c("ILC", "eTAC", "T_cell_gd", "T_cell_CD8", "Treg", "T_cell_gd", "MigDC")) %>%
  select(celltype_receptor, cytokine, log2FC) %>% 
  pivot_wider(names_from = cytokine, values_from = log2FC) %>%
  select(order(colnames(.))) %>%
  relocate(celltype_receptor)
  # select("cytokine", "Neutrophil", "Monocyte", "Macrophage", "cDC1", "cDC2", "MigDC", "pDC", "Treg", "NK_cell")

deg_spec_gene = as.data.frame(deg_spec_gene)
rownames(deg_spec_gene) =  deg_spec_gene$celltype_receptor
deg_spec_gene = deg_spec_gene[,-1]
deg_spec_gene[is.na(deg_spec_gene)] = 0
deg_spec_gene

Heatmap(deg_spec_gene, 
        rect_gp = gpar(col = "white", lwd = 2), 
        name = paste0("Log2\nFold\nchange"),
        row_title	= paste("Fc receptors in Immune Dictionary"),
        cluster_columns = F,
        cluster_rows = F,
        row_dend_width = unit(4, "mm"),
        column_dend_height = unit(4, "mm"),
        show_column_names = TRUE,
        # col = c("white", RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)),
        if (min(deg_spec_gene) < 0) {col  = circlize::colorRamp2(c(min(deg_spec_gene), 0, max(deg_spec_gene)), c("#40a8ed", "white", "#bf301d"))}
        else {col  = circlize::colorRamp2(c(0, max(deg_spec_gene)*0.9), c("white", "brown"))},
        row_names_side = "left", 
        row_dend_side = "right", 
        column_names_side = "top", 
        column_dend_side = "bottom",
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        row_dend_gp = gpar(col = "darkgrey"),
        column_dend_gp = gpar(col = "darkgrey"))