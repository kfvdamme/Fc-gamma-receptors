# Downloaded  seurat_v5_object.rds md5:100e6686ff4003b19d9807d2d582c169 from https://zenodo.org/records/15744698 on 23/08/2025
source(file.path("D:/Experiments/20230906 Tabula Sapiens and Muris/colors", 'Cell_pop_colors.R'))
library(Seurat)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)

seuratobj = readRDS("seurat_v5_object.rds")

Idents(seuratobj) = "minor_label"
DimPlot(seuratobj, reduction = "umap.rpca")
unique(seuratobj@meta.data$stim)


# Exploration ----

fcreceptors = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", 
                "FCER1G", "FCGRT")

FeaturePlot(seuratobj, features = fcreceptors, reduction = "umap.rpca", 
            order = F,
            max.cutoff = "q90") & NoAxes() & NoLegend() &
  scale_color_gradientn(colors = c("grey90", rev(magma(10))))

DotPlot(seuratobj, features = rev(fcreceptors), 
        cols= "RdYlBu",
        scale.by = "radius",
        scale = T) + 
  theme(axis.text.x=element_text(angle=90,vjust=0.5)) + 
  theme(axis.text=element_text(size=9)) +
  theme(axis.title=element_blank()) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.title=element_text(size=10)) +
  scale_size(breaks = c(0, 10, 20, 30)) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0)
  ) + 
  scale_y_discrete(position = "right") +
  coord_flip()

fcreceptors %in% VariableFeatures(seuratobj)

# Differential expression analysis ----

gene_data <- as_tibble(FetchData(seuratobj, vars = c(fcreceptors, "minor_label", "donor", "stim"))) %>%
  rename(celltype = minor_label, sample = donor)

summarise_receptor <- function(df, receptor) {
  df %>%
    select(celltype, sample, stim, !!sym(receptor)) %>%
    group_by(celltype, sample, stim) %>%
    filter(n() >= 3) %>%
    summarise(mean = mean(.data[[receptor]], na.rm = TRUE),
              expression_list = list(.data[[receptor]]),
              .groups = "drop") %>%
    mutate(fcreceptor = receptor)
}

data <- map_dfr(fcreceptors, ~ summarise_receptor(gene_data, .x))

results <- list()

for (fcr in unique(data$fcreceptor)) {
  
  data_fcr <- data %>%
    filter(fcreceptor == fcr)
  
  for (ct in unique(data$celltype)) {
    sub_data <- data_fcr %>%
      filter(celltype == ct)
    
    if (!any(sub_data$stim == "Baseline")) next
    
    control_expr <- sub_data$expression_list[sub_data$stim == "Baseline"][[1]]
    
    for (i in seq_along(sub_data$stim)) {
      if (sub_data$stim[i] != "Baseline") {
        treatment_expr <- sub_data$expression_list[[i]]
        
        if (length(treatment_expr) > 0 && length(control_expr) > 0) {
          mean_expr <- mean(c(treatment_expr, control_expr), na.rm = TRUE)
          p_value <- wilcox.test(treatment_expr, control_expr, 
                                 alternative = "two.sided", exact = FALSE)$p.value
          log2FC <- log2(mean(treatment_expr, na.rm = TRUE) / mean(control_expr, na.rm = TRUE))
          
          results[[paste(ct, sub_data$stim[i], fcr, sep = "_")]] <- list(
            celltype = ct, 
            stim = sub_data$stim[i], 
            log2FC = log2FC, 
            p_value = p_value, 
            fcreceptor = fcr,
            mean_expr = mean_expr
          )
        }
      }
    }
  }
}

results <- bind_rows(results) %>%
  group_by(fcreceptor) %>%
  filter(mean_expr > 0.1) %>% # filter results where expr is at least 0.1 in stimulated condition!
  mutate(
    FDR = p.adjust(p_value, method = "fdr") # fdr was calculated within each Fc receptor
  ) %>%
  ungroup() %>%
  filter(!is.infinite(log2FC), !is.na(log2FC))

write_csv(results, "fcreceptor_cytokine_influence_PBMCs.csv")

# Filter significant results and visualize -----

results_sign = results %>%
  filter(FDR < 0.05) %>%
  filter(abs(log2FC) > 1) %>%
  mutate(celltype_receptor = paste(fcreceptor, celltype)) 

sort(table(results_sign$celltype))

deg_spec_gene = results_sign %>% 
  filter(!fcreceptor %in% c("FCGR2C")) %>%
  select(celltype_receptor, stim, log2FC) %>% 
  pivot_wider(names_from = stim, values_from = log2FC) %>%
  select(order(colnames(.))) %>%
  select(celltype_receptor, LPS, IFNa, IL1B, IL2, IL4, IL10, TGFb1, TNFa, Cytostim, CD3, BCS)

deg_spec_gene = as.data.frame(deg_spec_gene)
rownames(deg_spec_gene) =  deg_spec_gene$celltype_receptor
deg_spec_gene = deg_spec_gene[,-1]
deg_spec_gene[is.na(deg_spec_gene)] = 0
deg_spec_gene

Heatmap(deg_spec_gene, 
        rect_gp = gpar(col = "white", lwd = 2), 
        name = paste0("Log2\nFold\nchange"),
        row_title	= paste("In vitro simulated PBMCs"),
        cluster_columns = F,
        cluster_rows = F,
        row_dend_width = unit(4, "mm"),
        column_dend_height = unit(4, "mm"),
        show_column_names = TRUE,
        # col = c("white", RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)),
        if (min(deg_spec_gene) < 0) {col  = circlize::colorRamp2(c(0.6*min(deg_spec_gene), 0, 0.8*max(deg_spec_gene)), c("#4a50a8", "white", "brown"))} # note that I trimmed coloration due to very strong downregulation of some genes
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
