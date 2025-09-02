library(ComplexHeatmap)
library(tidyverse)
fcreceptors_cytokines = readRDS("fcreceptors_cytokines_immune_dictionary.RDS")
fcreceptors = names(table(fcreceptors_cytokines$variable))

# Downloaded csv files for each population for all Fc receptors at immune-dictionary.org
# Reformatted and renamed csv files

populations = c("neutro", "mono", "mac", "cdc1", "cdc2", "migdc", "pdc", "b", "cd8", "cd4", "treg", "gd", "nk") # 13 populations

# i = "neutro"
fcreceptors_cytokines = lapply( # This reads the csv file of each immune population, and summarizes the mean expression for each Fc receptor
  populations,
  function(i) {
    csv = read.csv(paste0(i,".csv"), header = T, sep = ";", dec = ".")
    csv_sum = csv %>%
      group_by(celltype, variable, sample) %>%
      summarise(avg = mean(value)) %>%
      pivot_wider(names_from = sample,
                  values_from = avg)
  })

fcreceptors_cytokines = fcreceptors_cytokines %>% bind_rows()

## Start with PBS and average columns, and delta for each cytokine

i = "Fcgr1"
plotList <- lapply(
  fcreceptors,
  function(i) {
    fcreceptors_cytokines_spec = 
      fcreceptors_cytokines %>%
      filter(!celltype %in% c("T_cell_CD4", "T_cell_CD8", "T_cell_gd", "Treg", "MigDC", "pDC")) %>%
      # fcreceptors_cytokines_no_t_cells %>%
      filter(variable == i) %>%
      ungroup() %>%
      select(-c(variable))
    
    fcreceptors_cytokines_spec[is.na(fcreceptors_cytokines_spec)] = 0 # IL1a in cDC2s is missing!
    
    average_spec = 
      fcreceptors_cytokines_spec %>%
      select(-c(celltype, PBS)) %>%
      rowMeans()
    
    PBS_spec_df = cbind(as.numeric(replicate(87, fcreceptors_cytokines_spec$PBS)))
    
    differences_spec = fcreceptors_cytokines_spec %>% select(-c(celltype)) - PBS_spec_df
    
    differences_spec[c("celltype", "PBS", "average")] = c(fcreceptors_cytokines_spec$celltype, fcreceptors_cytokines_spec$PBS, average_spec)
    
    differences_spec = 
      differences_spec %>%
      relocate(celltype, PBS, average)
    
    differences_spec = as.data.frame(differences_spec)
    rownames(differences_spec) =  differences_spec$celltype
    differences_spec = differences_spec[,-1]
    differences_spec = mutate_all(differences_spec, function(x) as.numeric(as.character(x)))
    
    barplot_PBS = row_anno_barplot(matrix(differences_spec$PBS))
    barplot_avg = row_anno_barplot(matrix(differences_spec$average))
    row_ha = HeatmapAnnotation(PBS = barplot_PBS, μ = barplot_avg, which = "row", annotation_name_side	= "top", annotation_name_rot = 0, gap = unit(1, "mm"))
    plot = grid.grabExpr(draw(
      Heatmap(differences_spec[-c(1,2)], 
              left_annotation = row_ha,
              rect_gp = gpar(col = "white", lwd = 2), 
              name = paste0("Δ", i),
              row_title	= i,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              row_dend_width = unit(8, "mm"),
              column_dend_height = unit(8, "mm"),
              show_column_names = TRUE,
              # col = c("white", RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)),
              col  = circlize::colorRamp2(c(min(differences_spec)*0.9, 0, max(differences_spec)*0.9), c("blue", "white", "red")),
              row_names_side = "left", 
              row_dend_side = "right", 
              column_names_side = "top", 
              column_dend_side = "bottom",
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 9),
              row_names_gp = gpar(fontsize = 9),
              row_dend_gp = gpar(col = "darkgrey"),
              column_dend_gp = gpar(col = "darkgrey"))
    ))
    print(plot)
  }
)
allplots <- ggpubr::ggarrange(plotlist=plotList, ncol = 1, nrow = 5, align = "hv")
allplots # Saved 17 x 10