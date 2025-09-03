library(tidyverse)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)

yc_targets = readxl::read_xls("jem_20222052_tables3.xls", skip = 1)
yc_targets = as_tibble(yc_targets)
yc_targets

yc_targets_long = yc_targets %>%
  pivot_longer(
    cols = 3:255,
    names_to = "population_cytokine_replicate") %>%
  mutate(replicate = stringr::str_sub(population_cytokine_replicate, -1, -1)) %>%
  mutate(cytokine = sub('.*IL', '', population_cytokine_replicate)) %>%
  mutate(cytokine = sub("\\..*", "", cytokine)) %>%
  mutate(cytokine = paste0("IL", cytokine)) %>%
  mutate(cytokine = ifelse(
    cytokine != "IL2" & cytokine != "IL4" & cytokine != "IL7" & cytokine != "IL9" & cytokine != "IL15" & cytokine != "IL21",
    "PBS",
    cytokine)) %>%
  mutate(organ = ifelse(
    grepl("PC", population_cytokine_replicate, fixed = TRUE),
    "peritoneal Cavity",
    "spleen"
  )) %>%
  mutate(population = sub("\\..*", "", population_cytokine_replicate)) %>%
  select(GeneSymbol, population, cytokine, organ, `Cluster#`, replicate, value) %>%
  group_by(GeneSymbol, population, cytokine, organ, `Cluster#`) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>%
  filter(!is.na(value))

yc_targets_long$cytokine = factor(yc_targets_long$cytokine, 
                               levels = c("IL2", "IL4", "IL7", "IL9", "IL15", "IL21", "PBS"))
yc_targets_long$value = format(yc_targets_long$value, scientific = FALSE)

yc_targets_long %>%
  filter(cytokine != "PBS") %>%
  ggplot( aes(x=value, color=cytokine, fill=cytokine)) +
  geom_histogram(alpha=0.6, binwidth = 0.5) +
  viridis::scale_fill_viridis(discrete=TRUE) +
  viridis::scale_color_viridis(discrete=TRUE) +
  # ggthemes::theme_hc() +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  xlab("Log2-fold changes") +
  ylab("Number of genes") +
  facet_wrap(~cytokine, nrow = 2)

fcreceptors = c("Fcgr1", "Fcgr2b", "Fcgr3", "Fcgr4", "Fcer1g", "Fcgrt", "Trim21"
                , "Fcer1a", "Fcer2a", "Fcamr", "Faim3", "Pigr",
               "Tnf"
               )

yc_fcrs = yc_targets_long %>%
  filter(GeneSymbol %in% fcreceptors)

yc_fcrs_pop_wide = yc_fcrs %>%
  filter(cytokine != "PBS") %>%
  filter(population == "NK") %>%
  pivot_wider(names_from = "cytokine", values_from = "value") %>%
  ungroup() %>%
  select(GeneSymbol, IL2, IL4, IL7, IL15, IL21)
yc_fcrs_pop_wide = as.data.frame(yc_fcrs_pop_wide)
rownames(yc_fcrs_pop_wide) =  yc_fcrs_pop_wide$GeneSymbol
yc_fcrs_pop_wide = yc_fcrs_pop_wide[,-1]
# cytosig_wide[is.na(cytosig_wide)] = 0

mat = yc_fcrs_pop_wide %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

Heatmap(mat,
        rect_gp = gpar(col = "white", lwd = 2),
        name = "LFC",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_dend_width = unit(8, "mm"),
        column_dend_height = unit(5, "mm"),
        show_column_names = TRUE,
        col  = circlize::colorRamp2(c(min(mat), 0, max(mat)),
                                    c("darkblue", "white", "brown")),
        row_names_side = "left",
        row_dend_side = "right",
        column_title = "yc cytokine activity\non NK cells",
        column_names_side = "top",
        column_dend_side = "bottom",
        column_names_rot = 90,
        row_dend_gp = gpar(col = "darkgrey"),
        column_dend_gp = gpar(col = "darkgrey"))

i= "Fcgr3"
plotlist <- lapply(
  names(table(yc_fcrs$GeneSymbol)),
  function(i) {
    yc_fcrs_pop = yc_fcrs %>%
      filter(cytokine != "PBS") %>%
      filter(population == "NK") %>%
      pivot_wider(names_from = "cytokine", values_from = "value") %>%
      ungroup() %>%
      select(GeneSymbol, IL2, IL4, IL7, IL15, IL21) # no IL9 with effect on NK cells
    yc_fcrs_pop_transp = t(yc_fcrs_pop[-1])
    colnames(yc_fcrs_pop_transp) = pull(yc_fcrs_pop[,1])
    rn = rownames(yc_fcrs_pop_transp)
    yc_fcrs_pop_transp = as.data.frame(apply(yc_fcrs_pop_transp, 2, as.numeric))
    rownames(yc_fcrs_pop_transp) = rn
    yc_fcrs_pop_transp$x = rownames(yc_fcrs_pop_transp)
    yc_fcrs_pop_transp$y = as.numeric(scales::number(as.numeric(yc_fcrs_pop_transp[,i]), accuracy = 0.1))
    yc_fcrs_pop_transp$x = factor(yc_fcrs_pop_transp$x,
                                     levels = rev(c("IL2", "IL4", "IL7", "IL9", "IL15", "IL21")[c("IL2", "IL4", "IL7", "IL9", "IL15", "IL21") %in% yc_fcrs_pop_transp$x]))
    plot = yc_fcrs_pop_transp %>%
      ggplot(aes(x = x, y = y)) +
      geom_segment(aes(x = x,
                       xend = x,
                       y = 0, yend = y,
                       color = y), lwd = 1) +
                       # color = "gray60", lwd = 1) +
      geom_point(aes(color = y), size = 9, pch = 16) +
      geom_text(aes(label = y), color = "white", size = 3) +
      coord_flip() +
      ggthemes::theme_hc() +  
      ggtitle(paste0(i)) + 
      ylim(min(yc_fcrs_pop_transp[,1:3]*1.2), max(yc_fcrs_pop_transp[,1:4])*1.2) +
      theme(axis.ticks.y = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.title = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(colour = "grey60"),
            axis.ticks = element_line(colour = "grey60"),
            # text=element_text(family="sans"),
            # axis.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5)) +
      # labs(x = "NK cells", y = "log2-fold change\nafter yc cytokine") + 
      labs(x = NULL, y = NULL) + 
      scale_color_gradient2(limits = c(min(yc_fcrs_pop_transp[,1:3]*1.05), max(yc_fcrs_pop_transp[,1:4])*1.05),
                            low = "darkblue",
                            mid = "#fcfcb6",
                            high = "brown")
    print(plot)
  }
)
allplots = ggpubr::ggarrange(plotlist=plotlist, ncol = 4, nrow = 1, align = "hv")
allplots = ggpubr::annotate_figure(allplots, 
                                   top = text_grob("", color = "white", face = "bold", size = 4),
                                   # left = textGrob("transcription in NK cells", rot = 90, vjust = 1, gp = gpar(size = 10)),
                                   bottom = textGrob("log2 fold change in NK cells after cytokine administration", gp = gpar(size = 10)),
                                   fig.lab = "Reanalyzed from Baysoy et al., JEM 2023",
                                   fig.lab.face = "italic",
                                   fig.lab.pos = "top.right"
                                   )
                                  
allplots