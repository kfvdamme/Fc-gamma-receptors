devtools::load_all("D:/Git_repos/FlowSOM")
library(dplyr)
library(FlowSOM)
library(flowCore)
setwd("E:/Cytometry/KarelVanDamme/20241024 human FcyR spectral/")
source("helperFunctions.R")

wspFiles <- list.files(pattern =  ".wsp")
channels <- c("FcyRI" = "PE-CF594-A", # CD64
              "FcyRIIA" = "FITC-A", # CD32a
              "FcyRIIB" = "PE-Cy7-A", # CD32b
              "FcyRIII" = "BV570-A", # CD16
              "FcRn" = "AF647-A", #Synt-002
              "CD89" = "BUV805-A") 
#hardcoded because it is not stored in the FMOs
markers_of_interest <- names(channels)

files <- c(list.files(pattern = "FS1.*fcs", recursive = TRUE),
           list.files(pattern = "FMO.*fcs", recursive = TRUE))
file_to_wsp <- c("WLSM" = "Spectral_axis_IC")

wspFiles_match <- rep(NA, length(files))
names(wspFiles_match) <- files
for(file_type in names(file_to_wsp)){
  wspFiles_match[grep(file_type, files)] <- grep(file_to_wsp[file_type], 
                                                 wspFiles,
                                                 value = TRUE)
}


# Parse workspaces ----

medians <- list()
counts <- list()

for(file in files){
  print(file)
  tryCatch({
    data <- GetFlowJoLabels(file, 
                            wspFile = wspFiles_match[file], 
                            path = dirname(file),
                            group = "All Samples",
                            getData = TRUE)
    manual <- as.character(data$manual)
    manual <- gsub("DC1s", "DCs", 
                   gsub("DC2s", "DCs", manual))
    
    medians[[file]] <- sapply(channels, function(channel){
      tapply(exprs(data$flowFrame)[, channel], manual, median)
    })
  }, error = function(e) {print(paste("Failed for", file))})
  counts[[file]] <- table(manual)
  
}
saveRDS(medians, "medians_hFcyR in human blood - extracellular.RDS")
saveRDS(counts, "counts_hFcyR in human blood - extracellular.RDS")

# Or load previously parsed results ----

medians <- readRDS("medians_hFcyR in human blood - extracellular.RDS")
counts <- readRDS("counts_hFcyR in human blood - extracellular.RDS")

sort(unique(unlist(lapply(medians, rownames))))

clean_names <- function(x){ 
  x
}

medians <- lapply(medians, function(x){rownames(x) <- clean_names(rownames(x)); x})
counts <- lapply(counts, function(x){names(x) <- clean_names(names(x)); x})

# Plot ----

## Merge samples ----
files <- names(medians)
data <- data.frame(File = files,
                   Stain = gsub("([^ ]*) .*", "\\1", files),
                   Tissue = "Blood",
                   Sample_ID = gsub("... (FS1_)*(.*) WLSM.fcs", "\\2", files))

celltypes <- c("Neutrophils", "Eosinophils", "Basophils", 
               "CD14++ CD16-", "CD14+ CD16+", "CD14- CD16+", 
               "DCs", "pDCs",
               "B cells", "T cells", "NK cells")  
samples <- c("35", "36", "37", "38", "39", "40", 
             "FMO CD16", "FMO CD89",
             "FMO CD32b", "FMO CD32a", 
             "FMO CD64", "FMO FcRn")
tissue <- c("Blood")
col_description <- expand.grid(Sample = samples,
                               Tissue = tissue,
                               Marker = markers_of_interest)
col_description$id <- paste(col_description$Marker, col_description$Tissue, col_description$Sample,
                            sep = "_")

results <- matrix(NA, 
                  nrow = length(celltypes), 
                  ncol = nrow(col_description), 
                  dimnames = list(celltypes,col_description$id))
result_counts <- matrix(NA, nrow = length(celltypes), ncol = nrow(col_description), 
                        dimnames = list(celltypes,col_description$id))

for(sample in samples){
  for(tissue  in unique(data$Tissue)){
    selected_files <- data[data$Sample_ID == sample & data$Tissue == tissue, "File"]
    for(file in selected_files){
      selection <- medians[[file]]
      selected_celltypes <- celltypes[which(celltypes %in% rownames(selection))]
      results[selected_celltypes, 
              paste(colnames(selection), tissue, sample, sep = "_")] <- 
        selection[selected_celltypes, ]
      result_counts[selected_celltypes, 
                    paste(colnames(selection), tissue, sample, sep = "_")] <- 
        counts[[file]][selected_celltypes]
    }
  }
}



# Subtract FMOs

fmos <- grep("FMO", colnames(results), value = TRUE)
fmos <- c(grep("FcyRI_.*CD64", fmos, value = TRUE),
          grep("FcyRIIA_.*CD32a", fmos, value = TRUE),
          grep("FcyRIIB_.*CD32b", fmos, value = TRUE),
          grep("FcyRIII_.*CD16", fmos, value = TRUE),
          grep("FcRn_.*FcRn", fmos, value = TRUE),
          grep("CD89_.*CD89", fmos, value = TRUE))


#results["CD14+CD16+ mono's", fmos] <- results["CD14++CD16- mono's", fmos]
#results["CD14-CD16++ mono's", fmos] <- results["CD14++CD16- mono's", fmos]

results_diff <- results
for(fmo in fmos){
  cols_to_adapt <- grep(gsub("_FMO.*", "", fmo), colnames(results), value = TRUE)
  results_diff[,cols_to_adapt] <- results[,cols_to_adapt] - results[,rep(fmo, length(cols_to_adapt))]
}

results_diff

results_diff[results_diff < 0] <- 0

cols_to_plot <- grep("FMO", colnames(results_diff), invert = TRUE)
col_description <- col_description[cols_to_plot,]
to_plot <- results_diff[,cols_to_plot] #log10(result_counts)

max_value_per_marker <- sapply(markers_of_interest, function(marker){ 
  cols <- grep(paste0(marker,"_"), colnames(to_plot), value = TRUE)
  max(to_plot[, cols])
})

marker_per_col <- gsub("_.*", "", colnames(to_plot))
max_to_norm <- max_value_per_marker[marker_per_col]

to_plot <- t(apply(to_plot, 1, function(x) x / max_to_norm))
to_plot <- as.data.frame(to_plot) %>%
  tibble::rownames_to_column("id")

writexl::write_xlsx(list(medians = as.data.frame(results) %>% tibble::rownames_to_column("celltype"), 
                         fmo_corrected = as.data.frame(results_diff) %>% tibble::rownames_to_column("celltype"), 
                         normalized = to_plot,
                         counts = as.data.frame(result_counts) %>% tibble::rownames_to_column("celltype")),
                    "hFcyR in human.xlsx")

column_info <- data.frame(id = colnames(to_plot),
                          group = c("Marker", paste(as.character(col_description$Marker),
                                                    as.character(col_description$Tissue))),
                          name = c("", as.character(col_description$Sample)),
                          geom = c("text", rep("circle", nrow(col_description))),
                          palette = c(NA, as.character(col_description$Marker)),
                          legend = TRUE,
                          hjust = 1)
#column_info$name <- gsub("B FMO", "FMO", column_info$name)

column_groups <- data.frame(group = unique(column_info$group),
                            level1 = gsub(" .*", "", unique(column_info$group)),
                            level2 = gsub(".* ", "", unique(column_info$group)),
                            palette = gsub(" .*", "", unique(column_info$group)))
column_groups[1, "level2"] <- "Tissue"

p <- funkyheatmap::funky_heatmap(to_plot, 
                                 palettes = list("FcyRI" = colorRampPalette(c("white","#e00480"))(100),
                                                 "FcyRIIA" = colorRampPalette(c("white","#22a395"))(100),
                                                 "FcyRIIB" = colorRampPalette(c("white","#ff9a00"))(100),
                                                 "FcyRIII" =  colorRampPalette(c("white","#006131"))(100),
                                                 "FcRn" =  colorRampPalette(c("white","#ff5c27"))(100),
                                                 "CD89" =  colorRampPalette(c("white","darkblue"))(100)),
                                 column_info = column_info,
                                 column_groups = column_groups,
                                 scale_column = FALSE,
                                 col_annot_offset = 2,
                                 col_annot_angle = 90)
p

ggsave(filename = "hFcyR in human_heatmap.pdf",
       plot = p,
       width = 20, height = 5)


# summary plot with average over replicates ----

to_plot_summarised <- sapply(markers_of_interest, function(marker){ 
  cols <- grep(paste0(marker,"_"), colnames(to_plot), value = TRUE)
  m <- apply(to_plot[, cols], 1, mean)
  names(m) <- to_plot$id
  m
})
to_plot_summarised <- as.data.frame(to_plot_summarised) %>%
  tibble::rownames_to_column("id")

column_info <- data.frame(id = c(colnames(to_plot_summarised)),
                          group = c("Marker", paste(colnames(to_plot_summarised)[-1],
                                                    "Blood")),
                          name = c(""),
                          geom = c("text", rep("circle", ncol(to_plot_summarised[,-1]))),
                          palette = c(NA, colnames(to_plot_summarised)[-1]),
                          legend = TRUE,
                          hjust = 1)

column_groups <- data.frame(group = unique(column_info$group),
                            level1 = gsub(" .*", "", unique(column_info$group)),
                            level2 = gsub(".* ", "", unique(column_info$group)),
                            palette = gsub(" .*", "", unique(column_info$group)))
column_groups[1, "level2"] <- "Tissue"

p <- funkyheatmap::funky_heatmap(to_plot_summarised, 
                                 palettes = list("FcyRI" = colorRampPalette(c("white","#e00480"))(100),
                                                 "FcyRIIA" = colorRampPalette(c("white","#22a395"))(100),
                                                 "FcyRIIB" = colorRampPalette(c("white","#ff9a00"))(100),
                                                 "FcyRIII" =  colorRampPalette(c("white","#006131"))(100),
                                                 "FcRn" =  colorRampPalette(c("white","#ff5c27"))(100)),
                                 column_info = column_info,
                                 column_groups = column_groups,
                                 scale_column = FALSE,
                                 col_annot_offset = 2,
                                 col_annot_angle = 90)


ggsave(filename = "hFcyR in human_heatmap_summarised.pdf",
       plot = p,
       width = 5, height = 5)
