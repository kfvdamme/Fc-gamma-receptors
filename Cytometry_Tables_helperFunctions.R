
match_wsp_to_file <- function(files, wspFiles, file_to_wsp){
  wspFiles_match <- rep(NA, length(files))
  names(wspFiles_match) <- files
  for(file_type in names(file_to_wsp)){
    wspFiles_match[grep(file_type, files)] <- grep(file_to_wsp[file_type], wspFiles, value = TRUE)
  }
  return(wspFiles_match)
}

parse_workspaces <- function(files, wspFiles, channels, analysis){
  medians <- list()
  counts <- list()
  
  for(file in files){ 
    print(file)
    tryCatch({
      data <- GetFlowJoLabels(gsub("B09 FMO CD32b WLSM.fcs", "B09 FMO18 WLSM.fcs", 
                                   gsub("E09 FMO FcRn WLSM.fcs", "E09 FMO21 WLSM.fcs", file)), 
                              wspFile = wspFiles[file], 
                              path = dirname(file),
                              group = "All Samples",
                              getData = TRUE)
      manual <- as.character(data$manual)
      manual[manual == "res cDC1s"] <- "cDC1s"
      manual[manual == "Mig cDC1s"] <- "cDC1s"
      manual[manual == "res cDC2s"] <- "cDC2s"
      manual[manual == "Mig cDC2s"] <- "cDC2s"
      
      medians[[file]] <- sapply(channels, function(channel){
        tapply(exprs(data$flowFrame)[, channel], manual, median)
      })
    }, error = function(e) {print(paste("Failed for", file))})
    counts[[file]] <- table(manual)
    
  }
  saveRDS(medians, paste0("medians_", analysis, ".RDS"))
  saveRDS(counts, paste0("counts_", analysis, ".RDS"))
  message(paste0("Saved results in medians_", analysis, ".RDS and counts_", analysis, ".RDS"))
  
  return(list(medians = medians, 
              counts = counts))
}

medians_to_df <- function(markers_of_interest, celltypes, tissues, samples,
                          file_description, medians, counts){
  
  
  
  col_description <- expand.grid(Sample = samples,
                                 Tissue = tissues,
                                 Marker = markers_of_interest)
  col_description$id <- paste(col_description$Marker, 
                              col_description$Tissue, 
                              col_description$Sample,
                              sep = "_")
  
  
  results <- matrix(NA, nrow = length(celltypes), 
                    ncol = nrow(col_description), 
                    dimnames = list(celltypes,col_description$id))
  result_counts <- matrix(NA, 
                          nrow = length(celltypes), 
                          ncol = nrow(col_description), 
                          dimnames = list(celltypes,col_description$id))
  
  for(sample in samples){
    for(tissue  in tissues){
      
      selected_files <- file_description[file_description$Sample_ID == sample & 
                                           file_description$Tissue == tissue, 
                                         "File"]
      
      for(file in selected_files){
        selected_celltypes <- celltypes[which(celltypes %in% rownames(medians[[file]]))]
        if(length(selected_celltypes) > 0 )
        if(analysis == "mFcyR in wild-type mice"){
          if(length(grep("FcRn", basename(file))) > 0){
            selected_markers <- "FcRn"
          } else {
            selected_markers <- c("FcyRI", "FcyRIIB", "FcyRIII", "FcyRIV")
          }
        } else {
          selected_markers <- colnames(medians[[file]])
        }
        results[selected_celltypes, 
                paste(selected_markers, tissue, sample, sep = "_")] <- 
          medians[[file]][selected_celltypes, selected_markers]
        result_counts[selected_celltypes, 
                      paste(selected_markers, tissue, sample, sep = "_")] <- 
          counts[[file]][selected_celltypes]
      }
      
    }
  }
  return(list(results = results,
              result_counts = result_counts))
}

remove_populations <- function(results, populations_to_remove){
  for(tissue in names(populations_to_remove)){
    for(population in populations_to_remove[[tissue]]){
      results[grep(population, rownames(results)),
              grep(tissue, colnames(results))] <- NA
    }
  }
  return(results)
}

subtract_FMO <- function(results){
  fmos <- grep("FMO", colnames(results), value = TRUE)
  results_diff <- results
  for(fmo in fmos){
    cols_to_adapt <- grep(gsub("_[^_]*$", "", fmo), colnames(results), value = TRUE)
    results_diff[,cols_to_adapt] <- results[,cols_to_adapt] - results[,rep(fmo, length(cols_to_adapt))]
  }
  results_diff[results_diff < 0] <- 0
  results_diff[, grep("FMO", colnames(results), invert = TRUE)]
}

normalize_per_marker <- function(results, markers_of_interest){
  max_value_per_marker <- sapply(markers_of_interest, function(marker){ 
    cols <- grep(paste0(marker, "_"), colnames(results), value = TRUE)
    max(results[, cols], na.rm = TRUE)
  })
  
  marker_per_col <- gsub("_.*", "", colnames(results))
  max_to_norm <- max_value_per_marker[marker_per_col]
  
  results <- t(apply(results, 1, function(x) x / max_to_norm))
}

plot_heatmap <- function(to_plot){
  column_info <- data.frame(id = colnames(to_plot),
                            group = gsub("(.*)_([^_]*)$", "\\1", colnames(to_plot)),
                            name = gsub("(.*)_([^_]*)$", "\\2", colnames(to_plot)),
                            geom = c("text", rep("circle", ncol(to_plot)-1)),
                            palette = gsub("(.*)_(.*)_([^_]*)$", "\\1", colnames(to_plot)),
                            legend = TRUE,
                            hjust = 1)
  column_info[1, "palette"] <- NA
  
  column_groups <- data.frame(group = unique(column_info$group),
                              level1 = gsub("_.*", "", unique(column_info$group)),
                              level2 = gsub(".*_", "", unique(column_info$group)),
                              palette = gsub("_.*", "", unique(column_info$group)))
  column_groups[1, "level2"] <- "Tissue"
  
  p <- funkyheatmap::funky_heatmap(to_plot, 
                                   palettes = list("FcyRI" = colorRampPalette(c("white","#e00480"))(100),
                                                   "FcyRIIA" = colorRampPalette(c("white","#22a395"))(100),
                                                   "FcyRIIB" = colorRampPalette(c("white","#ff9a00"))(100),
                                                   "FcyRIII" =  colorRampPalette(c("white","#006131"))(100),
                                                   "FcyRIV" = colorRampPalette(c("white","#01bb00"))(100),
                                                   "FcRn" =  colorRampPalette(c("white","#ff5c27"))(100)),
                                   column_info = column_info,
                                   column_groups = column_groups,
                                   scale_column = FALSE,
                                   col_annot_offset = 2,
                                   col_annot_angle = 90)
  p
}

summarise <- function(to_plot, groups){
  to_plot_summarised <- sapply(groups, function(group){ 
    cols <- grep(paste0(group,"_"), colnames(to_plot), value = TRUE)
    m <- apply(to_plot[, cols], 1, mean, na.rm = TRUE)
    names(m) <- to_plot$id
    m
  })
  colnames(to_plot_summarised) <- paste0(colnames(to_plot_summarised), "_mean")
  to_plot_summarised <- as.data.frame(to_plot_summarised) %>%
    tibble::rownames_to_column("id")
}