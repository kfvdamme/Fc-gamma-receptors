library(Seurat)
library(anndata)
library(SeuratDisk)
library(tidyverse)

Convert("ALL_Tissue_Global_Clustering_Scanpy.h5ad", dest = "h5seurat", assay = 'RNA', overwrite = T)
seuratobj = LoadH5Seurat("ALL_Tissue_Global_Clustering_Scanpy.h5seurat") #, meta.data = FALSE, misc = FALSE)
saveRDS(seuratobj, "cynoatlas.rds")

# add metadata

metadata = read.delim("Metadata_ALL_Tissue_Global_Clustering.txt", header = TRUE)

## check similarity of metadata

head(metadata)
head(seuratobj@meta.data)
head(seuratobj@reductions$umap@cell.embeddings)
data.frame(metadata$percent_mito, seuratobj@meta.data$percent_mito)
# indeed similar
table(metadata$sing.celltype)
names(table(metadata$sing.celltype)) # 173 cell types in sing.celltype
seuratobj@meta.data$sing.celltype = metadata$sing.celltype
table(metadata$celltype)
names(table(metadata$celltype)) # 113 cell types in celltype
seuratobj@meta.data$celltype = metadata$celltype

p1 = DimPlot(seuratobj, group.by = "celltype", label = TRUE) + NoLegend()
p2 = DimPlot(seuratobj, group.by = "sing.celltype", label = TRUE) + NoLegend()
ggpubr::ggarrange(p1, p2)
# saveRDS(seuratobj, "cynoatlas.rds")

# Data exploration

Idents(seuratobj) = "celltype"
VlnPlot(seuratobj, "TRIM21", pt.size = 0) + NoLegend()

fc_receptors = c("FCGR1A", "FCGR2B", "FCGR3A", "FCGRT")
fc_receptor_data = FetchData(seuratobj, vars = c("celltype", fc_receptors))

sorted_celltypes_by_fc_receptors = fc_receptor_data %>%
  mutate(mean_fcrs = rowMeans(across(FCGR1A:FCGR3A))) %>%
  group_by(celltype) %>%
  summarize(across(FCGR1A:mean_fcrs, ~ mean(., na.rm = TRUE))) %>%
  arrange(desc(mean_fcrs))

DotPlot(seuratobj, fc_receptors)

sort(table(seuratobj@meta.data$celltype), decreasing = TRUE)
pop_per_organ = data.frame(unclass(table(seuratobj@meta.data$organ, seuratobj@meta.data$celltype)))
write.table(pop_per_organ, "pop_per_organ.csv")
names(sort(table(seuratobj@meta.data$celltype), decreasing = TRUE))
?write.csv

# Explore and rename annotation

Idents(seuratobj) <- "celltype"
seuratobj_subset = subset(seuratobj, downsample = 75)

cluster.markers <- lapply(names(table(seuratobj_subset[["celltype"]])),
                          function(x)FindMarkers(seuratobj_subset,ident.1 = x,min.pct = 0.30)) 

sapply(0:(length(cluster.markers)-1),function(x) {
  cluster.markers[[x+1]]$gene <<- rownames(cluster.markers[[x+1]])
  cluster.markers[[x+1]]$cluster <<- x})
as_tibble(do.call(rbind,cluster.markers)) %>% arrange(p_val_adj) -> cluster.markers
cluster.markers

top10 = cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

DoHeatmap(subset(ts_subset, downsample=20), features = top10$gene, assay = "RNA", size = 3,
          slot = "counts", label = TRUE, disp.min = -0.5, draw.lines = TRUE, lines.width = 1) + 
  NoLegend() +
  theme(text = element_text(size = 8)) +
  scale_fill_binned(type = "viridis", na.value = "white")


seuratobj = RenameIdents(object = seuratobj,
                         "Granule cell"                           = "parenchymal",                           
                         "Hepatocyte"                             = "parenchymal",                             
                         "Male reproductive epithelial cell"      = "epithelial", 
                         "Salivary acinar cell"                   = "epithelial", 
                         "Epididymis sterociliated cell"          = "epithelial", 
                         "Adipocyte progenitor/stromal cell"      = "stromal", 
                         "Smooth muscle cell"                     = "stromal", 
                         "Cortical excitatory neuron"             = "parenchymal", 
                         "Pancreas acinar cell"                   = "epithelial", 
                         "Endothelial cell"                       = "endothelial", 
                         "Adipocyte"                              = "parenchymal", 
                         "B cell"                                 = "B cells", 
                         "Type II skeletal muscle myonuclei"      = "stromal",
                         "Fallopian tube epithelial cell"         = "epithelial",         
                         "Mesothelial cell"                       = "epithelial", 
                         "T cell"                                 = "T cells", 
                         "Alveolar epithelial type 1 cell"        = "epithelial", 
                         "Gallbladder glandular cell"             = "epithelial", 
                         "Glomerulosa cell"                       = "parenchymal", 
                         "Stromal cell"                           = "Stromal", 
                         "Macrophage"                             = "macrophages", 
                         "Type IIa skeletal muscle myonuclei"     = "stromal", 
                         "Lactotrope"                             = "parenchymal", 
                         "Ascending loop of Henle cell"           = "epithelial", 
                         "Fasciculata cell"                       = "parenchymal", 
                         "Pinealocyte"                            = "parenchymal", 
                         "Proximal tubule cell"                   = "epithelial", 
                         "Somatotrope"                            = "parenchymal", 
                         "Bladder intermediate cell"              = "epithelial", 
                         "Granulosa cell"                         = "stromal", 
                         "Inhibitory neuron"                      = "parenchymal", 
                         "Uterus luminal epithelial cell"         = "epithelial", 
                         "Uterus stromal cell"                    = "stromal", 
                         "Salivary serous/acinar cell"            = "epithelial", 
                         "Oligodendrocyte"                        = "parenchymal", 
                         "Principal cell"                         = "epithelial", 
                         "Epidermal stromal cell"                 = "stromal", 
                         "Lingual epithelial cell"                = "epithelial", 
                         "Neuron"                                 = "parenchymal", 
                         "Thyroid follicular cell"                = "epithelial", 
                         "Astrocyte"                              = "parenchymal", 
                         "T cell/NKT cell"                        = "T cells/NKT cells", 
                         "Alveolar epithelial type 2 cell"        = "epithelial", 
                         "Male reproductive basal cell"           = "epithelial", 
                         "Gallbladder mucous cell"                = "epithelial", 
                         "Upper digestive tract epithelial cells" = "epithelial", 
                         "Type I skeletal muscle myonuclei"       = "stromal", 
                         "Respiratory ciliated cell"              = "epithelial", 
                         "Cardiomyocyte"                          = "stromal", 
                         "Tongue stromal cell"                    = "stromal", 
                         "Intercalated cell"                      = "epithelial", 
                         "Plasma B cell"                          = "B cells", 
                         "Corticotrope"                           = "parenchymal", 
                         "Female reproductive ciliated cell"      = "epithelial", 
                         "Retinal pigment epithelial cell"        = "epithelial", 
                         "Theca cell/granulosa"                   = "stromal", 
                         "Leyding/myoid cells"                    = "parenchymal", 
                         "Distal convoluted tubule cell"          = "epithelial", 
                         "Gonadotrope cell"                       = "parenchymal", 
                         "Goblet/serous cell"                     = "epithelial", 
                         "Retinal cell"                           = "parenchymal", 
                         "Liver endothelial cell 1"               = "endothelial", 
                         "Pancreas ductal and islet cell"         = "epithelial", 
                         "Epidermal endothelial cell"             = "endothelial", 
                         "Dscending loop of Henle cell"           = "epithelial", 
                         "Dopaminergic neuron"                    = "parenchymal", 
                         "Photoreceptor"                          = "parenchymal", 
                         "Erythroid cell"                         = "other", 
                         "Unknown cell"                           = "other", 
                         "Enteric glial cell"                     = "parenchymal", 
                         "Salivary intercalated cell"             = "epithelial", 
                         "Folliculostellate cell"                 = "stromal", 
                         "Keratinocyte "                          = "epithelial", 
                         "Adrenal endothelial cell"               = "endothelial", 
                         "Chromaffin cell"                        = "parenchymal", 
                         "Salivary goblet/mucous cell"            = "epithelial", 
                         "Bladder urothelial cell"                = "epithelial", 
                         "Pseudostratified columnar cell"         = "epithelial",         
                         "Molecular interneuron"                  = "parenchymal",  
                         "Monocyte"                               = "monocytes", 
                         "Dendritic cell"                         = "dendritic cells", 
                         "Liver immune cell"                      = "other", 
                         "Epidermal pericyte"                     = "stromal", 
                         "Lower digestive tract epithelial cell"  = "epithelial", 
                         "Epididymal epithelial cell"             = "epithelial", 
                         "OPC"                                    = "other", 
                         "Adrenal immune cell"                    = "other", 
                         "Melanocyte"                             = "parenchymal", 
                         "Connecting tubule cell"                 = "epithelial", 
                         "Immune mesothelial cell"                = "other", 
                         "Salivary stromal cell"                  = "stromal", 
                         "Digestive tract stromal cell"           = "stromal", 
                         "Podocyte"                               = "parenchymal", 
                         "Male reproductive smooth muscle cell"   = "stromal", 
                         "Digestive tract endothelial cell 1"     = "endothelial", 
                         "Liver endothelial cell 2"               = "endothelial", 
                         "Viseral adipose mesothelial cell"       = "epithelial", 
                         "Epidermal adipocyte "                   = "stromal", 
                         "Bile ductal cell"                       = "epithelial", 
                         "Digestive tract endothelial cell 2"     = "endothelial", 
                         "Salivary NKT "                          = "T cells/NKT cells", 
                         "Epidermal smooth muscle cell"           = "stromal", 
                         "Sertoli cell"                           = "epithelial", 
                         "Gallbladder endothelial cell"           = "endothelial", 
                         "Stellate cell"                          = "stromal", 
                         "Gallbladder smooth muscle cell"         = "stromal", 
                         "Thyroid endothelial cell"               = "endothelial", 
                         "Microglia"                              = "macrophages", 
                         "Thyroid stromal cell"                   = "stromal", 
                         "Epidydimal monocyte"                    = "monocytes", 
                         "Purkinje cell"                          = "parenchymal", 
                         "Pineal micoglia"                        = "macrophages", 
                         "Pineal oligodendrocyte"                 = "parenchymal")