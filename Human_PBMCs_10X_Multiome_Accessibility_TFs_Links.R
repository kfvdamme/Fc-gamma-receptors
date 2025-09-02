library(tidyverse)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(RColorBrewer)
library(patchwork)
library(rtracklayer)
source(file.path("D:/Experiments/20230906 Tabula Sapiens and Muris/colors", 'Cell_pop_colors.R'))
set.seed(1234)

pbmc = readRDS("scenicplus_clustered_consensusSIGNAC_obj.RDS")
pbmc$ATAC@fragments[[1]]@path = paste0(getwd(), 
                                       "/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

r2g_links <- read.table(
  'r2g.rho.bed', sep = '\t', skip = 1,
  col.names = c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp',
                'color', 'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName',
                'sourceStrand', 'targetChrom', 'targetStart', 'targetEnd', 'targetName',
                'targetStrand'))

r2g_links <- transform(r2g_links, start_ = pmin(sourceStart, targetEnd))
r2g_links <- transform(r2g_links, end_ = pmax(sourceStart, targetEnd))
r2g_links <- r2g_links[c('chrom', 'start_', 'end_', 'targetName', 'value', 'chrom', 'chromStart', 'chromEnd')]
colnames(r2g_links) <- c('seqnames', 'start', 'end', 'gene', 'score', 'chrom', 'chromStart', 'chromEnd')

eRegulon_regions <- read.table('eRegulons.bed', sep = '\t')
colnames(eRegulon_regions) <- c('seqnames', 'start', 'end', 'name')
gr_eRegulon_regions <- makeGRangesFromDataFrame(eRegulon_regions, keep.extra.columns = TRUE)

pbmc$GEX_celltype = factor(
  pbmc$GEX_celltype,
  levels = c(
    "CD14+_Monocytes",
    "FCGR3A+_Monocytes",
    "Dendritic_cells_1",
    "Dendritic_cells_2",
    "B_cells",
    "CD4_T_cells",
    "CD8_T_cells",
    "NK_cells"
  )
)

# For FCGR2A

region_to_plot <- "chr1-161475000-161575000" # limited genomic region
region_to_plot <- "chr1-161300000-161700000" # large genomic region
gene_to_plot <- "FCGR2A"
cell_line = "K562"

# For FCGRT

region_to_plot <- "chr19-49490000-49540000" # limited genomic region
region_to_plot <- "chr19-49450000-49670000" # large genomic region
gene_to_plot <- "FCGRT"
cell_line = "HepG2"

##### Plotting -----

r2g_links_gene <- subset(r2g_links, gene == gene_to_plot)
r2g_links_gr <- makeGRangesFromDataFrame(subset(r2g_links_gene, chromStart %in% eRegulon_regions$start)[c('seqnames', 'start', 'end', 'gene', 'score')], keep.extra.columns = TRUE)
Links(pbmc) <- r2g_links_gr

chip = vector()
for (file in list.files(paste0(getwd(), "/Encode_cell_lines/", cell_line))) {
  chip_add =  as_tibble(read.table(file.path(paste0("D:/Experiments/20230906 Tabula Sapiens and Muris/Scenic+ pbmcs/Encode_cell_lines/", cell_line), 
                                             file),
                                   sep = '\t')) %>% mutate(V4 = str_sub(file, end = -5))
  chip_add = chip_add[1:4]
  chip = rbind(chip, chip_add)
}
colnames(chip) = c('seqnames', 'start', 'end', 'name')
rm(chip_add)
chip_regions <- makeGRangesFromDataFrame(chip, keep.extra.columns = TRUE)
chip_region <- chip_regions[queryHits(GenomicRanges::findOverlaps(chip_regions, StringToGRanges(region_to_plot)))]
chip_region$name <- unname(sapply(chip_region$name, function(x) strsplit(x, "_")[[1]][1]))
chip_region <- data.frame(
  start = GenomicRanges::start(chip_region),
  end = GenomicRanges::start(chip_region),
  name = chip_region$name)
chip_region <- chip_region[order(chip_region$start, chip_region$name), ]
chip_region$id <- unname(sapply(chip_region$name, match, unique(chip_region$name)))
start <- strtoi(strsplit(region_to_plot, "-")[[1]][2])
end <- strtoi(strsplit(region_to_plot, "-")[[1]][3])
peak_plot <- ggplot(data = chip_region, mapping = aes(x = start + 250, y = id)) & 
  geom_point(size = 1) &
  xlim(start, end) & 
  geom_text(data = chip_region, mapping = aes(y = id, x = (max(chip_region$start) + 1000), label = name),
            hjust = 0, size = 3) & 
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_blank())

cov_plot <- CoveragePlot(
  object = pbmc,
  region = region_to_plot,
  group.by = 'GEX_celltype',
  annotation = FALSE,
  peaks = FALSE,
  links=FALSE,
  bigwig.scale = 'separate'
)  & scale_fill_manual(values=cell_colors)

gene_plot <- AnnotationPlot(
  object = pbmc,
  region = region_to_plot
)

link_plot <- LinkPlot(
  object = pbmc,
  region = region_to_plot,
  min.cutoff = -1) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Blues"), limits=c(0, NA)) 

p <- CombineTracks(
  plotlist = list(cov_plot, peak_plot, link_plot, gene_plot),
  #expression.plot = expr_plot,
  heights = c(5, 1, 1, 1) # ,
  # widths = c(10,4)
)
p
