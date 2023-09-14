### 0. Load packages -----------------------------------------------------------
set.seed(1234)

library(Seurat)
library(ggplot2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(venn)
library(geco.utils)
library(magrittr)
library(geco.visu)
library(viridis)
library(dittoSeq)
library(scuttle)
library(VennDiagram)
library(fgsea)
library(msigdbr)
source("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/GSEA_tools/Range_enrichment_SC.R")


### 1. Colors and directories  -------------------------------------------------

merged.directory <- paste0("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/merged/")
if(!dir.exists(merged.directory)){dir.create(merged.directory)}

my_scale_exp = scale_colour_gradientn(colors = c( "gray70", "#ffff00",  "red", "#cc0000"))

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Fibroblasts - HSCs - VSMCs" = "#8f35e5", 
                   "Biliary epithelial cells" = "#e53561", 
                   "HOX+ Progenitor" = "#ff9f2e",
                   "NK" = "pink",
                   "T cells" = "#36e22e",
                   "SAA2-CRP high hepatocytes/fibroblasts?" = "blue",
                   "Periportal bifurcations cells ?" = "forestgreen")

cols_status_cn_LOH <- c("cnLOH high confidence" = "red",
                        "cnLOH low confidence" = "orange",
                        "no cnLOH high confidence" = "darkblue",
                        "no cnLOH low confidence" = "blue",
                        "unknown" = "gray90")

cols_status_cn_LOH_BAF <- c("cnLOH" = "red", "no cnLOH" = "blue", "unknown" = "gray90")


final.cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "Endothelial cells" = "#3399FF", 
                   "Kupffer cells" = "#ff02f5", 
                   "Fibroblasts - HSCs" = "#66FF66", 
                   "Cholangiocytes" = "#e53561", 
                   "HOX+ Progenitor" = "#0019ef",
                   "Tumor cells conta." = "#8d00dc",
                   "VSMCs" = "yellow")


cols.numbat <- c("balanced" = "blue", "unknown" = "gray80", "cn-LOH" = "orange", "NA" = "gray80")

samp.cols <- c("CHC3559N" = "orange", "CHC3115N" = "red", "CHC4001N" = "#ff97f6", "CHC2996N" = "#99FFFF")

samp.cols2  <- c("#4001" = "darkblue", "#3559" = "#fbddf5", "#3115" = "#039a87", "#2996" = "#d0ff5d")



display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


### 2. Load seurat objects -----------------------------------------------------

CHC3115N <- readRDS(file="D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/CHC3115N/preprocess.QC/CHC3115NpostQC.rds")

CHC4001N <- readRDS(file="D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/CHC4001N/preprocess.QC/CHC4001NpostQC.rds")

CHC2996N <- readRDS(file="D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/CHC2996N/preprocess.QC/CHC2996NpostQC.rds")

CHC3559N <- readRDS(file="D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/CHC3559N/preprocess.QC/CHC3559NpostQC.rds")

length((colnames(CHC3115N)))
length((colnames(CHC4001N)))
length((colnames(CHC2996N)))
length((colnames(CHC3559N)))


### 3. CHC3115N  ---------------------------------------------------------------

# Add cn-LOH by cell info in meta data 
# snp method Théo 

recap_by_cell_3115 <- geco.load("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/BAF_single_cell/Results/recap_by_cell_CHC3115N.Rdata")

venn::venn(list(cells_SC_readcount = recap_by_cell_3115$ReadGroup, cells_seurat = colnames(CHC3115N)))

CHC3115N@meta.data <- CHC3115N@meta.data %>%
  mutate(ReadGroup = rownames(.)) %>% 
  dplyr::select(-intersect(colnames(recap_by_cell_3115)[-1], colnames(.))) %>% 
  left_join(recap_by_cell_3115) %>% 
  set_rownames(.$ReadGroup) %>% 
  mutate(status_cnLOH = ifelse(is.na(status_cnLOH), "unknown", status_cnLOH),
         status_cnLOH_BAF = ifelse(is.na(status_cnLOH_BAF), "unknown", status_cnLOH_BAF))

# numbat method
recap_by_cell_3115_numbat <- geco.load("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/BAF_single_cell/Results/recap_by_cell_numbat_CHC3115N.Rdata")

CHC3115N@meta.data <- CHC3115N@meta.data %>%
  dplyr::select(-intersect(colnames(recap_by_cell_3115_numbat)[-1], colnames(.))) %>% 
  left_join(recap_by_cell_3115_numbat) %>% 
  set_rownames(.$ReadGroup)

DefaultAssay(CHC3115N) <- "SCT"


# Rename Idents 
Idents(object = CHC3115N) <- CHC3115N@meta.data$cell.type

# Plot gene expression
FeaturePlot(CHC3115N, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC3115N, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC3115N, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC3115N, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp

DimPlot(CHC3115N, group.by = "cell.type", cols = cell.type.cols, reduction = "umap", raster=F, pt.size = 0.2)
DimPlot(CHC3115N, group.by = "status_cnLOH_BAF", cols = c("red", "blue", "gray80"), reduction = "umap", raster=F, pt.size = 0.1)
DimPlot(CHC3115N, group.by = "status_cnLOH", cols = cols_status_cn_LOH, reduction = "umap", raster=F, pt.size = 0.1)
DimPlot(CHC3115N, group.by = "status_cn_LOH_numbat", cols = cols.numbat, reduction = "umap", raster=F, pt.size = 0.1)


### 4. CHC4001N  ---------------------------------------------------------------

# Add cn-LOH by cell info in meta data

recap_by_cell_4001 <- geco.load("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/BAF_single_cell/Results/recap_by_cell_CHC4001N.Rdata")

venn::venn(list(cells_SC_readcount = recap_by_cell_4001$ReadGroup, cells_seurat = colnames(CHC4001N)))

CHC4001N@meta.data <- CHC4001N@meta.data %>%
  mutate(ReadGroup = rownames(.)) %>% 
  dplyr::select(-intersect(colnames(recap_by_cell_4001)[-1], colnames(.))) %>% 
  left_join(recap_by_cell_4001) %>% 
  set_rownames(.$ReadGroup) %>% 
  mutate(status_cnLOH = ifelse(is.na(status_cnLOH), "unknown", status_cnLOH),
         status_cnLOH_BAF = ifelse(is.na(status_cnLOH_BAF), "unknown", status_cnLOH_BAF))

DefaultAssay(CHC4001N) <- "SCT"

# Rename Idents 
Idents(object = CHC4001N) <- CHC4001N@meta.data$cell.type

# Plot gene expression
FeaturePlot(CHC4001N, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC4001N, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC4001N, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC4001N, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp

DimPlot(CHC4001N, group.by = "cell.type", cols = cell.type.cols, reduction = "umap", raster=F, pt.size = 0.2)
DimPlot(CHC4001N, group.by = "status_cnLOH_BAF", cols = c("orange", "blue", "gray80"), reduction = "umap", raster=F, pt.size = 0.1)
DimPlot(CHC4001N, group.by = "status_cnLOH", cols = cols_status_cn_LOH, reduction = "umap", raster=F, pt.size = 0.1)

### 5. CHC3559N  ---------------------------------------------------------------
# Add cn-LOH by cell info in meta data 

recap_by_cell_3559 <- geco.load("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/BAF_single_cell/Results/recap_by_cell_CHC3559N.Rdata")

venn::venn(list(cells_SC_readcount = recap_by_cell_3559$ReadGroup, cells_seurat = colnames(CHC3559N)))

CHC3559N@meta.data <- CHC3559N@meta.data %>%
  mutate(ReadGroup = rownames(.)) %>% 
  dplyr::select(-intersect(colnames(recap_by_cell_3559)[-1], colnames(.))) %>% 
  left_join(recap_by_cell_3559) %>% 
  set_rownames(.$ReadGroup) %>% 
  mutate(status_cnLOH = ifelse(is.na(status_cnLOH), "unknown", status_cnLOH),
         status_cnLOH_BAF = ifelse(is.na(status_cnLOH_BAF), "unknown", status_cnLOH_BAF))

DefaultAssay(CHC3559N) <- "SCT"

# Rename Idents 
Idents(object = CHC3559N) <- CHC3559N@meta.data$cell.type

# Plot gene expression
FeaturePlot(CHC3559N, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC3559N, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC3559N, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
FeaturePlot(CHC3559N, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp

DimPlot(CHC3559N, group.by = "cell.type", cols = cell.type.cols, reduction = "umap", raster=F, pt.size = 0.2)
DimPlot(CHC3559N, group.by = "status_cnLOH_BAF", cols = c("orange", "blue", "gray80"), reduction = "umap", raster=F, pt.size = 0.1)
DimPlot(CHC3559N, group.by = "status_cnLOH", cols = cols_status_cn_LOH, reduction = "umap", raster=F, pt.size = 0.1)


### 6. Merge -------------------------------------------------------------------
  # 6.1 Merge creation ---------------------------------------------------------

CHC2996N <- RenameCells(object = CHC2996N, add.cell.id = "CHC2996N")

CHC3115N <- RenameCells(object = CHC3115N, add.cell.id = "CHC3115N")

CHC4001N <- RenameCells(object = CHC4001N, add.cell.id = "CHC4001N")

CHC3559N <- RenameCells(object = CHC3559N, add.cell.id = "CHC3559N")

# Add missing info cnLOH

CHC2996N@meta.data <- CHC2996N@meta.data %>%
  mutate(status_cnLOH = "no cnLOH high confidence",
         status_cnLOH_BAF = "no cnLOH")

# merge

all.merged <- merge(CHC2996N, y = c(CHC3115N, CHC4001N, CHC3559N), 
                    merge.data = TRUE)  # keep merge.data = TRUE to keep previous normalization

# Normalization with SCTransform
obj.list <- SplitObject(all.merged, split.by = "orig.ident")

DefaultAssay(all.merged) <- "SCT"

obj.list <- lapply(X = obj.list, FUN = SCTransform)

# Create Variable features

my_integration_features <- SelectIntegrationFeatures(
  obj.list,
  nfeatures = 3000,
  verbose = TRUE)

VariableFeatures(all.merged) <- my_integration_features

# PCA on gene expression

all.merged <- RunPCA(all.merged, features = VariableFeatures(all.merged), npcs=50, verbose=T) # 50 components
ElbowPlot(all.merged, ndims=50)
all.merged <- RunUMAP(all.merged, dims=1:30, verbose=F)

# Clustering 

DefaultAssay(all.merged) <- "SCT"

all.merged <- FindNeighbors(all.merged, dims = 1:30, k.param = 15)
all.merged <- FindClusters(all.merged, resolution = 0.1)
all.merged$clusters_0_1 <- all.merged$seurat_clusters

all.merged <- FindClusters(all.merged, resolution = 0.2)
all.merged$clusters_0_2 <- all.merged$seurat_clusters


p <- DimPlot(all.merged, group.by = "clusters_0_2", reduction = "umap", raster=F, pt.size = 0.1) + 
  scale_color_viridis_d(option="turbo")
p
ggsave(paste0(merged.directory, "clusters_0_2.png"), plot = p, device = "png", width =11, height = 10, units = "cm", dpi = 320)


all.merged@meta.data <- all.merged@meta.data %>%
  mutate(final.cell.type = case_when(clusters_0_2 =="0" ~ "Hepatocytes",
                                   clusters_0_2 == "1" ~ "Hepatocytes",
                                   clusters_0_2 == "2" ~ "VSMCs",
                                   clusters_0_2 == "3" ~ "Hepatocytes",
                                   clusters_0_2 == "4" ~ "Endothelial cells",
                                   clusters_0_2 == "5" ~ "Hepatocytes",
                                   clusters_0_2 == "6" ~ "Kupffer cells",
                                   clusters_0_2 == "7" ~ "Cholangiocytes",
                                   clusters_0_2 == "8" ~ "Fibroblasts - HSCs",
                                   clusters_0_2 == "9" ~ "Hepatocytes",
                                   clusters_0_2 == "10" ~ "HOX+ Progenitor",
                                   clusters_0_2 == "11" ~ "Tumor cells conta."))


all.merged@meta.data <- all.merged@meta.data %>%
  mutate(orig.ident = case_when(orig.ident =="CHC4001N" ~ "#4001",
                                orig.ident == "CHC3559N" ~ "#3559",
                                orig.ident == "CHC2996N" ~ "#2996", 
                                orig.ident == "CHC3115N" ~ "#3115"))

saveRDS(all.merged, file=paste0(merged.directory, "all.merged.rds"))

# Load directly merged object saved 
all.merged <- readRDS("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/merged/all.merged.rds")

# Look at QCs 

all.merged$nFeature_RNA_log <- log10(all.merged$nFeature_RNA)

p<-VlnPlot(
  object = all.merged,
  features = c("nCount_RNA_log","nFeature_RNA_log", "percent.mt"),
  group.by = "orig.ident",
  cols=samp.cols2,
  ncol = 3,
  pt.size = 0
) & NoLegend()


p


length(colnames(all.merged))

median(all.merged@meta.data$nCount_RNA)

  # 6.2 Analysis with all cells ------------------------------------------------

# Create source file umap figure 5

umap_coordinates <- all.merged[["umap"]]@cell.embeddings
clust <- all.merged[[c("final.cell.type", "status_cnLOH_BAF", "BAF")]]
Source_figure_5 = cbind(clust, umap_coordinates)

write.table(Source_figure_5, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/Source_figure_5.txt", sep="\t")


# UMAP with cell types and orig ident annotations 

p <- DimPlot(all.merged, group.by = "final.cell.type", reduction = "umap", raster=F, pt.size = 0.1, cols = final.cell.type.cols)  + ggtitle("") + 
  theme(axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=18))


p
ggsave(paste0(merged.directory, "final.cell.type.png"), plot = p, device = "png", width =17, height = 10, units = "cm", dpi = 320)

p <- DimPlot(all.merged, group.by = "cell.type", cols = cell.type.cols, reduction = "umap", raster=F, pt.size = 0.1)
p
ggsave(paste0(merged.directory, "cell.type.from.each.sample.png"), plot = p, device = "png", width =17, height = 10, units = "cm", dpi = 320)

p <- DimPlot(all.merged, group.by = "orig.ident", cols = samp.cols2, reduction = "umap", raster=F, pt.size = 0.1) + ggtitle("") + 
  theme(axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=18))


p
ggsave(paste0(merged.directory, "orig.ident2.png"), plot = p, device = "png", width =13, height = 10, units = "cm", dpi = 320)

p<-DimPlot(all.merged, group.by = "status_cnLOH_BAF", cols = c("red", "blue", "gray90"), 
           reduction = "umap", raster=F, pt.size = 0.01, shuffle = T) + ggtitle("") + 
  theme(axis.line = element_line(linewidth = 1),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=18))


p
ggsave(paste0(merged.directory, "status_cnLOH_BAF.png"), plot = p, device = "png", width =15, height = 10, units = "cm", dpi = 320)

# Plot gene expression
p <- FeaturePlot(all.merged, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "IGF2.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "H19.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="LGR5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "LGR5.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "HAL.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "CYP2E1.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "CD163.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "FLT1.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "KRT19.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "KRT7.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="HOXB4",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "HOXB4.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

p <- FeaturePlot(all.merged, reduction="umap", features="CD274",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(merged.directory, "CD274.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)


# percentage cell type per patient 

all.merged.sce <- as.SingleCellExperiment(all.merged)
all.merged.sce@colData$final.cell.type= factor(all.merged.sce@colData$final.cell.type, levels = c("Tumor cells conta.", "Hepatocytes", 
                                                                                                  "Cholangiocytes", 
                                                                                                  "Endothelial cells",
                                                                                                  "VSMCs",
                                                                                                  "Fibroblasts - HSCs", 
                                                                                                  "Kupffer cells", 
                                                                                                  "HOX+ Progenitor"))


p <- dittoBarPlot(all.merged.sce, var = "final.cell.type", group.by = "orig.ident") +
  scale_fill_manual(values = final.cell.type.cols) + ggtitle("") + labs(y="Proportion of each cell type")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=18),
        axis.line = element_line(linewidth = 1),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=12))

p


ggsave(paste0(merged.directory, "prop.cell.type.png"), plot = p, device = "png", width =10, height = 10, units = "cm", dpi = 320)

# Probability cnLOH per cell type in mosaic samples

all.mosaic <- subset(all.merged, subset=orig.ident %in% c("#3115", "#3559", "#4001"))

all.mosaic@meta.data$final.cell.type=factor(all.mosaic@meta.data$final.cell.type, levels=c("Hepatocytes", 
                                                                                           "Cholangiocytes", 
                                                                                           "Tumor cells conta.",
                                                                                           "Endothelial cells",
                                                                                           "VSMCs",
                                                                                           "Fibroblasts - HSCs", 
                                                                                           "Kupffer cells", 
                                                                                           "HOX+ Progenitor"))



p<-VlnPlot(
  object = all.mosaic,
  group.by = "final.cell.type",
  features = c("BAF"),
  cols = final.cell.type.cols,
  pt.size = 0,adjust = 0.9) + 
  ylab("B-allele frequency (BAF)") +
  theme(
    axis.text.x = element_text(size=16, angle=45),
    axis.line = element_line(linewidth = 1),
    axis.text.y = element_text(size=16),
    axis.title.y = element_text(size=16, margin =margin(t = 0, r = 8, b = 0, l = 0)),
    axis.title.x = element_blank()) + labs(title="") & NoLegend() 

p

ggsave(paste0(merged.directory, "BAF.per.cell.type.png"), plot = p, device = "png", width =15, height = 10, units = "cm", dpi = 320)


p<-RidgePlot(
  object = all.mosaic,
  group.by = "final.cell.type",
  features = c("BAF"),
  cols = final.cell.type.cols) + 
  ylab("B-allele frequency (BAF)") +
  theme(
    axis.text.x = element_text(size=12, angle=45),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=12, margin =margin(t = 0, r = 8, b = 0, l = 0))) + labs(title="") & NoLegend() 

p


# Table with markers of each cell type 

all.merged = PrepSCTFindMarkers(all.merged)

Idents(all.merged) = all.merged@meta.data$final.cell.type
sample.markers <- FindAllMarkers(object = all.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1, group.by="final.cell.type")
write.table(sample.markers, quote = F,row.names = FALSE, paste(merged.directory, "_allmarkers.txt",sep=""),  sep="\t")

# Plot markers of each cell type 
Genes_of_interest <- c("CYP2B6",
                       "APOB",
                       "BICC1",
                        "KRT7", 
                       "KRT19", 
                       "FLT1", 
                       "RELN", 
                       "VWF",
                       "COL1A1", 
                       "COL6A3", 
                       "COLEC11",
                       "HGF", 
                       "ACTA2",
                       "MYH11", 
                       "TAGLN",
                       "CD163",
                       "CD80",
                       "MKI67", 
                       "TOP2A", 
                       "AFP",
                       "HOXC6",
                       "HOXA3", 
                       "HOXB3")



p <- dittoDotPlot(all.merged.sce,
             vars=Genes_of_interest,
             group.by= "final.cell.type", y.reorder = c(8, 1,7,5,6,  4, 3, 2), scale=T) +
  ylab("") +
  theme(
    axis.text.x = element_text(size=10, angle=90),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=12, margin =margin(t = 0, r = 8, b = 0, l = 0)),
    axis.title.x = element_blank()) + labs(title="") & NoLegend() 

p

ggsave(paste0(merged.directory, "marker.per.cell.type.png"), plot = p, device = "png", width =13, height = 20, units = "cm", dpi = 320)


  # 6.3 Analysis with hepatocytes only -----------------------------------------
    # 6.3.1 clustering ---------------------------------------------------------

# subset hepatocytes

merged_hepatocytes <- subset(all.merged, subset=final.cell.type %in% c("Hepatocytes"))
ElbowPlot(merged_hepatocytes)
merged_hepatocytes <- RunPCA(merged_hepatocytes, features = VariableFeatures(merged_hepatocytes), npcs=50, verbose=T) 
merged_hepatocytes <- RunUMAP(merged_hepatocytes, dims=1:30, verbose=F)
merged_hepatocytes <- FindNeighbors(merged_hepatocytes, dims = 1:30)
merged_hepatocytes <- FindClusters(merged_hepatocytes, resolution = 0.05)
merged_hepatocytes$clusters_0_05 <- merged_hepatocytes$seurat_clusters

# Read hepatocytes subset object 
merged_hepatocytes <- readRDS(file="D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/merged/merged_hepatocytes.rds")
length(colnames(merged_hepatocytes))

    # 6.3.2 visualize genes of interest and annotations ------------------------
# Plot gene expression and groups cnLOH
p <- DimPlot(merged_hepatocytes, group.by = "clusters_0_05", reduction = "umap", raster=F, pt.size = 0.1, cols =c("darkorange", 
                                                                                                     "cyan",

                                                                                                                                                                                                        "pink")) 
p
ggsave(paste0(merged.directory, "clusters_0_05_merged_hepatocytes.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)

p <-DimPlot(merged_hepatocytes, group.by = "final.cell.type", cols = final.cell.type.cols, 
            reduction = "umap", raster=F, pt.size = 0.1) + ggtitle("")
p
ggsave(paste0(merged.directory, "final.cell.type_merged_hepatocytes.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)


p <-DimPlot(merged_hepatocytes, group.by = "orig.ident", cols = samp.cols2, 
            reduction = "umap", raster=F, pt.size = 0.1) + ggtitle("") + 
  theme(axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=18))
p
ggsave(paste0(merged.directory, "orig.ident_merged_hepatocytes.png"), plot = p, device = "png", width =13, height = 10, units = "cm", dpi = 320)


p<- DimPlot(merged_hepatocytes, group.by = "status_cnLOH_BAF", cols = c("red", "blue", "gray90"), 
            reduction = "umap", raster=F, pt.size = 0.1) + ggtitle("") + 
  theme(axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=18))
  
p
ggsave(paste0(merged.directory, "status_cnLOH_BAF_merged_hepatocytes.png"), plot = p, device = "png", width =14, height = 10, units = "cm", dpi = 320)

p<- FeaturePlot(merged_hepatocytes, reduction="umap", features="IGF2",  pt.size = 0.1) + 
  my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18)) 
p
ggsave(paste0(merged.directory, "IGF2_merged_hepatocytes.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)

p<-FeaturePlot(merged_hepatocytes, reduction="umap", features="H19",  pt.size = 0.1) + 
  my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18)) 
p
ggsave(paste0(merged.directory, "H19_merged_hepatocytes.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)

p<-FeaturePlot(merged_hepatocytes, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18)) 
p
ggsave(paste0(merged.directory, "CYP2E1_merged_hepatocytes.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)

p<-FeaturePlot(merged_hepatocytes, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18)) 
p
ggsave(paste0(merged.directory, "HAL_merged_hepatocytes.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)

saveRDS(merged_hepatocytes, file=paste0(merged.directory, "merged_hepatocytes.rds"))

    # 6.3.3 Differential expression analysis------------------------------------

# Directly load merged hepatocytes object
merged_hepatocytes <- readRDS("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/merged/merged_hepatocytes.rds")

# Diff exp snRNAseq

merged_hepatocytes = PrepSCTFindMarkers(merged_hepatocytes)

mos_vs_non_mos <- FindMarkers(object = merged_hepatocytes, only.pos = TRUE,logfc.threshold = 0.1, 
                               ident.1 = "cnLOH", ident.2 = "no cnLOH",  group.by="status_cnLOH_BAF")

mos_vs_non_mos
dim(mos_vs_non_mos)

non_mos_vs_mos <- FindMarkers(object = merged_hepatocytes, only.pos = TRUE, logfc.threshold = 0.1, 
                              ident.1 = "no cnLOH", ident.2 = "cnLOH",  group.by="status_cnLOH_BAF")

non_mos_vs_mos
dim(non_mos_vs_mos)

DE.mos.snRNAseq <- FindMarkers(object = merged_hepatocytes, logfc.threshold = 0.1, 
                              ident.1 = "cnLOH", ident.2 = "no cnLOH",  group.by="status_cnLOH_BAF")

DE.mos.snRNAseq
dim(DE.mos.snRNAseq)

#write.table(DE.mos.snRNAseq, quote = F,row.names = TRUE, paste(merged.directory, "DE.mos.snRNAseq.txt",sep=""),  sep="\t")

de_snRNAseq.pos <- mos_vs_non_mos %>%
  dplyr::filter(p_val_adj <=0.01)
dim(de_snRNAseq.pos)

de_snRNAseq.pos = rownames(de_snRNAseq.pos)

de_snRNAseq.neg <- non_mos_vs_mos %>%
  dplyr::filter(p_val_adj <=0.01)

de_snRNAseq.neg = rownames(de_snRNAseq.neg)

# Positive enrichment
de_snRNAseq.pos <- de_snRNAseq.pos %>% intersect(rownames(all.merged))
all.merged$de_snRNAseq.pos <- apply(GetAssayData(object = all.merged, slot = "data")[de_snRNAseq.pos, ], 2, mean) 

p<-FeaturePlot(all.merged, reduction="umap", features="de_snRNAseq.pos",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18))  + ggtitle("")
p
ggsave(paste0(merged.directory, "all.merged_de_snRNAseq.pos.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)


# Load signature from bulkRNAseq and spatial transcriptomics 
# Differentially expressed genes in RNAseq  
de_RNAseq = read.delim("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/Limma/All_NT_HB_without_3370_age_adjusted/Mosaic_vs_Non_Mos_HB_without_CHC3370N_all_ages_age_adjusted.txt", sep= "\t")

de_RNAseq.pos = de_RNAseq %>%
  dplyr::filter(P.Value <=0.01,
         logFC>0)
de_RNAseq.pos = de_RNAseq.pos$Gene

de_RNAseq.neg = de_RNAseq %>%
  dplyr::filter(P.Value <=0.01,
         logFC<0)
de_RNAseq.neg = de_RNAseq.neg$Gene

# Negative enrichment
de_RNAseq.neg <- de_RNAseq.neg %>% intersect(rownames(all.merged))
all.merged$de_RNAseq.neg <- apply(GetAssayData(object = all.merged, slot = "data")[de_RNAseq.neg, ], 2, mean) 

# Positive enrichment
de_RNAseq.pos <- de_RNAseq.pos %>% intersect(rownames(all.merged))
all.merged$de_RNAseq.pos <- apply(GetAssayData(object = all.merged, slot = "data")[de_RNAseq.pos, ], 2, mean) 

p<-FeaturePlot(all.merged, reduction="umap", features="de_RNAseq.pos",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18))  + ggtitle("")
p
ggsave(paste0(merged.directory, "all.merged_de_bulkRNAseq.pos.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)


# Negative enrichment
de_RNAseq.neg <- de_RNAseq.neg %>% intersect(rownames(all.merged))
all.merged$de_RNAseq.neg <- apply(GetAssayData(object = all.merged, slot = "data")[de_RNAseq.neg, ], 2, mean) 

p<-FeaturePlot(all.merged, reduction="umap", features="de_RNAseq.neg",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18))  + ggtitle("")
p
ggsave(paste0(merged.directory, "all.merged_de_bulkRNAseq.neg.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)


# differentially expressed genes in visium 

de_spatial <- read.delim("D:/Dropbox/11p15.5 mosaicism/Visium/GSEA/Mosaic_vs_non_mosaic_all_patients.txt", sep="\t")

de_spatial.pos = de_spatial %>%
  dplyr::filter(p_val_adj <=0.01,
                avg_log2FC>0)
de_spatial.pos = rownames(de_spatial.pos)

de_spatial.neg = de_spatial %>%
  dplyr::filter(p_val_adj <=0.01,
                avg_log2FC<0)
de_spatial.neg = rownames(de_spatial.neg)

# Negative enrichment
de_spatial.neg <- de_spatial.neg %>% intersect(rownames(all.merged))
all.merged$de_spatial.neg <- apply(GetAssayData(object = all.merged, slot = "data")[de_spatial.neg, ], 2, mean) 

# Positive enrichment
de_spatial.pos <- de_spatial.pos %>% intersect(rownames(all.merged))
all.merged$de_spatial.pos <- apply(GetAssayData(object = all.merged, slot = "data")[de_spatial.pos, ], 2, mean) 

p<-FeaturePlot(all.merged, reduction="umap", features="de_spatial.pos",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18))  + ggtitle("")
p
ggsave(paste0(merged.directory, "all.merged_de_spatial.pos.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)

p<-FeaturePlot(all.merged, reduction="umap", features="de_spatial.neg",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18))  + ggtitle("")
p 
ggsave(paste0(merged.directory, "all.merged_de_spatial.neg.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)



# Venn diagram 

x <- list(
  Spatial = de_spatial.pos,
  snRNAseq = de_snRNAseq.pos,
  Bulk.RNAseq = de_RNAseq.pos
)


display_venn(
  x, 
  category.names = c("Spatial" , "snRNAseq", "Bulk.RNAseq"),
  fill = c("Yellow", "#009E73", "#cfe2f3"),
  lwd = 2,
  lty = 'blank',
  cex = 2,
  fontface = "italic",
  # Set names
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  margin=0.2
)

common.genes.up <- intersect(de_spatial.pos, de_snRNAseq.pos)
length(common.genes.up)

intersect(common.genes.up, de_RNAseq.pos)

# Source Figure 5e-f

write.table(de_spatial.pos, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/de_spatial.pos.txt", sep="\t")
write.table(de_snRNAseq.pos, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/de_snRNAseq.pos.txt", sep="\t")
write.table(de_RNAseq.pos, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/de_RNAseq.pos.txt", sep="\t")
write.table(common.genes.up, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/common_snRNAseq_spatial.txt", sep="\t")


# Plot 60 genes signature on UMAP
all.merged$mos.sig <- apply(GetAssayData(object = all.merged, slot = "data")[common.genes.up, ], 2, mean) 

p<-FeaturePlot(all.merged, reduction="umap", features="mos.sig",  pt.size = 0.1) + my_scale_exp + 
  theme(
    axis.line = element_line(linewidth = 1),
    title = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title = element_text(size=18))  + ggtitle("")
p

ggsave(paste0(merged.directory, "mos.sig.spatial.snRNAseq.png"), plot = p, device = "png", width =12, height = 10, units = "cm", dpi = 320)


#common.genes.down <- intersect(de_spatial.neg, de_snRNAseq.neg)
#common.genes <- c(common.genes.up, common.genes.down)

# Expression per group cn_LOH 
merged_hepatocytes.sce <- as.SingleCellExperiment(merged_hepatocytes)

celltype_mean <- aggregateAcrossCells(merged_hepatocytes.sce,  
                                      ids = merged_hepatocytes$status_cnLOH_BAF, 
                                      statistics = "mean",
                                      use.assay.type = "counts", 
                                      subset.row = common.genes.up)



p <- dittoHeatmap(celltype_mean,
             assay = "counts", cluster_cols = TRUE, 
             scale = "row",  # "none" no scaling otherwise "row"
             heatmap.colors = viridis(n=10),
             fontsize = 9, 
           #  cutree_rows = 2 ,
             treeheight_row = 10,
             treeheight_col = 10,
             annot.by = c("status_cnLOH_BAF" , "ncells") ,
             annotation_colors = list(status_cnLOH_BAF = cols_status_cn_LOH_BAF,
                                      ncells = plasma(100))
) 

p
ggsave(paste0(merged.directory, "heatmap.venn.png"), plot = p, device = "png", width =12, height =17, units = "cm", dpi = 320)


# Bcat target genes 

p<-VlnPlot(
  object = all.merged,
  features = c("TBX3"),
  group.by = "status_cnLOH_BAF", log = T,
  cols= cols_status_cn_LOH_BAF,
  ncol = 1,
  pt.size = 0
) + theme(axis.title.x = element_blank()) & NoLegend()


p

p<-VlnPlot(
  object = all.merged,
  features = c("GLUL"),
  group.by = "status_cnLOH_BAF", log = T,
  cols= cols_status_cn_LOH_BAF,
  ncol = 1,
  pt.size = 0
) + theme(axis.title.x = element_blank()) & NoLegend()


p

p<-VlnPlot(
  object = all.merged,
  features = c("LGR5"),
  group.by = "status_cnLOH_BAF", log = T,
  cols= cols_status_cn_LOH_BAF,
  ncol = 1,
  pt.size = 0
) + theme(axis.title.x = element_blank()) & NoLegend()


p

p<-VlnPlot(
  object = all.merged,
  features = c("AXIN2"),
  group.by = "status_cnLOH_BAF", log = T,
  cols= cols_status_cn_LOH_BAF,
  ncol = 1,
  pt.size = 0
) + theme(axis.title.x = element_blank()) & NoLegend()


p



    # 6.3.4  Gene set enrichment -----------------------------------------------

toppfun <- read.delim("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/merged/Toppfun_58_markers_enrichment.txt", sep="\t")
toppfun = toppfun %>%
  dplyr::filter(Category %in% c("Coexpression", "Disease", "GO: Biological Process", "GO: Cellular Component", "GO: Molecular Function",
                                "Human Phenotype", "Pathway")) %>%
  arrange(q.value.Bonferroni)

pathways <- c("DESCARTES_FETAL_LIVER_HEPATOBLASTS", 
              "DESCARTES_ORGANOGENESIS_HEPATOCYTES", 
              "blood microparticle", 
              "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
              "lipid metabolic process", 
              "steroid metabolic process",
              "HALLMARK_COAGULATION", 
              "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES", 
              "response to xenobiotic stimulus")

Enrich.plot = toppfun %>%
  dplyr::filter(Name %in% pathways)  %>%
  mutate(log.Bonferroni = -log(q.value.Bonferroni)) %>%
  arrange(log(log.Bonferroni))

Enrich.plot = Enrich.plot %>%
  dplyr::filter(!duplicated(Name))

Enrich.plot$Name = factor(Enrich.plot$Name, levels=Enrich.plot$Name)



ggplot(Enrich.plot, aes(x=log.Bonferroni, y=Name, fill=log.Bonferroni, size=log.Bonferroni)) +
  geom_point(alpha=0.8, shape=19, color="red") +
  scale_size("log.Bonferroni") + 
  theme_classic() + 
  theme(#axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=18),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(1,1,1,1, "cm")) 

# Clean Seurat object for publication ------------------------------------------

all.merged <- readRDS("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/merged/all.merged.rds")


all.merged@meta.data$SCT_snn_res.0.1 <- NULL
all.merged@meta.data$seurat_clusters <- NULL
all.merged@meta.data$clusters_0_1 <- NULL
all.merged@meta.data$cell.type.old <- NULL
all.merged@meta.data$scDblFinder.class <- NULL
all.merged@meta.data$clusters_0_1_postQC <- NULL
all.merged@meta.data$clusters_0_2_postQC <- NULL
all.merged@meta.data$SCT_snn_res.0.05 <- NULL
all.merged@meta.data$clusters_0_05_postQC <- NULL
all.merged@meta.data$SCT_snn_res.0.3 <- NULL
all.merged@meta.data$clusters_0_3_postQC <- NULL
all.merged@meta.data$cell.type <- NULL
all.merged@meta.data$Phase <- NULL
all.merged@meta.data$SCT_snn_res.0.4 <- NULL
all.merged@meta.data$clusters_0_4_postQC <- NULL
all.merged@meta.data$ReadGroup <- NULL
all.merged@meta.data$pval_if_no_cnLOH <- NULL
all.merged@meta.data$qval_if_no_cnLOH <- NULL
all.merged@meta.data$status_cn_LOH_numbat <- NULL
all.merged@meta.data$proba_cn_LOH_numbat <- NULL
all.merged@meta.data$pval_if_cnLOH <- NULL
all.merged@meta.data$clusters_0_2 <- NULL
all.merged@meta.data$clusters_0_3 <- NULL
all.merged@meta.data$SCT_snn_res.0.5 <- NULL
all.merged@meta.data$clusters_0_5_postQC <- NULL
all.merged@meta.data$SCT_snn_res.0.6 <- NULL
all.merged@meta.data$clusters_0_6_postQC <- NULL
all.merged@meta.data$SCT_snn_res.0.8 <- NULL
all.merged@meta.data$clusters_0_8_postQC <- NULL
all.merged@meta.data$SCT_snn_res.0.9 <- NULL
all.merged@meta.data$clusters_0_9_postQC <- NULL
all.merged@meta.data$total <- NULL
all.merged@meta.data$paternal_reads <- NULL
all.merged@meta.data$maternal_reads <- NULL

names(all.merged@meta.data)

all.merged@meta.data$cell.population <- all.merged@meta.data$final.cell.type

table(all.merged@meta.data$ce)

saveRDS(all.merged, file=paste0(merged.directory, "snRNAseq_merged_seurat.rds"))
