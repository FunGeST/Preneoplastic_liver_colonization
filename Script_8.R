### 1. Load packages -----------------------------------------------------------

set.seed(1234)

library(Seurat)
library(ggplot2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(viridis)
library(Signac)
library(SoupX)

### 1. CHC2996N ----------------------------------------------------------------

sample_name <- "CHC2996N"
print(sample_name)

# Create useful directories paths-----------------------------------------------

data_dir <- paste0("//10.93.23.19/hepato partage/GENOMIC_DATA/scRNAseq/Project_HB_tumors/cellranger_outputs/PJ2302056_cellranger_v5/", sample_name, "/filtered_feature_bc_matrix/")
samp.directory <- paste0("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/", sample_name, "/")

output.directory <- paste0(samp.directory, "preprocess.QC/")
if(!dir.exists(output.directory)){dir.create(output.directory)}

plots.dir.preQC <- paste0(output.directory, "plots/preQC/")
if(!dir.exists(plots.dir.preQC)){dir.create(plots.dir.preQC)}

plots.dir.postQC <- paste0(output.directory, "plots/postQC/")
if(!dir.exists(plots.dir.postQC)){dir.create(plots.dir.postQC)}

# load the RNA data-------------------------------------------------------------
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz , atac_fragments.tsv.gz.tbi 
data <- Read10X(data.dir = data_dir) # Load expression data 

# get gene annotations for hg38-------------------------------------------------
hg38 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(hg38) <- "hg38"
hg38 <- renameSeqlevels(hg38, mapSeqlevels(seqlevels(hg38), "UCSC"))

# create a Seurat object containing the RNA data--------------------------------
sample <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = sample_name
)


# UMAP and clustering  ---------------------------------------------------------
# colors and parameters --------------------------------------------------------
my_scale_exp = scale_colour_gradientn(colors = c( "gray70", "#ffff00",  "red", "#cc0000"))

# Normalization with SCTransform -----------------------------------------------
DefaultAssay(sample) <- "RNA"
sample <- SCTransform(sample, verbose = TRUE, variable.features.n=3000) # sctransform = par d?faut 3000 nmostvar et pas 2000 comme lognormalize
sample <- ScaleData(sample, assay = "RNA")
DefaultAssay(sample) <- "SCT"

# PCA on gene expression--------------------------------------------------------
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=50) # 50 components

# Clustering --------------------------------------------------------------------
nbcomp=30
sample <- RunUMAP(sample, dims=1:nbcomp, verbose=F)
resolution_for_clusters <- 0.1  # change resolution to change nb of clusters determined

sample <- FindNeighbors(sample, dims = 1:30)
sample <- FindClusters(sample, resolution = 0.1)
sample$clusters_0_1 <- sample$seurat_clusters

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "clusters_0_1", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_1_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

# add mitochondrial percentage -------------------------------------------------
DefaultAssay(sample) <- "RNA"
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern="^MT-") 


# visualize QC in UMAP ---------------------------------------------------------
sample$nCount_RNA_log <- log10(sample$nCount_RNA)
sample$percent.mt_log <- log2(sample$percent.mt)

p<-FeaturePlot(sample, reduction="umap", features="percent.mt_log",  pt.size = 0.1) + my_scale_exp & NoLegend()
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_percent_mt_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nCount_RNA_log",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nCount_RNA_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Show specific markers to annotate cell types ---------------------------------

DefaultAssay(sample) <- "SCT"

p<-FeaturePlot(sample, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CD163_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CD163_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT19_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT7_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CYP2E1_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HAL_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HOXB3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HOXB3_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_VIM_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Annotate cell types ----------------------------------------------------------

sample@meta.data <- sample@meta.data%>%
  mutate(cell.type.old = case_when(clusters_0_1 =="0" ~ "High ambient RNA and MT genes",
                               clusters_0_1 == "1" ~ "Hepatocytes",
                               clusters_0_1 == "2" ~ "SAA2-CRP high hepatocytes/fibroblasts?",
                               clusters_0_1 == "3" ~ "Endothelial cells",
                               clusters_0_1 == "4" ~ "Macrophages",
                               clusters_0_1 == "5" ~ "Biliary epithelial cells",
                               clusters_0_1 == "6" ~ "HOX+ Progenitor",
                               clusters_0_1 == "7" ~ "Tumor cells"))

# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "High ambient RNA and MT genes" = "#36e22e", 
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages" = "#357ee5", 
                   "Biliary epithelial cells" = "#e53561", 
                   "SAA2-CRP high hepatocytes/fibroblasts?" = "blue",
                   "HOX+ Progenitor" = "#ff9f2e",
                   "Tumor cells" = "#0c0994")

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "cell.type.old", reduction = "umap", cols = cell.type.cols) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_1_cell_type_old_preQC.png"), plot = p, device = "png", width = 17, height = 10, units = "cm", dpi = 320)


# Differentially expressed genes  ----------------------------------------------

sample.markers <- FindAllMarkers(object = sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
png(paste(plots.dir.preQC, sample_name, "_DoHeatmap_top10markers_preQC.png",sep=""),width=22,height=15,units="cm",res=1000)
print(DoHeatmap(sample, features = top10$gene,cells=colnames(sample), group.colors = cell.type.cols, group.by = "cell.type.old", size = 2) + 
        scale_fill_viridis(option="mako") +
        theme(text = element_text(size = 6)) +
        NoLegend())
dev.off()

write.table(sample.markers, quote = F,row.names = FALSE, paste(plots.dir.preQC, sample_name, "_allmarkers_preQC.txt",sep=""),  sep="\t")

# Perform quality controls  ---------------------------------------------------

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("percent.mt"),
  cols = cell.type.cols,
  pt.size = 0
) + geom_hline(yintercept=5)
p

ggsave(paste0(plots.dir.preQC, sample_name, "_MT_percentage_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("percent.mt"),log = T,
  cols = cell.type.cols,
  pt.size = 0
) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_MT_percentage_log_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nCount_RNA_log"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_nCount_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nFeature_RNA"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0)  & NoLegend()

p

ggsave(paste0(plots.dir.preQC, sample_name, "_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


#  Filter out low quality hepatocytes ------------------------------------------

sample.postqc <- subset(sample, subset = cell.type.old =="High ambient RNA and MT genes", invert=TRUE)

sample.postqc <- subset(
  x = sample.postqc,
  subset=nFeature_RNA>1000 & percent.mt<5)

print(length(colnames(sample)))
print(length(colnames(sample.postqc)))

# Remove doublets predicted with scrublet---------------------------------------

tmp_sobj <- Seurat::CreateSeuratObject(counts=Seurat::GetAssayData(sample.postqc,assay='RNA',slot="counts"))
tmp_sobj@assays$RNA@data <- Seurat::GetAssayData(sample.postqc,assay='RNA',slot="data")
tmp_sobj$orig.ident <- sample.postqc$orig.ident
sce <- Seurat::as.SingleCellExperiment(tmp_sobj)

sample.postqc$scDblFinder.class <- Seurat::as.Seurat(scDblFinder::scDblFinder(sce,samples=tmp_sobj$orig.ident))$scDblFinder.class
sample.postqc$scDblFinder.class <- unname(sample.postqc$scDblFinder.class == "doublet") # turn into TRUE/FALSE

p<-DimPlot(sample.postqc, group.by = "scDblFinder.class", cols = c("gray80","red"), reduction = "umap", raster=F) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "scDblFinder.class_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

print(length(colnames(sample.postqc)))

sample.postqc <- subset(
  x = sample.postqc,
  subset = scDblFinder.class=="FALSE")

print(length(colnames(sample.postqc)))

# Clustering --------------------------------------------------------------------
DefaultAssay(sample.postqc)<- "SCT"

sample.postqc <- RunPCA(sample.postqc, features = VariableFeatures(sample.postqc), npcs=30) 
sample.postqc <- RunUMAP(sample.postqc, dims=1:30, verbose=T)
sample.postqc <- FindNeighbors(sample.postqc, dims = 1:30)

sample.postqc <- FindClusters(sample.postqc, resolution = 0.1)
sample.postqc$clusters_0_1_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.2)
sample.postqc$clusters_0_2_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.05)
sample.postqc$clusters_0_05_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.3)
sample.postqc$clusters_0_3_postQC <- sample.postqc$seurat_clusters

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample.postqc, group.by = "clusters_0_05_postQC", reduction = "umap", pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_05_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320,)

p<-DimPlot(sample.postqc, group.by = "clusters_0_1_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_1_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320,)

p<-DimPlot(sample.postqc, group.by = "clusters_0_2_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_2_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_3_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_3_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="percent.mt",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_percent_mt_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_nFeature_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# visualize specific gene expression -------------------------------------------
DefaultAssay(sample.postqc) <- "SCT"

# endothelial cells-------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD34",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD34_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VWF",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VWF_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANGPT2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANGPT2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="PTGDS",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PTGDS_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FLT1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGFBP3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGFBP3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CDH5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CDH5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC14A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC14A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="SLCO2A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SLCO2A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="RELN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_RELN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC4M",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC4M_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC1B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC1B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FCN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FCN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OIT3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OIT3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hepatic stellate cells -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COLEC11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COLEC11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LRAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LRAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NGFR",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NGFR_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# 11p-related genes ------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGF2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_H19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# progenitor markers -----------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GPC3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GPC3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="DLK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DLK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="EPCAM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_EPCAM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KIT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KIT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AFP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AFP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LIN28B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LIN28B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="THY1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_THY1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hematopoiesis ----------------------------------------------------------------
p<-FeaturePlot(sample.postqc, reduction="umap", features="HBB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HBM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Fibroblasts ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL6A3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL6A3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OGN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OGN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="POSTN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_POSTN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VIM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FN1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FN1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Macrophages-------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD163_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD68",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD68_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD80",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD80_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD74",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD74_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HLA-DRA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HLA-DRA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LYZ",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LYZ_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# B lymphocytes ----------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="JCHAIN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_JCHAIN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="IGKC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGKC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD79A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD79A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MZB1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MZB1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# T lymphocytes ----------------------------------------------------------------


p<-FeaturePlot(sample.postqc, reduction="umap", features="CCL5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CCL5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMK",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMK_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD3D",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD3D_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NKG7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NKG7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# NK cells ---------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GNLY",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GNLY_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD247",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD247_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD160",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD160_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NCR1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NCR1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# VSMC -------------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TAGLN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TAGLN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MYH11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MYH11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ACTA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ACTA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Lipogenesis ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOC1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOC1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ALB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ALB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="F2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_F2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Periportal markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="PCK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PCK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CPS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CPS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HAL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Perivenous markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP2E1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="GSTM3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GSTM3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NOTUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NOTUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FASN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FASN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="DPP4",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DPP4_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Wnt pathway activation -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TBX3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TBX3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AXIN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AXIN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LGR5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LGR5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GLUL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GLUL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Biliary epithelial cells------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SOX9",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SOX9_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SPP1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SPP1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------


p<-FeaturePlot(sample.postqc, reduction="umap", features="ITGB2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ITGB2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SELP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SELP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Annotate cell types ----------------------------------------------------------

sample.postqc@meta.data <- sample.postqc@meta.data%>%
  mutate(cell.type = case_when(clusters_0_2_postQC =="0" ~ "SAA2-CRP high hepatocytes/fibroblasts?",
                               clusters_0_2_postQC == "1" ~ "Hepatocytes",
                               clusters_0_2_postQC == "2" ~ "Hepatocytes",
                               clusters_0_2_postQC == "3" ~ "Endothelial cells",
                               clusters_0_2_postQC == "4" ~ "Biliary epithelial cells",
                               clusters_0_2_postQC == "5" ~ "Macrophages - Kupffer cells",
                               clusters_0_2_postQC == "6" ~ "Tumor cells",
                               clusters_0_2_postQC == "7" ~ "HOX+ Progenitor"))


# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Biliary epithelial cells" = "#e53561", 
                   "SAA2-CRP high hepatocytes/fibroblasts?" = "blue",
                   "HOX+ Progenitor" = "#ff9f2e",
                   "Tumor cells" = "#0c0994")



# Differentially expressed genes  ----------------------------------------------

Idents(sample.postqc) <- sample.postqc@meta.data$cell.type

sample.markers <- FindAllMarkers(object = sample.postqc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1, group.by="cell.type")
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

png(paste(plots.dir.postQC, sample_name, "_DoHeatmap_top10markers_postQC.png",sep=""),width=22,height=15,units="cm",res=1000)
print(DoHeatmap(sample.postqc, features = top10$gene, cells=colnames(sample.postqc), group.colors = cell.type.cols, group.by = "cell.type", size = 2) + 
        scale_fill_viridis(option="mako") +
        theme(text = element_text(size = 6)) +
        NoLegend())
dev.off()

write.table(sample.markers, quote = F,row.names = FALSE, paste(plots.dir.postQC, sample_name, "_allmarkers_postQC.txt",sep=""),  sep="\t")


DefaultAssay(sample.postqc) <- "SCT"

# visualize QC after filtering--------------------------------------------------

p<-VlnPlot(
  object = sample.postqc,
  features = c("nCount_RNA","nFeature_RNA", "percent.mt"),
  group.by = "cell.type",
  cols=cell.type.cols,
  ncol = 3,
  pt.size = 0
) & NoLegend()


p
ggsave(paste0(plots.dir.postQC, sample_name, "_postQC_plots_postQC.png"), plot = p, device = "png", width = 22, height = 14, units = "cm", dpi = 320)


# Cell cycle scoring------------------------------------------------------------

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

DefaultAssay(sample.postqc) <- "RNA"
sample.postqc <- NormalizeData(sample.postqc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

length(s.genes)
setdiff(s.genes, rownames(sample.postqc))
s.genes <- s.genes[which(s.genes %in% rownames(sample.postqc))]
length(s.genes)

length(g2m.genes)
setdiff(g2m.genes, rownames(sample.postqc))
g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(sample.postqc))]
length(g2m.genes)

sample.postqc <- CellCycleScoring(sample.postqc, s.features = s.genes, g2m.features = g2m.genes)

# view cell cycle scores and phase assignments----------------------------------
head(sample.postqc[[]])

p<-FeaturePlot(sample.postqc, reduction="umap", features="G2M.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_G2M.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="S.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_S.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-DimPlot(sample.postqc, group.by = "cell.type", reduction = "umap", cols = cell.type.cols, pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_cell_type.png"), plot = p, device = "png", width = 18, height = 10, units = "cm", dpi = 320,)


saveRDS(sample.postqc, file=paste0(output.directory, sample_name, "postQC.rds"))
saveRDS(sample, file=paste0(output.directory, sample_name, "preQC.rds"))


### 2. CHC3559N ----------------------------------------------------------------

sample_name <- "CHC3559N"
print(sample_name)

# Create useful directories paths-----------------------------------------------

data_dir <- paste0("//10.93.23.19/hepato partage/GENOMIC_DATA/scRNAseq/Project_HB_tumors/cellranger_outputs/PJ2302056_cellranger_v5/", sample_name, "/filtered_feature_bc_matrix/")
samp.directory <- paste0("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/", sample_name, "/")

output.directory <- paste0(samp.directory, "preprocess.QC/")
if(!dir.exists(output.directory)){dir.create(output.directory)}

plots.dir.preQC <- paste0(output.directory, "plots/preQC/")
if(!dir.exists(plots.dir.preQC)){dir.create(plots.dir.preQC)}

plots.dir.postQC <- paste0(output.directory, "plots/postQC/")
if(!dir.exists(plots.dir.postQC)){dir.create(plots.dir.postQC)}

# load the RNA data-------------------------------------------------------------
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz , atac_fragments.tsv.gz.tbi 
data <- Read10X(data.dir = data_dir) # Load expression data 

# get gene annotations for hg38-------------------------------------------------
hg38 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(hg38) <- "hg38"
hg38 <- renameSeqlevels(hg38, mapSeqlevels(seqlevels(hg38), "UCSC"))

# create a Seurat object containing the RNA data--------------------------------
sample <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = sample_name
)

# UMAP and clustering  ---------------------------------------------------------
# colors and parameters --------------------------------------------------------
my_scale_exp = scale_colour_gradientn(colors = c( "gray70", "#ffff00",  "red", "#cc0000"))

# Normalization with SCTransform -----------------------------------------------
DefaultAssay(sample) <- "RNA"
sample <- SCTransform(sample, verbose = TRUE, variable.features.n=3000) # sctransform = par d?faut 3000 nmostvar et pas 2000 comme lognormalize
sample <- ScaleData(sample, assay = "RNA")
DefaultAssay(sample) <- "SCT"

# PCA on gene expression--------------------------------------------------------
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=50) # 50 components

# Clustering --------------------------------------------------------------------
nbcomp=30
sample <- RunUMAP(sample, dims=1:nbcomp, verbose=F)

sample <- FindNeighbors(sample, dims = 1:30, k.param = 20)
sample <- FindClusters(sample, resolution = 0.1)
sample$clusters_0_1 <- sample$seurat_clusters

sample <- FindClusters(sample, resolution = 0.2)
sample$clusters_0_2 <- sample$seurat_clusters

sample <- FindClusters(sample, resolution = 0.3)
sample$clusters_0_3 <- sample$seurat_clusters

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "clusters_0_1", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_1_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample, group.by = "clusters_0_2", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_2_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample, group.by = "clusters_0_3", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_3_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

# add mitochondrial percentage -------------------------------------------------
DefaultAssay(sample) <- "RNA"
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern="^MT-") 

# visualize QC in UMAP ---------------------------------------------------------
sample$nCount_RNA_log <- log10(sample$nCount_RNA)
sample$percent.mt_log <- log2(sample$percent.mt)

p<-FeaturePlot(sample, reduction="umap", features="percent.mt_log",  pt.size = 0.1) + my_scale_exp & NoLegend()
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_percent_mt_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nCount_RNA_log",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nCount_RNA_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Show specific markers to annotate cell types ---------------------------------

DefaultAssay(sample) <- "SCT"

p<-FeaturePlot(sample, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CD163_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_FLT1_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT19_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT7_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CYP2E1_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HAL_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HOXB3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HOXB3_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_VIM_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Annotate cell types ----------------------------------------------------------

sample@meta.data <- sample@meta.data%>%
  mutate(cell.type.old = case_when(clusters_0_3 =="0" ~ "High ambient RNA and MT genes",
                               clusters_0_3 == "1" ~ "Periportal bifurcations cells ?",
                               clusters_0_3 == "2" ~ "Hepatocytes",
                               clusters_0_3 == "3" ~ "Endothelial cells",
                               clusters_0_3 == "4" ~ "Macrophages",
                               clusters_0_3 == "5" ~ "Fibroblasts",
                               clusters_0_3 == "6" ~ "HOX+ Progenitor",
                               clusters_0_3 == "7" ~ "Biliary epithelial cells"))

# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "High ambient RNA and MT genes" = "#36e22e",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages" = "#357ee5", 
                   "Fibroblasts" = "#8f35e5", 
                   "Biliary epithelial cells" = "#e53561", 
                   "HOX+ Progenitor" = "#ff9f2e", 
                   "Periportal bifurcations cells ?" = "forestgreen")

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "cell.type.old", reduction = "umap", cols = cell.type.cols) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_3_cell_type_old_preQC.png"), plot = p, device = "png", width = 17, height = 10, units = "cm", dpi = 320)


# Differentially expressed genes  ----------------------------------------------


sample.markers <- FindAllMarkers(object = sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

png(paste(plots.dir.preQC, sample_name, "_DoHeatmap_top10markers_preQC.png",sep=""),width=22,height=15,units="cm",res=1000)
print(DoHeatmap(sample, features = top10$gene,cells=colnames(sample), group.colors = cell.type.cols, group.by = "cell.type.old", size = 2) + 
        scale_fill_viridis(option="mako") +
        theme(text = element_text(size = 6)) +
        NoLegend())
dev.off()

write.table(sample.markers, quote = F,row.names = FALSE, paste(plots.dir.preQC, sample_name, "_allmarkers_preQC.txt",sep=""),  sep="\t")


# Perform quality controls  ---------------------------------------------------

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("percent.mt"),
  cols = cell.type.cols,
  pt.size = 0
) + geom_hline(yintercept=5)
p

ggsave(paste0(plots.dir.preQC, sample_name, "_MT_percentage_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("percent.mt"),log = T,
  cols = cell.type.cols,
  pt.size = 0
) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_MT_percentage_log_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nCount_RNA_log"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_nCount_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nFeature_RNA"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0)  & NoLegend()

p

ggsave(paste0(plots.dir.preQC, sample_name, "_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


# Filter out low quality hepatocytes ------------------------------------------

sample.postqc <- subset(sample, subset = cell.type.old =="High ambient RNA and MT genes", invert=TRUE)

sample.postqc <- subset(
  x = sample.postqc,
  subset=nFeature_RNA>1000 & percent.mt<5)

print(length(colnames(sample)))
print(length(colnames(sample.postqc)))

# Remove doublets predicted with scrublet---------------------------------------

tmp_sobj <- Seurat::CreateSeuratObject(counts=Seurat::GetAssayData(sample.postqc,assay='RNA',slot="counts"))
tmp_sobj@assays$RNA@data <- Seurat::GetAssayData(sample.postqc,assay='RNA',slot="data")
tmp_sobj$orig.ident <- sample.postqc$orig.ident
sce <- Seurat::as.SingleCellExperiment(tmp_sobj)

sample.postqc$scDblFinder.class <- Seurat::as.Seurat(scDblFinder::scDblFinder(sce,samples=tmp_sobj$orig.ident))$scDblFinder.class
sample.postqc$scDblFinder.class <- unname(sample.postqc$scDblFinder.class == "doublet") # turn into TRUE/FALSE

p<-DimPlot(sample.postqc, group.by = "scDblFinder.class", cols = c("gray80","red"), reduction = "umap", raster=F) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "scDblFinder.class_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

print(length(colnames(sample.postqc)))

sample.postqc <- subset(
  x = sample.postqc,
  subset = scDblFinder.class=="FALSE")

print(length(colnames(sample.postqc)))

# Clustering --------------------------------------------------------------------

ElbowPlot(sample.postqc, ndims=15)
sample.postqc <- RunPCA(sample.postqc, features = VariableFeatures(sample.postqc), npcs=30) 
sample.postqc <- RunUMAP(sample.postqc, dims=1:20, verbose=T, reduction = "pca") 
sample.postqc <- FindNeighbors(sample.postqc, dims = 1:10, k.param = 10) 

sample.postqc <- FindClusters(sample.postqc, resolution = 0.1)
sample.postqc$clusters_0_1_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.2)
sample.postqc$clusters_0_2_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.3)
sample.postqc$clusters_0_3_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.5)
sample.postqc$clusters_0_5_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.6)
sample.postqc$clusters_0_6_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.8)
sample.postqc$clusters_0_8_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.9)
sample.postqc$clusters_0_9_postQC <- sample.postqc$seurat_clusters

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample.postqc, group.by = "clusters_0_1_postQC", reduction = "umap", pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_1_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320,)

p<-DimPlot(sample.postqc, group.by = "clusters_0_2_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_2_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_3_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_3_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_5_postQC", reduction = "umap", pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_5_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_6_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_6_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_8_postQC", reduction = "umap",  pt.size = 0.1) 
p + scale_color_viridis_d(option="turbo")
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_8_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)


p<-DimPlot(sample.postqc, group.by = "clusters_0_9_postQC", reduction = "umap",  pt.size = 0.1) 
p + scale_color_viridis_d(option="turbo")

ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_9_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="percent.mt",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_percent_mt_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_nFeature_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# visualize specific gene expression -------------------------------------------
DefaultAssay(sample.postqc) <- "SCT"

# endothelial cells-------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD34",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD34_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VWF",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VWF_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANGPT2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANGPT2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="PTGDS",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PTGDS_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FLT1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGFBP3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGFBP3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CDH5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CDH5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC14A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC14A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="SLCO2A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SLCO2A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="RELN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_RELN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC4M",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC4M_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC1B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC1B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FCN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FCN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OIT3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OIT3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hepatic stellate cells -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COLEC11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COLEC11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LRAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LRAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NGFR",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NGFR_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HGF",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HGF_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# 11p-related genes ------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGF2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_H19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# progenitor markers -----------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GPC3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GPC3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="DLK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DLK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="EPCAM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_EPCAM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KIT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KIT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AFP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AFP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LIN28B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LIN28B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="THY1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_THY1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hematopoiesis ----------------------------------------------------------------
p<-FeaturePlot(sample.postqc, reduction="umap", features="HBB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HBM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Fibroblasts ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL6A3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL6A3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OGN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OGN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="POSTN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_POSTN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VIM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FN1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FN1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Macrophages-------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD163_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD68",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD68_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD80",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD80_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD74",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD74_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HLA-DRA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HLA-DRA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LYZ",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LYZ_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# B lymphocytes ----------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="JCHAIN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_JCHAIN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="IGKC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGKC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD79A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD79A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MZB1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MZB1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# T lymphocytes ----------------------------------------------------------------


p<-FeaturePlot(sample.postqc, reduction="umap", features="CCL5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CCL5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMK",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMK_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD3D",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD3D_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NKG7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NKG7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# NK cells ---------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GNLY",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GNLY_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD247",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD247_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD160",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD160_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NCR1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NCR1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# VSMC -------------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TAGLN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TAGLN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MYH11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MYH11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)
KAZNmacrop
p<-FeaturePlot(sample.postqc, reduction="umap", features="ACTA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ACTA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Lipogenesis ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOC1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOC1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ALB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ALB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="F2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_F2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Periportal markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="PCK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PCK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CPS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CPS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HAL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Perivenous markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP2E1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="GSTM3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GSTM3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NOTUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NOTUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FASN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FASN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="DPP4",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DPP4_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Wnt pathway activation -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TBX3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TBX3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AXIN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AXIN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LGR5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LGR5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GLUL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GLUL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Biliary epithelial cells------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SOX9",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SOX9_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SPP1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SPP1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------


p<-FeaturePlot(sample.postqc, reduction="umap", features="ITGB2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ITGB2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SELP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SELP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Annotate cell types ----------------------------------------------------------

BEC <- sample.postqc[["umap"]]@cell.embeddings[sample.postqc[["umap"]]@cell.embeddings[,1] > -6.5 & 
                                          sample.postqc[["umap"]]@cell.embeddings[,1] < -3 &
                                          sample.postqc[["umap"]]@cell.embeddings[,2] < 1,]

BEC = rownames(BEC)

sample.postqc@meta.data <- sample.postqc@meta.data %>%
  mutate(cell.type = case_when(clusters_0_5_postQC =="0" ~ "Periportal bifurcations cells ?",
                               clusters_0_5_postQC == "1" ~ "Periportal bifurcations cells ?",
                               clusters_0_5_postQC == "2" ~ "Hepatocytes",
                               clusters_0_5_postQC == "3" ~ "Hepatocytes",
                               colnames(sample.postqc) %in% BEC ~ "Biliary epithelial cells",
                               clusters_0_5_postQC == "4" ~ "Endothelial cells",
                               clusters_0_5_postQC == "5" ~ "Hepatocytes",
                               clusters_0_5_postQC == "6" ~ "Macrophages - Kupffer cells",
                               clusters_0_5_postQC == "7" ~ "Fibroblasts - HSCs - VSMCs",
                               clusters_0_5_postQC == "8" ~ "HOX+ Progenitor"))


DefaultAssay(sample.postqc) <- "SCT"


cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Fibroblasts - HSCs - VSMCs" = "#8f35e5",  
                   "HOX+ Progenitor" = "#ff9f2e",
                   "Biliary epithelial cells" = "#e53561", 
                   "Periportal bifurcations cells ?" = "forestgreen")


# Differentially expressed genes  ----------------------------------------------

Idents(sample.postqc) <- sample.postqc@meta.data$cell.type

sample.markers <- FindAllMarkers(object = sample.postqc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1, group.by="cell.type")
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

png(paste(plots.dir.postQC, sample_name, "_DoHeatmap_top10markers_postQC.png",sep=""),width=22,height=15,units="cm",res=1000)
print(DoHeatmap(sample.postqc, features = top10$gene, cells=colnames(sample.postqc), group.colors = cell.type.cols, group.by = "cell.type", size = 2) + 
        scale_fill_viridis(option="mako") +
        theme(text = element_text(size = 6)) +
        NoLegend())
dev.off()

write.table(sample.markers, quote = F,row.names = FALSE, paste(plots.dir.postQC, sample_name, "_allmarkers_postQC.txt",sep=""),  sep="\t")


DefaultAssay(sample.postqc) <- "SCT"

# visualize QC after filtering--------------------------------------------------

p<-VlnPlot(
  object = sample.postqc,
  features = c("nCount_RNA","nFeature_RNA", "percent.mt"),
  group.by = "cell.type",
  cols=cell.type.cols,
  ncol = 3,
  pt.size = 0
) & NoLegend()


p
ggsave(paste0(plots.dir.postQC, sample_name, "_postQC_plots.png"), plot = p, device = "png", width = 22, height = 14, units = "cm", dpi = 320)


# Cell cycle scoring------------------------------------------------------------

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

DefaultAssay(sample.postqc) <- "RNA"
sample.postqc <- NormalizeData(sample.postqc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

length(s.genes)
setdiff(s.genes, rownames(sample.postqc))
s.genes <- s.genes[which(s.genes %in% rownames(sample.postqc))]
length(s.genes)

length(g2m.genes)
setdiff(g2m.genes, rownames(sample.postqc))
g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(sample.postqc))]
length(g2m.genes)

sample.postqc <- CellCycleScoring(sample.postqc, s.features = s.genes, g2m.features = g2m.genes)

# view cell cycle scores and phase assignments
head(sample.postqc[[]])

p<-FeaturePlot(sample.postqc, reduction="umap", features="G2M.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_G2M.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="S.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_S.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-DimPlot(sample.postqc, group.by = "cell.type", reduction = "umap", cols = cell.type.cols, pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_cell_type.png"), plot = p, device = "png", width = 18, height = 10, units = "cm", dpi = 320,)


saveRDS(sample.postqc, file=paste0(output.directory, sample_name, "postQC.rds"))
saveRDS(sample, file=paste0(output.directory, sample_name, "preQC.rds"))


### 3. CHC4001N ----------------------------------------------------------------

sample_name <- "CHC4001N"
print(sample_name)

# Create useful directories paths-----------------------------------------------

data_dir <- paste0("//10.93.23.19/hepato partage/GENOMIC_DATA/scRNAseq/Project_HB_tumors/cellranger_outputs/PJ2302056_cellranger_v5/", sample_name, "/filtered_feature_bc_matrix/")
samp.directory <- paste0("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/", sample_name, "/")

output.directory <- paste0(samp.directory, "preprocess.QC/")
if(!dir.exists(output.directory)){dir.create(output.directory)}

plots.dir.preQC <- paste0(output.directory, "plots/preQC/")
if(!dir.exists(plots.dir.preQC)){dir.create(plots.dir.preQC)}

plots.dir.postQC <- paste0(output.directory, "plots/postQC/")
if(!dir.exists(plots.dir.postQC)){dir.create(plots.dir.postQC)}

# load the RNA data-------------------------------------------------------------
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz , atac_fragments.tsv.gz.tbi 
data <- Read10X(data.dir = data_dir) # Load expression data 

# get gene annotations for hg38-------------------------------------------------
hg38 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(hg38) <- "hg38"
hg38 <- renameSeqlevels(hg38, mapSeqlevels(seqlevels(hg38), "UCSC"))

# create a Seurat object containing the RNA data--------------------------------
sample <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = sample_name
)

# UMAP and clustering  ---------------------------------------------------------
# colors and parameters --------------------------------------------------------
my_scale_exp = scale_colour_gradientn(colors = c( "gray70", "#ffff00",  "red", "#cc0000"))

# Normalization with SCTransform -----------------------------------------------
DefaultAssay(sample) <- "RNA"
sample <- SCTransform(sample, verbose = TRUE, variable.features.n=3000) # sctransform = par d?faut 3000 nmostvar et pas 2000 comme lognormalize
sample <- ScaleData(sample, assay = "RNA")
DefaultAssay(sample) <- "SCT"

# PCA on gene expression--------------------------------------------------------
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=50) # 50 components

# Clustering --------------------------------------------------------------------
nbcomp=30
sample <- RunUMAP(sample, dims=1:nbcomp, verbose=F)
resolution_for_clusters <- 0.1  # change resolution to change nb of clusters determined

sample <- FindNeighbors(sample, dims = 1:30, k.param = 15)
sample <- FindClusters(sample, resolution = 0.1)
sample$clusters_0_1 <- sample$seurat_clusters

sample <- FindClusters(sample, resolution = 0.2)
sample$clusters_0_2 <- sample$seurat_clusters

sample <- FindClusters(sample, resolution = 0.3)
sample$clusters_0_3 <- sample$seurat_clusters

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "clusters_0_1", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_1_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample, group.by = "clusters_0_2", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_2_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample, group.by = "clusters_0_3", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_3_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

# add mitochondrial percentage -------------------------------------------------
DefaultAssay(sample) <- "RNA"
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern="^MT-") 


# visualize QC in UMAP ---------------------------------------------------------
sample$nCount_RNA_log <- log10(sample$nCount_RNA)
sample$percent.mt_log <- log2(sample$percent.mt)

p<-FeaturePlot(sample, reduction="umap", features="percent.mt_log",  pt.size = 0.1) + my_scale_exp & NoLegend()
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_percent_mt_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nCount_RNA_log",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nCount_RNA_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Show specific markers to annotate cell types ---------------------------------

DefaultAssay(sample) <- "SCT"

p<-FeaturePlot(sample, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CD163_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_FLT1_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT19_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT7_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CYP2E1_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HAL_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HOXB3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HOXB3_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_VIM_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

FeaturePlot(sample, reduction="umap", features="HGF",  pt.size = 0.1) + my_scale_exp

# Annotate cell types ----------------------------------------------------------



sample@meta.data <- sample@meta.data%>%
  mutate(cell.type.old = case_when(clusters_0_2 =="0" ~ "Hepatocytes",
                                   clusters_0_2 == "1" ~ "Periportal bifurcations cells ?",
                                   clusters_0_2 == "2" ~ "Hepatocytes",
                                   clusters_0_2 == "3" ~ "High ambient RNA and MT genes",
                                   clusters_0_2 == "4" ~ "Endothelial cells",
                                   clusters_0_2 == "5" ~ "HOX+ Progenitor",
                                   clusters_0_2 == "6" ~ "Macrophages - Kupffer cells",
                                   clusters_0_2 == "7" ~ "Fibroblasts - HSCs"))

# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "High ambient RNA and MT genes" = "#36e22e", 
                   "Periportal bifurcations cells ?" = "forestgreen",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Fibroblasts - HSCs" = "#8f35e5",
                   "Biliary epithelial cells" = "#e53561", 
                   "HOX+ Progenitor" = "#ff9f2e")

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "cell.type.old", reduction = "umap", cols = cell.type.cols) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_3_cell_type_old_preQC.png"), plot = p, device = "png", width = 17, height = 10, units = "cm", dpi = 320)


# Differentially expressed genes  ----------------------------------------------


sample.markers <- FindAllMarkers(object = sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

png(paste(plots.dir.preQC, sample_name, "_DoHeatmap_top10markers_preQC.png",sep=""),width=22,height=15,units="cm",res=1000)
print(DoHeatmap(sample, features = top10$gene,cells=colnames(sample), group.colors = cell.type.cols, group.by = "cell.type.old", size = 2) + 
        scale_fill_viridis(option="mako") +
        theme(text = element_text(size = 6)) +
        NoLegend())
dev.off()

write.table(sample.markers, quote = F,row.names = FALSE, paste(plots.dir.preQC, sample_name, "_allmarkers_preQC.txt",sep=""),  sep="\t")


# Perform quality controls  ---------------------------------------------------

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("percent.mt"),log = T,
  cols = cell.type.cols,
  pt.size = 0
) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_MT_percentage_log_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nCount_RNA_log"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_nCount_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nFeature_RNA"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0)  & NoLegend()

p

ggsave(paste0(plots.dir.preQC, sample_name, "_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


#  Filter out low quality hepatocytes ------------------------------------------

sample.postqc <- subset(sample, subset = cell.type.old =="High ambient RNA and MT genes", invert=TRUE)

sample.postqc <- subset(
  x = sample.postqc,
  subset=nFeature_RNA>1000 & percent.mt<5)

print(length(colnames(sample)))
print(length(colnames(sample.postqc)))

# Remove doublets predicted with scrublet---------------------------------------

tmp_sobj <- Seurat::CreateSeuratObject(counts=Seurat::GetAssayData(sample.postqc,assay='RNA',slot="counts"))
tmp_sobj@assays$RNA@data <- Seurat::GetAssayData(sample.postqc,assay='RNA',slot="data")
tmp_sobj$orig.ident <- sample.postqc$orig.ident
sce <- Seurat::as.SingleCellExperiment(tmp_sobj)

sample.postqc$scDblFinder.class <- Seurat::as.Seurat(scDblFinder::scDblFinder(sce,samples=tmp_sobj$orig.ident))$scDblFinder.class
sample.postqc$scDblFinder.class <- unname(sample.postqc$scDblFinder.class == "doublet") # turn into TRUE/FALSE

p<-DimPlot(sample.postqc, group.by = "scDblFinder.class", cols = c("gray80","red"), reduction = "umap", raster=F) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "scDblFinder.class_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

print(length(colnames(sample.postqc)))

sample.postqc <- subset(
  x = sample.postqc,
  subset = scDblFinder.class=="FALSE")

print(length(colnames(sample.postqc)))

# Clustering --------------------------------------------------------------------

ElbowPlot(sample.postqc, ndims=15)
sample.postqc <- RunPCA(sample.postqc, features = VariableFeatures(sample.postqc), npcs=20) 
sample.postqc <- RunUMAP(sample.postqc, dims=1:20, verbose=T, reduction = "pca")
sample.postqc <- FindNeighbors(sample.postqc, dims = 1:20, k.param = 15)

sample.postqc <- FindClusters(sample.postqc, resolution = 0.1)
sample.postqc$clusters_0_1_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.2)
sample.postqc$clusters_0_2_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.3)
sample.postqc$clusters_0_3_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.5)
sample.postqc$clusters_0_5_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.6)
sample.postqc$clusters_0_6_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.8)
sample.postqc$clusters_0_8_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.9)
sample.postqc$clusters_0_9_postQC <- sample.postqc$seurat_clusters


# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample.postqc, group.by = "clusters_0_1_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_1_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320,)

p<-DimPlot(sample.postqc, group.by = "clusters_0_2_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_2_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_3_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_3_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_5_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_5_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_6_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_6_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_8_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_8_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)


p<-DimPlot(sample.postqc, group.by = "clusters_0_9_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_9_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="percent.mt",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_percent_mt_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_nFeature_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# visualize specific gene expression -------------------------------------------
DefaultAssay(sample.postqc) <- "SCT"

# endothelial cells-------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD34",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD34_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VWF",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VWF_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANGPT2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANGPT2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="PTGDS",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PTGDS_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FLT1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGFBP3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGFBP3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CDH5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CDH5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC14A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC14A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="SLCO2A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SLCO2A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="RELN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_RELN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC4M",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC4M_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC1B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC1B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FCN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FCN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OIT3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OIT3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hepatic stellate cells -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COLEC11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COLEC11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LRAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LRAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NGFR",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NGFR_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# 11p-related genes ------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGF2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_H19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# progenitor markers -----------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GPC3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GPC3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="DLK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DLK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="EPCAM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_EPCAM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KIT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KIT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AFP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AFP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LIN28B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LIN28B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="THY1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_THY1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hematopoiesis ----------------------------------------------------------------
p<-FeaturePlot(sample.postqc, reduction="umap", features="HBB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HBM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Fibroblasts ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL6A3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL6A3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OGN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OGN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="POSTN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_POSTN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VIM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FN1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FN1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Macrophages-------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD163_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD68",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD68_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD80",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD80_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD74",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD74_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HLA-DRA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HLA-DRA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LYZ",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LYZ_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# B lymphocytes ----------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="JCHAIN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_JCHAIN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="IGKC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGKC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD79A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD79A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MZB1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MZB1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# T lymphocytes ----------------------------------------------------------------


p<-FeaturePlot(sample.postqc, reduction="umap", features="CCL5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CCL5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMK",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMK_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD3D",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD3D_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NKG7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NKG7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# NK cells ---------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GNLY",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GNLY_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD247",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD247_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD160",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD160_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NCR1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NCR1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# VSMC -------------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TAGLN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TAGLN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MYH11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MYH11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ACTA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ACTA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Lipogenesis ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOC1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOC1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ALB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ALB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="F2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_F2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Periportal markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="PCK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PCK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CPS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CPS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HAL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Perivenous markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP2E1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="GSTM3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GSTM3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NOTUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NOTUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FASN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FASN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="DPP4",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DPP4_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Wnt pathway activation -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TBX3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TBX3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AXIN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AXIN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LGR5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LGR5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GLUL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GLUL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Biliary epithelial cells------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SOX9",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SOX9_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SPP1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SPP1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------


p<-FeaturePlot(sample.postqc, reduction="umap", features="ITGB2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ITGB2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SELP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SELP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)



# Annotate cell types ----------------------------------------------------------

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample.postqc, group.by = "cell.type.old", reduction = "umap", cols = cell.type.cols) 
p

sample.postqc@meta.data <- sample.postqc@meta.data%>%
  mutate(cell.type = case_when(clusters_0_9_postQC =="0" ~ "Hepatocytes",
                               clusters_0_9_postQC == "1" ~ "Hepatocytes",
                               clusters_0_9_postQC == "2" ~ "Hepatocytes",
                               clusters_0_9_postQC == "3" ~ "Endothelial cells",
                               clusters_0_9_postQC == "4" ~ "Hepatocytes",
                               clusters_0_9_postQC == "5" ~ "Fibroblasts - HSCs",
                               clusters_0_9_postQC == "6" ~ "HOX+ Progenitor",
                               clusters_0_9_postQC == "7" ~ "Macrophages - Kupffer cells",))

# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Fibroblasts - HSCs" = "#8f35e5",
                   "Biliary epithelial cells" = "#e53561", 
                   "HOX+ Progenitor" = "#ff9f2e")


DefaultAssay(sample.postqc) <- "SCT"

# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "High ambient RNA and MT genes" = "#36e22e", 
                   "Periportal bifurcations cells ?" = "forestgreen",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Fibroblasts - HSCs" = "#8f35e5",
                   "Biliary epithelial cells" = "#e53561", 
                   "HOX+ Progenitor" = "#ff9f2e")

# visualize QC after filtering--------------------------------------------------

p<-VlnPlot(
  object = sample.postqc,
  features = c("nCount_RNA","nFeature_RNA", "percent.mt"),
  group.by = "cell.type",
  cols=cell.type.cols,
  ncol = 3,
  pt.size = 0
) & NoLegend()


p
ggsave(paste0(plots.dir.postQC, sample_name, "_postQC_plots_postQC.png"), plot = p, device = "png", width = 22, height = 14, units = "cm", dpi = 320)


# Cell cycle scoring------------------------------------------------------------

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

DefaultAssay(sample.postqc) <- "RNA"
sample.postqc <- NormalizeData(sample.postqc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

length(s.genes)
setdiff(s.genes, rownames(sample.postqc))
s.genes <- s.genes[which(s.genes %in% rownames(sample.postqc))]
length(s.genes)

length(g2m.genes)
setdiff(g2m.genes, rownames(sample.postqc))
g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(sample.postqc))]
length(g2m.genes)

sample.postqc <- CellCycleScoring(sample.postqc, s.features = s.genes, g2m.features = g2m.genes)

# view cell cycle scores and phase assignments
head(sample.postqc[[]])

p<-FeaturePlot(sample.postqc, reduction="umap", features="G2M.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_G2M.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="S.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_S.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-DimPlot(sample.postqc, group.by = "cell.type", reduction = "umap", cols = cell.type.cols, pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_cell_type.png"), plot = p, device = "png", width = 15, height = 10, units = "cm", dpi = 320,)

saveRDS(sample.postqc, file=paste0(output.directory, sample_name, "postQC.rds"))
saveRDS(sample, file=paste0(output.directory, sample_name, "preQC.rds"))

### 4. CHC3115N ----------------------------------------------------------------

sample_name <- "CHC3115N"
print(sample_name)

# Create useful directories paths-----------------------------------------------

data_dir <- paste0("//10.93.23.19/hepato partage/GENOMIC_DATA/scRNAseq/Project_HB_tumors/cellranger_outputs/PJ2302056_cellranger_v5/", sample_name, "/filtered_feature_bc_matrix/")

samp.directory <- paste0("D:/Dropbox/11p15.5 mosaicism/sn-RNAseq/samp_2996_3115_3559_4001/results/", sample_name, "/")

output.directory <- paste0(samp.directory, "preprocess.QC/")
if(!dir.exists(output.directory)){dir.create(output.directory)}

plots.dir.preQC <- paste0(output.directory, "plots/preQC/")
if(!dir.exists(plots.dir.preQC)){dir.create(plots.dir.preQC)}

plots.dir.postQC <- paste0(output.directory, "plots/postQC/")
if(!dir.exists(plots.dir.postQC)){dir.create(plots.dir.postQC)}

# load the RNA data-------------------------------------------------------------
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz , atac_fragments.tsv.gz.tbi 
data <- Read10X(data.dir = data_dir) # Load expression data 

# get gene annotations for hg38-------------------------------------------------
hg38 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(hg38) <- "hg38"
hg38 <- renameSeqlevels(hg38, mapSeqlevels(seqlevels(hg38), "UCSC"))

# create a Seurat object containing the RNA data--------------------------------
sample <- CreateSeuratObject(
  counts = data,
  assay = "RNA",
  project = sample_name
)


# UMAP and clustering  ---------------------------------------------------------
# colors and parameters --------------------------------------------------------
my_scale_exp = scale_colour_gradientn(colors = c( "gray70", "#ffff00",  "red", "#cc0000"))


# Normalization with SCTransform -----------------------------------------------
DefaultAssay(sample) <- "RNA"
sample <- SCTransform(sample, verbose = TRUE, variable.features.n=3000) # sctransform = par d?faut 3000 nmostvar et pas 2000 comme lognormalize
sample <- ScaleData(sample, assay = "RNA")
DefaultAssay(sample) <- "SCT"

# PCA on gene expression--------------------------------------------------------
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs=50) # 50 components

# Clustering --------------------------------------------------------------------
sample <- RunUMAP(sample, dims=1:30, verbose=F)
sample <- FindNeighbors(sample, dims = 1:30)
sample <- FindClusters(sample, resolution = 0.1)
sample$clusters_0_1 <- sample$seurat_clusters

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "clusters_0_1", reduction = "umap") 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_1_preQC.png"), plot = p, device = "png", width = 12, height = 10, units = "cm", dpi = 320)

# add mitochondrial percentage -------------------------------------------------
DefaultAssay(sample) <- "RNA"
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern="^MT-") 


# visualize QC in UMAP ---------------------------------------------------------
sample$nCount_RNA_log <- log10(sample$nCount_RNA)
sample$percent.mt_log <- log2(sample$percent.mt)

p<-FeaturePlot(sample, reduction="umap", features="percent.mt_log",  pt.size = 0.1) + my_scale_exp & NoLegend()
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_percent_mt_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="nCount_RNA_log",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_nCount_RNA_log_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Show specific markers to annotate cell types ---------------------------------

DefaultAssay(sample) <- "SCT"

p<-FeaturePlot(sample, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CD163_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CD163_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT19_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_KRT7_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_CYP2E1_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HAL_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="HOXB3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_HOXB3_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_VIM_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample, reduction="umap", features="COL6A3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.preQC, sample_name, "_UMAP_COL6A3_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Annotate cell types ----------------------------------------------------------

sample@meta.data <- sample@meta.data%>%
  mutate(cell.type.old = case_when(clusters_0_1 =="0" ~ "Hepatocytes",
                                   clusters_0_1 == "1" ~ "High ambient RNA and MT genes",
                                   clusters_0_1 == "2" ~ "Endothelial cells",
                                   clusters_0_1 == "3" ~ "Macrophages - Kupffer cells",
                                   clusters_0_1 == "4" ~ "Fibroblasts - HSCs - VSMCs",
                                   clusters_0_1 == "5" ~ "Biliary epithelial cells",
                                   clusters_0_1 == "6" ~ "HOX+ Progenitor"))

# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "High ambient RNA and MT genes" = "#36e22e", 
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Fibroblasts - HSCs - VSMCs" = "#8f35e5", 
                   "Biliary epithelial cells" = "#e53561", 
                   "HOX+ Progenitor" = "#ff9f2e")

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample, group.by = "cell.type.old", reduction = "umap", cols = cell.type.cols) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "_clusters_0_1_cell_type_old_preQC.png"), plot = p, device = "png", width = 17, height = 10, units = "cm", dpi = 320)


# Differentially expressed genes  ----------------------------------------------


sample.markers <- FindAllMarkers(object = sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

png(paste(plots.dir.preQC, sample_name, "_DoHeatmap_top10markers_preQC.png",sep=""),width=22,height=15,units="cm",res=1000)
print(DoHeatmap(sample, features = top10$gene,cells=colnames(sample), group.colors = cell.type.cols, group.by = "cell.type.old", size = 2) + 
        scale_fill_viridis(option="mako") +
        theme(text = element_text(size = 6)) +
        NoLegend())
dev.off()

write.table(sample.markers, quote = F,row.names = FALSE, paste(plots.dir.preQC, sample_name, "_allmarkers_preQC.txt",sep=""),  sep="\t")

# Perform quality controls  ---------------------------------------------------

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("percent.mt"),
  cols = cell.type.cols,
  pt.size = 0
) + geom_hline(yintercept=5) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_MT_percentage_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("percent.mt"),log = T,
  cols = cell.type.cols,
  pt.size = 0
) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_MT_percentage_log_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nCount_RNA_log"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0) & NoLegend()
p

ggsave(paste0(plots.dir.preQC, sample_name, "_nCount_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)

p<-VlnPlot(
  object = sample,
  group.by = "cell.type.old",
  features = c("nFeature_RNA"),
  cols = cell.type.cols,
  ncol = 1,
  pt.size = 0)  & NoLegend()

p

ggsave(paste0(plots.dir.preQC, sample_name, "_nFeature_RNA_preQC.png"), plot = p, device = "png", width = 10, height = 15, units = "cm", dpi = 320)


#  Filter out low quality hepatocytes ------------------------------------------

sample.postqc <- subset(sample, subset = cell.type.old =="High ambient RNA and MT genes", invert=TRUE)

sample.postqc <- subset(
  x = sample.postqc,
  subset=nFeature_RNA>1000 & percent.mt<5)

print(length(colnames(sample)))
print(length(colnames(sample.postqc)))

# Remove doublets predicted with scrublet---------------------------------------

tmp_sobj <- Seurat::CreateSeuratObject(counts=Seurat::GetAssayData(sample.postqc,assay='RNA',slot="counts"))
tmp_sobj@assays$RNA@data <- Seurat::GetAssayData(sample.postqc,assay='RNA',slot="data")
tmp_sobj$orig.ident <- sample.postqc$orig.ident
sce <- Seurat::as.SingleCellExperiment(tmp_sobj)

sample.postqc$scDblFinder.class <- Seurat::as.Seurat(scDblFinder::scDblFinder(sce,samples=tmp_sobj$orig.ident))$scDblFinder.class
sample.postqc$scDblFinder.class <- unname(sample.postqc$scDblFinder.class == "doublet") # turn into TRUE/FALSE

p<-DimPlot(sample.postqc, group.by = "scDblFinder.class", cols = c("gray80","red"), reduction = "umap", raster=F) 
p
ggsave(paste0(plots.dir.preQC, sample_name, "scDblFinder.class_preQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

print(length(colnames(sample.postqc)))

sample.postqc <- subset(
  x = sample.postqc,
  subset = scDblFinder.class=="FALSE")

print(length(colnames(sample.postqc)))

# Clustering --------------------------------------------------------------------

sample.postqc <- RunPCA(sample.postqc, features = VariableFeatures(sample.postqc), npcs=20) 
sample.postqc <- RunUMAP(sample.postqc, dims=1:20, verbose=T)
sample.postqc <- FindNeighbors(sample.postqc, dims = 1:20)

sample.postqc <- FindClusters(sample.postqc, resolution = 0.1)
sample.postqc$clusters_0_1_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.2)
sample.postqc$clusters_0_2_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.05)
sample.postqc$clusters_0_05_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.3)
sample.postqc$clusters_0_3_postQC <- sample.postqc$seurat_clusters

sample.postqc <- FindClusters(sample.postqc, resolution = 0.4)
sample.postqc$clusters_0_4_postQC <- sample.postqc$seurat_clusters

# Plot UMAP --------------------------------------------------------------------

p<-DimPlot(sample.postqc, group.by = "clusters_0_05_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_05_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320,)

p<-DimPlot(sample.postqc, group.by = "clusters_0_1_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_1_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320,)

p<-DimPlot(sample.postqc, group.by = "clusters_0_2_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_2_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)

p<-DimPlot(sample.postqc, group.by = "clusters_0_3_postQC", reduction = "umap",  pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_clusters_0_3_postQC.png"), plot = p, device = "png", width = 10, height = 10, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="percent.mt",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_percent_mt_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="nFeature_RNA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_nFeature_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# visualize specific gene expression -------------------------------------------
DefaultAssay(sample.postqc) <- "SCT"

# endothelial cells-------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD34",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD34_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VWF",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VWF_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANGPT2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANGPT2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="PTGDS",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PTGDS_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FLT1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FLT1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGFBP3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGFBP3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CDH5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CDH5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC14A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC14A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="SLCO2A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SLCO2A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="RELN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_RELN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC4M",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC4M_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CLEC1B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CLEC1B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FCN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FCN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OIT3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OIT3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hepatic stellate cells -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COLEC11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COLEC11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LRAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LRAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NGFR",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NGFR_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# 11p-related genes ------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="IGF2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGF2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="H19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_H19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# HOX+ Progenitor markers -----------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GPC3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GPC3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="DLK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DLK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="EPCAM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_EPCAM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KIT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KIT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AFP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AFP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LIN28B",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LIN28B_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="THY1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_THY1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Hematopoiesis ----------------------------------------------------------------
p<-FeaturePlot(sample.postqc, reduction="umap", features="HBB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ANK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ANK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HBM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HBM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Fibroblasts ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="COL6A3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_COL6A3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OGN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OGN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="POSTN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_POSTN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="VIM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_VIM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FN1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FN1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Macrophages-------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD163",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD163_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD68",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD68_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD80",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD80_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="C1QC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_C1QC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD74",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD74_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="HLA-DRA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HLA-DRA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LYZ",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LYZ_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# B lymphocytes ----------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="JCHAIN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_JCHAIN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="IGKC",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_IGKC_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD79A",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD79A_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MZB1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MZB1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# T lymphocytes ----------------------------------------------------------------


p<-FeaturePlot(sample.postqc, reduction="umap", features="CCL5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CCL5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMA",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMA_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GZMK",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GZMK_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD3D",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD3D_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NKG7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NKG7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# NK cells ---------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="GNLY",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GNLY_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD247",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD247_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CD160",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CD160_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NCR1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NCR1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# VSMC -------------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TAGLN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TAGLN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="MYH11",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_MYH11_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ACTA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ACTA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Lipogenesis ------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOC1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOC1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ALB",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ALB_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="APOA2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_APOA2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="F2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_F2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Periportal markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="PCK1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_PCK1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="CPS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CPS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASS1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASS1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="ASL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ASL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="HAL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_HAL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


# Perivenous markers------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP2E1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP2E1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="OAT",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_OAT_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="GSTM3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GSTM3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="NOTUM",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_NOTUM_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="CYP1A1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_CYP1A1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="FASN",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_FASN_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-FeaturePlot(sample.postqc, reduction="umap", features="DPP4",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_DPP4_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Wnt pathway activation -------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="TBX3",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_TBX3_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="AXIN2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_AXIN2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="LGR5",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_LGR5_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="GLUL",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_GLUL_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Biliary epithelial cells------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT7",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT7_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="KRT19",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_KRT19_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SOX9",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SOX9_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SPP1",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SPP1_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

# Platelets --------------------------------------------------------------------

p<-FeaturePlot(sample.postqc, reduction="umap", features="ITGB2",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_ITGB2_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="SELP",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_SELP_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)



# Annotate cell types ----------------------------------------------------------

sample.postqc@meta.data <- sample.postqc@meta.data %>%
  mutate(cell.type = case_when(clusters_0_05_postQC =="0" ~ "Hepatocytes",
                               clusters_0_05_postQC == "1" ~ "Endothelial cells",
                               clusters_0_05_postQC == "2" ~ "Biliary epithelial cells",
                               clusters_0_05_postQC == "3" ~ "Macrophages - Kupffer cells",
                               clusters_0_05_postQC == "4" ~ "Fibroblasts - HSCs - VSMCs",
                               clusters_0_05_postQC == "5" ~ "HOX+ Progenitor"))

# Cell type colors -------------------------------------------------------------

cell.type.cols = c("Hepatocytes" = "#07dcf3",
                   "Endothelial cells" = "#f3ef07", 
                   "Macrophages - Kupffer cells" = "#357ee5", 
                   "Fibroblasts - HSCs - VSMCs" = "#8f35e5", 
                   "Biliary epithelial cells" = "#e53561", 
                   "HOX+ Progenitor" = "#ff9f2e")


# Differentially expressed genes  ----------------------------------------------

Idents(sample.postqc) <- sample.postqc@meta.data$cell.type

sample.markers <- FindAllMarkers(object = sample.postqc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1, group.by="cell.type")
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

png(paste(plots.dir.postQC, sample_name, "_DoHeatmap_top10markers_postQC.png",sep=""),width=22,height=15,units="cm",res=1000)
print(DoHeatmap(sample.postqc, features = top10$gene, cells=colnames(sample.postqc), group.colors = cell.type.cols, group.by = "cell.type", size = 2) + 
        scale_fill_viridis(option="mako") +
        theme(text = element_text(size = 6)) +
        NoLegend())
dev.off()

write.table(sample.markers, quote = F,row.names = FALSE, paste(plots.dir.postQC, sample_name, "_allmarkers_postQC.txt",sep=""),  sep="\t")


DefaultAssay(sample.postqc) <- "SCT"

# Transcriptomic program cluster "Hematopoietic stem cells" --------------------

library(dittoSeq)
library(scuttle)

sample.markers <- FindMarkers(object = sample.postqc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1, group.by="cell.type", ident.1 = "HOX+ Progenitor")

sample.sce <- as.SingleCellExperiment(sample.postqc)

Genes_of_interest <- sample.markers[1:60,]
Genes_of_interest = rownames(Genes_of_interest)

celltype_mean <- aggregateAcrossCells(sample.sce,  
                                      ids = sample.sce$cell.type, 
                                      statistics = "mean",
                                      use.assay.type = "counts", 
                                      subset.row = Genes_of_interest)

dittoHeatmap(celltype_mean,
             assay = "counts", cluster_cols = TRUE, 
             scale = "row",  # "none" no scaling otherwise "row"
             heatmap.colors = viridis(100), 
             fontsize = 8,
             annot.by = c("cell.type", "ncells") ,
             annotation_colors = list(cell.type = cell.type.cols,  
                                      ncells = plasma(100))
)

`%nin%` = Negate(`%in%`)
Hox_genes = rownames(sample.sce)[grepl("HOX", rownames(sample.sce))]
Hox_genes = Hox_genes[Hox_genes %nin% c("SHOX", "RHOXF1-AS1", "RHOXF1")]
Hox_genes = Hox_genes[!grepl("AS", Hox_genes)]

celltype_mean <- aggregateAcrossCells(sample.sce,  
                                      ids = sample.sce$cell.type, 
                                      statistics = "mean",
                                      use.assay.type = "counts", 
                                      subset.row = Hox_genes)

dittoHeatmap(celltype_mean,
             assay = "counts", cluster_cols = TRUE, 
             scale = "row",  # "none" no scaling otherwise "row"
             heatmap.colors = viridis(100), 
             fontsize = 8,
             annot.by = c("cell.type", "ncells") ,
             annotation_colors = list(cell.type = cell.type.cols,  
                                      ncells = plasma(100))
)


HSC_markers = c("SPINK2", "CD34", "CD27", "KIT")

celltype_mean <- aggregateAcrossCells(sample.sce,  
                                      ids = sample.sce$cell.type, 
                                      statistics = "mean",
                                      use.assay.type = "counts", 
                                      subset.row = HSC_markers)

dittoHeatmap(celltype_mean,
             assay = "counts", cluster_cols = TRUE, 
             scale = "row",  # "none" no scaling otherwise "row"
             heatmap.colors = viridis(100), 
             fontsize = 8,
             annot.by = c("cell.type", "ncells") ,
             annotation_colors = list(cell.type = cell.type.cols,  
                                      ncells = plasma(100))
)


# visualize QC after filtering--------------------------------------------------

p<-VlnPlot(
  object = sample.postqc,
  features = c("nCount_RNA","nFeature_RNA", "percent.mt"),
  group.by = "cell.type",
  cols=cell.type.cols,
  ncol = 3,
  pt.size = 0
) & NoLegend()


p
ggsave(paste0(plots.dir.postQC, sample_name, "_postQC_plots_postQC.png"), plot = p, device = "png", width = 22, height = 14, units = "cm", dpi = 320)


# Cell cycle scoring------------------------------------------------------------

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

DefaultAssay(sample.postqc) <- "RNA"
sample.postqc <- NormalizeData(sample.postqc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

length(s.genes)
setdiff(s.genes, rownames(sample.postqc))
s.genes <- s.genes[which(s.genes %in% rownames(sample.postqc))]
length(s.genes)

length(g2m.genes)
setdiff(g2m.genes, rownames(sample.postqc))
g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(sample.postqc))]
length(g2m.genes)

sample.postqc <- CellCycleScoring(sample.postqc, s.features = s.genes, g2m.features = g2m.genes)

# view cell cycle scores and phase assignments
head(sample.postqc[[]])

p<-FeaturePlot(sample.postqc, reduction="umap", features="G2M.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_G2M.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)

p<-FeaturePlot(sample.postqc, reduction="umap", features="S.Score",  pt.size = 0.1) + my_scale_exp
p
ggsave(paste0(plots.dir.postQC, sample_name, "_UMAP_S.Score_postQC.png"), plot = p, device = "png", width = 10, height = 8, units = "cm", dpi = 320)


p<-DimPlot(sample.postqc, group.by = "cell.type", reduction = "umap", cols = cell.type.cols, pt.size = 0.1) 
p
ggsave(paste0(plots.dir.postQC, sample_name, "UMAP_cell_type.png"), plot = p, device = "png", width = 15, height = 10, units = "cm", dpi = 320,)


saveRDS(sample.postqc, file=paste0(output.directory, sample_name, "postQC.rds"))
saveRDS(sample, file=paste0(output.directory, sample_name, "preQC.rds"))

