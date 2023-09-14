library(ggplot2)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
library(devtools)
library(pheatmap)
library(viridis)
library(ggpubr)
library(clustree)
library(harmony)

### 1. SEURAT 4001 NT ----------------------------------------------------------
  # 1.1 Create image and Seurat object ----------------------------------------------------
image.dir = "//10.93.23.19/HEPATO PARTAGE/GENOMIC_DATA/Visium/Visum_results/Visium_1_PJ2108330/4001NT/for_R/"
data.dir = "//10.93.23.19/HEPATO PARTAGE/GENOMIC_DATA/Visium/Visum_results/Visium_1_PJ2108330/4001NT/for_R/"

obj_4001NT = Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Samp_4001NT",
  filter.matrix = TRUE,
  to.upper = FALSE #,
#  image = Img_4001NT
)

  # 1.2 QC -------------------------------------------------------------------------------
VlnPlot(obj_4001NT, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(obj_4001NT, features = "nCount_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")

VlnPlot(obj_4001NT, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(obj_4001NT, features = "nFeature_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")

FeatureScatter(obj_4001NT, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

SpatialFeaturePlot(obj_4001NT, features = c("nCount_Spatial", "nFeature_Spatial"))

  # 1.3 Normalization --------------------------------------------------------------------

obj_4001NT = SCTransform(obj_4001NT, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE)

  # 1.4 Find variable features------------------------------------------------------------
obj_4001NT_var <- FindVariableFeatures(obj_4001NT, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj_4001NT_var), 10)
top20 <- head(VariableFeatures(obj_4001NT_var), 20)
top10

  # 1.5 Clustering -------------------------------------------------------------
    # 1.5.1 PCA analysis------------------------------------------------------------------ 
obj_4001NT  <- RunPCA(obj_4001NT , assay = "SCT", verbose = FALSE, features = VariableFeatures(object = obj_4001NT_var))
DimPlot(obj_4001NT, reduction = "pca")
DimHeatmap(obj_4001NT, dims = 1:8, balanced = TRUE, cells =500)
ElbowPlot(obj_4001NT) 

obj_4001NT  <- FindNeighbors(obj_4001NT , reduction = "pca", dims = 1:20)
obj_4001NT  <- FindClusters(obj_4001NT , verbose = FALSE, resolution = 0.25) # Choose resolution to change cluster number
head(Idents(obj_4001NT), 5)

Spatial_cols = c("#FFD54F","indianred", "#1A237E", "#880E4F", "#F57F17", "#64B5F6", "red")

DimPlot(obj_4001NT, reduction = "pca", label = FALSE, 
        cols = Spatial_cols, 
        pt.size=1)

SpatialDimPlot(obj_4001NT, label = FALSE, pt.size.factor = 1.6, image.alpha = 0.5) + scale_fill_manual(values = Spatial_cols)

    # 1.5.2 UMAP --------------------------------------------------------------------------

obj_4001NT <- RunUMAP(obj_4001NT , reduction = "pca", dims = 1:20)

DimPlot(obj_4001NT, reduction = "umap", label = FALSE, cols = Spatial_cols, pt.size=1)
SpatialDimPlot(obj_4001NT, label = FALSE, label.color = "black", label.size = 3) + scale_fill_manual(values = Spatial_cols)

new.cluster.ids <- c("Non_mosaic_hepatocytes", "Mosaic_hepatocytes", "Portal_tracts", "Tumor_capsule", "Tumor", "Tumor2",
                     "cluster_6", "B_Lymphocytes")
names(new.cluster.ids) <- levels(obj_4001NT)
obj_4001NT <- RenameIdents(obj_4001NT, new.cluster.ids)

#Highlight clusters
SpatialDimPlot(obj_4001NT, cells.highlight = CellsByIdentities(object = obj_4001NT, 
                                                               idents = c("Non_mosaic_hepatocytes", "Mosaic_hepatocytes", "Portal_tracts", "Tumor_capsule", "Tumor", "Tumor2",
                                                                          "cluster_6", "B_Lymphocytes")),
               facet.highlight = TRUE, ncol = 4, pt.size.factor = 1.6)

SpatialDimPlot(obj_4001NT, label = FALSE, label.size = 3, pt.size.factor = 1.6,label.color = "black", image.alpha = 0.5) + scale_fill_manual(values = Spatial_cols)


  # 1.6 Find differentially expressed genes ----------------------------------------------
    # 1.6.1 Find DE markers --------------------------------------------------------------
de_markers <- FindMarkers(obj_4001NT, ident.1 = "Mosaic_hepatocytes", ident.2 = "Non_mosaic_hepatocytes") 

#resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/Visium/Results/");if(!file.exists(resdir))	dir.create(resdir) 
#write.table(as.data.frame(de_markers), file= file.path(resdir,paste0("Mosaic_vs_non_mosaic_hepatocytes", ".txt")), sep="\t")
de.4001.pos = de_markers %>%
  filter(p_val_adj<=0.01,
         avg_log2FC>0)

de.4001.neg = de_markers %>%
  filter(p_val_adj<=0.01,
         avg_log2FC<0)


SpatialFeaturePlot(object = obj_4001NT, features = rownames(de_markers)[1:5], 
                   alpha = c(0.3, 1), ncol = 5, pt.size.factor = 1.6)


    # 1.6.2 Represent top deregulated genes ------------------------------------------------

VlnPlot(obj_4001NT, features = "THY1", cols = Spatial_cols, pt.size=0.1, log = FALSE)
RidgePlot(obj_4001NT, features = "THY1", cols = Spatial_cols, log = FALSE)
DotPlot(obj_4001NT, features = c("IGF2"), cols = c("white", "darkblue"))
FeaturePlot(obj_4001NT, features = c("SDS", "GLUL"),
            #cols = rev(brewer.pal(11,"Spectral")),
            cols = c("lightgrey", "#E8EAF6", "#9FA8DA", "#3F51B5",  "darkblue"),
            pt.size = 0.8, ncol = 2)
SpatialFeaturePlot(obj_4001NT, features = "P2RY1",
                   #alpha = c(0.1,0.6), 
                   interactive = FALSE, image.alpha = 0.6, ncol =1)

    # 1.6.3 Heatmap with top differentially expressed genes-------------------------------

samp_4001_markers <- FindAllMarkers(obj_4001NT, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
samp_4001_markers = samp_4001_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = p_val) %>%
  top_n(n = 20, wt = p_val) -> top10

DoHeatmap(obj_4001NT, features = top10$gene, group.colors = Spatial_cols) + NoLegend()
DoHeatmap(obj_4001NT, features = as.character(c("CPS1","PCK1", "CYP2E1", "ALB", "IGF2", "ADH1A", "GPC3", "VIM")),
          group.colors = Spatial_cols, label = FALSE) + NoLegend()

    # 1.6.4 Apply different genes signatures -------------------------------------
# load tables
library(readxl)
Genes_list_dir <- "D:/Dropbox/11p15.5 mosaicism/Visium/Genes_list_ref"
genes_for_MCP <- read_xlsx(file.path(Genes_list_dir, "20170721_genesMCPcounter_orientationFonctionelle.xlsx"),
                           sheet = "new.list") %>% as.data.frame()
genes_Danaher <- read_xlsx(file.path(Genes_list_dir, "Danaher_2017.xlsx"),
                           sheet = "S4. Selected markers") %>% 
  rename(HUGO.symbols = Gene, Cell.population = `Cell Type`)

genes_Massalha <- read_xlsx(file.path(Genes_list_dir, "Massalha_supp2.xlsx") 
                    #choose all genes or genes represented in fig1
                    #,sheet = "Genes_represented_fig1"
                    ) %>% as.data.frame() %>% 
  rename(HUGO.symbols = markers, Cell.population = `cell_type`)

genes_dobie <- read_xlsx(file.path(Genes_list_dir, "dobie_supp.xlsx"), sheet = "markers_derived_from_fig1") %>% as.data.frame() %>% 
  rename(HUGO.symbols = gene, Cell.population = `cluster`)
# pas d'intersect avec FB et HSC et les rownames de 4001_NT


genes_payen <- read_xlsx(file.path(Genes_list_dir, "payen_Supp.xlsx")) %>% as.data.frame() 

# Create_signatures -------------

# Mosaic signature from 4001 DE genes positive enrichment -----------
de.4001.pos = rownames(de.4001.pos)

obj_4001NT$de.4001.pos <- apply(GetAssayData(object = obj_4001NT, slot = "data")[de.4001.pos, ], 2, mean) 

SpatialFeaturePlot(obj_4001NT, features = "de.4001.pos",
                   interactive = FALSE, pt.size.factor = 1.6, ncol =1, image.alpha = 0.6,alpha = c(0.1, 1))

# Mosaic signature from 4001 DE genes negative enrichment -----------
de.4001.neg = rownames(de.4001.neg)

obj_4001NT$de.4001.neg <- apply(GetAssayData(object = obj_4001NT, slot = "data")[de.4001.neg, ], 2, mean) 

SpatialFeaturePlot(obj_4001NT, features = "de.4001.neg",
                   interactive = FALSE, pt.size.factor = 1.6, ncol =1, image.alpha = 0.6,alpha = c(0.1, 1))

# MCP counter.-------------------
T_cells <- genes_for_MCP %>% filter(Cell.population == "T cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
B_cells <- genes_for_MCP %>% filter(Cell.population == "B cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
endo <- genes_for_MCP %>% filter(Cell.population == "Endothelial cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
NK <- genes_for_MCP %>% filter(Cell.population == "NK cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Cytotoxic_cells <- genes_for_MCP %>% filter(Cell.population == "Cytotoxic") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
fibroblasts <- genes_for_MCP %>% filter(Cell.population == "Fibroblasts") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))


obj_4001NT$Cytotoxic_cells <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Cytotoxic_cells, ], 2, mean) 
obj_4001NT$T_cell_score <- apply(GetAssayData(object = obj_4001NT, slot = "data")[T_cells, ], 2, mean) 
obj_4001NT$B_cell_score <- apply(GetAssayData(object = obj_4001NT, slot = "data")[B_cells, ], 2, mean) 
obj_4001NT$endo_score <- apply(GetAssayData(object = obj_4001NT, slot = "data")[endo, ], 2, mean) 
obj_4001NT$fibroblasts <- apply(GetAssayData(object = obj_4001NT, slot = "data")[fibroblasts, ], 2, mean) 
obj_4001NT$NK <- apply(GetAssayData(object = obj_4001NT, slot = "data")[NK, ], 2, mean) 


# Massalha -------------------
B_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "B cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
CAFs_Massalha <- genes_Massalha %>% filter(Cell.population == "CAFs_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
cDC1_Massalha <- genes_Massalha %>% filter(Cell.population == "cDC1_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
cDC2_Massalha <- genes_Massalha %>% filter(Cell.population == "cDC2_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Kupffer_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "Kupffer cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
LSEC_Massalha <- genes_Massalha %>% filter(Cell.population == "LSEC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
LVEC_Massalha <- genes_Massalha %>% filter(Cell.population == "LVEC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
LVECt_Massalha <- genes_Massalha %>% filter(Cell.population == "LVECt_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Pericytes_Massalha<- genes_Massalha %>% filter(Cell.population == "Pericytes_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
SAMs_Massalha <- genes_Massalha %>% filter(Cell.population == "SAMs_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Stellate_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "Stellate cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
T_cells_Massalha<- genes_Massalha %>% filter(Cell.population == "T cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
TM1_Massalha <- genes_Massalha %>% filter(Cell.population == "TM1_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
vSMC_Massalha <- genes_Massalha %>% filter(Cell.population == "vSMC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))

obj_4001NT$B_cells_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[B_cells_Massalha, ], 2, mean) 
obj_4001NT$CAFs_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[CAFs_Massalha, ], 2, mean) 
obj_4001NT$cDC1_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[cDC1_Massalha, ], 2, mean) 
obj_4001NT$cDC2_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[cDC2_Massalha, ], 2, mean) 
obj_4001NT$Kupffer_cells_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Kupffer_cells_Massalha, ], 2, mean) 
obj_4001NT$LSEC_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LSEC_Massalha, ], 2, mean) 
obj_4001NT$LVEC_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LVEC_Massalha, ], 2, mean) 
obj_4001NT$LVECt_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LVECt_Massalha, ], 2, mean) 
obj_4001NT$Pericytes_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Pericytes_Massalha, ], 2, mean) 
obj_4001NT$SAMs_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[SAMs_Massalha, ], 2, mean) 
obj_4001NT$Stellate_cells_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Stellate_cells_Massalha, ], 2, mean) 
obj_4001NT$T_cells_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[T_cells_Massalha, ], 2, mean) 
obj_4001NT$TM1_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[TM1_Massalha, ], 2, mean) 
obj_4001NT$vSMC_Massalha <- apply(GetAssayData(object = obj_4001NT, slot = "data")[vSMC_Massalha, ], 2, mean) 

# Dobie  -------------------
FB_dobie<- genes_dobie %>% filter(Cell.population == "FB") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
VSMC_dobie<- genes_dobie %>% filter(Cell.population == "VSMC") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
HSC_dobie<- genes_dobie %>% filter(Cell.population == "HSC") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))

obj_4001NT$VSMC_dobie <- apply(GetAssayData(object = obj_4001NT, slot = "data")[VSMC_dobie, ], 2, mean)
obj_4001NT$FB_dobie <- apply(GetAssayData(object = obj_4001NT, slot = "data")[FB_dobie, ], 2, mean)
obj_4001NT$HSC_dobie <- apply(GetAssayData(object = obj_4001NT, slot = "data")[HSC_dobie, ], 2, mean)

# Payen  -------------------

VSMC_payen <- genes_payen %>% filter(Cell.population == "VSMC score") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
HSC_payen <- genes_payen %>% filter(Cell.population == "HSC score") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
HSC1_payen <- genes_payen %>% filter(Cell.population == "HSC1 signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
HSC2_payen <- genes_payen %>% filter(Cell.population == "HSC2 signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Xenobiotic_payen <- genes_payen %>% filter(Cell.population == "Xenobiotic metabolism") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Fatty_acid_payen <- genes_payen %>% filter(Cell.population == "Fatty acid biosynthesis") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Retinoic_acid_payen <- genes_payen %>% filter(Cell.population == "Retinoid metabolic process") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
GAG_metabolism_payen <- genes_payen %>% filter(Cell.population == "GAG metabolism") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Elastic_fibers_payen <- genes_payen %>% filter(Cell.population == "Elastic fiber") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Secretion_payen <- genes_payen %>% filter(Cell.population == "Secretion") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Kupffer_payen <- genes_payen %>% filter(Cell.population == "Kupffer signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))


obj_4001NT$VSMC_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[VSMC_payen, ], 2, mean) 
obj_4001NT$HSC_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[HSC_payen, ], 2, mean) 
obj_4001NT$HSC1_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[HSC1_payen, ], 2, mean) 
obj_4001NT$HSC2_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[HSC2_payen, ], 2, mean) 
obj_4001NT$Xenobiotic_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Xenobiotic_payen, ], 2, mean) 
obj_4001NT$Fatty_acid_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Fatty_acid_payen, ], 2, mean) 
obj_4001NT$Retinoic_acid_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Retinoic_acid_payen, ], 2, mean) 
obj_4001NT$GAG_metabolism_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[GAG_metabolism_payen, ], 2, mean) 
obj_4001NT$Elastic_fibers_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Elastic_fibers_payen, ], 2, mean) 
obj_4001NT$Secretion_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Secretion_payen, ], 2, mean) 
obj_4001NT$Kupffer_payen <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Kupffer_payen, ], 2, mean) 

# My own signatures---------
my.sig = read_xlsx(file.path(Genes_list_dir, "my_gene_sets.xlsx"),
          sheet = "GeneSet") %>% as.data.frame()

HSC <- my.sig %>% filter(GeneSet == "HSC") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Proliferation <- my.sig %>% filter(GeneSet == "Proliferation") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
VSMC <- my.sig %>% filter(GeneSet == "VSMC") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Iron <- my.sig %>% filter(GeneSet == "Iron ion homeostasis") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
LVEC <- my.sig %>% filter(GeneSet == "LVEC") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
LSEC <- my.sig %>% filter(GeneSet == "LSEC") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
LT <- my.sig %>% filter(GeneSet == "LT") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
NK <- my.sig %>% filter(GeneSet == "NK") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Angiogenesis <- my.sig %>% filter(GeneSet == "Angiogenesis") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
DC <- my.sig %>% filter(GeneSet == "DC") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
CAF <- my.sig %>% filter(GeneSet == "CAF") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Perivenous <- my.sig %>% filter(GeneSet == "Perivenous") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Periportal <- my.sig %>% filter(GeneSet == "Periportal") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Hematopoiesis <- my.sig %>% filter(GeneSet == "Hematopoiesis") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Progenitor <- my.sig %>% filter(GeneSet == "Progenitor markers") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Stem <- my.sig %>% filter(GeneSet == "Hepatic stem cell") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
Cholangiocytes <- my.sig %>% filter(GeneSet == "Cholangiocytes") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
LB <- my.sig %>% filter(GeneSet == "LB") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))
platelets <- my.sig %>% filter(GeneSet == "Activated_platelets") %>% pull(Gene) %>% intersect(rownames(obj_4001NT))


obj_4001NT$Proliferation <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Proliferation, ], 2, mean) 
obj_4001NT$HSC <- apply(GetAssayData(object = obj_4001NT, slot = "data")[HSC, ], 2, mean) 
obj_4001NT$VSMC <- apply(GetAssayData(object = obj_4001NT, slot = "data")[VSMC, ], 2, mean) 
obj_4001NT$Iron <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Iron, ], 2, mean) 
obj_4001NT$LVEC <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LVEC, ], 2, mean) 
obj_4001NT$LSEC <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LSEC, ], 2, mean) 
obj_4001NT$LT <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LT, ], 2, mean) 
obj_4001NT$NK <- apply(GetAssayData(object = obj_4001NT, slot = "data")[NK, ], 2, mean) 
obj_4001NT$Angiogenesis <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Angiogenesis, ], 2, mean) 
obj_4001NT$DC <- apply(GetAssayData(object = obj_4001NT, slot = "data")[DC, ], 2, mean) 
obj_4001NT$CAF <- apply(GetAssayData(object = obj_4001NT, slot = "data")[CAF, ], 2, mean) 
obj_4001NT$Perivenous <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Perivenous, ], 2, mean) 
obj_4001NT$Periportal <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Periportal, ], 2, mean) 
obj_4001NT$Hematopoiesis <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Hematopoiesis, ], 2, mean) 
obj_4001NT$Progenitor <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Progenitor, ], 2, mean) 
obj_4001NT$Stem <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Stem, ], 2, mean) 
obj_4001NT$Cholangiocytes <- apply(GetAssayData(object = obj_4001NT, slot = "data")[Cholangiocytes, ], 2, mean) 
obj_4001NT$LB <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LB, ], 2, mean) 
obj_4001NT$platelets <- apply(GetAssayData(object = obj_4001NT, slot = "data")[LB, ], 2, mean) 

plot = SpatialFeaturePlot(obj_4001NT, alpha = 1, features =  "Perivenous", pt.size.factor = 1.6, image.alpha = 0) +  
  scale_fill_viridis(option="plasma", guide = guide_colourbar(title.theme = element_blank(),
                                                              label.theme = element_text(size=20),
                                                              barwidth = 10,
                                                              barheight = 2)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) & NoGrid() & NoAxes() 

#ggsave("/Users/jill/Desktop/plot_test.png", plot = plot, device = "png", width = 15, height = 15, units = "cm", dpi = 320)


### 2. SEURAT 3559 NT --------------------------------------------------------
  # 2.1 Create image and Seurat object ------------------------------------------
image.dir = "//10.93.23.19/HEPATO PARTAGE/GENOMIC_DATA/Visium/Visum_results/Visium_1_PJ2108330/3559NT/For_R/"
data.dir = "//10.93.23.19/HEPATO PARTAGE/GENOMIC_DATA/Visium/Visum_results/Visium_1_PJ2108330/3559NT/For_R/"


obj_3559NT = Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Samp_3559NT",
  filter.matrix = TRUE,
  to.upper = FALSE
)

  # 2.2 QC ---------------------------------------------------------------------
VlnPlot(obj_3559NT, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(obj_3559NT, features = "nCount_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")

VlnPlot(obj_3559NT, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(obj_3559NT, features = "nFeature_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")

FeatureScatter(obj_3559NT, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
SpatialFeaturePlot(obj_3559NT, features = c("nCount_Spatial", "nFeature_Spatial"))


  # 2.3 Normalization ----------------------------------------------------------

obj_3559NT = SCTransform(obj_3559NT, assay = "Spatial", verbose = FALSE)


  # 2.4 apply gene signatures --------------------------------------------------
#Mosaic signature from 4001 DE genes positive enrichment -----------

obj_3559NT$de.4001.pos <- apply(GetAssayData(object = obj_3559NT, slot = "data")[de.4001.pos, ], 2, mean) 

SpatialFeaturePlot(obj_3559NT, features = "de.4001.pos",
                   interactive = FALSE, pt.size.factor = 1.6, ncol =1, image.alpha = 0.6,alpha = c(0.1, 1))

#Mosaic signature from 4001 DE genes negative enrichment -----------
obj_3559NT$de.4001.neg <- apply(GetAssayData(object = obj_3559NT, slot = "data")[de.4001.neg, ], 2, mean) 

SpatialFeaturePlot(obj_3559NT, features = "de.4001.neg",
                   interactive = FALSE, pt.size.factor = 1.6, ncol =1, image.alpha = 0.6,alpha = c(0.1, 1))

# MCP counter.-------------------
T_cells <- genes_for_MCP %>% filter(Cell.population == "T cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
B_cells <- genes_for_MCP %>% filter(Cell.population == "B cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
endo <- genes_for_MCP %>% filter(Cell.population == "Endothelial cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
NK <- genes_for_MCP %>% filter(Cell.population == "NK cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
Cytotoxic_cells <- genes_for_MCP %>% filter(Cell.population == "Cytotoxic") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))
fibroblasts <- genes_for_MCP %>% filter(Cell.population == "Fibroblasts") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_4001NT))


obj_3559NT$Cytotoxic_cells <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Cytotoxic_cells, ], 2, mean) 
obj_3559NT$T_cell_score <- apply(GetAssayData(object = obj_3559NT, slot = "data")[T_cells, ], 2, mean) 
obj_3559NT$B_cell_score <- apply(GetAssayData(object = obj_3559NT, slot = "data")[B_cells, ], 2, mean) 
obj_3559NT$endo_score <- apply(GetAssayData(object = obj_3559NT, slot = "data")[endo, ], 2, mean) 
obj_3559NT$fibroblasts <- apply(GetAssayData(object = obj_3559NT, slot = "data")[fibroblasts, ], 2, mean) 
obj_3559NT$NK <- apply(GetAssayData(object = obj_3559NT, slot = "data")[NK, ], 2, mean) 


# Massalha -------------------
B_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "B cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
CAFs_Massalha <- genes_Massalha %>% filter(Cell.population == "CAFs_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
cDC1_Massalha <- genes_Massalha %>% filter(Cell.population == "cDC1_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
cDC2_Massalha <- genes_Massalha %>% filter(Cell.population == "cDC2_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Kupffer_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "Kupffer cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
LSEC_Massalha <- genes_Massalha %>% filter(Cell.population == "LSEC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
LVEC_Massalha <- genes_Massalha %>% filter(Cell.population == "LVEC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
LVECt_Massalha <- genes_Massalha %>% filter(Cell.population == "LVECt_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Pericytes_Massalha<- genes_Massalha %>% filter(Cell.population == "Pericytes_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
SAMs_Massalha <- genes_Massalha %>% filter(Cell.population == "SAMs_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Stellate_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "Stellate cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
T_cells_Massalha<- genes_Massalha %>% filter(Cell.population == "T cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
TM1_Massalha <- genes_Massalha %>% filter(Cell.population == "TM1_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
vSMC_Massalha <- genes_Massalha %>% filter(Cell.population == "vSMC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))

obj_3559NT$B_cells_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[B_cells_Massalha, ], 2, mean) 
obj_3559NT$CAFs_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[CAFs_Massalha, ], 2, mean) 
obj_3559NT$cDC1_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[cDC1_Massalha, ], 2, mean) 
obj_3559NT$cDC2_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[cDC2_Massalha, ], 2, mean) 
obj_3559NT$Kupffer_cells_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Kupffer_cells_Massalha, ], 2, mean) 
obj_3559NT$LSEC_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LSEC_Massalha, ], 2, mean) 
obj_3559NT$LVEC_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LVEC_Massalha, ], 2, mean) 
obj_3559NT$LVECt_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LVECt_Massalha, ], 2, mean) 
obj_3559NT$Pericytes_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Pericytes_Massalha, ], 2, mean) 
obj_3559NT$SAMs_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[SAMs_Massalha, ], 2, mean) 
obj_3559NT$Stellate_cells_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Stellate_cells_Massalha, ], 2, mean) 
obj_3559NT$T_cells_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[T_cells_Massalha, ], 2, mean) 
obj_3559NT$TM1_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[TM1_Massalha, ], 2, mean) 
obj_3559NT$vSMC_Massalha <- apply(GetAssayData(object = obj_3559NT, slot = "data")[vSMC_Massalha, ], 2, mean) 

# Dobie  -------------------
FB_dobie<- genes_dobie %>% filter(Cell.population == "FB") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
VSMC_dobie<- genes_dobie %>% filter(Cell.population == "VSMC") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
HSC_dobie<- genes_dobie %>% filter(Cell.population == "HSC") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))

obj_3559NT$VSMC_dobie <- apply(GetAssayData(object = obj_3559NT, slot = "data")[VSMC_dobie, ], 2, mean)
obj_3559NT$FB_dobie <- apply(GetAssayData(object = obj_3559NT, slot = "data")[FB_dobie, ], 2, mean)
obj_3559NT$HSC_dobie <- apply(GetAssayData(object = obj_3559NT, slot = "data")[HSC_dobie, ], 2, mean)


# Payen  -------------------

VSMC_payen <- genes_payen %>% filter(Cell.population == "VSMC score") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
HSC_payen <- genes_payen %>% filter(Cell.population == "HSC score") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
HSC1_payen <- genes_payen %>% filter(Cell.population == "HSC1 signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
HSC2_payen <- genes_payen %>% filter(Cell.population == "HSC2 signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Xenobiotic_payen <- genes_payen %>% filter(Cell.population == "Xenobiotic metabolism") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Fatty_acid_payen <- genes_payen %>% filter(Cell.population == "Fatty acid biosynthesis") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Retinoic_acid_payen <- genes_payen %>% filter(Cell.population == "Retinoid metabolic process") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
GAG_metabolism_payen <- genes_payen %>% filter(Cell.population == "GAG metabolism") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Elastic_fibers_payen <- genes_payen %>% filter(Cell.population == "Elastic fiber") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Secretion_payen <- genes_payen %>% filter(Cell.population == "Secretion") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))
Kupffer_payen <- genes_payen %>% filter(Cell.population == "Kupffer signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3559NT))


obj_3559NT$VSMC_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[VSMC_payen, ], 2, mean) 
obj_3559NT$HSC_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[HSC_payen, ], 2, mean) 
obj_3559NT$HSC1_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[HSC1_payen, ], 2, mean) 
obj_3559NT$HSC2_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[HSC2_payen, ], 2, mean) 
obj_3559NT$Xenobiotic_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Xenobiotic_payen, ], 2, mean) 
obj_3559NT$Fatty_acid_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Fatty_acid_payen, ], 2, mean) 
obj_3559NT$Retinoic_acid_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Retinoic_acid_payen, ], 2, mean) 
obj_3559NT$GAG_metabolism_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[GAG_metabolism_payen, ], 2, mean) 
obj_3559NT$Elastic_fibers_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Elastic_fibers_payen, ], 2, mean) 
obj_3559NT$Secretion_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Secretion_payen, ], 2, mean) 
obj_3559NT$Kupffer_payen <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Kupffer_payen, ], 2, mean) 


# My own signatures---------
HSC <- my.sig %>% filter(GeneSet == "HSC") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Proliferation <- my.sig %>% filter(GeneSet == "Proliferation") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
VSMC <- my.sig %>% filter(GeneSet == "VSMC") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Iron <- my.sig %>% filter(GeneSet == "Iron ion homeostasis") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
LVEC <- my.sig %>% filter(GeneSet == "LVEC") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
LSEC <- my.sig %>% filter(GeneSet == "LSEC") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
LT <- my.sig %>% filter(GeneSet == "LT") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
NK <- my.sig %>% filter(GeneSet == "NK") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Angiogenesis <- my.sig %>% filter(GeneSet == "Angiogenesis") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
DC <- my.sig %>% filter(GeneSet == "DC") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
CAF <- my.sig %>% filter(GeneSet == "CAF") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Perivenous <- my.sig %>% filter(GeneSet == "Perivenous") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Periportal <- my.sig %>% filter(GeneSet == "Periportal") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Hematopoiesis <- my.sig %>% filter(GeneSet == "Hematopoiesis") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Progenitor <- my.sig %>% filter(GeneSet == "Progenitor markers") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Stem <- my.sig %>% filter(GeneSet == "Hepatic stem cell") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
Cholangiocytes <- my.sig %>% filter(GeneSet == "Cholangiocytes") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
LB <- my.sig %>% filter(GeneSet == "LB") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))
platelets <- my.sig %>% filter(GeneSet == "Activated_platelets") %>% pull(Gene) %>% intersect(rownames(obj_3559NT))


obj_3559NT$Proliferation <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Proliferation, ], 2, mean) 
obj_3559NT$HSC <- apply(GetAssayData(object = obj_3559NT, slot = "data")[HSC, ], 2, mean) 
obj_3559NT$VSMC <- apply(GetAssayData(object = obj_3559NT, slot = "data")[VSMC, ], 2, mean) 
obj_3559NT$Iron <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Iron, ], 2, mean) 
obj_3559NT$LVEC <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LVEC, ], 2, mean) 
obj_3559NT$LSEC <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LSEC, ], 2, mean) 
obj_3559NT$LT <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LT, ], 2, mean) 
obj_3559NT$NK <- apply(GetAssayData(object = obj_3559NT, slot = "data")[NK, ], 2, mean) 
obj_3559NT$Angiogenesis <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Angiogenesis, ], 2, mean) 
obj_3559NT$DC <- apply(GetAssayData(object = obj_3559NT, slot = "data")[DC, ], 2, mean) 
obj_3559NT$CAF <- apply(GetAssayData(object = obj_3559NT, slot = "data")[CAF, ], 2, mean) 
obj_3559NT$Perivenous <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Perivenous, ], 2, mean) 
obj_3559NT$Periportal <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Periportal, ], 2, mean) 
obj_3559NT$Hematopoiesis <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Hematopoiesis, ], 2, mean) 
obj_3559NT$Progenitor <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Progenitor, ], 2, mean) 
obj_3559NT$Stem <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Stem, ], 2, mean) 
obj_3559NT$Cholangiocytes <- apply(GetAssayData(object = obj_3559NT, slot = "data")[Cholangiocytes, ], 2, mean) 
obj_3559NT$LB <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LB, ], 2, mean) 
obj_3559NT$platelets <- apply(GetAssayData(object = obj_3559NT, slot = "data")[LB, ], 2, mean) 


plot = SpatialFeaturePlot(obj_3559NT, alpha = 1, features =  "Perivenous", pt.size.factor = 1.6, image.alpha = 0) +  
  scale_fill_viridis(option="plasma", guide = guide_colourbar(title.theme = element_blank(),
                                                              label.theme = element_text(size=20),
                                                              barwidth = 10,
                                                              barheight = 2)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) & NoGrid() & NoAxes() 

#ggsave("/Users/jill/Desktop/plot_test.png", plot = plot, device = "png", width = 15, height = 15, units = "cm", dpi = 320)

  # 2.5 Clustering -------------------------------------------------------------
    # 2.5.1 PCA analysis------------------------------------------------------------------ 
obj_3559NT_var <- FindVariableFeatures(obj_3559NT, selection.method = "vst", nfeatures = 5000)

obj_3559NT  <- RunPCA(obj_3559NT , assay = "SCT", verbose = FALSE, features = VariableFeatures(object = obj_3559NT_var))
DimPlot(obj_3559NT, reduction = "pca")
DimHeatmap(obj_3559NT, dims = 1:8, balanced = TRUE, cells =500)
ElbowPlot(obj_3559NT) 

obj_3559NT  <- FindNeighbors(obj_3559NT , reduction = "pca", dims = 1:20)
obj_3559NT  <- FindClusters(obj_3559NT , verbose = FALSE, resolution = 0.3) # Choose resolution to change cluster number

obj_3559NT <- RunUMAP(obj_3559NT , reduction = "pca", dims = 1:15)

Spatial_cols = viridis(6, direction = -1)
DimPlot(obj_3559NT, reduction = "pca", label = FALSE, pt.size=1) + scale_fill_manual(values = Spatial_cols)
SpatialDimPlot(obj_3559NT, label = FALSE, label.size = 3, pt.size.factor = 1.6, image.alpha = 0.5) + scale_fill_manual(values = Spatial_cols)


### 3. SEURAT 3115 NT ----------------------------------------------------------
# 3.1 Create image and Seurat object ----------------------------------------------------
image.dir = "//10.93.23.19/HEPATO PARTAGE/GENOMIC_DATA/Visium/Visum_results/Visium_3_PJ2203114/3115-N/for_R/"
data.dir = "//10.93.23.19/HEPATO PARTAGE/GENOMIC_DATA/Visium/Visum_results/Visium_3_PJ2203114/3115-N/for_R/"

obj_3115NT = Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Samp_3115NT",
  filter.matrix = TRUE,
  to.upper = FALSE)

# 3.2 QC -------------------------------------------------------------------------------
VlnPlot(obj_3115NT, features = "nCount_Spatial", pt.size = 0.1) + NoLegend() #+ scale_y_continuous(breaks = c(1000, 2000, 3000, 4000, 5000))
SpatialFeaturePlot(obj_3115NT, features = "nCount_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")

VlnPlot(obj_3115NT, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(obj_3115NT, features = "nFeature_Spatial", pt.size.factor = 1.6) + theme(legend.position = "right")

FeatureScatter(obj_3115NT, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

SpatialFeaturePlot(obj_3115NT, features = c("nCount_Spatial", "nFeature_Spatial"))

obj_3115NT = obj_3115NT[, obj_3115NT$nCount_Spatial > 1500]

# 3.3 Normalization --------------------------------------------------------------------

obj_3115NT = SCTransform(obj_3115NT, assay = "Spatial", verbose = FALSE, return.only.var.genes = FALSE)

# 3.4 Find variable features------------------------------------------------------------
obj_3115NT_var <- FindVariableFeatures(obj_3115NT, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj_3115NT_var), 10)
top20 <- head(VariableFeatures(obj_3115NT_var), 20)

top10

# 3.5 Clustering -------------------------------------------------------------
# 3.5.1 PCA analysis------------------------------------------------------------------ 
obj_3115NT  <- RunPCA(obj_3115NT , assay = "SCT", verbose = FALSE, features = VariableFeatures(object = obj_3115NT_var))
DimPlot(obj_3115NT, reduction = "pca")
DimHeatmap(obj_3115NT, dims = 1:8, balanced = TRUE, cells =500)
ElbowPlot(obj_3115NT) 

obj_3115NT  <- FindNeighbors(obj_3115NT , reduction = "pca", dims = 1:20)
obj_3115NT  <- FindClusters(obj_3115NT , verbose = FALSE, resolution = 0.2) # Choose resolution to change cluster number
head(Idents(obj_3115NT), 5)

Spatial_cols = c("indianred","#FFD54F", "blue", "red")
DimPlot(obj_3115NT, reduction = "pca", label = FALSE, cols = Spatial_cols, pt.size=1)
SpatialDimPlot(obj_3115NT, label = FALSE, label.size = 3,  pt.size.factor = 1.6, image.alpha = 0.5) + scale_fill_manual(values=Spatial_cols)

# 3.5.2 UMAP --------------------------------------------------------------------------

obj_3115NT <- RunUMAP(obj_3115NT , reduction = "pca", dims = 1:20)
obj_3115NT
DimPlot(obj_3115NT, reduction = "umap", label = FALSE, cols = Spatial_cols, pt.size=1)
SpatialDimPlot(obj_3115NT, label = FALSE, label.color = "black", label.size = 3) + scale_fill_manual(values=Spatial_cols)


new.cluster.ids <- c("Mosaic_hepatocytes", "Non_mosaic_hepatocytes","Portal_tracts", "effet.de.bord")
names(new.cluster.ids) <- levels(obj_3115NT)
obj_3115NT <- RenameIdents(obj_3115NT, new.cluster.ids)

#Highlight clusters

SpatialDimPlot(obj_3115NT, label = FALSE, label.size = 3,  pt.size.factor = 1.6,label.color = "black", image.alpha = 0.5) + scale_fill_manual(values=Spatial_cols)


# 3.6 Find differentially expressed genes ----------------------------------------------
# 3.6.1 Find DE markers --------------------------------------------------------------
de_markers <- FindMarkers(obj_3115NT, ident.1 = "Mosaic_hepatocytes", ident.2 = "Non_mosaic_hepatocytes", slot = "counts") 

#resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/Visium/Results/");if(!file.exists(resdir))	dir.create(resdir) 
#write.table(as.data.frame(de_markers), file= file.path(resdir,paste0("Mosaic_vs_non_mosaic_hepatocytes", ".txt")), sep="\t")
de.3115.pos = de_markers %>%
  filter(p_val_adj<=0.01,
         avg_log2FC>0)

de.3115.neg = de_markers %>%
  filter(p_val_adj<=0.01,
         avg_log2FC<0)


SpatialFeaturePlot(obj_3115NT, alpha = 1, features =  "GPC3", pt.size.factor = 1.6, image.alpha = 0) +  scale_fill_viridis(option="plasma")  & NoGrid() & NoAxes() 

# 3.6.2 Represent top deregulated genes ------------------------------------------------

VlnPlot(obj_3115NT, features = "THY1", cols = Spatial_cols, pt.size=0.1, log = FALSE)
RidgePlot(obj_3115NT, features = "THY1", cols = Spatial_cols, log = FALSE)
DotPlot(obj_3115NT, features = c("IGF2"), cols = c("white", "darkblue"))
FeaturePlot(obj_3115NT, features = c("SDS", "GLUL"),
            #cols = rev(brewer.pal(11,"Spectral")),
            cols = c("lightgrey", "#E8EAF6", "#9FA8DA", "#3F51B5",  "darkblue"),
            pt.size = 0.8, ncol = 2)
SpatialFeaturePlot(obj_3115NT, features = c("SDS","GLUL", "SLC1A2"),
                   #alpha = c(0.1,0.6), 
                   interactive = FALSE, image.alpha = 0.6, ncol =3)

# 3.6.3 Heatmap with top differentially expressed genes-------------------------------

samp_3115_markers <- FindAllMarkers(obj_3115NT, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
samp_3115_markers = samp_3115_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = p_val) %>%
  top_n(n = 20, wt = p_val) -> top10

DoHeatmap(obj_3115NT, features = top10$gene, group.colors = Spatial_cols) + NoLegend()
DoHeatmap(obj_3115NT, features = as.character(c("CPS1","PCK1", "CYP2E1", "ALB", "IGF2", "ADH1A", "GPC3", "VIM")),
          group.colors = Spatial_cols, label = FALSE) + NoLegend()

# 3.6.4 Apply different genes signatures -------------------------------------
# load tables
library(readxl)
Genes_list_dir <- "D:/Dropbox/11p15.5 mosaicism/Visium/Genes_list_ref"
genes_for_MCP <- read_xlsx(file.path(Genes_list_dir, "20170721_genesMCPcounter_orientationFonctionelle.xlsx"),
                           sheet = "new.list") %>% as.data.frame()
genes_Danaher <- read_xlsx(file.path(Genes_list_dir, "Danaher_2017.xlsx"),
                           sheet = "S4. Selected markers") %>% 
  rename(HUGO.symbols = Gene, Cell.population = `Cell Type`)

genes_Massalha <- read_xlsx(file.path(Genes_list_dir, "Massalha_supp2.xlsx") 
                            #choose all genes or genes represented in fig1
                            #,sheet = "Genes_represented_fig1"
) %>% as.data.frame() %>% 
  rename(HUGO.symbols = markers, Cell.population = `cell_type`)

genes_dobie <- read_xlsx(file.path(Genes_list_dir, "dobie_supp.xlsx"), sheet = "markers_derived_from_fig1") %>% as.data.frame() %>% 
  rename(HUGO.symbols = gene, Cell.population = `cluster`)
# pas d'intersect avec FB et HSC et les rownames de 3115_NT


genes_payen <- read_xlsx(file.path(Genes_list_dir, "payen_Supp.xlsx")) %>% as.data.frame() 

# Create_signatures -------------

# Mosaic signature from 3115 DE genes positive enrichment -----------
de.3115.pos = rownames(de.3115.pos)

obj_3115NT$de.3115.pos <- apply(GetAssayData(object = obj_3115NT, slot = "data")[de.3115.pos, ], 2, mean) 

SpatialFeaturePlot(obj_3115NT, features = "de.3115.pos",
                   interactive = FALSE, pt.size.factor = 1.6, ncol =1, image.alpha = 0.6,alpha = c(0.1, 1))

# Mosaic signature from 3115 DE genes negative enrichment -----------
de.3115.neg = rownames(de.3115.neg)

obj_3115NT$de.3115.neg <- apply(GetAssayData(object = obj_3115NT, slot = "data")[de.3115.neg, ], 2, mean) 

SpatialFeaturePlot(obj_3115NT, features = "de.3115.neg",
                   interactive = FALSE, pt.size.factor = 1.6, ncol =1, image.alpha = 0.6,alpha = c(0.1, 1))


# MCP counter.-------------------
T_cells <- genes_for_MCP %>% filter(Cell.population == "T cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
B_cells <- genes_for_MCP %>% filter(Cell.population == "B cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
endo <- genes_for_MCP %>% filter(Cell.population == "Endothelial cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
NK <- genes_for_MCP %>% filter(Cell.population == "NK cells") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Cytotoxic_cells <- genes_for_MCP %>% filter(Cell.population == "Cytotoxic") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
fibroblasts <- genes_for_MCP %>% filter(Cell.population == "Fibroblasts") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))


obj_3115NT$Cytotoxic_cells <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Cytotoxic_cells, ], 2, mean) 
obj_3115NT$T_cell_score <- apply(GetAssayData(object = obj_3115NT, slot = "data")[T_cells, ], 2, mean) 
obj_3115NT$B_cell_score <- apply(GetAssayData(object = obj_3115NT, slot = "data")[B_cells, ], 2, mean) 
obj_3115NT$endo_score <- apply(GetAssayData(object = obj_3115NT, slot = "data")[endo, ], 2, mean) 
obj_3115NT$fibroblasts <- apply(GetAssayData(object = obj_3115NT, slot = "data")[fibroblasts, ], 2, mean) 
obj_3115NT$NK <- apply(GetAssayData(object = obj_3115NT, slot = "data")[NK, ], 2, mean) 


# Massalha -------------------
B_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "B cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
CAFs_Massalha <- genes_Massalha %>% filter(Cell.population == "CAFs_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
cDC1_Massalha <- genes_Massalha %>% filter(Cell.population == "cDC1_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
cDC2_Massalha <- genes_Massalha %>% filter(Cell.population == "cDC2_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Kupffer_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "Kupffer cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
LSEC_Massalha <- genes_Massalha %>% filter(Cell.population == "LSEC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
LVEC_Massalha <- genes_Massalha %>% filter(Cell.population == "LVEC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
LVECt_Massalha <- genes_Massalha %>% filter(Cell.population == "LVECt_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Pericytes_Massalha<- genes_Massalha %>% filter(Cell.population == "Pericytes_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
SAMs_Massalha <- genes_Massalha %>% filter(Cell.population == "SAMs_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Stellate_cells_Massalha <- genes_Massalha %>% filter(Cell.population == "Stellate cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
T_cells_Massalha<- genes_Massalha %>% filter(Cell.population == "T cells_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
TM1_Massalha <- genes_Massalha %>% filter(Cell.population == "TM1_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
vSMC_Massalha <- genes_Massalha %>% filter(Cell.population == "vSMC_Massalha") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))

obj_3115NT$B_cells_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[B_cells_Massalha, ], 2, mean) 
obj_3115NT$CAFs_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[CAFs_Massalha, ], 2, mean) 
obj_3115NT$cDC1_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[cDC1_Massalha, ], 2, mean) 
obj_3115NT$cDC2_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[cDC2_Massalha, ], 2, mean) 
obj_3115NT$Kupffer_cells_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Kupffer_cells_Massalha, ], 2, mean) 
obj_3115NT$LSEC_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LSEC_Massalha, ], 2, mean) 
obj_3115NT$LVEC_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LVEC_Massalha, ], 2, mean) 
obj_3115NT$LVECt_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LVECt_Massalha, ], 2, mean) 
obj_3115NT$Pericytes_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Pericytes_Massalha, ], 2, mean) 
obj_3115NT$SAMs_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[SAMs_Massalha, ], 2, mean) 
obj_3115NT$Stellate_cells_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Stellate_cells_Massalha, ], 2, mean) 
obj_3115NT$T_cells_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[T_cells_Massalha, ], 2, mean) 
obj_3115NT$TM1_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[TM1_Massalha, ], 2, mean) 
obj_3115NT$vSMC_Massalha <- apply(GetAssayData(object = obj_3115NT, slot = "data")[vSMC_Massalha, ], 2, mean) 

# Dobie  -------------------
FB_dobie<- genes_dobie %>% filter(Cell.population == "FB") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
VSMC_dobie<- genes_dobie %>% filter(Cell.population == "VSMC") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
HSC_dobie<- genes_dobie %>% filter(Cell.population == "HSC") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))

obj_3115NT$VSMC_dobie <- apply(GetAssayData(object = obj_3115NT, slot = "data")[VSMC_dobie, ], 2, mean)
obj_3115NT$FB_dobie <- apply(GetAssayData(object = obj_3115NT, slot = "data")[FB_dobie, ], 2, mean)
obj_3115NT$HSC_dobie <- apply(GetAssayData(object = obj_3115NT, slot = "data")[HSC_dobie, ], 2, mean)

# Payen  -------------------

VSMC_payen <- genes_payen %>% filter(Cell.population == "VSMC score") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
HSC_payen <- genes_payen %>% filter(Cell.population == "HSC score") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
HSC1_payen <- genes_payen %>% filter(Cell.population == "HSC1 signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
HSC2_payen <- genes_payen %>% filter(Cell.population == "HSC2 signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Xenobiotic_payen <- genes_payen %>% filter(Cell.population == "Xenobiotic metabolism") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Fatty_acid_payen <- genes_payen %>% filter(Cell.population == "Fatty acid biosynthesis") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Retinoic_acid_payen <- genes_payen %>% filter(Cell.population == "Retinoid metabolic process") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
GAG_metabolism_payen <- genes_payen %>% filter(Cell.population == "GAG metabolism") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Elastic_fibers_payen <- genes_payen %>% filter(Cell.population == "Elastic fiber") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Secretion_payen <- genes_payen %>% filter(Cell.population == "Secretion") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))
Kupffer_payen <- genes_payen %>% filter(Cell.population == "Kupffer signature") %>% pull(HUGO.symbols) %>% intersect(rownames(obj_3115NT))

obj_3115NT$VSMC_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[VSMC_payen, ], 2, mean) 
obj_3115NT$HSC_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[HSC_payen, ], 2, mean) 
obj_3115NT$HSC1_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[HSC1_payen, ], 2, mean) 
obj_3115NT$HSC2_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[HSC2_payen, ], 2, mean) 
obj_3115NT$Xenobiotic_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Xenobiotic_payen, ], 2, mean) 
obj_3115NT$Fatty_acid_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Fatty_acid_payen, ], 2, mean) 
obj_3115NT$Retinoic_acid_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Retinoic_acid_payen, ], 2, mean) 
obj_3115NT$GAG_metabolism_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[GAG_metabolism_payen, ], 2, mean) 
obj_3115NT$Elastic_fibers_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Elastic_fibers_payen, ], 2, mean) 
obj_3115NT$Secretion_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Secretion_payen, ], 2, mean) 
obj_3115NT$Kupffer_payen <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Kupffer_payen, ], 2, mean) 


# My own signatures---------
my.sig = read_xlsx(file.path(Genes_list_dir, "my_gene_sets.xlsx"),
                   sheet = "GeneSet") %>% as.data.frame()

HSC <- my.sig %>% filter(GeneSet == "HSC") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Proliferation <- my.sig %>% filter(GeneSet == "Proliferation") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
VSMC <- my.sig %>% filter(GeneSet == "VSMC") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Iron <- my.sig %>% filter(GeneSet == "Iron ion homeostasis") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
LVEC <- my.sig %>% filter(GeneSet == "LVEC") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
LSEC <- my.sig %>% filter(GeneSet == "LSEC") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
LT <- my.sig %>% filter(GeneSet == "LT") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
NK <- my.sig %>% filter(GeneSet == "NK") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Angiogenesis <- my.sig %>% filter(GeneSet == "Angiogenesis") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
DC <- my.sig %>% filter(GeneSet == "DC") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
CAF <- my.sig %>% filter(GeneSet == "CAF") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Perivenous <- my.sig %>% filter(GeneSet == "Perivenous") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Periportal <- my.sig %>% filter(GeneSet == "Periportal") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Hematopoiesis <- my.sig %>% filter(GeneSet == "Hematopoiesis") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Progenitor <- my.sig %>% filter(GeneSet == "Progenitor markers") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Stem <- my.sig %>% filter(GeneSet == "Hepatic stem cell") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
Cholangiocytes <- my.sig %>% filter(GeneSet == "Cholangiocytes") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
LB <- my.sig %>% filter(GeneSet == "LB") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))
platelets <- my.sig %>% filter(GeneSet == "Activated_platelets") %>% pull(Gene) %>% intersect(rownames(obj_3115NT))



obj_3115NT$Proliferation <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Proliferation, ], 2, mean) 
obj_3115NT$HSC <- apply(GetAssayData(object = obj_3115NT, slot = "data")[HSC, ], 2, mean) 
obj_3115NT$VSMC <- apply(GetAssayData(object = obj_3115NT, slot = "data")[VSMC, ], 2, mean) 
obj_3115NT$Iron <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Iron, ], 2, mean) 
obj_3115NT$LVEC <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LVEC, ], 2, mean) 
obj_3115NT$LSEC <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LSEC, ], 2, mean) 
obj_3115NT$LT <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LT, ], 2, mean) 
obj_3115NT$NK <- apply(GetAssayData(object = obj_3115NT, slot = "data")[NK, ], 2, mean) 
obj_3115NT$Angiogenesis <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Angiogenesis, ], 2, mean) 
obj_3115NT$DC <- apply(GetAssayData(object = obj_3115NT, slot = "data")[DC, ], 2, mean) 
obj_3115NT$CAF <- apply(GetAssayData(object = obj_3115NT, slot = "data")[CAF, ], 2, mean) 
obj_3115NT$Perivenous <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Perivenous, ], 2, mean) 
obj_3115NT$Periportal <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Periportal, ], 2, mean) 
obj_3115NT$Hematopoiesis <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Hematopoiesis, ], 2, mean) 
obj_3115NT$Progenitor <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Progenitor, ], 2, mean) 
obj_3115NT$Stem <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Stem, ], 2, mean) 
obj_3115NT$Cholangiocytes <- apply(GetAssayData(object = obj_3115NT, slot = "data")[Cholangiocytes, ], 2, mean) 
obj_3115NT$LB <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LB, ], 2, mean) 
obj_3115NT$platelets <- apply(GetAssayData(object = obj_3115NT, slot = "data")[LB, ], 2, mean) 

plot = SpatialFeaturePlot(obj_3115NT, alpha = 1, features =  "Hematopoiesis", pt.size.factor = 2, image.alpha = 0) +  
  scale_fill_viridis(option="plasma", guide = guide_colourbar(title.theme = element_blank(),
                                                              label.theme = element_text(size=20) ,
                                                              barwidth = 10,
                                                              barheight = 2)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) & NoGrid() & NoAxes() 

#ggsave("/Users/jill/Desktop/plot_test.png", plot = plot, device = "png", width = 15, height = 15, units = "cm", dpi = 320)

###. 4. Merged analysis 4001NT, 3559NT and 3115NT-------------------------------
  # 4.1 Create merge data set --------------------------------------------------

obj_3559NT$orig.ident = "3559"
obj_4001NT$orig.ident = "4001"
obj_3115NT$orig.ident = "3115"
mos.merge <- merge(obj_4001NT, c(obj_3559NT, obj_3115NT))

mos.merge <- readRDS("D:/Dropbox/11p15.5 mosaicism/Visium/spatial_transcriptomics_merged_seurat.rds")


  # 4.2 Clustering -------------------------------------------------------------
DefaultAssay(mos.merge) <- "SCT"
VariableFeatures(mos.merge)<-c(VariableFeatures(obj_4001NT), VariableFeatures(obj_3559NT), VariableFeatures(obj_3115NT))
mos.merge <- RunPCA(mos.merge, verbose = FALSE, assay = "SCT", features = VariableFeatures(object = mos.merge))
mos.merge <- FindNeighbors(mos.merge, dims = 1:40)
mos.merge <- FindClusters(mos.merge, verbose = TRUE, resolution = 0.4)
mos.merge <- RunUMAP(mos.merge, dims = 1:40)


DimPlot(mos.merge, reduction = "umap",
        cols =c("0" = "#99FFFF", "1" = "orange", "2" ="#3399FF" , 
                "3" = "#990000","4" = "#FF3300",  "5" = "#1f1aed",
                "6" = "yellow","7" = "#8d00dc", "8" = "#66FF66"), label.size = 10)

SpatialDimPlot(mos.merge, cols = c("0" = "#99FFFF", "1" = "orange", "2" ="#3399FF" , 
                                   "3" = "#990000","4" = "#FF3300",  "5" = "#1f1aed",
                                   "6" = "yellow","7" = "#8d00dc", "8" = "#66FF66"), 
               alpha = 1, pt.size.factor = 1.6) 

plot1 <- SpatialDimPlot(mos.merge, cols = c("0" = "#99FFFF", "1" = "orange", "2" ="#3399FF" , 
                                           "3" = "#990000","4" = "#FF3300",  "5" = "#1f1aed",
                                           "6" = "yellow","7" = "#8d00dc", "8" = "#66FF66"), 
                       alpha = 1, pt.size.factor = 1.6, images="Samp_4001NT") +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot2 <- SpatialDimPlot(mos.merge, cols = c("0" = "#99FFFF", "1" = "orange", "2" ="#3399FF" , 
                                            "3" = "#990000","4" = "#FF3300",  "5" = "#1f1aed",
                                            "6" = "yellow","7" = "#8d00dc", "8" = "#66FF66"), 
                        alpha = 1, pt.size.factor = 1.6, images="Samp_3559NT")  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot3 <- SpatialDimPlot(mos.merge, cols = c("0" = "#99FFFF", "1" = "orange", "2" ="#3399FF" , 
                                            "3" = "#990000","4" = "#FF3300",  "5" = "#1f1aed",
                                            "6" = "yellow","7" = "#8d00dc", "8" = "#66FF66"), 
                        alpha = 1, pt.size.factor = 1.6, images="Samp_3115NT")  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot <- ggarrange(plot1, plot2, plot3, common.legend = FALSE, legend = "bottom", ncol=3) 
  
ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/visum_clust_res0.4.png",
       width = 28, height = 15, units = "cm")

DimPlot(mos.merge, reduction = "umap", group.by = "orig.ident",
        cols = c("4001" = "darkblue", "3559" = "#fbddf5", "3115" = "#039a87")) & NoGrid() & NoAxes() #& NoLegend()

DimPlot(mos.merge, reduction = "umap", group.by = "SCT_snn_res.0.4",
        cols = c("0" = "#99FFFF", "1" = "orange", "2" ="#3399FF" , 
                 "3" = "#990000","4" = "#FF3300",  "5" = "#1f1aed",
                 "6" = "yellow","7" = "#8d00dc", "8" = "#66FF66")) & NoGrid() & NoAxes() #& NoLegend()

# Resolution 0.3

mos.merge <- FindClusters(mos.merge, verbose = TRUE, resolution = 0.3)
mos.merge <- RunUMAP(mos.merge, dims = 1:40)

viridis.col <-viridis(direction = -1, n=9, option = "turbo")
viridis.col <- c("0" = "#7A0403FF",
                 "1" = "#D23105FF",
                 "2" = "#FB8022FF",
                 "3" = "#EDD03AFF",
                 "4" = "#A2FC3CFF",
                 "5" = "#31F299FF",
                 "6" = "#28BBECFF",
                 "7" = "#466BE3FF",
                 "8" = "#30123BFF")

plot1 <- SpatialDimPlot(mos.merge,
                        alpha = 1,  pt.size.factor = 1.6, images="Samp_4001NT", cols = viridis.col) +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot2 <- SpatialDimPlot(mos.merge,
                        alpha = 1,pt.size.factor = 1.6, images="Samp_3559NT", cols = viridis.col)  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot3 <- SpatialDimPlot(mos.merge,
                        alpha = 1, pt.size.factor = 1.6, images="Samp_3115NT", cols = viridis.col)  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot <- ggarrange(plot1, plot2, plot3, common.legend = FALSE, legend = "bottom", ncol=3) 

ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/visum_clust_res0.3.png",
       width = 25, height = 15, units = "cm")


DimPlot(mos.merge, reduction = "umap", group.by = "orig.ident",
        cols = c("4001" = "darkblue", "3559" = "#fbddf5", "3115" = "#039a87")) & NoGrid() & NoAxes() #& NoLegend()

DimPlot(mos.merge, reduction = "umap", group.by = "SCT_snn_res.0.3", cols = viridis.col) & NoGrid() & NoAxes() #& NoLegend()

# Resolution 0.5

mos.merge <- FindClusters(mos.merge, verbose = TRUE, resolution = 0.5)
mos.merge <- RunUMAP(mos.merge, dims = 1:40)

viridis.col <-viridis(direction = -1, n=12, option = "turbo")
viridis.col <- c("0" = "#7A0403FF",
                 "1" = "#BE2102FF",
                 "2" = "#EA4F0DFF",
                 "3" = "#FE922AFF",
                 "4" = "#F1CA3AFF",
                 "5" = "#C1F334FF",
                 "6" = "#7DFF56FF",
                 "7" = "#29EFA2FF",
                 "8" = "#1FC8DEFF",
                 "9" = "#4490FEFF",
                 "10" = "#4454C4FF",
                 "11" = "#30123BFF")

plot1 <- SpatialDimPlot(mos.merge,
                        alpha = 1, pt.size.factor = 1.6, images="Samp_4001NT", cols = viridis.col) +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot2 <- SpatialDimPlot(mos.merge,
                        alpha = 1, pt.size.factor = 1.6, images="Samp_3559NT", cols = viridis.col)  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot3 <- SpatialDimPlot(mos.merge,
                        alpha = 1, pt.size.factor = 1.6, images="Samp_3115NT", cols = viridis.col)  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot <- ggarrange(plot1, plot2, plot3, common.legend = FALSE, legend = "bottom", ncol=3) 

ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/visum_clust_res0.5.png",
       width = 25, height = 15, units = "cm")


DimPlot(mos.merge, reduction = "umap", group.by = "orig.ident",
        cols = c("4001" = "darkblue", "3559" = "#fbddf5", "3115" = "#039a87")) & NoGrid() & NoAxes() #& NoLegend()

DimPlot(mos.merge, reduction = "umap", group.by = "SCT_snn_res.0.5", cols = viridis.col) & NoGrid() & NoAxes() #& NoLegend()


# Build clustering tree 

res <- seq(0.1, 1, 0.1)
for (i in res) { 
  mos.merge <- FindClusters(mos.merge, resolution = i)
}

mos.merge@meta.data$SCT_snn_res.0.25 <- NULL

cols =c("0" = "#99FFFF", "1" = "orange", "2" ="#3399FF" , 
        "3" = "#990000","4" = "#FF3300",  "5" = "#1f1aed",
        "6" = "yellow","7" = "#8d00dc", "8" = "#66FF66")

plot <- clustree(mos.merge, assay ="SCT", node_alpha = 1) + 
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10))
ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/clustree_colored_by_resolution.png", 
       width = 15, height = 17, units = "cm")

plot <- clustree(mos.merge, assay ="SCT", node_colour = "cluster", node_alpha =1) + scale_color_manual(values = cols) + 
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10))

plot
ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/clustree_colored_by_cluster.png", 
       width = 15, height = 15, units = "cm")

# Run Harmony

mos.merge.harmony <- RunHarmony(mos.merge, reduction = "pca", assay.use = "SCT",
                                group.by.vars = "orig.ident")

mos.merge.harmony <- FindNeighbors(mos.merge.harmony, dims = 1:40, reduction = "harmony")
mos.merge.harmony <- FindClusters(mos.merge.harmony, verbose = TRUE, resolution = 0.3)
mos.merge.harmony <- RunUMAP(mos.merge.harmony, dims = 1:40, reduction = "harmony")

cols =c("0" = "#990000", "1" = "#99FFFF", "2" ="#3399FF" , 
        "3" = "#FF3300","4" = "#1f1aed",  "5" = "orange",
        "6" = "#8d00dc","7" = "yellow", "8" = "#66FF66")


plot1 <- SpatialDimPlot(mos.merge.harmony,
                        alpha = 1, pt.size.factor = 1.6, images="Samp_4001NT", cols = cols) +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot2 <- SpatialDimPlot(mos.merge.harmony,
                        alpha = 1, pt.size.factor = 1.6, images="Samp_3559NT", cols = cols)  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot3 <- SpatialDimPlot(mos.merge.harmony,
                        alpha = 1, pt.size.factor = 1.6, images="Samp_3115NT", cols = cols)  +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))

plot <- ggarrange(plot1, plot2, plot3, common.legend = FALSE, legend = "bottom", ncol=3) 


ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/visum_clust_res0.3_harmony.png",
       width = 25, height = 15, units = "cm")

plot <- DimPlot(mos.merge.harmony, reduction = "umap", group.by = "orig.ident",
        cols = c("4001" = "darkblue", "3559" = "#fbddf5", "3115" = "#039a87")) +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8)))& NoGrid() & NoAxes() #& NoLegend()

ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/visum_UMAP_res0.3_harmony_orig_ident.png",
       width = 15, height = 10, units = "cm")

plot <- DimPlot(mos.merge.harmony, reduction = "umap", group.by = "SCT_snn_res.0.3",
        cols = cols) +
  theme(legend.text = element_text(size=20)) +
  guides(fill=guide_legend(override.aes = list(size=8))) & NoGrid() & NoAxes() #& NoLegend()

ggsave(plot = plot, filename = "D:/Dropbox/11p15.5 mosaicism/Visium/visum_UMAP_res0.3_harmony_groups.png",
       width = 15, height = 10, units = "cm")


DoHeatmap(mos.merge.harmony, features = c("IGF2", "IGFBP1","IGFBP2","IGFBP3", 
                                  "APOC1", "APOC2", "ALB",
                                  "C3",  "F2", "TF",
                                  "VIM", "ACTA2", 
                                  "COL1A1", "THY1",  "ELN", 
                                  "TAGLN", "CNN1",  "PLN", 
                                  "MUC5B", "MUC6", 
                                  "GPC3", "GLUL","DLK1","AFP", "DUSP9","AXIN2", "HBB", "HBA1", "HBD"), group.colors = cols, group.bar.height = 0.03, size = 8) + NoLegend() + 
  scale_fill_viridis(option="mako") +
  theme(text = element_text(size = 16))



  # 4.3 Assign new cluster names------------------------------------------------
new.cluster.ids <- c("Non_mosaic_hepatocytes", 
                    "Mosaic_hepatocytes_3559",
                    "Portal_system", 
                    "Mosaic_hepatocytes_4001", 
                    "Mosaic_hepatocytes_3115", 
                    "Capsule_3559",
                    "Large_vessel", 
                    "Tumor",
                    "Capsule_4001")

old.cluster.ids = levels(mos.merge$seurat_clusters)


mos.merge@meta.data$new.cluster.ids <- plyr::mapvalues(x = mos.merge@meta.data$seurat_clusters, from = old.cluster.ids, to = new.cluster.ids)

names(new.cluster.ids) <- levels(mos.merge)
mos.merge <- RenameIdents(mos.merge, new.cluster.ids)

#mos.merge <- readRDS(file = "D:/Dropbox/11p15.5 mosaicism/Visium/merge_3115_3559_4001.rds")

  # 4.4 Visualize clusters with new names --------------------------------------
cols = c("Non_mosaic_hepatocytes" = "#99FFFF", 
         "Mosaic_hepatocytes_3559" = "orange",
         "Portal_system"="#3399FF", 
         "Mosaic_hepatocytes_4001"="#990000", 
         "Mosaic_hepatocytes_3115" = "#FF3300", 
         "Capsule_3559" = "#1f1aed",
         "Large_vessel"="yellow", 
         "Tumor"="#8d00dc",
         "Capsule_4001"="#66FF66")

SpatialDimPlot(mos.merge, cols = cols, alpha = 1, pt.size.factor = 1.6, image.alpha = 0)   & NoGrid() & NoAxes() & NoLegend()

my_levels = c("Non_mosaic_hepatocytes", 
              "Mosaic_hepatocytes_3559",
              "Mosaic_hepatocytes_3115", 
              "Mosaic_hepatocytes_4001", 
              "Portal_system",
              "Capsule_3559",
              "Capsule_4001",
              "Large_vessel", 
              "Tumor")

levels(mos.merge) = my_levels

DimPlot(mos.merge, cols = cols, label.size = 10)   & NoGrid() & NoAxes() #& NoLegend()

# Create source file umap figure 4c-d

Source_figure_4c <- mos.merge[["umap"]]@cell.embeddings
Source_figure_4cbis <- mos.merge[["orig.ident"]]

Source_figure_4c_3115 <- GetTissueCoordinates(mos.merge, image = "Samp_3115NT")
Source_figure_4c_3559 <- GetTissueCoordinates(mos.merge, image = "Samp_3559NT")
Source_figure_4c_4001 <- GetTissueCoordinates(mos.merge, image = "Samp_4001NT")
Source_figure_4ctiers <- rbind(Source_figure_4c_3115, Source_figure_4c_3559)
Source_figure_4ctiers <- rbind(Source_figure_4ctiers, Source_figure_4c_4001)
order <- rownames(Source_figure_4c) 
Source_figure_4ctiers = Source_figure_4ctiers[order(match(rownames(Source_figure_4ctiers), order)), ,drop = FALSE]

Source_figure_4c_all <- cbind(Source_figure_4c, Source_figure_4cbis)
Source_figure_4c_all <- cbind(Source_figure_4c_all, Source_figure_4ctiers)
write.table(Source_figure_4c_all, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/Source_figure_4c.txt", sep="\t")


# plot supp data signatures


p1 = SpatialFeaturePlot(obj_4001NT, alpha = 1, features =  "Stem", pt.size.factor = 1.8, image.alpha = 0) +  scale_fill_viridis(option="plasma") +
  theme(text = element_text(size = 8)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) & NoGrid() & NoAxes() 


plot <- plot_grid(p1, p2, p3, p4, p5, p6, p7,p8,p9,p10,p11,p12,nrow = 4, align = "hv")
#ggsave("/Users/jill/Desktop/plot_test.png", plot = plot, device = "png", width = 15, height = 28, units = "cm", dpi = 320)



  # 4.5 Apply different genes signatures ---------------------------------------

# Create features variable to look at 
features_Massalha = c("Stellate_cells_Massalha", "CAFs_Massalha","VSMC_Massalha","LVEC_Massalha",
                      "LSEC_Massalha", "LVECt_Massalha", "Pericytes_Massalha")

features_enrichment = c("de.4001.neg", "de.4001.pos","de.mos.sn.RNAseq.neg",
                        "de.mos.sn.RNAseq.pos", "de_RNAseq.pos","de_RNAseq.neg")[6]



# Do Violin plot with stack data to visualize multiple markers
# Use "fill.by" parameter to color by cell identity and not slide identity
VlnPlot(mos.merge, features = c("IGF2", "VIM"), 
        pt.size = 0, fill.by = "ident", cols = cols, stack = FALSE, flip = TRUE)


FeaturePlot(mos.merge, features = features_enrichment,
            #cols = rev(brewer.pal(11,"Spectral")),
            cols = c("lightgrey", "#E8EAF6", "#9FA8DA", "#3F51B5",  "darkblue"),
            pt.size = 0.8)

  # 4.6 Differential gene expression analysis for all 3 patients ---------------

resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/Visium/GSEA");if(!file.exists(resdir))	dir.create(resdir)

mos.merge = PrepSCTFindMarkers(mos.merge, assay = "SCT")
de_markers <- FindMarkers(mos.merge, ident.1 = c("Mosaic_hepatocytes_3559",
                                                      "Mosaic_hepatocytes_3115", 
                                                      "Mosaic_hepatocytes_4001"), ident.2 = "Non_mosaic_hepatocytes",
                              logfc.threshold = 0.1) 

write.table(de_markers,file.path(resdir,paste0("Mosaic_vs_non_mosaic_all_patients.txt")), sep="\t")

  # 4.7 GSEA enrichment visium -------------------------------------------------
    # 4.7.1 Merged GSEA with 4001 and 3559--------------------------------------

library(Biobase)
library(geco.NGS)
library(geco.supervised)
library(genefilter)
library(fgsea)
library(msigdbr)
library(geco.RNAseq)
source("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/GSEA_tools/Range_enrichment_SC.R")
library(geco.visu)
library(edgeR)

load("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/GSEA_tools/MSIG_v6.1.RData")
m_hum = msigdbr(species = "Homo sapiens")
#m_hum = m_hum %>%
#  filter(gs_cat %in% c("H","C6", "C7", "C8")) # Choose the categories you're interested in

m_list_hum = m_hum %>% split(x = .$gene_symbol, f = .$gs_name)


resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/Visium/GSEA");if(!file.exists(resdir))	dir.create(resdir)

## CREATE RANK
## Compute rank for infinite values (pval=0)
de_markers[which(de_markers$p_val_adj==0),"p_val_adj"]=1e-301 

de_markers$RANK<- de_markers$avg_log2FC * (-log10(de_markers$p_val_adj))
de_markers<- de_markers[order(de_markers$RANK,decreasing = T),]
de_markers<- de_markers[!is.na(de_markers$RANK),]

RANK_LIMMA<- de_markers$RANK
names.list = rownames(de_markers)
names(RANK_LIMMA)<- names.list

### RUN GSEA
fgseaRes_LIMMA <- fgsea(pathways = m_list_hum ,
                        stats = RANK_LIMMA,
                        minSize=25,
                        maxSize=200,
                        nperm=1000)


### LOOK at the RESULTS
head(fgseaRes_LIMMA[order(ES,decreasing = T), ])

RES_LIMMA<- as.data.frame(fgseaRes_LIMMA)
RES_LIMMA$leadingEdge <- vapply(RES_LIMMA$leadingEdge, paste, collapse = ", ", character(1L))
resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/Visium/GSEA");if(!file.exists(resdir))	dir.create(resdir)
#write.table(RES_LIMMA,file.path(resdir,paste0("FGGSEA_Mosaic_vs_non_mosaic_all_patients.txt")), sep="\t")

RANK2_LIMMA<- de_markers$RANK
names(RANK2_LIMMA)<- rownames(de_markers)

 
    # 4.7.2 Make graph enrichment plot -----------------------------------------

RES_LIMMA = RES_LIMMA %>% arrange(NES)
RES_LIMMA$pathway = factor(RES_LIMMA$pathway, levels = RES_LIMMA$pathway)
  
### MAKE FINAL GRAPH 
Enrich.plot = RES_LIMMA %>%
  filter(pathway %in% c("DESCARTES_FETAL_LIVER_HEPATOBLASTS",
                        "GOBP_PLATELET_DEGRANULATION",
                        "HALLMARK_COAGULATION",
                        "REACTOME_HEMOSTASIS",
                        "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
                        "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                       "REACTOME_TRANSPORT_OF_SMALL_MOLECULES",
                        "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
                        "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
                       "GOBP_LIPID_CATABOLIC_PROCESS"))


Enrich.plot$pathway = factor(Enrich.plot$pathway, levels = Enrich.plot$pathway)
ggplot(Enrich.plot, aes(x=pathway, y=ES, fill =pval)) + 
  geom_bar(color="black", position="stack", stat="identity", width = 0.8, size=1)+
  scale_fill_gradient(low = "red", 
                       high = "blue", space = "Lab",limits=c(0,0.04),breaks = c(0, 0.01, 0.02, 0.03, 0.04)) +
  scale_y_continuous(limits = c(0,1.1),breaks = c(0,1))+
  coord_flip()+
  theme_bw()+
  theme(legend.key = element_blank(),
        panel.border = element_rect(size = 1, fill = NA, linetype = "solid"),
        panel.background = element_blank(),
             legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0.1, "cm"),
  #      axis.text=element_blank()
  )

  # 4.8 Markers of each group -------------------------------------------------
#Supervised clustering
cols = c("Non_mosaic_hepatocytes" = "#99FFFF", 
         "Mosaic_hepatocytes_3559" = "orange",
         "Portal_system"="#3399FF", 
         "Mosaic_hepatocytes_4001"="#990000", 
         "Mosaic_hepatocytes_3115" = "#FF3300", 
         "Capsule_3559" = "#1f1aed",
         "Large_vessel"="yellow", 
         "Tumor"="#8d00dc",
         "Capsule_4001"="#66FF66")

my_levels = c("Non_mosaic_hepatocytes", 
              "Mosaic_hepatocytes_4001", 
              "Mosaic_hepatocytes_3559",
              "Mosaic_hepatocytes_3115", 
              "Portal_system",
              "Capsule_3559",
              "Capsule_4001",
              "Large_vessel", 
              "Tumor")

levels(mos.merge) = my_levels


DoHeatmap(mos.merge, features = c("IGF2", "IGFBP1","IGFBP2","IGFBP3", 
                                  "APOC1", "APOC2", "ALB",
                                  "C3",  "F2", "TF",
                                  "VIM", "ACTA2", 
                                  "COL1A1", "THY1",  "ELN", 
                                  "TAGLN", "CNN1",  "PLN", 
                                  "MUC5B", "MUC6", 
                                  "GPC3", "GLUL","DLK1","AFP", "DUSP9","AXIN2"), group.colors = cols, group.bar.height = 0.03, size = 1) + NoLegend() + 
  scale_fill_viridis(option="mako") +
  theme(text = element_text(size = 16))

genes <- c("IGF2", "IGFBP1","IGFBP2","IGFBP3", 
           "APOC1", "APOC2", "ALB",
           "C3",  "F2", "TF",
           "VIM", "ACTA2", 
           "COL1A1", "THY1",  "ELN", 
           "TAGLN", "CNN1",  "PLN", 
           "MUC5B", "MUC6", 
           "GPC3", "GLUL","DLK1","AFP", "DUSP9","AXIN2")

Source_figure_4e <- mos.merge@assays$SCT@scale.data[genes,]
write.table(Source_figure_4e, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/Source_figure_4e.txt", sep="\t")


#Clustering platelets markers
platelets_genes = my.sig %>% filter(GeneSet=="Activated_platelets")
platelets_genes =platelets_genes$Gene
my_genes = my.sig$Gene
DoHeatmap(mos.merge, features = platelets_genes, group.colors = cols, group.bar.height = 0.02, size = 5)  + 
  scale_fill_gradient2(low = "darkblue", mid="white", high="yellow")+
  theme(text = element_text(size = 15))

#saveRDS(mos.merge, file = "D:/Dropbox/11p15.5 mosaicism/Visium/merge_3115_3559_4001.rds")

  # 4.9 Progenitor/stem cell signature and Wnt/betacat--------------------------

cols = c("Non_mosaic_hepatocytes" = "#99FFFF", 
         "Mosaic_hepatocytes_3559" = "orange",
         "Portal_system"="#3399FF", 
         "Mosaic_hepatocytes_4001"="#990000", 
         "Mosaic_hepatocytes_3115" = "#FF3300", 
         "Capsule_3559" = "#1f1aed",
         "Large_vessel"="yellow", 
         "Tumor"="#8d00dc",
         "Capsule_4001"="#66FF66")

my_levels = c("Non_mosaic_hepatocytes", 
              "Mosaic_hepatocytes_4001", 
              "Mosaic_hepatocytes_3559",
              "Mosaic_hepatocytes_3115", 
              "Portal_system",
              "Capsule_3559",
              "Capsule_4001",
              "Large_vessel", 
              "Tumor")

levels(mos.merge) = my_levels

DefaultAssay(mos.merge)

# Create source data file supp figure 13
genes <- c("GPC3", "DLK1","EPCAM","KIT", 
           "AFP", "LIN28B", "THY1", "CD34", "TBX3", "GLUL","AXIN2","LGR5", 
           "LEF1")
Source_supp_figure_13 <- mos.merge@assays$SCT@scale.data[genes,]
write.table(Source_supp_figure_13, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/Source_supp_figure_13.txt", sep="\t")

DotPlot(mos.merge, features = c("GPC3", "DLK1","EPCAM","KIT", 
                                  "AFP", "LIN28B", "THY1", "CD34", "TBX3", "GLUL","AXIN2","LGR5", 
                                "LEF1"), cols=c("#fffb1c", "#ff1c1c")) +  #NoLegend() +
  theme(text = element_blank())


SCpubr::do_BoxPlot(sample = subset(mos.merge, new.cluster.ids %in% c("Non_mosaic_hepatocytes", "Mosaic_hepatocytes_4001", 
                                                                    "Mosaic_hepatocytes_3559", "Mosaic_hepatocytes_3115")),
                        feature = "GPC3", colors.use = cols,
                   use_test = TRUE,
                   comparisons = list(c("Non_mosaic_hepatocytes", "Mosaic_hepatocytes_4001"),
                                      c("Non_mosaic_hepatocytes", "Mosaic_hepatocytes_3559"),
                                      c("Non_mosaic_hepatocytes", "Mosaic_hepatocytes_3115")),
                   map_signif_level = FALSE, font.size=14, plot.title = NULL, 
                   xlab = "", test="wilcox.test", outlier.color = "gray80") 

# Source file supp figures 14-15
Source_CAF <- mos.merge[["CAF"]]
Source_HSCs <- mos.merge[["HSC"]]
Source_VSMC <- mos.merge[["VSMC"]]
Source_Hematopoiesis <- mos.merge[["Hematopoiesis"]]
Source_Cholangiocytes <- mos.merge[["Cholangiocytes"]]
Source_LVEC <- mos.merge[["LVEC"]]
Source_Periportal <- mos.merge[["Periportal"]]
Source_Perivenous <- mos.merge[["Perivenous"]]

genes <- c("VIM", "FAP", "ACTA2", "HGF", "LRAT", "RGS5", "COLEC11", "NGFR", "TAGLN", "CNN1", "PLN", "HBB", "GYPE",
           "KRT17", "SPP1", "CD34", "VWF", "CPE", "SLCO2A1", "CLEC14A", "TGM2", "CPS1", "ASS1", "PCK1", "CYP1A2", "CYP2E1")
  
  
Source_supp_figure_14 <- t(mos.merge@assays$SCT@scale.data[genes,])

order <- rownames(Source_supp_figure_14) 
coordinates = Source_figure_4ctiers[order(match(rownames(Source_figure_4ctiers), order)), ,drop = FALSE]

Source_supp_figure_14_all <- cbind(Source_supp_figure_14, coordinates)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_CAF)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_HSCs)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_VSMC)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_Hematopoiesis)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_Cholangiocytes)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_LVEC)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_Periportal)
Source_supp_figure_14_all <- cbind(Source_supp_figure_14_all, Source_Perivenous)

write.table(Source_supp_figure_14_all, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/Source_supp_figure_14_15.txt", sep="\t")


# Clean merged object for publication ------------------------------------------
mos.merge <- readRDS(file = "D:/Dropbox/11p15.5 mosaicism/Visium/merge_3115_3559_4001.rds")

mos.merge@meta.data$de.mos.sn.RNAseq.pos <- NULL
mos.merge@meta.data$de.mos.sn.RNAseq.neg <- NULL
mos.merge@meta.data$de_RNAseq.neg <- NULL
mos.merge@meta.data$de_RNAseq.pos <- NULL
mos.merge@meta.data$de.4001.pos <- NULL
mos.merge@meta.data$de.4001.neg <- NULL
mos.merge@meta.data$Cytotoxic_cells <- NULL
mos.merge@meta.data$T_cell_score <- NULL
mos.merge@meta.data$B_cell_score <- NULL
mos.merge@meta.data$endo_score <- NULL
mos.merge@meta.data$fibroblasts <- NULL
mos.merge@meta.data$NK <- NULL
mos.merge@meta.data$B_cells_Massalha <- NULL
mos.merge@meta.data$CAFs_Massalha <- NULL
mos.merge@meta.data$cDC1_Massalha <- NULL
mos.merge@meta.data$cDC2_Massalha <- NULL
mos.merge@meta.data$Kupffer_cells_Massalha <- NULL
mos.merge@meta.data$LVEC_Massalha <- NULL
mos.merge@meta.data$LVECt_Massalha <- NULL
mos.merge@meta.data$LSEC_Massalha <- NULL
mos.merge@meta.data$Pericytes_Massalha <- NULL
mos.merge@meta.data$SAMs_Massalha <- NULL
mos.merge@meta.data$Stellate_cells_Massalha <- NULL
mos.merge@meta.data$T_cells_Massalha <- NULL
mos.merge@meta.data$VSMC_dobie <- NULL
mos.merge@meta.data$HSC_payen <- NULL
mos.merge@meta.data$vSMC_Massalha <- NULL                   
mos.merge@meta.data$FB_dobie <- NULL
mos.merge@meta.data$HSC_dobie <- NULL
mos.merge@meta.data$TM1_Massalha <- NULL
mos.merge@meta.data$HSC1_payen <- NULL
mos.merge@meta.data$HSC2_payen <- NULL
mos.merge@meta.data$Xenobiotic_payen <- NULL
mos.merge@meta.data$Fatty_acid_payen <- NULL
mos.merge@meta.data$Retinoic_acid_payen <- NULL
mos.merge@meta.data$GAG_metabolism_payen <- NULL
mos.merge@meta.data$Elastic_fibers_payen <- NULL
mos.merge@meta.data$Secretion_payen <- NULL
mos.merge@meta.data$Kupffer_payen <- NULL
mos.merge@meta.data$de.3115.pos <- NULL
mos.merge@meta.data$de.3115.neg <- NULL
mos.merge@meta.data$SCT_snn_res.0.25 <- NULL
mos.merge@meta.data$SCT_snn_res.0.2 <- NULL
mos.merge@meta.data$VSMC_payen <- NULL
mos.merge@meta.data$Proliferation <- NULL
mos.merge@meta.data$Iron <- NULL
mos.merge@meta.data$seurat_clusters <- NULL
mos.merge@meta.data$LSEC <- NULL
mos.merge@meta.data$LT <- NULL
mos.merge@meta.data$Angiogenesis <- NULL
mos.merge@meta.data$DC <- NULL
mos.merge@meta.data$platelets <- NULL
mos.merge@meta.data$Progenitor <- NULL
mos.merge@meta.data$Stem <- NULL
mos.merge@meta.data$LB <- NULL

names(mos.merge@meta.data)

mos.merge@meta.data$cell.population <- mos.merge@meta.data$new.cluster.ids

table(mos.merge@meta.data$cell.population)

saveRDS(mos.merge, file = "D:/Dropbox/11p15.5 mosaicism/Visium/spatial_transcriptomics_merged_seurat.rds")

