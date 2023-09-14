### 0. Libraries and functions --------
library(geco.utils)
library(geco.visu)
library(ConsensusClusterPlus)
library(geco.unsupervised)
library(seriation)
library(scatterplot3d)
library(ggplot2)
library(magrittr)
library(clinfun)
library(DescTools)
library(Rmisc)
library(devtools)
library(RColorBrewer)
library(survival)
library(rms) # for Kapplan Meier
library(ggrepel)
library(tidyr)
library(scales)
library(matrixStats)
library(Hmisc)
library(ggforce)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(GSVA)
library(dplyr)
library(viridis)
library(gplots)
library(svglite)
library(car)

resdir <- "D:/Dropbox/11p15.5 mosaicism"
if(!file.exists(resdir))  dir.create(resdir)
source("D:/Dropbox/11p15.5 mosaicism not shared/color_extract.R", chdir = TRUE)

ConvertColors <- function(vec, old, new, all_other = NA) {
  if (length(old) != length(new)) stop("replacement is not the same length as original")
  if (!is.na(all_other)) vec[!vec %in% old] <- all_other
  for (i in 1:length(old)) {
    vec[vec %in% old[i]] <- new[i]
  }
  return(vec)
}

### 1. Load data -------------------------------------------------------------------------
annot <- geco.load("D:/Dropbox/11p15.5 mosaicism/Capsule_temporelle_Integrated_table/2022_05_04_integrated_pediatric_table.Rdata")
annot = as.data.frame(annot)
dim(annot)

annot. = annot %>%
  dplyr::rename(GPC3.Mut = GPC3) %>%
  dplyr::rename(CDKN1C.Mut = CDKN1C)


  # 1.1 Add RNAseq expression levels--------------------------------------------------------

gene_exp <-geco.load("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/mos_paper_by_gene_name/exp_190s.Rdata")

gene_exp = as.data.frame(t(gene_exp))
gene_exp$CHCID = rownames(gene_exp)
dim(gene_exp)
Periportal =c("HAL","SDS","PCK1","ASL","ASS1","CPS1")
Perivenous=c("GLUL","LGR5","AXIN2","CYP2E1","CYP1A2","OAT")

gene_exp = gene_exp[,c("IGF1", "IGF2", "H19", "IGF1R", "IGF2BP1", "IGF2BP3",  "CDKN1C", "KCNQ1", "KCNQ1OT1", Perivenous, Periportal, "CDH1", "TBX3", "APC", "ADH1A","HNF4A",  "ADH1B", "UGT2B7","IDH1", "BAAT",
                       "ADH1A","HNF4A",  "ADH1B", "UGT2B7","IDH1", "BAAT","CDH5", "CD160","CD3D", "CD8A",  "CD19", "HLA-DOA", "CD1A", "CD80", "CXCR1", "SPP1","KRT17",
                       "TWIST1", "VIM", "VCAN","MMP2", "MMP7", "KRT23","SOX9","PROM1", "KRT19", "THY1", "CD34", "KIT", "EPCAM","COL1A1", "COL1A2", "COL5A1", "COL5A2", "FAP",
                       "GPC3", "DLK1", "AFP","LIN28B","CCNA2", "CDC20", "BUB1","ANK1", "HBM", "HBG2", "HBE1","GYPE", "APOC1", "ALB", "SERPINE1", "MAPK8", "TBX3", "GLUL", "AXIN2", "LGR5", "LEF1")]
gene_exp$CHCID = rownames(gene_exp)


# Create annot. with RNAseq exp expression data

annot. = left_join(annot., gene_exp, by = "CHCID")
table(annot.$Histological.Diagnosis)

  # 1.2 Add IGF2 promoter use---------------------------------------------------

Prom <- read.delim("D:/Dropbox/11p15.5 mosaicism/IGF2_quantif_transcripts/df_ratio.txt", sep="\t", as.is=T, na.strings=c("",".","NA","na","#N/D","<NA>"))
dim(Prom)
rownames(Prom) = Prom$Sample
Prom = Prom %>%
  select(!c("Age", "Diag", "G1G6", "exp_IGF2")) %>%
  dplyr::rename(CHCID = Sample)

annot. = left_join(annot., Prom, by = "CHCID")

  # 1.3 Add IC1/ IC2 and allelic discrimination assay results ------------------
annot. = annot. %>%
  mutate(IC1_bval_merged = coalesce(IC1_b_value_MLPA, IC1_b_value),
         IC2_bval_merged = coalesce(IC2_b_value_MLPA, IC2_b_value)) 

Allelic.d <- read.delim("D:/Dropbox/11p15.5 mosaicism/MS-MLPA_allelic_discrimination/MS-MLPA_allelic_discrimination_recap.txt", sep="\t", as.is=T, na.strings=c("",".","NA","na","#N/D","<NA>"))
rownames(Allelic.d) = Allelic.d$CHCID

Allelic.d = Allelic.d %>%
  select(c(grep("rs", colnames(Allelic.d))), "CHCID")

annot. = left_join(annot., Allelic.d, by = "CHCID")


  # 1.4 Add methylation promoter data ------------------------------------------
meth_prom <- geco.load("D:/Dropbox/11p15.5 mosaicism/Promoter_methylation/mean_meth_in_promoters.RData")
rownames(meth_prom) = meth_prom$CHCID
annot. = left_join(annot., meth_prom, by = "CHCID")


  # 1.5  Filter annot ----------------------------------------------------------

`%nin%` = Negate(`%in%`)
annot. = annot. %>%
  mutate(MethRatio = IC1_bval_merged/IC2_bval_merged,
         Mosaic = ifelse(is.na(Mosaic), "wt", Mosaic))

annot.NT = annot. %>%
  filter(Papier_mosaiques_2021=="yes",
         grepl("N|S|FF", CHCID)) %>%
  mutate(Mosaic_alteration = case_when(CHCID=="CHC3168N" ~ "GAIN",
                                       CHCID %in% c("CHC3180N", "CHC3996N", "CHC2914N") ~ "LOM_IC2",
                                       CHCID %in% c("CHC4078N", "CHC3383N", "CHC3115N", "CHC3559N", "CHC3549N", "CHC3377N", "CHC4001N", "CHC3131N", "CHC3608N") ~ "cn-LOH",
                                       CHCID=="CHC3370N" ~ "wt",
                                       Mosaic=="wt" ~ "wt"),
         Mosaic_per_sample = case_when(CHCID=="CHC3370N" ~ "wt",
                                       Mosaic=="wt" ~ "wt",
                                       Mosaic=="Mosaic" ~ "Mosaic"),
         Age.surg.years = case_when(grepl("FF", CHCID) ~ (as.numeric(sub("SA,", "", sub("Foie fetal ", "", Miror.Bloc)))-40)/12.25,
                                    !is.na(Age.at.surgery.recode.months) ~ as.numeric(Age.at.surgery.recode.months)/12.25,
                                    T ~ NA_real_)) %>%
  as.data.frame()  

dim(annot.NT)
names = annot.NT$CHCID
table(annot.NT$Mosaic_per_sample)

# Create annot.NT2 with RNAseq data 

annot.NT2= annot.NT %>%
  filter(transcriptomic_mos=="yes",
         RNAseq=="yes") %>%
  filter(CHCID !="CHC2944N") %>% # bad quality
  filter(CHCID != "CHC3370N") # non mosaic and non normal

NDUFA11_manual_correct = c("CHC3101N", "CHC3090N")
NDUFB9_manual_correct = "CHC3369N"
TJP2_manual_correct = "CHC3092N"
ABCB11_manual_correct ="CHC3563N"
FAH_manual_correct = c("CHC3389N", "CHC3380N")
APC_manual_correct = c("CHC2967N", "CHC3538S", "CHC3970N")
GPC3_manual_correct = "CHC4008N"
HNF1A_manual_correct = c("CHC05237T", "CHC522N")
BRCA1_manual_correct = "CHC2964N"
BRCA2_manual_correct = c("CHC3017N", "CHC2824N")
PAX5_manual_correct ="CHC3369N"
AXIN1_manual_correct ="CHC3934N"

annot.NT2 = annot.NT2 %>% 
  mutate(NDUFA11.germ = ifelse(CHCID %in% NDUFA11_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         NDUFB9.germ = ifelse(CHCID %in% NDUFB9_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         TJP2.germ = ifelse(CHCID %in% TJP2_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         ABCB11.germ = ifelse(CHCID %in% ABCB11_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         FAH.germ = ifelse(CHCID %in% FAH_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         APC.germ = ifelse(CHCID %in% APC_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         GPC3.germ = ifelse(CHCID %in% GPC3_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         HNF1A.germ = ifelse(CHCID %in% HNF1A_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         BRCA1.germ = ifelse(CHCID %in% BRCA1_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         PAX5.germ = ifelse(CHCID %in% PAX5_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         BRCA2.germ = ifelse(CHCID %in% BRCA2_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)),
         AXIN1.germ = ifelse(CHCID %in% AXIN1_manual_correct, "germline", ifelse(WXS_serie=="yes", "wt", NA)))


cols.NT = c("FF" = "#004040", "HB" = "#CDADFF", "HCC" = "darkred", "FLC" = "#673AB7", "HCA" = "#81C784")

### 2. Perivenous Periportal correlations supp figure 11a -----------------------
  # 2.1 Prepare table ----------------------------------------------------------

annot.NT2$Age.surg.years = as.numeric(annot.NT2$Age.surg.years)
annot.PV.PP = annot.NT2 %>%
  mutate(PP_score = rowMeans(annot.NT2[,which(colnames(annot.NT2) %in% Periportal)], na.rm = T),
         PV_score = rowMeans(annot.NT2[,which(colnames(annot.NT2) %in% Perivenous)], na.rm = T)) %>%
  filter(Histological.Diagnosis=="HB") %>%
  filter(Mosaic_per_sample=="wt")

Periportal_col = "#1976D2"
Perivenous_col = "#AEDDD1"

dim(annot.PV.PP)

  # 2.2 Correlations periportal/perivenous with IGF2 expression supp figure 11a--

ggplot(annot.PV.PP, aes(x = IGF2, y = ASL)) +
  geom_point(shape = 20, size = 6, colour = Periportal_col) +
  geom_smooth(colour = "black", fill = "white", alpha =0.1, method = 'lm', size=1.5) +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank())

cor.test(annot.PV.PP$IGF2,annot.PV.PP$CYP1A2, method="pearson")

Source_supp_figure_11a <- annot.PV.PP %>% select(CHCID, IGF2, ASL, PCK1, SDS)
#write.table(Source_supp_figure_11a, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_11a.txt", sep="\t")

### 3. Supervised Clustering figure 1d -----------------------------------------
  # 3.1 Choose Plot specific genes ---------------------------------------------
genes_to_plot_list <- list(IGF2_related = c("IGF2", "H19","IGF1R", "IGFBP5","IGFBP6", "IGFBP3", "IGF1"),
                           Erythro = c("ANK1", "HBM", "HBG1"),
                           prolif = c("MKI67", "CDC20", "BUB1"),
                           Hepatoblast = c("GPC3", "DLK1","EPCAM", "KIT", "AFP","LIN28B"),
                           Stem = c("THY1", "CD34"),
                           CAFs = c("VIM", "FAP", "ACTA2"),
                           VSMCs = c("CNN1", "TAGLN"),
                           Col = c("COL1A1", "COL1A2"),
                           EMT = c("MMP2", "MMP9", "KRT23"),
                           diff.metabolism = c("HNF4A", "ADH1B", "UGT2B7"),
                           wnt.targets = c("TBX3", "GLUL", "AXIN2", "LGR5"))

  # 3.2 Table filtering ------------------------------------------------------------------
annot.NT2  = annot.NT2 %>%
  filter(Histological.Diagnosis %in% c("FF", "HB")) %>%
  mutate(Hist.diag = case_when(Mosaic_per_sample=="Mosaic" ~ "Mosaic",
                               Histological.Diagnosis=="FF" ~ "FF",
                               Mosaic_per_sample=="wt" ~ "HB")) %>%
  arrange((Age.surg.years)) %>%
  arrange(factor(Mosaic_per_sample, levels = c("Mosaic", "wt")))  %>%
  arrange(factor(Histological.Diagnosis, levels = c("FF", "HB"))) 

annot.NT2$CHCID = factor(annot.NT2$CHCID, levels = annot.NT2$CHCID)

supervised_order = annot.NT2$CHCID

exp <- geco.load("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/mos_paper_by_gene_name/exp_190s.Rdata")

col_vec = viridis(100)
col_vec = viridis(100, option = "plasma") # viridis, plasma, magma, inferno, cividis
center_input = c(T, F)[1] 
change_range = c(T, F)[1]

input = exp[unlist(genes_to_plot_list),] # Select gene of interest
input = input[,colnames(input) %in% supervised_order] # Semect samples of interest 
input = input[,match(annot.NT2$CHCID, colnames(input))] #Reorder input table with samples ordered


if (center_input) { # centering
  input = input - apply(input, 1, . %>% mean(., na.rm = T))
}

if (change_range) { # changing range
  for(i in 1:nrow(input)){ # change the range of each row (gene)
    input[i, ] = geco.changeRange(input[i, ], newmin = 0, newmax = 1)
  }
}

Source_figure_1d <- input
write.table(Source_figure_1d, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure1d.txt", sep="\t")

  # 3.3 plot each pathway ------------------------------------------------------ 
par(mar=c(0.5, 0.5, 0.5, 0.5), las=0)
p <- heatmap.2(as.matrix(input),
          density.info="none", 
          trace="none", 
          Colv=FALSE, 
          Rowv=FALSE, 
          col = col_vec,
          dendrogram="none", 
          # now the separations:
          rowsep=c(7, 10, 13,19,21,24,26,28,31,34),
          sepwidth=c(0.1,0.1),
          sepcolor = "white",
          labRow = T,
          labCol = F,
          srtRow = NULL,
          key=FALSE, 
          symkey=FALSE)

##ggsave("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/RNAseq.clustering.png", plot = p, device = "png", width = 10, height = 14, units = "cm", dpi = 320)


  # 3.4 Create annotations -----------------------------------------------------
# Create function 
color_annot_generated <- function(annot_table,
                                  col_interet = colnames(annot_table), 
                                  plotLegend = T, 
                                  plotLegendFile = NULL, 
                                  maxnumcateg =2,
                                  computer = c("Mac", "PC")[1],
                                  color.matrix = "D:/Dropbox/11p15.5 mosaicism not shared/Color_matrix_JP27012020.txt"){
  
  Color_matrix = data.frame(fread(color.matrix,
                                  sep = "\t", header= TRUE))
  
  qual_col_mat = factoall(Color_matrix)
  for (j in 1:ncol(annot_table)) annot_table[which(annot_table[, j] ==""),j] <- NA
  
  
  annot_table = annot_table[,which(is.element(colnames(annot_table), col_interet))]
  
  
  
  
  annot_colMat = unique(qual_col_mat[,c("Name_annot", "Type")])
  old_table = intersect(annot_colMat[,1], colnames(annot_table))
  vect_old_type = annot_colMat[match(old_table, annot_colMat[,"Name_annot"]),"Type"]
  
  
  new_annot = setdiff(colnames(annot_table), annot_colMat[,"Name_annot"])
  
  if(length(new_annot) > 0){		
    new_mat_annot = annot_table[,new_annot]
    if(length(new_annot) ==1){	
      anotype <- rep("categ", 1)
      names(anotype) <- colnames(new_mat_annot)
      classes <- class(new_mat_annot)
      nmodal <- length(unique(setdiff(new_mat_annot, NA)))
      
      
      
    }else{
      anotype <- rep("categ", ncol(new_mat_annot))
      names(anotype) <- colnames(new_mat_annot)
      classes <- sapply(1:ncol(new_mat_annot), function(j) class(new_mat_annot[,j]))
      nmodal <- sapply(1:ncol(new_mat_annot), function(j) length(unique(setdiff(new_mat_annot[,j], NA))))
      
    }
    
    anotype[which(classes %in% c("integer", "numeric") & nmodal > maxnumcateg)] <- "quantit"
    anotype[which(nmodal == 2)] <- "binary"
    
    vect_old_new = c(rep("old", length(old_table)), rep("new", length(new_annot)))
    
    tot_mat = data.frame(c(old_table, new_annot), c(vect_old_type, anotype), vect_old_new)
  }else{
    tot_mat = data.frame(old_table, vect_old_type, rep("old", length(old_table)))
  }
  
  colnames(tot_mat) = c("Name_annotation", "Type_annotation", "Old_new_annotation")
  tot_mat = tot_mat[match(colnames(annot_table), tot_mat[,"Name_annotation"]),]
  
  anocol <- annot_table
  
  if (plotLegend) 
    pdf(plotLegendFile)
  
  categCol <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", 
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", 
                "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca", 
                "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", 
                "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776")
  k = 1
  for (j in which(tot_mat[,"Type_annotation"] == "categ")) {
    
    tmp_annot = colnames(anocol)[j]
    tmp <- as.factor(anocol[, j])
    classes <- as.character(levels(tmp))
    
    if(length(classes) > 0){
      if(tot_mat[j,"Old_new_annotation"] == "old"){
        color_tmp <- qual_col_mat[which(qual_col_mat[,"Name_annot"] == tmp_annot),]
        if(length(classes) > 1){
          fill <- as.character(color_tmp[match(classes, color_tmp[,"Value_annot"]),"Color"])
        }else{
          fill <- as.character(color_tmp["Color"])
        }
        levels(tmp) <- fill
        anocol[, j] <- as.character(tmp)
        
      }else{
        ncat <- length(levels(tmp))
        if (k + ncat > length(categCol)) 
          categCol <- c(categCol, categCol)
        levels(tmp) <- categCol[k:(k + ncat - 1)]
        fill <- as.character(levels(tmp))
        anocol[, j] <- as.character(tmp)
        k <- k + ncat
      }
      
      if (plotLegend) {
        par(mar = c(0, 0, 0, 0))
        plot(-10, axes = F, xlim = c(0, 5), ylim = c(0, 5), 
             xlab = "", ylab = "")
        legend(1, 5, legend = classes, fill = fill, title = colnames(anocol)[j], 
               xjust = 0.5, yjust = 1)
      }
    }
    
    
  }
  
  memcol <- c()
  
  for (j in which(tot_mat[,"Type_annotation"] == "binary")) {
    tmp <- as.factor(anocol[,j])
    
    if(tot_mat[j,"Old_new_annotation"] == "old"){
      tmp_annot = colnames(anocol)[j]
      
      classes <- as.character(levels(tmp))
      color_tmp <- qual_col_mat[which(qual_col_mat[,"Name_annot"] == tmp_annot),]
      if(length(classes) > 1){
        fill <- as.character(color_tmp[match(classes, color_tmp[,"Value_annot"]),"Color"])
      }
      levels(tmp) <- fill
      
    }else{
      new <- setdiff(anocol[, j], c(NA, names(memcol)))
      if (length(new) == 2) {
        memcol <- c(memcol, c("#E7E8E2", "#1B1E26"))
        if(setequal(new, c("NM", "M"))){
          names(memcol)[(length(memcol) - 1):length(memcol)] <- c("NM", "M")
        }else if(setequal(new, c("0", "1"))){
          names(memcol)[(length(memcol) - 1):length(memcol)] <- c("0", "1")
        }else if(setequal(new, c("no", "yes"))){
          names(memcol)[(length(memcol) - 1):length(memcol)] <- c("no", "yes")
        }else{
          names(memcol)[(length(memcol) - 1):length(memcol)] <- sort(new)
        }
      }else if (length(new) == 1){
        memcol <- c(memcol, "#E7E8E2")
        if(new == "NM"){
          names(memcol)[length(memcol)] <- "NM"
        }
      }
      classes <- intersect(names(memcol), annot_table[, j])
      fill <- memcol[is.element(names(memcol), classes)]
      
      levels(tmp) <- fill[match(levels(tmp), names(fill))]
      
    }
    anocol[, j] <- as.character(tmp)
    if (plotLegend) {
      par(mar = c(0, 0, 0, 0))
      plot(-10, axes = F, xlim = c(0, 5), ylim = c(0, 5), 
           xlab = "", ylab = "")
      
      legend(1, 5, legend = classes, fill = fill, title = colnames(anocol)[j], 
             xjust = 0.5, yjust = 1)
    }
  }  	
  
  quantitCol <- c("#57001F", "midnightblue", "firebrick4", "deeppink4")
  
  k <- 1
  
  for (j in which(tot_mat[,"Type_annotation"] == "quantit")) {
    
    if(tot_mat[j,"Old_new_annotation"] == "old"){
      col_quant_tmp = as.character(qual_col_mat[match(tot_mat[j,"Name_annotation"], qual_col_mat[,"Name_annot"]),"Color"])
    }else{
      col_quant_tmp = quantitCol[k]
      if (k < length(quantitCol)) {
        k <- k + 1
      }
      else {
        k <- 1
      }
    }
    # changer la selection des couleur quantitative pour les déja enregistrées
    
    colrange <- colorRampPalette(c("white", col_quant_tmp))(17)
    
    anocol[, j] <- colrange[round(geco.changeRange(as.numeric(anocol[,j], newmin = 1, newmax = 50)))]
    
    if (plotLegend) {
      par(mar = c(8, 2, 5, 1))
      lims <- seq(-1, 1, length.out = 200)
      image(matrix(lims, nc = 1), col = colrange, axes = F, xlab = colnames(anocol)[j])
    }   
  }    
  
  
  if (plotLegend) 
    dev.off()
  for (j in 1:ncol(anocol)) anocol[which(is.na(anocol[, j])), j] <- "white"
  return(as.matrix(anocol))
  
}

# Generate annotations ---------------------------------------------------------

annotCol = c("Mosaic", "Mosaic_alteration")
cols <-c("Mosaic" = "Indianred", "wt" = "white", "GAIN" = "#FFC107", "LOM_IC2" = "#EF9A9A","cn-LOH" = "#B71C1C")

anocol = color_annot_generated(annot.NT2[, annotCol], col_interet = annotCol, plotLegend = T, computer = "PC",
                               color.matrix = "D:/Dropbox/11p15.5 mosaicism not shared/Color_matrix_JP27012020.txt")

pdf(file.path(resdir, "Oncoprint.pdf"), width = 5.04, height = 4)
geco.imageCol(anocol, xlab.cex = 0, ylab.cex = 0,drawLines = c("h"))
dev.off()

Source_figure_1dbis <- annot.NT2[, c("CHCID", "Mosaic_alteration", "Hist.diag", "Age.surg.years")]
write.table(Source_figure_1dbis, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_1dbis.txt", sep="\t")

  # 3.5 Add Age at surgery  - figure 1d annotation -----------------------------
cols = c("FF" = "#004040", "Mosaic"="indianred", "HB" = "#CDADFF")

ggplot(annot.NT2, aes(x=CHCID, Age.surg.years), y=Age.surg.years, color=Hist.diag) + 
  geom_point(aes(color=Hist.diag), size = 3)+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =0.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank(),
       axis.text.y=element_blank(),
        axis.text.x = element_blank() )

  # 3.6 Correlations Wnt/bcat target genes and IGF2 expression - supp figure 8 -

annot.bcat.targets = annot.NT2

cols = c("FF" = "#004040", "Mosaic"="indianred", "HB" = "#CDADFF")
dim(annot.bcat.targets)

Source_supp_figure_8 <- annot.bcat.targets[, c("CHCID", "Hist.diag", "TBX3", "GLUL", "AXIN2", "LGR5")]
write.table(Source_supp_figure_8, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_8.txt", sep="\t")

ggplot(annot.bcat.targets, aes(x = IGF2, y = TBX3, colour = Hist.diag)) +
  geom_point(shape = 20, size = 6) +
  scale_color_manual(values = cols) + 
  geom_smooth(colour = "black", fill = "white", alpha =0.1, method = 'lm', size=1.5) +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black", size =.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank())

cor.test(annot.bcat.targets$IGF2,annot.bcat.targets$TBX3, method="pearson")

ggplot(annot.bcat.targets, aes(x = IGF2, y = GLUL, colour = Hist.diag)) +
  geom_point(shape = 20, size = 6) +
  scale_color_manual(values = cols) + 
  geom_smooth(colour = "black", fill = "white", alpha =0.1, method = 'lm', size=1.5) +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black", size =.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank())

cor.test(annot.bcat.targets$IGF2,annot.bcat.targets$GLUL, method="pearson")

ggplot(annot.bcat.targets, aes(x = IGF2, y = AXIN2, colour = Hist.diag)) +
  geom_point(shape = 20, size = 6) +
  scale_color_manual(values = cols) + 
  geom_smooth(colour = "black", fill = "white", alpha =0.1, method = 'lm', size=1.5) +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black", size =.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank())

cor.test(annot.bcat.targets$IGF2,annot.bcat.targets$AXIN2, method="pearson")

ggplot(annot.bcat.targets, aes(x = IGF2, y = LGR5, colour = Hist.diag)) +
  geom_point(shape = 20, size = 6) +
  scale_color_manual(values = cols) + 
  geom_smooth(colour = "black", fill = "white", alpha =0.1, method = 'lm', size=1.5) +
  theme_bw() +
  labs(x="IGF2 exp", y="LGR5 exp")+
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(size=15, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text.y = element_text(size=15, margin = margin(t = 0, r = 0, b = 0, l = 10)),,
        axis.title = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black", size =.5),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank())

cor.test(annot.bcat.targets$IGF2,annot.bcat.targets$LGR5, method="pearson")

### 4. Compare gene expression in the different groups -------------------------
  # 4.1 IGF2 expression in transcriptomic groups -------------------------------
dim(annot.NT2)
cols = c("FF" = "#004040", "Mosaic"="indianred", "HB" = "#CDADFF")

annot.NT2$Hist.diag = factor(annot.NT2$Hist.diag, levels = c("FF", "Mosaic", "HB"))
par(mar=c(4, 4, 4, 4), mgp=c(1, 1, 0), las=0) 
boxplot(annot.NT2$IGF2~ annot.NT2$Hist.diag ,col=cols, ylab = c(""), xaxt = "n",  ann=F, ylim=c(min(annot.NT2$IGF2-1), 1+max(annot.NT2$IGF2)))

#Multiple linear regression model to adjust for age

annot.stat = annot.NT2
annot.stat$Hist.diag = factor(annot.stat$Hist.diag, levels = c("Mosaic", "FF", "HB"))
model2 <- lm(IGF2 ~ Age_at_surg_AM + Hist.diag, data = annot.stat)
summary(model2)

  # 4.2 Wnt targets expression in ML, NML and FL - supp figure 8a --------------

annot.NT2$Hist.diag = factor(annot.NT2$Hist.diag, levels = c("FF", "Mosaic", "HB"))
par(mar=c(4, 4, 4, 4), mgp=c(1, 1, 0), las=0) 
boxplot(annot.NT2$TBX3~ annot.NT2$Hist.diag ,col=cols, ylab = c(""), xaxt = "n",  ann=F, ylim=c(min(annot.NT2$TBX3-1), 1+max(annot.NT2$TBX3)))

annot.NT2$Hist.diag = factor(annot.NT2$Hist.diag, levels = c("FF", "Mosaic", "HB"))
par(mar=c(4, 4, 4, 4), mgp=c(1, 1, 0), las=0) 
boxplot(annot.NT2$GLUL~ annot.NT2$Hist.diag ,col=cols, ylab = c(""), xaxt = "n",  ann=F, ylim=c(min(annot.NT2$GLUL-1), 1+max(annot.NT2$GLUL)))

annot.NT2$Hist.diag = factor(annot.NT2$Hist.diag, levels = c("FF", "Mosaic", "HB"))
par(mar=c(4, 4, 4, 4), mgp=c(1, 1, 0), las=0) 
boxplot(annot.NT2$AXIN2~ annot.NT2$Hist.diag ,col=cols, ylab = c(""), xaxt = "n",  ann=F, ylim=c(min(annot.NT2$AXIN2-1), 1+max(annot.NT2$AXIN2)))

annot.NT2$Hist.diag = factor(annot.NT2$Hist.diag, levels = c("FF", "Mosaic", "HB"))
par(mar=c(4, 4, 4, 4), mgp=c(1, 1, 0), las=0) 
boxplot(annot.NT2$LGR5~ annot.NT2$Hist.diag ,col=cols, ylab = c(""), xaxt = "n",  ann=F, ylim=c(min(annot.NT2$LGR5-1), 1+max(annot.NT2$LGR5)))

annot.NT2$Hist.diag = factor(annot.NT2$Hist.diag, levels = c("FF", "Mosaic", "HB"))
par(mar=c(4, 4, 4, 4), mgp=c(1, 1, 0), las=0) 
boxplot(annot.NT2$LEF1.y ~ annot.NT2$Hist.diag ,col=cols, ylab = c(""), xaxt = "n",  ann=F, ylim=c(min(annot.NT2$LEF1.y-1), 1+max(annot.NT2$LEF1.y)))

  # 4.3 PV and PP scores in mosaic and non mosaic livers - supp figure 11b-------
annot.stat = annot.NT2 %>%
  mutate(PP_score = rowMeans(annot.NT2[,which(colnames(annot.NT2) %in% Periportal)], na.rm = T),
         PV_score = rowMeans(annot.NT2[,which(colnames(annot.NT2) %in% Perivenous)], na.rm = T),
         Hist.diag.with.age = case_when(Hist.diag=="HB" & Age.surg.years<=3.3 ~ "HB_young",
                                         Hist.diag=="HB" & Age.surg.years>3.3 ~ "HB_old", 
                                         Hist.diag=="Mosaic" ~ "Mosaic", 
                                         Hist.diag=="FF" ~ "FF")) 

# Statistical test
annot.stat2 = annot.stat %>%
  filter(Hist.diag.with.age %in% c("Mosaic", "HB_old"))
wilcox.test(annot.stat2$PP_score ~ annot.stat2$Hist.diag.with.age)

annot.stat = annot.stat %>%
  select(c("CHCID", "Hist.diag.with.age", "PV_score", "PP_score", "Age.surg.years")) %>%
  pivot_longer(cols=c("PP_score", "PV_score"), names_to = "PPPV", values_to = "score")

Periportal_col = "#1976D2"
Perivenous_col = "#AEDDD1"

annot.stat$Hist.diag.with.age = factor(annot.stat$Hist.diag.with.age, levels = c("FF", "Mosaic","HB_young", "HB_old"))
annot.stat$PPPV = factor(annot.stat$PPPV, levels = c("PV_score", "PP_score"))

Source_supp_figure_11b <- annot.stat 
write.table(Source_supp_figure_11b, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_11b.txt", sep="\t")


ggplot(annot.stat, aes(x=Hist.diag.with.age, y=score)) +
  geom_boxplot(aes(fill=PPPV), color = "black", width=4, outlier.alpha = .3) +
  scale_fill_manual(values = c(Perivenous_col, Periportal_col)) +
  scale_y_continuous("Score", limits=c(6,20))+
  scale_x_discrete("") +
  theme_bw() +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=15))
  

### 5. Differential gene expression analysis -----------------------------------
  # 5.1  Mosaic livers vs non mosaic livers-------------------------------------
    # 5.1.1 Load packages ------------------------------------------------------
library(readxl)
library(Biobase)
library(limma)
library(geco.NGS)
library(geco.supervised)
library(genefilter)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(geco.RNAseq)
source("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/GSEA_tools/Range_enrichment_SC.R")
library(geco.visu)
library(edgeR)

    # 5.1.2 load data ----------------------------------------------------------

exp <-geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/exp_190s.Rdata")
fpkm <- geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/fpkm_classical_190s.Rdata")
gene_table <-  geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/gene_table_filtered_190s.Rdata")
count_matrix <- geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/count_matrix_190s.Rdata")

all(rownames(fpkm)==gene_table$gene_name)

load("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/GSEA_tools/MSIG_v6.1.RData")

    # 5.1.3 Select samples of interest -----------------------------------------
annot.GSEA = annot.NT2 %>%
  filter(Histological.Diagnosis =="HB")

annot.GSEA$CHCID = as.character(annot.GSEA$CHCID)
exp = exp[,annot.GSEA$CHCID]
count_matrix = count_matrix[,annot.GSEA$CHCID]
dim(count_matrix)

    # 5.1.4 Select protein coding genes ----------------------------------------

gene_table<-gene_table[which(gene_table$gene_type=="protein_coding"),] %T>%
{print(dim(.))}

    # 5.1.5 Limma differential gene expression  --------------------------------

Counts=count_matrix
samples = colnames(Counts)
gene_names = rownames(Counts)
y <- DGEList(counts = Counts, genes = as.vector(gene_names))
A <- rowSums(y$counts)
isexpr <- A > 30
hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y)

y <- calcNormFactors(y)

annot.GSEA$Mosaic_per_sample= factor(annot.GSEA$Mosaic_per_sample,levels = c("wt", "Mosaic"))  # Choose Mosaic for mosaic per patient or Mosaic per sample

# Correct for age-------------------------------------------------------------------------
design <- model.matrix(~0+Mosaic_per_sample+Age_at_surg_AM, data = data.frame(annot.GSEA))
v <- voom(y, design, plot = F)

cont.matrix <- makeContrasts( 
  mos_vs_HB =  Mosaic_per_sampleMosaic - Mosaic_per_samplewt,  #  for mosaic vs NT HB : mos_vs_wt =  Mosaic_per_sampleMosaic - Mosaic_per_samplewt,
  levels=design)
fit <- lmFit(v, design) 
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.eb   <- eBayes(fit.cont)
results_tot = topTable(fit.eb, adjust="BH", n=nrow(y$counts))
all_resullts_limma = results_tot
resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/Limma/");if(!file.exists(resdir))	dir.create(resdir)
all_resullts_limma_row = cbind(Gene = rownames(all_resullts_limma), all_resullts_limma)
#write.table(all_resullts_limma_row, file = paste(resdir, "Mosaic_vs_Non_Mos_HB_without_CHC3370N_all_ages_age_adjusted.txt", sep = ""), sep = "\t", col.names = T, row.names = F)

dim(all_resullts_limma)

  # 5.2  Fetal livers vs mosaic livers------------------------------------------
    # 5.2.1 load data --------------------------------------------------------------

exp <-geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/exp_190s.Rdata")
fpkm <- geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/fpkm_classical_190s.Rdata")
gene_table <-  geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/gene_table_filtered_190s.Rdata")
count_matrix <- geco.load("//10.93.23.19/HEPATO PARTAGE/GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/count_matrix_190s.Rdata")

load("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/GSEA_tools/MSIG_v6.1.RData")

    # 5.2.2 Select samples of interest ---------------------------------------------
annot.GSEA = annot.NT2 %>%
  filter(Hist.diag %in% c("FF", "Mosaic"))

annot.GSEA$CHCID = as.character(annot.GSEA$CHCID)
exp = exp[,annot.GSEA$CHCID]
count_matrix = count_matrix[,annot.GSEA$CHCID]
dim(count_matrix)

    # 5.2.3 Select protein coding genes ----------------------------------------

gene_table<-gene_table[which(gene_table$gene_type=="protein_coding"),] %T>%
  {print(dim(.))}

    # 5.2.4 Limma differential gene expression  --------------------------------

Counts=count_matrix
samples = colnames(Counts)
gene_names = rownames(Counts)
y <- DGEList(counts = Counts, genes = as.vector(gene_names))
A <- rowSums(y$counts)
isexpr <- A > 30
hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, , keep.lib.size = FALSE]
dim(y)

y <- calcNormFactors(y)

annot.GSEA$Hist.diag= factor(annot.GSEA$Hist.diag,levels = c("Mosaic", "FF"))  


# Correct for age-------------------------------------------------------------------------
design <- model.matrix(~0+Hist.diag+Age_at_surg_AM, data = data.frame(annot.GSEA))
v <- voom(y, design, plot = F)

cont.matrix <- makeContrasts( 
  mos_vs_FF =   Hist.diagFF - Hist.diagMosaic ,  
  levels=design)
fit <- lmFit(v, design) 
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.eb   <- eBayes(fit.cont)
results_tot = topTable(fit.eb, adjust="BH", n=nrow(y$counts))
all_resullts_limma = results_tot
resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/Limma/");if(!file.exists(resdir))	dir.create(resdir)
#all_resullts_limma_row = cbind(Gene = rownames(all_resullts_limma), all_resullts_limma)
write.table(all_resullts_limma_row, file = paste(resdir, "FL_vs_mosaic_livers_without_CHC3370N_all_ages_age_adjusted.txt", sep = ""), sep = "\t", col.names = T, row.names = F)

dim(all_resullts_limma)

### 6. Enrichment LIMMA --------------------------------------------------------
  # 6.1 Load genesets database--------------------------------------------------
load("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/GSEA_tools/MSIG_v6.1.RData")
m_hum = msigdbr(species = "Homo sapiens")
#m_hum = m_hum %>%
#  filter(gs_cat %in% c("H","C6", "C7", "C8")) # Choose the categories you're interested in

m_list_hum = m_hum %>% split(x = .$gene_symbol, f = .$gs_name)


  # 6.2 Create rank ------------------------------------------------------------
all_resullts_limma$RANK<- all_resullts_limma$logFC * (-log10(all_resullts_limma$P.Value))
all_resullts_limma<- all_resullts_limma[order(all_resullts_limma$RANK,decreasing = T),]
all_resullts_limma<- all_resullts_limma[!is.na(all_resullts_limma$RANK),]

RANK_LIMMA<- all_resullts_limma$RANK
names.list = rownames(all_resullts_limma)
names(RANK_LIMMA)<- names.list

  # 6.3 Run fgsea---------------------------------------------------------------
fgseaRes_LIMMA <- fgsea(pathways = m_list_hum ,
                        stats = RANK_LIMMA,
                        minSize=25,
                        maxSize=200,
                        nperm=1000)


# Look at the results
head(fgseaRes_LIMMA[order(ES,decreasing = T), ])

RES_LIMMA<- as.data.frame(fgseaRes_LIMMA)
RES_LIMMA$leadingEdge <- vapply(RES_LIMMA$leadingEdge, paste, collapse = ", ", character(1L))
resdir <- paste0("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/Limma/");if(!file.exists(resdir))	dir.create(resdir)
write.table(RES_LIMMA,file.path(resdir,paste0("FGGSEA_LIMMA_Mosaic_vs_non_mos_HB_without_3370.txt")), sep="\t")

RANK2_LIMMA<- all_resullts_limma$RANK
names(RANK2_LIMMA)<- rownames(all_resullts_limma)


  # 6.4 Plot enriched pathways -------------------------------------------------
 
Choosen_pathways = c("REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
                     "HALLMARK_ANGIOGENESIS", 
                     "REACTOME_COLLAGEN_FORMATION",
                     "HP_ARTERIOVENOUS_MALFORMATION",
                     "TSENG_ADIPOGENIC_POTENTIAL_UP",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                     "DESCARTES_MAIN_FETAL_SMOOTH_MUSCLE_CELLS",
                     "DESCARTES_FETAL_HEART_STROMAL_CELLS",
                     "GOBP_PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
                     "GOBP_PLATELET_ACTIVATION",
                     "DESCARTES_FETAL_LIVER_STELLATE_CELLS",
                     "GOBP_DRUG_CATABOLIC_PROCESS",
                     "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
                     "GOBP_GLUTAMATE_METABOLIC_PROCESS",
                     "KEGG_RETINOL_METABOLISM",
                     "GOBP_SERINE_FAMILY_AMINO_ACID_METABOLIC_PROCESS")

Enrich.plot = RES_LIMMA %>%
  filter(pathway %in% Choosen_pathways) %>%
  arrange(ES)

Enrich.plot$pathway = factor(Enrich.plot$pathway, levels = Enrich.plot$pathway)


ggplot(Enrich.plot, aes(x=pathway, y=ES)) + 
  geom_bar(aes(fill=pval), position="stack", stat="identity", width = 0.8, size=1)+
  scale_fill_gradient(low = "red", 
                      high = "blue", space = "Lab" ) +
  scale_y_continuous(limits = c(-1.1,1.1),breaks = c(-1,-0.5,0, 0.5, 1))+
  coord_flip()+
  theme_bw()+
  theme(legend.key = element_blank(),
        panel.border = element_rect(size = 1, fill = NA, linetype = "solid"),
        panel.background = element_blank(),
      #      legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0.1, "cm"),
      axis.text=element_blank()
  )

