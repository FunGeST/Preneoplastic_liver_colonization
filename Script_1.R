### 0. Libraries and functions -------------------------------------------------
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
library(ggrepel)
library(tidyr)
library(scales)
library(matrixStats)
library(Hmisc)
library(ggforce)
library(ggpubr)
library(plyr)
library(corrplot)
library(ComplexHeatmap)
library(GSVA)
library(dplyr)
library(viridis)
library(ggridges)
library(readxl)
library(data.table)


resdir <- "D:/Dropbox/11p15.5 mosaicism"
if(!file.exists(resdir))  dir.create(resdir)

ConvertColors <- function(vec, old, new, all_other = NA) {
  if (length(old) != length(new)) stop("replacement is not the same length as original")
  if (!is.na(all_other)) vec[!vec %in% old] <- all_other
  for (i in 1:length(old)) {
    vec[vec %in% old[i]] <- new[i]
  }
  return(vec)
}

### 1. Load data ---------------------------------------------------------------
annot <- read_xlsx(file.path("D:/Dropbox/11p15.5 mosaicism/Capsule_temporelle_Integrated_table/", "Pediatric_table.xlsx")) %>% as.data.frame()
annot. = as.data.frame(annot)

dim(annot.)
table(annot.$Type)

annot. = annot. %>%
  dplyr::rename(GPC3.Mut = GPC3) %>%
  dplyr::rename(CDKN1C.Mut = CDKN1C)

  # 1.1 Add RNAseq expression levels--------------------------------------------

gene_exp <-geco.load("D:/Dropbox/11p15.5 mosaicism/RNAseq_clust/mos_paper_by_gene_name/exp_190s.Rdata")

gene_exp = as.data.frame(t(gene_exp))
gene_exp$CHCID = rownames(gene_exp)
dim(gene_exp)

gene_exp = gene_exp[,c("IGF1", "IGF2", "H19", "IGF1R", "IGF2BP1", "IGF2BP3",  
                       "CDKN1C", "KCNQ1", "KCNQ1OT1","ASCL2","AFP", "DLK1", 
                       "GPC3", "LIN28B", "EPCAM", "KIT","CD34", "THY1", "APC", 
                       "TBX3", "GLUL", "LEF1", "AXIN2", "LGR5")]
gene_exp$CHCID = rownames(gene_exp)


# Create annot. with RNAseq exp expression data

annot. = left_join(annot., gene_exp, by = "CHCID")

  # 1.2 Add IGF2 promoter use---------------------------------------------------

# Add proportion
Prom <- read.delim("D:/Dropbox/11p15.5 mosaicism/IGF2_quantif_transcripts/df_ratio.txt", sep="\t", as.is=T, na.strings=c("",".","NA","na","#N/D","<NA>"))
dim(Prom)
rownames(Prom) = Prom$Sample
Prom = Prom %>%
  select(!c("Age", "Diag", "G1G6", "exp_IGF2")) %>%
  dplyr::rename(CHCID = Sample)

annot. = left_join(annot., Prom, by = "CHCID")

#Add reads quantification 
Prom.reads <- read.csv("D:/Dropbox/11p15.5 mosaicism/IGF2_quantif_transcripts/isoform_estimate_2021.txt", sep=",", as.is=T, na.strings=c("",".","NA","na","#N/D","<NA>"))
rownames(Prom.reads) = Prom.reads$Sample
Prom.reads= rename(Prom.reads, p0p1_reads = prom_0_p1)
Prom.reads= rename(Prom.reads, p1_reads = prom_1)
Prom.reads= rename(Prom.reads, p2_reads = prom_2)
Prom.reads= rename(Prom.reads, p3_reads = prom_3)
Prom.reads= rename(Prom.reads, p4_reads = prom_4)

Prom.reads = Prom.reads %>%
  dplyr::rename(CHCID = Sample)

annot. = left_join(annot., Prom.reads, by = "CHCID")

  # 1.3  Add IC1/ IC2 and allelic discrimination assay results ---------------------------
annot. = annot. %>%
mutate(IC1_bval_merged = coalesce(IC1_b_value_MLPA, IC1_b_value),
       IC2_bval_merged = coalesce(IC2_b_value_MLPA, IC2_b_value)) 

Allelic.d <- read.delim("D:/Dropbox/11p15.5 mosaicism/MS-MLPA_allelic_discrimination/MS-MLPA_allelic_discrimination_recap.txt", sep="\t", as.is=T, na.strings=c("",".","NA","na","#N/D","<NA>"))
rownames(Allelic.d) = Allelic.d$CHCID

Allelic.d = Allelic.d %>%
  select(c(grep("rs", colnames(Allelic.d))), "CHCID")

annot. = left_join(annot., Allelic.d, by = "CHCID")


  # 1.4 Add methylation promoter data ------------------------------------------------------
meth_prom <- geco.load("D:/Dropbox/11p15.5 mosaicism/Promoter_methylation/mean_meth_in_promoters.RData")
rownames(meth_prom) = meth_prom$CHCID
annot. = left_join(annot., meth_prom, by = "CHCID")

annot.$fetal_mean_meth = apply(annot.[,c("mean_P2", "mean_P3", "mean_P4")], 1, mean)


  # 1.5 Add new annotations ----------------------------------------------------

`%nin%` = Negate(`%in%`)
annot. = annot. %>%
  mutate(MethRatio = IC1_bval_merged/IC2_bval_merged,
         Mosaic = ifelse(is.na(Mosaic), "wt", Mosaic))
  

# Add germline mutations
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

annot. = annot. %>% 
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



# Create Mosaic alteration annotation
annot.[annot.$CHCID=="CHC3538S", "Age.at.surgery.recode.months"] <- 126.2
annot.[annot.$CHCID=="CHC3122S", "Age.at.surgery.recode.months"] <- 190.9

annot. = annot. %>%
  mutate(Mosaic_alteration = case_when(CHCID =="CHC3168N" ~ "GAIN",
                                       CHCID %in% c("CHC3180N", "CHC3996N", "CHC2914N") ~ "LOM_IC2",
                                       CHCID %in% c("CHC4078N", "CHC3383N", "CHC3115N", "CHC3559N", "CHC3549N", "CHC3377N", "CHC4001N", "CHC3131N", "CHC3608N") ~ "cn-LOH",
                                       CHCID =="CHC3370N" ~ "wt",
                                       Mosaic =="wt" ~ "wt"),
         Cell_fraction = case_when(CHCID=="CHC3383N" ~ 37,
                                   CHCID=="CHC3377N" ~ 15,
                                   CHCID=="CHC3549N" ~ 6,
                                   CHCID=="CHC4001N" ~ 6,
                                   CHCID=="CHC3559N" ~ 30,
                                   CHCID=="CHC3115N" ~ 58,
                                   CHCID=="CHC3131N" ~ 3,
                                   CHCID=="CHC3168N" ~ 10,
                                   CHCID=="CHC4078N" ~ 33,
                                   CHCID=="CHC3608N" ~ 48,
                                   CHCID=="CHC3180N" ~ 26,
                                   CHCID=="CHC3996N" ~ 20,
                                   CHCID=="CHC2914N" ~ 56),
          Mosaic.per.sample = case_when(CHCID=="CHC3370N" ~ "wt",
                                        Mosaic=="wt" ~ "wt",
                                        Mosaic=="Mosaic" ~ "Mosaic"),
         Hist.diag= case_when(CHCID=="CHC3370N" ~ "HB",
                              Mosaic =="Mosaic" ~ "Mosaic",
                              Histological.Diagnosis=="HCC" ~ "HCC",
                              Histological.Diagnosis=="HB" ~ "HB",
                              Histological.Diagnosis=="FLC" ~ "FLC",
                              Histological.Diagnosis=="HCA" ~ "HCA",
                              Histological.Diagnosis=="FF" ~ "FF"),
          Age.surg.years = case_when(grepl("FF", CHCID) ~ (as.numeric(sub("SA,", "", sub("Foie fetal ", "", Miror.Bloc)))-40)/12.25,
                                  !is.na(Age.at.surgery.recode.months) ~ as.numeric(Age.at.surgery.recode.months)/12.25,
                                  T ~ NA_real_), 
         Histological.Diagnosis2 = case_when(CHCID=="CHC3370N" ~ "HB",
                                             Mosaic_alteration=="cn-LOH" ~ "cn-LOH",
                                             Mosaic_alteration=="GAIN" ~ "GAIN",
                                             Mosaic_alteration=="LOM_IC2" ~ "LOM_IC2",
                                             Histological.Diagnosis=="HCC" ~ "HCC",
                                             Histological.Diagnosis=="HB" ~ "HB",
                                             Histological.Diagnosis=="FLC" ~ "FLC",
                                             Histological.Diagnosis=="HCA" ~ "HCA",
                                             Histological.Diagnosis=="FF" ~ "FF"),
         Hist.diag.with.age = case_when(Hist.diag =="Mosaic" ~ "Mosaic",
                                        Hist.diag =="HCC" ~ "HCC",
                                        Hist.diag =="FLC" ~ "FLC",
                                        Hist.diag =="HCA" ~ "HCA",
                                        Hist.diag =="FF" ~ "FF",
                                        Hist.diag =="HB" & Age.surg.years<=3.3 ~ "HB_young",
                                        Hist.diag =="HB" & Age.surg.years>3.3 ~ "HB_old"),
         adult.prom=((p0_p1_ratio + p1_ratio)*100),
         fetal.prom=((p2_ratio + p3_ratio + p4_ratio)*100))
         

#Manual correction of 11p15 status from mosaic patients and series
annot.[which(annot.$CHCID=="CHC05303T"),"To_keep_for_survival_analysis_worst"] ="yes"
annot.[which(annot.$CHCID=="CHC3126T"), "To_keep_for_survival_analysis_worst"] = "yes"
annot.[which(annot.$CHCID=="CHC4087T"), "To_keep_for_survival_analysis_worst"] = "yes"
annot.[which(annot.$CHCID=="CHC2914T"), "status_11p15_ext"] = "LOM_IC2"
annot.[which(annot.$CHCID=="CHC3169T"), "status_11p15_ext"] = "GAIN"
annot.[which(annot.$CHCID=="CHC3168T"), "status_11p15_ext"] = "GAIN"

  # 1.6 Create annot.NT table --------------------------------------------------
annot.all <- annot. %>%
  filter(Papier_mosaiques_2021=="yes") %>%
  as.data.frame()  

#write.table(annot.all, "D:/Dropbox/11p15.5 mosaicism/annot.all.txt", sep="\t")

annot.NT = annot.all %>%
  filter(grepl("N|S|FF", CHCID)) %>%
  as.data.frame()  

dim(annot.NT)
names = annot.NT$CHCID
table(annot.NT$Mosaic)
table(annot.NT$Histological.Diagnosis)

#write.table(annot.NT, "D:/Dropbox/11p15.5 mosaicism/annot.NT.txt", sep="\t")

  # 1.7 Create annot.NT2 table with RNAseq data  -------------------------------

annot.NT2= annot.NT %>%
  filter(transcriptomic_mos=="yes") 

cols.NT = c("FF" = "#004040", "HB" = "#CDADFF", "HCC" = "darkred", "FLC" = "#673AB7", "HCA" = "#81C784")

#write.table(annot.NT2, "D:/Dropbox/11p15.5 mosaicism/annot.NT2.txt", sep="\t")

### 2. Age at surgery NT plot - Figure 1a --------------------------------------
  # 2.1 Plot ages with all diagnoses--------------------------------------------
cols.NT2 = c("FF"="#004040","HB"= "#CDADFF", "HCC"="#039BE5",  "FLC"= "#673AB7", "HCA"= "#81C784", "GAIN" = "#FFC107", "LOM_IC2" = "#EF9A9A","cn-LOH" = "#B71C1C")
table(annot.NT$Histological.Diagnosis)

annot.NT$Histological.Diagnosis = factor(annot.NT$Histological.Diagnosis, levels = c("FF","HCA", "FLC", "HCC", "HB"))
annot.NT$Age.surg.years = as.numeric(annot.NT$Age.surg.years)

table(annot.NT$Histological.Diagnosis)

age.plot <- ggplot(annot.NT, aes(x = Age.surg.years ,  y = Histological.Diagnosis, color = Histological.Diagnosis2, fill = Histological.Diagnosis2)) +
  geom_segment(aes(x =Age.surg.years , xend = Age.surg.years+0.05 , y = Histological.Diagnosis, yend=Histological.Diagnosis,  linewidth=2)) +
  scale_color_manual(values = cols.NT2) +
  theme_classic() +
  theme(legend.key=element_blank(),
        legend.position = "none",
        axis.title=element_blank())

age.plot

Source_figure_1a <- annot.NT %>% select(CHCID, Histological.Diagnosis, Histological.Diagnosis2, Age.surg.years)
table(Source_figure_1a$Histological.Diagnosis)
#write.table(Source_figure_1a, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Souce_Figure1a.txt", sep="\t")
#ggsave("/Users/jill/Desktop/age.plot.png", plot = age.plot, device = "png", width = 11, height = 6, units = "cm", dpi = 320)

  # 2.2 Plot zoom with only HB mos and non mos ---------------------------------
Median_HB = median(annot.NT[which(annot.NT$Hist.diag=="HB"),"Age.surg.years"])
Median_Mos = median(annot.NT[which(annot.NT$Hist.diag=="Mosaic"),"Age.surg.years"])

table(annot.NT$Hist.diag)

age.mos <- ggplot(subset(annot.NT, Hist.diag %in% c("Mosaic","HB")), aes(Age.surg.years ,  y = Hist.diag, color = Hist.diag, fill = Hist.diag)) +
  geom_density_ridges2(rel_min_height = .01, size = 0.5, 
                       bandwidth=1.5, alpha=0.8, scale=3) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(values = c("#CDADFF", "#B71C1C")) +
  scale_fill_manual(values = c("#CDADFF", "#B71C1C")) +
  coord_cartesian(clip = "off") +
  geom_segment(aes(x = Median_HB, y=1, xend=Median_HB, yend=2.8), linetype="dashed", 
               color = "#8f79b2", size=1) +
  geom_segment(aes(x = Median_Mos, y=2, xend=Median_Mos, yend=3.8), linetype="dashed", 
               color = "#B71C1C", size=1) +
  theme_ridges(center = TRUE) +
  theme(legend.key=element_blank(),
        legend.position = "none",
        axis.title=element_blank(),
        axis.text.y = element_blank())

#ggsave("/Users/jill/Desktop/age.mos.png", plot = age.mos, device = "png", width = 13, height = 6, units = "cm", dpi = 320)

test.data = annot.NT %>%
  filter(Hist.diag %in% c("HB", "Mosaic"))
t.test(test.data$Age.surg.years ~ test.data$Hist.diag, alternative="two.sided")

Source_figure_1abis <- annot.NT %>% subset(Hist.diag %in% c("Mosaic","HB")) %>% select("CHCID", "Hist.diag", "Age.surg.years")
table(Source_figure_1abis$Hist.diag)
#write.table(Source_figure_1abis, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_1abis.txt", sep="\t")

### 3. IC1/IC2 methylation values in non tumor samples - Figure 1b -------------

annot.methratio = annot.NT %>%
  filter(!is.na(MethRatio)) %>%
  arrange(MethRatio) 
annot.methratio$CHCID = factor(annot.methratio$CHCID, levels = annot.methratio$CHCID)

table(annot.methratio$Hist.diag)

cols = c("wt" = "Gray92", "GAIN" = "#FFC107", "LOM_IC2" = "#EF9A9A","cn-LOH" = "#B71C1C")

rat<- ggplot(annot.methratio, aes(CHCID, MethRatio, y=MethRatio, color=Mosaic_alteration)) + 
  geom_point(shape = 21,size=3, alpha=0.8, aes(color=Mosaic_alteration, fill=Mosaic_alteration))  +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits=c(0.8,3.5))+
  coord_cartesian(clip = 'off') + # To avoid clipped points 
  theme(legend.key = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", linewidth =0.3),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.ticks.length.x = unit(0.05, "cm"),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        plot.margin = margin(t = 40,  # Top margin
                             r = 40,  # Right margin
                             b = 40,  # Bottom margin
                             l = 40)) # Left margin)

rat

Source_figure_1b <- annot.methratio %>% select("CHCID", "Mosaic_alteration", "MethRatio")
table(Source_figure_1b$Mosaic_alteration)
#write.table(Source_figure_1b, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_1b.txt", sep="\t")

 #ggsave("/Users/jill/Desktop/rat.png", plot = rat, device = "png", width = 9, height = 10, units = "cm", dpi = 320)
 
### 4. IGF2 promoter usage  Supplementary Figure 7 -----------------------------

  # 4.1 Promoter value per diagnosis - Supp figure 7b---------------------------
cols = c("FF"="#004040","HB"= "#CDADFF", "HCC"="#039BE5",  "FLC"= "#673AB7", "HCA"= "#81C784", "GAIN" = "#FFC107", "LOM_IC2" = "#EF9A9A","cn-LOH" = "#B71C1C")

annot.prom.mos = annot.NT2 %>%
  arrange(Age.surg.years) %>%
  arrange(factor(Hist.diag, levels=c( "FF", "Mosaic", "HB", "HCC", "FLC", "HCA"))) %>%
  as.data.frame() 

annot.prom.mos$CHCID = factor(annot.prom.mos$CHCID , levels = annot.prom.mos$CHCID)

ggplot(annot.prom.mos, aes(x=CHCID, y=fetal.prom, color=Histological.Diagnosis2)) + 
  geom_point(shape = 21,size=3, alpha=0.8, aes(color=Histological.Diagnosis2, fill=Histological.Diagnosis2))  +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  theme_classic() + 
  theme(
        axis.text = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))

# adult promoters
ggplot(annot.prom.mos, aes(x=CHCID, y=adult.prom, color=Histological.Diagnosis2)) + 
  geom_point(shape = 21,size=3, alpha=0.8, aes(color=Histological.Diagnosis2, fill=Histological.Diagnosis2))  +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  theme_classic() + 
  theme(
    axis.text = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))


# Age at diagnosis
ggplot(annot.prom.mos, aes(x=CHCID, y=Age.surg.years, color=Histological.Diagnosis2)) + 
   geom_point(shape = 21,size=3, alpha=0.8, aes(color=Histological.Diagnosis2, fill=Histological.Diagnosis2))  +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() + 
  theme(
       axis.text = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, color = "#666666"))
 
# Plot IGF2 expression 
ggplot(annot.prom.mos, aes(x=CHCID, y=IGF2, color=Histological.Diagnosis2)) + 
  geom_point(shape = 21,size=3, alpha=0.8, aes(color=Histological.Diagnosis2, fill=Histological.Diagnosis2))  +
  scale_y_continuous(limits = c(10,20), breaks=c(10,15,20)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() + 
  theme(
    axis.text = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))

# Plot H19 expression 
ggplot(annot.prom.mos, aes(x=CHCID, y=H19, color=Histological.Diagnosis2)) + 
  geom_point(shape = 21,size=3, alpha=0.8, aes(color=Histological.Diagnosis2, fill=Histological.Diagnosis2))  +
  scale_y_continuous(limits = c(12,20), breaks =c(15,20)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() + 
  theme(
    axis.text = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))

# IC1/ IC2 

ggplot(annot.prom.mos, aes(x=CHCID, y=MethRatio, color=Histological.Diagnosis2)) + 
  geom_point(shape = 21,size=3, alpha=0.8, aes(color=Histological.Diagnosis2, fill=Histological.Diagnosis2))  +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(-1,4), breaks = c(0,2,4)) +
  theme_classic() + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))

Source_supp_figure_7b <- annot.prom.mos %>% select("CHCID",
                                                   "Age.surg.years", 
                                                   "Histological.Diagnosis2", 
                                                   "MethRatio",
                                                   "IGF2",
                                                   "H19",
                                                   "fetal.prom", 
                                                   "adult.prom")
table(Source_supp_figure_7b$Histological.Diagnosis2)
#write.table(Source_supp_figure_7b, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_7b.txt", sep="\t")

  # 4.2 Compare fetal promoter distribution - Supp Figure 7c -------------------
loess.distrib = annot.prom.mos %>%
  filter(Hist.diag %in% c("HB", "Mosaic"),
         Age.surg.years<=3.3)  %>%
  arrange(Age.surg.years)

dim(loess.distrib)
cols = c("HB"= "#CDADFF", "Mosaic" = "indianred")
table(loess.distrib$Hist.diag)

# Loess regression fetal promoters

ggplot(loess.distrib, aes(x=Age.surg.years, y=fetal.prom, color=Hist.diag)) + 
  geom_smooth(method = "loess", aes(fill = Hist.diag, color = Hist.diag), size = 3, se = T, level = 0.65, alpha = 0.2) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(-30,130), breaks = c(0,50, 100)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))


# Wilcoxon test fetal promoters 

ggplot(loess.distrib, aes(x=Hist.diag, y=fetal.prom, color=Hist.diag)) + 
  geom_violin(aes(fill = Hist.diag, color = Hist.diag), alpha = 0.7, trim = F) +
  stat_summary(fun.data=mean_sdl, size = 1,
               geom="pointrange", color="white") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
#  scale_y_continuous(limits = c(-0.2,0.5), breaks = c(0,0.2, 0.4)) +
  theme_classic() +
  theme(
       axis.text = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))

wilcox.test(loess.distrib$fetal.prom ~ loess.distrib$Hist.diag)

Source_supp_figure_7c <- loess.distrib %>% select("CHCID",
                                                   "Hist.diag", 
                                                   "fetal.prom")
table(Source_supp_figure_7c$Hist.diag)
#write.table(Source_supp_figure_7c, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_7c.txt", sep="\t")

  # 4.3 Methylation promoters - Supp Figure 7a ---------------------------------

`%nin%` = Negate(`%in%`)
cols.meth = c("FF" = "#004040", "HB" = "#CDADFF", "Mosaic"="indianred")

meth.data = annot.NT %>% 
  select(c("CHCID", "Hist.diag", "Histological.Diagnosis", "Mosaic.per.sample", "mean_P1", "mean_P2", "mean_P3", "mean_P4", "Age.surg.years", "Type")) %>%
  filter(!is.na(mean_P2))  %>% #Filter NA for P2 because 2 NA for Promoter P1
  gather(., key ="promoter", value="value", c("mean_P1", "mean_P2", "mean_P3", "mean_P4")) %>%
  arrange(Age.surg.years) %>%
  filter(Type %in% c("FF", "N")) %>%
  arrange(factor(Histological.Diagnosis, levels = c("FF", "HB")))  %>%  # Don't filter to generate source file but filter to generate graph
 filter(promoter=="mean_P4") 


Source_supp_figure_7a <- meth.data %>% select("CHCID", "Hist.diag", "promoter", "Age.surg.years")
table(Source_supp_figure_7a$Hist.diag)
#write.table(Source_supp_figure_7a, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_7a.txt", sep="\t")


meth.data$CHCID = factor(meth.data$CHCID, levels = meth.data$CHCID)

# Plot promoter methylation

ggplot(meth.data, aes(x=CHCID, y=value,fill = Hist.diag)) + 
  geom_bar(stat="identity", color = "black")+
  scale_fill_manual(values = cols.meth) +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =1),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank(),
        #    axis.text.y=element_blank(),
        axis.text.x=element_blank())

# Plot age 

ggplot(meth.data, aes(x=CHCID, Age.surg.years), y=Age.surg.years, color=Hist.diag) + 
  geom_point(shape = 21,size=10, alpha=0.8, aes(color=Hist.diag, fill=Hist.diag))  +
  scale_color_manual(values = cols.meth) +
  scale_fill_manual(values = cols.meth) +
  scale_y_continuous(limits=c(-5,20))+
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =2),
        axis.ticks.length = unit(0.3, "cm"),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x = element_blank() )

### 5. Compare BWS and mosaic frequencies - Figure 1c --------------------------

annot.BWS = annot.NT %>%
  filter(Histological.Diagnosis=="HB")  # Restrict to HB

  # 5.1 Calculate frequencies --------------------------------------------------

loc11p.table = as.data.frame(table(annot.BWS$Mosaic, annot.BWS$Mosaic_alteration))
colnames(loc11p.table) = c("Diag", "loc11p", "nb")

loc11p.table = loc11p.table %>%
  filter(Diag=="Mosaic",
         loc11p!="wt")

nb_Mosaic = as.numeric(loc11p.table %>% 
                         filter (Diag == "Mosaic") %>%
                         summarise(nb = sum(nb)))

nb_Mosaic
  
  # 5.2 Calculate frequencies and add BWS Brioude data -------------------------

loc11p.table = loc11p.table %>%
  mutate(freq = case_when(Diag=="Mosaic" ~ nb/nb_Mosaic*100)) %>%
  add_row(Diag = "BWS_Brioude_2018", loc11p="LOM_IC2", freq = 50) %>%
  add_row(Diag = "BWS_Brioude_2018", loc11p="GOM_IC1", freq = 5) %>%
  add_row(Diag = "BWS_Brioude_2018", loc11p="cn-LOH", freq = 20) %>%
  add_row(Diag = "BWS_Brioude_2018", loc11p="DEL-LOH", freq = 2) %>%
  add_row(Diag = "BWS_Brioude_2018", loc11p="CDKN1C_mut", freq = 5) %>%
  add_row(Diag = "BWS_Brioude_2018", loc11p="GAIN", freq = 4) %>%
  add_row(Diag = "BWS_Brioude_2018", loc11p="wt", freq = 14)


  # 5.3 Make graph frequencies of alterations - Figure 1c ----------------------
new_ord <-c("wt", "GAIN",
            "CDKN1C_mut", 
            "LOM_IC2", "GOM_IC1", "DEL-LOH", "cn-LOH")
cols <-c("cn-LOH" = "#B71C1C", "GOM_IC1" = "RoyalBlue", "LOM_IC2" = "#EF9A9A", 
         "CDKN1C_mut" = "lightblue", 
         "wt" = "white", "GAIN" = "#FFC107", "DEL-LOH" = "Blue")
loc11p.table$loc11p = factor(loc11p.table$loc11p, levels = new_ord)

loc11p.table$Diag = factor(loc11p.table$Diag, levels = c("BWS_Brioude_2018","Mosaic"))

ggplot(loc11p.table, aes(x=Diag, y = freq)) +
 # coord_flip() +
  geom_col(aes(fill = loc11p), position = position_stack(), width = 0.9, colour="black", size=1)  +
  scale_fill_manual(values = cols)  +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() + 
  theme(
       axis.text = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, color = "#666666"))

Mos_var = c(0,1,0,2,0,0,9) # wt, GAIN, CDKN1C, LOM IC2, GOM IC1, DEL-LOH, CN-LOH
res <- chisq.test(Mos_var, p = c(0.14, 0.04, 0.05, 0.5, 0.05, 0.02, 0.2))
res

Source_figure_1c <- loc11p.table %>% select("Diag", "loc11p", "freq")
head(Source_figure_1c)
#write.table(Source_figure_1c, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_1c.txt", sep="\t")


### 6. Tumors analysis - Figure 6, supp figure 18 ------------------------------
  # 6.1 Create table with tumors -----------------------------------------------
annot.asso = annot. %>% 
  filter(Papier_mosaiques_2021=="yes") %>%
  mutate(Hist.diag = case_when(Mosaic == "Mosaic" ~ "Mosaic",
                               Histological.Diagnosis == "HB" ~ "HB",
                               Histological.Diagnosis == "HCC" ~ "HCC",
                               Histological.Diagnosis == "FLC" ~ "FLC",
                               Histological.Diagnosis == "HCA" ~ "HCA",
                               Histological.Diagnosis == "HB.META" ~ "HB",
                               Histological.Diagnosis == "HB.RELAPSE" ~ "HB",
                               Histological.Diagnosis == "HCA.RELAPSE" ~ "HCA",
                               Histological.Diagnosis == "FF" ~ "FF")) %>%
  as.data.frame()  
 


dim(annot.asso)
table(annot.asso$Histological.Diagnosis)
table(annot.asso$Type)
names = annot.asso$CHCID

cols.T = c("HB" = "#CDADFF", "HB.MR" = "#FF1B91","HCC" = "#039BE5", "FLC" = "#673AB7", "HCA" = "#81C784", "Mosaic" = "indianred")

  # 6.2 Calculation of bcat and 11p frequency-----------------------------------

annot.all.T = annot.asso %>%
  filter(grepl("T", CHCID))

annot.all.T$Histological.Diagnosis

Recap_by_patient <- function(df) {
  df %>%
    count(Diag, alteration) %>%
    group_by(Diag) %>%
    summarise(nmut = sum(n[alteration]),
              ntot = sum(n),
              prop_percent = 100*sum(n[alteration])/sum(n))
}


loc11p.calc = annot.all.T %>%
  mutate(Histological.Diagnosis = gsub("\\.RELAPSE", "", gsub("\\.META", "", Histological.Diagnosis))) %>%
  filter(!is.na(status_11p15_ext)) %>%
  group_by(Patient.identification) %>%
  summarise(Diag = dplyr::first(Histological.Diagnosis),
            alteration = any(status_11p15_ext != "wt")) %>%
  Recap_by_patient()

bcat.calc = annot.all.T %>%
  mutate(Histological.Diagnosis = gsub("\\.RELAPSE", "", gsub("\\.META", "", Histological.Diagnosis))) %>%
  filter(!is.na(CTNNB1_simple)) %>%
  group_by(Patient.identification) %>%
  summarise(Diag = dplyr::first(Histological.Diagnosis),
            alteration = any(CTNNB1_simple != "wt")) %>%
  Recap_by_patient()

apc.calc = annot.all.T %>%
  mutate(Histological.Diagnosis = gsub("\\.RELAPSE", "", gsub("\\.META", "", Histological.Diagnosis))) %>%
  filter(!is.na(APC.x)) %>%
  group_by(Patient.identification) %>%
  summarise(Diag = dplyr::first(Histological.Diagnosis),
            alteration = any(APC.x != "wt")) %>%
  Recap_by_patient()


axin1.calc = annot.all.T %>%
  mutate(Histological.Diagnosis = gsub("\\.RELAPSE", "", gsub("\\.META", "", Histological.Diagnosis))) %>%
  filter(!is.na(AXIN1)) %>%
  group_by(Patient.identification) %>%
  summarise(Diag = dplyr::first(Histological.Diagnosis),
            alteration = any(AXIN1 != "wt")) %>%
  Recap_by_patient()


# Count total series of HB 
count.HB = annot.all.T  %>% mutate(Histological.Diagnosis = gsub("\\.RELAPSE", "", gsub("\\.META", "", Histological.Diagnosis))) %>% filter(Histological.Diagnosis=="HB")
length(unique(count.HB$Patient.identification))
dim(count.HB)

# Add 11p15.5 annotation per patient 
annot.all.T = annot.all.T %>%
  group_by(Patient.identification) %>% 
  mutate(status_11p_by_patient = case_when(any(status_11p15_ext=="cn-LOH") ~ "cn-LOH",
                                           any(status_11p15_ext=="DEL-LOH") ~ "DEL-LOH",
                                           any(status_11p15_ext=="CDKN1C_mut") ~ "CDKN1C_mut",
                                           any(status_11p15_ext=="GOM_IC1") ~ "GOM_IC1",
                                           any(status_11p15_ext=="LOM_IC2") ~ "LOM_IC2",
                                           any(status_11p15_ext=="GAIN") ~ "GAIN",
                                           any(status_11p15_ext=="wt") ~ "wt",
                                           is.na(status_11p15_ext) ~ NA_character_)) %>%
  ungroup()

  # 6.3 Plot tumor promoter proportions by diag - Figure 6c --------------------

annot.all.T = annot.all.T %>%
  mutate(Histological.Diagnosis = case_when(Mosaic == "Mosaic" ~ "Mosaic",
                                            Histological.Diagnosis == "HB" ~ "HB",
                                            Histological.Diagnosis == "HCC" ~ "HCC",
                                            Histological.Diagnosis == "FLC" ~ "FLC",
                                            Histological.Diagnosis == "HCA" ~ "HCA",
                                            Histological.Diagnosis == "HB.META" ~ "HB.MR",
                                            Histological.Diagnosis == "HB.RELAPSE" ~ "HB.MR",
                                            Histological.Diagnosis == "HCA.RELAPSE" ~ "HCA")) 

#write.table(annot.all.T, "D:/Dropbox/11p15.5 mosaicism/annot.all.T.txt", sep="\t")

# Count diag with prom data
count. = annot.all.T  %>% filter(!is.na(fetal.prom))
table(count.$Histological.Diagnosis)


Annot.prop= annot.all.T %>%
  filter(!is.na(fetal.prom)) %>%
  gather(., key ="Promoter", value="Proportion",  c("fetal.prom", "adult.prom")) %>% 
  select(c("CHCID", "Promoter", "Proportion", "Histological.Diagnosis")) %>%
  group_by(Histological.Diagnosis, Promoter) %>%
  summarise_at(vars(Proportion), list(Proportion = mean)) %>%
  as.data.frame()

Source_figure_6cbis <- annot.all.T %>%
  filter(!is.na(fetal.prom)) %>%
  gather(., key ="Promoter", value="Proportion",  c("fetal.prom", "adult.prom")) %>% 
  select(c("CHCID", "Promoter", "Proportion", "Histological.Diagnosis")) 
table(Source_figure_6cbis$Histological.Diagnosis)
#write.table(Source_figure_6cbis, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_6cbis.txt", sep="\t")

Annot.prop$Promoter = factor(Annot.prop$Promoter, levels = c("adult.prom", "fetal.prom"))
Annot.prop$Histological.Diagnosis = factor(Annot.prop$Histological.Diagnosis, levels = c("Mosaic", "HB", "HB.MR",  "HCC", "FLC", "HCA"))
prom <- ggplot(Annot.prop, aes(fill=Promoter, y=Proportion, x=Histological.Diagnosis)) + 
  geom_bar(position="stack", stat="identity", color = "black", width=0.8, size = 1.5)+
  scale_fill_manual(values = c("#FFEB3B", "#0319F4")) +
  theme_classic() +
  theme(legend.key = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =1.5),
        axis.text = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(size=1.5),
        panel.background = element_blank())

prom

#ggsave("/Users/jill/Desktop/promoters.png", plot = prom, device = "png", width = 15, height = 15, units = "cm", dpi = 320)

Source_figure_6c <- Annot.prop 
head(Source_figure_6c)
#write.table(Source_figure_6c, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_6c.txt", sep="\t")

  # 6.4 Area plot B-cat figure 6a ----------------------------------------------

annot.T.HB = annot.all.T %>%
  filter(Hist.diag %in% c("HB", "Mosaic")) %>%
  mutate(Age_simple = ceiling(Age.surg.years)) %>%
  filter(To_keep_for_survival_analysis_worst %in% c("yes", "yes?"))

dim(annot.T.HB)

#write.table(annot.T.HB, "D:/Dropbox/11p15.5 mosaicism/annot.T.txt", sep="\t")

viridis(n = 4)

cols = c("wt" = "#FDE725FF", "Missense" = "#35B779FF", "Small_deletion" = "#31688EFF", "Large_deletion" = "#440154FF")


    # 6.4.1 all HB figure 6a ---------------------------------------------------
annot.T.HB.nona = annot.T.HB %>%
  filter(!is.na(CTNNB1_simple))

table(annot.T.HB.nona$CTNNB1_simple)
dim(annot.T.HB.nona)

Source_figure_6a <- annot.T.HB.nona %>% select(c("CHCID", "Histological.Diagnosis", "Age_simple", "CTNNB1_simple")) 
head(Source_figure_6a)
#write.table(Source_figure_6a, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_6a.txt", sep="\t")


test.nonmos = as.data.frame(table(annot.T.HB.nona$Age_simple, annot.T.HB.nona$CTNNB1_simple))
colnames(test.nonmos) = c("Age_simple", "CTNNB1", "Freq")
test.nonmos$Freq = as.numeric(test.nonmos$Freq)
test.nonmos$Age_simple = as.numeric(test.nonmos$Age_simple)
    
    # 6.4.2 mosaic HB figure 6a ------------------------------------------------

annot.T2.mosaic = annot.T.HB.nona %>%
  filter(Histological.Diagnosis=="Mosaic")
test.mos = as.data.frame(table(annot.T2.mosaic$Age_simple, annot.T2.mosaic$CTNNB1_simple))

colnames(test.mos) = c("Age_simple", "CTNNB1", "Freq")
test.mos$Freq = as.numeric(test.mos$Freq)
test.mos$Age_simple = as.numeric(test.mos$Age_simple)

    # 6.4.3 Plot all figure 6a -------------------------------------------------

plot.bcat <- ggplot() + 
  geom_area(data = test.nonmos, alpha=1 , size=1,  aes(x=as.numeric(Age_simple), y=as.numeric(Freq), fill=CTNNB1, color=CTNNB1)) +
  scale_x_continuous(limits = c(0,18), breaks = c(0, 5, 10, 15))+
  scale_y_continuous(limits = c(0,25), breaks = c(0, 12.5, 25))+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_area(data = test.mos,  alpha=1 , size=1,  aes(x=as.numeric(Age_simple), y=as.numeric(Freq), fill=CTNNB1), color ="black",  linetype = "dashed") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =2),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size=2))

#ggsave("/Users/jill/Desktop/plot.bcat.png", plot = plot.bcat, device = "png", width = 22, height = 15, units = "cm", dpi = 320)
plot.bcat

  # 6.5 Area plot 11p-----------------------------------------------------------
    # 6.5.1 All HB  ------------------------------------------------------------

annot.T.HB.nona = annot.T.HB %>%
  filter(!is.na(status_11p_by_patient))

test.nonmos = as.data.frame(table(annot.T.HB.nona$Age_simple, annot.T.HB.nona$status_11p_by_patient))
colnames(test.nonmos) = c("Age_simple", "status_11p_by_patient", "Freq")
test.nonmos$status_11p_by_patient = factor(test.nonmos$status_11p_by_patient, levels = c("cn-LOH", "DEL-LOH", "GOM_IC1", "LOM_IC2", "CDKN1C_mut", "GAIN","wt"))

Source_figure_6abis <- annot.T.HB.nona %>% select(c("CHCID", "Histological.Diagnosis", "Age_simple", "status_11p_by_patient")) 
head(Source_figure_6abis)
#write.table(Source_figure_6abis, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_6abis.txt", sep="\t")


    # 6.5.2 Mosaic HB ----------------------------------------------------------

annot.T2.mosaic = annot.T.HB.nona %>%
  filter(Histological.Diagnosis=="Mosaic")
test.mos = as.data.frame(table(annot.T2.mosaic$Age_simple, annot.T2.mosaic$status_11p_by_patient))

colnames(test.mos) = c("Age_simple", "status_11p_by_patient", "Freq")

    # 6.5.3 Plot all------------------------------------------------------------

table(annot.T.HB.nona$status_11p_by_patient)

dim(annot.T.HB.nona)

cols.11p = c("CDKN1C_mut" = "lightblue", "cn-LOH" = "#B71C1C", "DEL-LOH" = "Blue", "GOM_IC1" = "Royalblue", "LOM_IC2" = "#EF9A9A", "GAIN" = "orange", "wt" = "white")

plot.11p <- ggplot() + 
  geom_area(data = test.nonmos, alpha=1 , size=.5,  aes(x=as.numeric(Age_simple), y=as.numeric(Freq), fill=status_11p_by_patient, color=status_11p_by_patient)) +
  scale_x_continuous(limits = c(0,18), breaks = c(0,5,10,15))+
  scale_y_continuous(limits = c(0,25), breaks = c(0, 12.5, 25))+
  scale_color_manual(values = cols.11p) +
  scale_fill_manual(values = cols.11p) +
  geom_area(data = test.mos,  alpha=1 , size=1,  aes(x=as.numeric(Age_simple), y=as.numeric(Freq), fill=status_11p_by_patient),  color ="black",linetype = "dashed") +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =2),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(size=2))


#ggsave("/Users/jill/Desktop/plot.11p.png", plot = plot.11p, device = "png", width = 22, height = 15, units = "cm", dpi = 320)

  # 6.6 Barplot Bcat------------------------------------------------------------

annot.. = annot.T.HB %>%
  mutate(Mosaic = case_when(Mosaic=="Mosaic" ~  "Mosaic",
                             Age.surg.years<=3.3 ~ "young",
                             Age.surg.years>3.3 ~ "old")) %>%
  filter(!is.na(CTNNB1_simple)) %>%
  group_by(Mosaic, CTNNB1_simple) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.. %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot..$CTNNB1_simple = factor(annot..$CTNNB1_simple, levels  =c("Large_deletion","Small_deletion","Missense", "wt"))
annot..$Mosaic = factor(annot..$Mosaic, levels  =c("Mosaic", "young", "old"))

annot..$count = as.numeric(annot..$count)

cols = c("wt" = "#FDE725FF", "Missense" = "#35B779FF", "Small_deletion" = "#31688EFF", "Large_deletion" = "#440154FF")

barplot.bcat <- ggplot(annot.., aes(x=Mosaic, y=perc,fill = CTNNB1_simple)) + 
  geom_bar(stat="identity", col = "black", size=2)+
  scale_fill_manual(values = cols) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =2),
        axis.ticks = element_line(size=2),
        axis.ticks.length = unit(0.3, "cm"),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x = element_blank())

# ggsave("/Users/jill/Desktop/barplot.bcat.png", plot = barplot.bcat, device = "png", width = 13, height = 20, units = "cm", dpi = 320)

annot.stat = annot..%>%
  pivot_wider(names_from = Mosaic,  
              values_from=count) %>%
  select(-perc) %>%
  group_by(CTNNB1_simple) %>%
  summarise(Mosaic = sum(Mosaic, na.rm = T),
            old = sum(old, na.rm = T),
            young = sum(young, na.rm = T))

names = annot.stat$CTNNB1_simple

annot.stat = annot.stat %>%
  select(c(Mosaic, old)) # choose groups to compare

rownames(annot.stat) = names

chisq.test(annot.stat)

  # 6.7 Barplot 11p15-----------------------------------------------------------

annot.. = annot.T.HB %>%
  mutate(Mosaic = case_when(Mosaic=="Mosaic" ~  "Mosaic",
                            Age.surg.years<=3.3 ~ "young",
                            Age.surg.years>3.3 ~ "old")) %>%
  filter(!is.na(status_11p_by_patient)) %>%
  group_by(Mosaic, status_11p_by_patient) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.. %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot..$status_11p_by_patient = factor(annot..$status_11p_by_patient, levels  =c("wt", "GAIN", "CDKN1C_mut", "LOM_IC2", "GOM_IC1", "DEL-LOH", "cn-LOH"))
annot..$Mosaic = factor(annot..$Mosaic, levels  =c("Mosaic", "young", "old"))

cols.11p = c("CDKN1C_mut" = "lightblue", "cn-LOH" = "#B71C1C", "DEL-LOH" = "Blue", "GOM_IC1" = "Royalblue", "LOM_IC2" = "#EF9A9A", "GAIN" = "orange", "wt" = "white")

barplot.11p <- ggplot(annot.., aes(x=Mosaic, y=perc,fill = status_11p_by_patient)) + 
  geom_bar(stat="identity", col = "black", size=2)+
  scale_fill_manual(values = cols.11p) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_line(colour = "black", size =2),
        axis.ticks = element_line(size=2),
        axis.ticks.length = unit(0.3, "cm"),
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x = element_blank())

#ggsave("/Users/jill/Desktop/barplot.11p.png", plot = barplot.11p, device = "png", width = 13, height = 20, units = "cm", dpi = 320)


annot.stat = annot..%>%
  pivot_wider(names_from = Mosaic,  
              values_from=count) %>%
  select(-perc) %>%
  group_by(status_11p_by_patient) %>%
  summarise(Mosaic = sum(Mosaic, na.rm = T),
            old = sum(old, na.rm = T),
            young = sum(young, na.rm = T))

names = annot.stat$status_11p_by_patient

annot.stat = annot.stat %>%
 filter(status_11p_by_patient !="CDKN1C_mut")  %>% # for comparision between Mosaic and old remove CDKN1C 0 for both old and mosaic
# filter(status_11p_by_patient != "LOM_IC2" & status_11p_by_patient!="GAIN") %>% #for comparision between young and old remove LOM IC2 and GAIN both 0
  select(c(Mosaic, old))  # choose groups to compare

rownames(annot.stat) = names

chisq.test(annot.stat)

  # 6.8 Create oncoprint -------------------------------------------------------
    # 6.8.1 Create function ----------------------------------------------------

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

    # 6.8.2 Create annotation for 11p timing -----------------------------------

timing.11p <- read_xlsx("D:/Dropbox/11p15.5 mosaicism/Timing_11p/Timing_11p_from_SuppFig5_Cancer_Discovery.xlsx") %>% as.data.frame()

timing.11p$Patient.identification = as.character(timing.11p$Patient.identification)

annot.T.HB$Patient.identification = as.character(annot.T.HB$Patient.identification)
annot.T.HB = left_join(annot.T.HB, timing.11p, by="Patient.identification")

annot.cooc = annot.T.HB %>%
  mutate(CTNNB1_simple=ifelse(is.na(CTNNB1_simple), "NA_gray", CTNNB1_simple), # Change NA in NA_gray to have NA in gray
         status_11p_by_patient=ifelse(is.na(status_11p_by_patient), "NA_gray", status_11p_by_patient),
         timing_11p = case_when(Mosaic=="Mosaic" ~ "Mosaic",
                                timing_11p=="trunk" ~  "NA_gray",
                                timing_11p=="wt" ~  "wt",
                                timing_11p=="private" ~  "private",
                                status_11p_by_patient =="wt" ~ "wt",
                                is.na(timing_11p) ~ "NA_gray")) %>%
  arrange(status_11p_by_patient) %>%
  arrange(factor(CTNNB1_simple, levels=c("wt", "Missense", "Small_deletion", "Large_deletion")))%>%
  arrange(Age.surg.years) %>%
  arrange(Mosaic) %>%
  as.data.frame()

dim(annot.cooc)
annot.cooc$timing_11p

annot.cooc[which(annot.cooc$CHCID=="CHC3082T"),"timing_11p"] ="private"  # Add 2 private because 11p alteration in one but not in another sample but no tree
annot.cooc[which(annot.cooc$CHCID=="CHC3575T"),"timing_11p"] ="private"  # Add 2 private because 11p alteration in one but not in another sample but no tree

Multiple.samples <- annot.cooc %>% select(c("CHCID", "Patient.identification",  "timing_11p", "Timing_alteration", "status_11p_by_patient", "Age.surg.years"))

Source_figure_6b <- Multiple.samples
head(Source_figure_6b)
#write.table(Source_figure_6b, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_6b.txt", sep="\t")

write.table(Multiple.samples, "D:/Dropbox/11p15.5 mosaicism/Timing_11p/multiple.samples.txt", sep="\t")
                        
    # 6.8.3 Heatmap one sample per patient--------------------------------------

annotCol = c("Age.surg.years", "timing_11p",  "status_11p_by_patient","CTNNB1_simple")

anocol = color_annot_generated(annot.cooc[, annotCol], col_interet = annotCol, plotLegend = T, computer = "PC",
                               color.matrix = "D:/Dropbox/11p15.5 mosaicism not shared/Color_matrix_JP27012020.txt")

geco.imageCol(anocol, xlab.cex = 0, ylab.cex = 0,drawLines = c("h"))

  # 6.9 Co-occurence between bcat deletion  and GOM IC1 ------------------------

# one sample per patient all del
annot.cooc = annot.cooc %>%
  mutate(CTNNB1_del_yes_no = ifelse(CTNNB1_simple %in% c("Large_deletion", "Small_deletion"), "yes", ifelse(status_11p15_ext=="NA_gray", "NA_gray", "no")),
         Epi_IC1_yes_no = ifelse(status_11p15_ext=="GOM_IC1", "yes", ifelse(status_11p15_ext=="NA_gray", "NA_gray", "no")))

test = annot.cooc %>%
  filter(Epi_IC1_yes_no %in% c("yes", "no"),
         CTNNB1_del_yes_no %in% c("yes", "no")) %>%
  group_by(Epi_IC1_yes_no, CTNNB1_del_yes_no) %>% 
  summarise(count = n())


test = test %>%
  pivot_wider(names_from = CTNNB1_del_yes_no, values_from = count, names_prefix = "CTNNB1_del_")

names = test$Epi_IC1_yes_no
test = test[,2:3]
rownames(test) = names
fisher.test(test)

  # 6.10.IGF2 expression in different groups - supp figure 18 ------------------

annot.all.T = annot.all.T %>%
  mutate(Diag_with_11p = case_when(Histological.Diagnosis=="Mosaic" ~ "Mosaic",
                                   Histological.Diagnosis=="HB" & status_11p15_ext!="wt" ~ "HBalt",
                                   Histological.Diagnosis=="HB" & status_11p15_ext=="wt" ~ "HBwt",
                                   Histological.Diagnosis=="HB.MR" & status_11p15_ext!="wt" ~ "HB.MRalt",
                                   Histological.Diagnosis=="HB.MR" & status_11p15_ext=="wt" ~ "HB.MRwt",
                                   Histological.Diagnosis=="HCC" & status_11p15_ext!="wt" ~ "HCCalt",
                                   Histological.Diagnosis=="HCC" & status_11p15_ext=="wt" ~ "HCCwt",
                                   Histological.Diagnosis=="FLC" & status_11p15_ext!="wt" ~ "FLCalt",
                                   Histological.Diagnosis=="FLC" & status_11p15_ext=="wt" ~ "FLCwt",
                                   Histological.Diagnosis=="HCA" & status_11p15_ext!="wt" ~ "HCAalt",
                                   Histological.Diagnosis=="HCA" & status_11p15_ext=="wt" ~ "HCAwt"))

cols.all = c("HBalt" = "#CDADFF", 
             "HBwt" = "#CDADFF", 
             "HB.MRalt" = "#FF1B91",
             "HB.MRwt" = "#FF1B91",
             "HCCalt" = "#039BE5", 
             "HCCwt" = "#039BE5", 
             "FLCalt" = "#673AB7", 
             "FLCwt" = "#673AB7", 
             "HCAalt" = "#81C784", 
             "HCAwt" = "#81C784", 
             "Mosaic" = "indianred")

annot.all.T$Diag_with_11p = factor(annot.all.T$Diag_with_11p, levels = c("Mosaic","HBwt","HBalt","HB.MRwt","HB.MRalt","HCCwt","HCCalt","FLCwt","FLCalt","HCAwt","HCAalt"))

annot.all.T.nona = annot.all.T %>%
  filter(!is.na(Diag_with_11p),
         !is.na(IGF2))
table(annot.all.T.nona$Diag_with_11p)

Source_supp_figure_18b <- annot.all.T.nona %>% select(c("CHCID", "Diag_with_11p", "IGF2", "H19", "CDKN1C")) 
table(Source_supp_figure_18b$Diag_with_11p)
#write.table(Source_supp_figure_18b, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_18b.txt", sep="\t")

# IGF2 expression 

IGF2 <- ggplot(annot.all.T.nona, aes(x = Diag_with_11p, y=IGF2, fill=Diag_with_11p)) + 
  geom_violin(width=1, alpha=.8, size=0.7,trim = F) +
  geom_boxplot(width=.3, cex=.7) +
  scale_color_manual(values = cols.all) +
  scale_y_continuous(limits=c(5,23))+
  scale_fill_manual(values = cols.all) +
  theme_bw() +
  theme(legend.position="none",
    plot.title = element_text(size=11),
    axis.text.y = element_text(size=20),
    axis.text.x = element_blank(),
    axis.line = element_line(colour = "black", size =1),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) 

IGF2

#ggsave("/Users/jill/Desktop/IGF2.png", plot = IGF2, device = "png", width = 20, height =7, units = "cm", dpi = 320)


annot.all.test = annot.all.T.nona %>%
  filter(Diag_with_11p %in% c("HCAwt", "HCAalt"))
wilcox.test(annot.all.test$IGF2 ~ annot.all.test$Diag_with_11p)

# H19 expression 

 H19 <-ggplot(annot.all.T.nona, aes(x = Diag_with_11p, y=H19,fill=Diag_with_11p)) + 
  geom_violin(width=1, alpha=.8, size=0.7,trim = F) +
  geom_boxplot(width=.3, cex=.7) +
  scale_color_manual(values = cols.all) +
  scale_y_continuous(limits=c(5,25) )+
  scale_fill_manual(values = cols.all) +
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(size=11),
        axis.text.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black", size =1),
        axis.ticks.length = unit(0.3, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
 
H19

#ggsave("/Users/jill/Desktop/H19.png", plot = H19, device = "png", width = 20, height = 7, units = "cm", dpi = 320)

annot.all.test = annot.all.T.nona %>%
  filter(Diag_with_11p %in% c("HCAwt", "HCAalt"))
wilcox.test(annot.all.test$H19 ~ annot.all.test$Diag_with_11p)

# CDKN1C expression 

CDKN1C <- ggplot(annot.all.T.nona, aes(x = Diag_with_11p, y=CDKN1C,fill=Diag_with_11p)) + 
  geom_violin(alpha=.8, size=0.7, trim = F) +
  geom_boxplot(width=.2, cex=.7) +
  scale_color_manual(values = cols.all) +
  scale_y_continuous(limits=c(0,15) )+
  scale_fill_manual(values = cols.all) +
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(size=11),
        axis.text.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black", size =1),
        axis.ticks.length = unit(0.3, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
CDKN1C

#ggsave("/Users/jill/Desktop/CDKN1C.png", plot = CDKN1C, device = "png", width = 20, height = 7, units = "cm", dpi = 320)

annot.all.test = annot.all.T.nona %>%
  filter(Diag_with_11p %in% c("HCAwt", "HCAalt"))
wilcox.test(annot.all.test$CDKN1C ~ annot.all.test$Diag_with_11p)

  # 6.11 Boxplot RRBS methylation values ---------------------------------------

annot.methyl = annot.asso %>%
  filter(!is.na(fetal_mean_meth)) %>% # remove NAs
  mutate(Histological.Diagnosis = case_when(Mosaic == "Mosaic" ~ "Mosaic",
                                            Histological.Diagnosis == "FF" ~ "FF",
                                            Histological.Diagnosis == "HB" ~ "HB",
                                            Histological.Diagnosis == "HCC" ~ "HCC",
                                            Histological.Diagnosis == "FLC" ~ "FLC",
                                            Histological.Diagnosis == "HCA" ~ "HCA",
                                            Histological.Diagnosis == "HB.META" ~ "HB.MR",
                                            Histological.Diagnosis == "HB.RELAPSE" ~ "HB.MR",
                                            Histological.Diagnosis == "HCA.RELAPSE" ~ "HCA"),
          Diag_for_meth = case_when(grepl("N", CHCID) & Histological.Diagnosis=="Mosaic" ~ "NMosaic",
                                   grepl("T", CHCID) & Histological.Diagnosis=="Mosaic" ~ "TMosaic", 
                                   grepl("N", CHCID) & Histological.Diagnosis!="Mosaic" ~ "NTL",
                                   Histological.Diagnosis =="FF" ~ "FF",
                                   grepl("T", CHCID) & Histological.Diagnosis=="HB" ~ "HB", 
                                   grepl("T", CHCID) & Histological.Diagnosis=="FLC" ~ "FLC", 
                                   grepl("T", CHCID) & Histological.Diagnosis=="HCC" ~ "HCC", 
                                   grepl("T", CHCID) & Histological.Diagnosis=="HB.MR" ~ "HB.MR"))

annot.methyl$Diag_for_meth =  factor(annot.methyl$Diag_for_meth, levels =c("FF","NTL","NMosaic", "TMosaic", "HB", "HB.MR", "HCC","FLC"))

Source_supp_figure_18a <- annot.methyl %>% select(c("CHCID", "Diag_for_meth", "fetal_mean_meth", "mean_P1")) 
table(Source_supp_figure_18a$Diag_for_meth)
#write.table(Source_supp_figure_18a, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_18a.txt", sep="\t")

cols = c("FF" = "#004040", "NTL"="#B2DFDB", "NMosaic"="indianred","TMosaic"="darkred","HB"= "#CDADFF",  "HB.MR" = "#CDADFF", "HCC"="#039BE5",  "FLC"= "#673AB7")

# Fetal promoters

fetal.methyl <- ggplot(annot.methyl, aes(x = Diag_for_meth, y=fetal_mean_meth ,fill=Diag_for_meth)) + 
  geom_violin(alpha=.8, size=0.7, trim = F,width=2) +
  geom_boxplot(width=.2, cex=.7) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(size=11),
        axis.text.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black", size =1),
        axis.ticks.length = unit(0.3, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

fetal.methyl
#ggsave("/Users/jill/Desktop/fetal.methyl.png", plot = fetal.methyl, device = "png", width = 20, height = 10, units = "cm", dpi = 320)

stat.test = annot.methyl %>%
  filter(Diag_for_meth %in% c("FF", "NTL"))
wilcox.test(stat.test$fetal_mean_meth ~ stat.test$Diag_for_meth)

# Adult promoters

P1_methyl <- ggplot(annot.methyl, aes(x = Diag_for_meth, y=mean_P1 ,fill=Diag_for_meth)) + 
  geom_violin(alpha=.8, size=0.7, trim = F,width=1.8) +
  geom_boxplot(width=.2, cex=.7) +
  scale_y_continuous(limits=c(-0.5, 1.8)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(size=11),
        axis.text.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black", size =1),
        axis.ticks.length = unit(0.3, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

stat.test = annot.methyl %>%
  filter(Diag_for_meth %in% c("NTL", "FF"))
wilcox.test(stat.test$mean_P1 ~ stat.test$Diag_for_meth)

#ggsave("/Users/jill/Desktop/P1_methyl.png", plot = P1_methyl, device = "png", width = 20, height = 10, units = "cm", dpi = 320)
