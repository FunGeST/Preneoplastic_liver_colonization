# libraries -- -- ----------
library(ComplexHeatmap)
library(ggplot2)
library(glue)
library(magrittr)
library(dplyr)

# functions -- -- ---------------

LoadRdata <- function(file_name) {
  if (file.exists(file_name)) return(eval(parse(text = load(file_name))))
  cat(paste("error: file", file_name, "not found"))
  NULL
}

FormatPValue <- function(pval) {
  if (pval >= 0.01) return(round(pval, 3))
  if (pval >= 0.001) return(round(pval, 4))
  return(formatC(pval, format = "e", digits = 2))
}

PlotPng <- function(plot_call,
                    filename, 
                    width = 300, 
                    height = 300, 
                    res = 100,
                    dirname = res_dir,
                    also_print = F) {
  if(also_print) print(plot_call)
  png(file.path(dirname, paste0(glue(filename), ".png")), width = width, height = height, res = res)
  print(plot_call)
  dev.off()
}

# directories -- -- ----------
server_dir <- c("//10.93.23.19", "/Volumes")[1]
hepato_dir <- file.path(server_dir, "hepato partage")
dropbox_dir <- "D:/TZH/Dropbox"
GEPELIN_dropbox_dir <- file.path(dropbox_dir, "GEPELIN_GENOMIC_ANALYSES")
mosaic_dropbox_dir <- file.path(dropbox_dir, "11p15.5 mosaicism")
data_dir <- file.path(GEPELIN_dropbox_dir, "Data")
manu_dir <- file.path(mosaic_dropbox_dir, "MANUSCRIPT")
manu_data_dir <- file.path(manu_dir, "Supp_tables")
main_dir <- file.path(manu_dir, "Nat_com_revisions/new.analyses/TZH_results")
res_dir <- file.path(main_dir, "Classifications")


# load data -- -- -----------
# load ped annot -----------
ped_annot_all <- grep("integrated_pediatric_table.Rdata", list.files(data_dir), value = T) %T>% {print(.)} %>%  
  first() %>% 
  file.path(data_dir, .) %>% 
  LoadRdata()
  
paper_annot <- readxl::read_xlsx(file.path(manu_data_dir, "annot.xlsx"))

HB_tumors <- paper_annot %>% 
  filter(grepl("T", CHCID), grepl("HB", Histological.Diagnosis)) %>% 
  pull(CHCID)

ped_annot_all %>% 
  filter(CHCID %in% HB_tumors) %>% 
  count(status_11p15_ext)


# colors -- -- ---------
yes_no_cols <- c(no = "#8F8F8F", yes = "#B3381D")
yes_no_cols_grey <- c(no = "#E0E0E0", yes = "#B3381D")
wt_alt_cols_grey <- c(wt = "#E0E0E0", alt = "#B3381D")
H_LP_M_cols <- c("H" = "#DE81C0", "LP" = "#500787", "M" = "#A84325")
H_LP_M_limits <- names(H_LP_M_cols)
Nagae_cols <- c("HEP" = "#E4A4C0", "PRO" = "#8D73FF", "MES" = "#A8735C")
Nagae_limits <- names(Nagae_cols)
Cairo_cols <- c("C1" = "forestgreen", "C2" = "red")
Cairo_limits <- names(Cairo_cols)
Hooks_cols <- c("C1" = "forestgreen", "C2A" = "red", "C2B" = "orange")
Hooks_limits <- names(Hooks_cols)
alt_11p_cols <- c("wt"= "white", "GAIN" = "#FFC407", "CDKN1C_mut" = "#C1F0FF", 
                  "LOM_IC2" = "#EF9A9A", "GOM_IC1" = "#4169E1",
                  "DEL-LOH" = "#0000EF", "cn-LOH" = "#B71C1C")
alt_11p_labels <- c("wt", "Pat. dup.", "CDKN1C mut.", "Epimutation IC2", "Epimutation IC1", "LOH", "cn-LOH")

# Analysis -- -- -- -- -- -- -- -- -- -- ------------
# status 11p15 vs different classifications -- -- -----------------
# Hirsch vs 11p15 alt ext ----------
# HB_tumors_Hirsch <- ped_annot_all %>% 
#   filter(CHCID %in% HB_tumors, !is.na(H_LP_M_RNA_fluidigm), !is.na(status_11p15_ext)) %>% 
#   mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", H_LP_M_RNA_fluidigm))

HB_tumors_Hirsch <- paper_annot %>% 
  left_join(select(ped_annot_all, CHCID, H_LP_M_RNA_fluidigm)) %>% 
  filter(CHCID %in% HB_tumors, !is.na(H_LP_M_RNA_fluidigm), status_11p15_ext != "NA") %>% 
  mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", H_LP_M_RNA_fluidigm))

HB_tumors_Hirsch %>% count(status_11p15_ext)

gg_Hirsch_color_alt <- HB_tumors_Hirsch %>% 
  ggplot(aes(x = H_LP_M_RNA_fluidigm, fill = factor(status_11p15_ext, levels = names(alt_11p_cols)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = names(alt_11p_cols), values = alt_11p_cols, labels = alt_11p_labels, "11p15.5 locus") +
  scale_y_continuous(labels = scales::percent, "", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = H_LP_M_limits, "Hirsch 2021 Classification",
                   labels = paste(H_LP_M_limits, table(HB_tumors_Hirsch$H_LP_M_RNA_fluidigm)[H_LP_M_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

chisq.test(table(HB_tumors_Hirsch$status_11p15_ext, HB_tumors_Hirsch$H_LP_M_RNA_fluidigm))
gg_Hirsch_color_alt
gg_Hirsch_color_alt %>% PlotPng("Alt_11p15_detail_vs_Hirsch_2021", width = 400)



gg_Hirsch_color_group <- HB_tumors_Hirsch %>% 
  ggplot(aes(x = H_LP_M_RNA_fluidigm, fill = factor(status_11p15_ext_yes_no, levels = c("no", H_LP_M_limits)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = c("no", H_LP_M_limits), values = c(no = "white", H_LP_M_cols), guide = "none") +
  scale_y_continuous(labels = scales::percent, "Proportion of 11p15.5 alterations", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = H_LP_M_limits, "Hirsch 2021 Classification",
                   labels = paste(H_LP_M_limits, table(HB_tumors_Hirsch$H_LP_M_RNA_fluidigm)[H_LP_M_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) 

chisq.test(table(HB_tumors_Hirsch$status_11p15_ext == "wt", HB_tumors_Hirsch$H_LP_M_RNA_fluidigm))
gg_Hirsch_color_group
gg_Hirsch_color_group %>% PlotPng("Alt_11p15_vs_Hirsch_2021")

# Nagae vs 11p15 alt ext ----------
# HB_tumors_Nagae <- ped_annot_all %>% 
#   filter(CHCID %in% HB_tumors, !is.na(Nagae_RNA_fluidigm), !is.na(status_11p15_ext)) %>% 
#   mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", Nagae_RNA_fluidigm))

HB_tumors_Nagae <- paper_annot %>% 
  left_join(select(ped_annot_all, CHCID, Nagae_RNA_fluidigm)) %>% 
  filter(CHCID %in% HB_tumors, !is.na(Nagae_RNA_fluidigm), status_11p15_ext != "NA") %>% 
  mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", Nagae_RNA_fluidigm))

gg_Nagae_color_alt <- HB_tumors_Nagae %>% 
  ggplot(aes(x = Nagae_RNA_fluidigm, fill = factor(status_11p15_ext, levels = names(alt_11p_cols)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = names(alt_11p_cols), values = alt_11p_cols, labels = alt_11p_labels, "11p15.5 locus") +
  scale_y_continuous(labels = scales::percent, "", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = Nagae_limits, "Nagae 2021 Classification",
                   labels = paste(Nagae_limits, table(HB_tumors_Nagae$Nagae_RNA_fluidigm)[Nagae_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

chisq.test(table(HB_tumors_Nagae$status_11p15_ext, HB_tumors_Nagae$Nagae_RNA_fluidigm))
gg_Nagae_color_alt
gg_Nagae_color_alt %>% PlotPng("Alt_11p15_detail_vs_Nagae_2021", width = 400)

gg_Nagae_color_group <- HB_tumors_Nagae %>% 
  ggplot(aes(x = Nagae_RNA_fluidigm, fill = factor(status_11p15_ext_yes_no, levels = c("no", Nagae_limits)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = c("no", Nagae_limits), values = c(no = "white", Nagae_cols), guide = "none") +
  scale_y_continuous(labels = scales::percent, "Proportion of 11p15.5 alterations", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = Nagae_limits, "Nagae 2021 Classification",
                   labels = paste(Nagae_limits, table(HB_tumors_Nagae$Nagae_RNA_fluidigm)[Nagae_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) 

chisq.test(table(HB_tumors_Nagae$status_11p15_ext == "wt", HB_tumors_Nagae$Nagae_RNA_fluidigm))
gg_Nagae_color_group
gg_Nagae_color_group %>% PlotPng("Alt_11p15_vs_Nagae_2021")


# Hooks vs 11p15 alt ext ----------
# HB_tumors_Hooks <- ped_annot_all %>% 
#   filter(CHCID %in% HB_tumors, !is.na(Hooks_RNA_fluidigm), !is.na(status_11p15_ext)) %>% 
#   mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", Hooks_RNA_fluidigm))

HB_tumors_Hooks <- paper_annot %>% 
  left_join(select(ped_annot_all, CHCID, Hooks_RNA_fluidigm)) %>% 
  filter(CHCID %in% HB_tumors, !is.na(Hooks_RNA_fluidigm), status_11p15_ext != "NA") %>% 
  mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", Hooks_RNA_fluidigm))


gg_Hooks_color_alt <- HB_tumors_Hooks %>% 
  ggplot(aes(x = Hooks_RNA_fluidigm, fill = factor(status_11p15_ext, levels = names(alt_11p_cols)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = names(alt_11p_cols), values = alt_11p_cols, labels = alt_11p_labels, "11p15.5 locus") +
  scale_y_continuous(labels = scales::percent, "", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = Hooks_limits, "Hooks 2018 Classification",
                   labels = paste(Hooks_limits, table(HB_tumors_Hooks$Hooks_RNA_fluidigm)[Hooks_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

chisq.test(table(HB_tumors_Hooks$status_11p15_ext, HB_tumors_Hooks$Hooks_RNA_fluidigm))
gg_Hooks_color_alt
gg_Hooks_color_alt %>% PlotPng("Alt_11p15_detail_vs_Hooks_2018", width = 400)

gg_Hooks_color_group <- HB_tumors_Hooks %>% 
  ggplot(aes(x = Hooks_RNA_fluidigm, fill = factor(status_11p15_ext_yes_no, levels = c("no", Hooks_limits)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = c("no", Hooks_limits), values = c(no = "white", Hooks_cols), guide = "none") +
  scale_y_continuous(labels = scales::percent, "Proportion of 11p15.5 alterations", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = Hooks_limits, "Hooks 2018 Classification",
                   labels = paste(Hooks_limits, table(HB_tumors_Hooks$Hooks_RNA_fluidigm)[Hooks_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) 

chisq.test(table(HB_tumors_Hooks$status_11p15_ext == "wt", HB_tumors_Hooks$Hooks_RNA_fluidigm))
gg_Hooks_color_group
gg_Hooks_color_group %>% PlotPng("Alt_11p15_vs_Hooks_2018")

# Cairo vs 11p15 alt ext ----------
# HB_tumors_Cairo <- ped_annot_all %>% 
#   filter(CHCID %in% HB_tumors, !is.na(Cairo_RNA_fluidigm), !is.na(status_11p15_ext)) %>% 
#   mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", Cairo_RNA_fluidigm))

HB_tumors_Cairo <- paper_annot %>% 
  left_join(select(ped_annot_all, CHCID, Cairo_RNA_fluidigm)) %>% 
  filter(CHCID %in% HB_tumors, !is.na(Cairo_RNA_fluidigm), status_11p15_ext != "NA") %>% 
  mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "no", Cairo_RNA_fluidigm))


gg_Cairo_color_alt <- HB_tumors_Cairo %>% 
  ggplot(aes(x = Cairo_RNA_fluidigm, fill = factor(status_11p15_ext, levels = names(alt_11p_cols)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = names(alt_11p_cols), values = alt_11p_cols, labels = alt_11p_labels, "11p15.5 locus") +
  scale_y_continuous(labels = scales::percent, "", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = Cairo_limits, "Cairo 2008 Classification",
                   labels = paste(Cairo_limits, table(HB_tumors_Cairo$Cairo_RNA_fluidigm)[Cairo_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

chisq.test(table(HB_tumors_Cairo$status_11p15_ext, HB_tumors_Cairo$Cairo_RNA_fluidigm))
gg_Cairo_color_alt
gg_Cairo_color_alt %>% PlotPng("Alt_11p15_detail_vs_Cairo_2008", width = 400)

gg_Cairo_color_group <- HB_tumors_Cairo %>% 
  ggplot(aes(x = Cairo_RNA_fluidigm, fill = factor(status_11p15_ext_yes_no, levels = c("no", Cairo_limits)))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(limits = c("no", Cairo_limits), values = c(no = "white", Cairo_cols), guide = "none") +
  scale_y_continuous(labels = scales::percent, "Proportion of 11p15.5 alterations", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = Cairo_limits, "Cairo 2008 Classification",
                   labels = paste(Cairo_limits, table(HB_tumors_Cairo$Cairo_RNA_fluidigm)[Cairo_limits],
                                  sep = "\nn=")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) 

chisq.test(table(HB_tumors_Cairo$status_11p15_ext == "wt", HB_tumors_Cairo$Cairo_RNA_fluidigm))
gg_Cairo_color_group
gg_Cairo_color_group %>% PlotPng("Alt_11p15_vs_Cairo_2008")


################## NOT USED IN THE PAPER ###############################
# 
# # load RNAseq data ----------
# exp <- LoadRdata(file.path(hepato_dir, "GEPELIN/RNAseq/Expression_matrix/mos_paper_by_gene_name/exp_190s.Rdata"))
# 
# # gene expression among different classifications -- -- -----------------
# HB_tumors_Hirsch_RNAseq <- paper_annot %>% 
#   left_join(select(ped_annot_all, CHCID, H_LP_M_RNA_fluidigm)) %>% 
#   filter(CHCID %in% HB_tumors, !is.na(H_LP_M_RNA_fluidigm), status_11p15_ext != "NA", CHCID %in% colnames(exp)) %>% 
#   mutate(status_11p15_ext_yes_no = ifelse(status_11p15_ext == "wt", "wt", "alt")) %>% 
#   left_join(data.frame(CHCID = colnames(exp), IGF2_exp = exp["IGF2", ])) %>% 
#   group_by(H_LP_M_RNA_fluidigm) %>% 
#   mutate(for_facet = glue("{H_LP_M_RNA_fluidigm} (n={n()})"))
# 
# HB_tumors_Hirsch_RNAseq %>% 
#   arrange(IGF2_exp)
# 
# gg_Hirsch_IGF2 <- HB_tumors_Hirsch_RNAseq %>% 
#   ggplot(aes(x = status_11p15_ext_yes_no, y = IGF2_exp,
#              fill = status_11p15_ext_yes_no)) +
#   geom_boxplot(outlier.shape = NA) +
#   scale_fill_manual(values = wt_alt_cols_grey, "11p15.5 locus", guide = "none") +
#   geom_point(position = position_jitter(width = .1)) +
#   facet_wrap(~ for_facet) +
#   scale_y_continuous("IGF2 expression") +
#   scale_x_discrete("11p15.5 locus", limits = names(wt_alt_cols_grey)) +
#   theme_classic() +
#   theme(axis.text = element_text(color = "black"))
# 
# gg_Hirsch_IGF2
# gg_Hirsch_IGF2 %>% PlotPng("IGF2_exp_vs_Hirsch_2021", width = 400)
# 
