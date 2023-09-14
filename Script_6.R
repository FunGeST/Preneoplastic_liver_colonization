### 0. Libraries and functions -------------------------------------------------
library(ggplot2)
library(magrittr)
library(tidyr)
library(ggrepel)
library(scales)
library(matrixStats)
library(Hmisc)
library(ggforce)
library(ggpubr)
library(GSVA)
library(plyr)
library(dplyr)
library(rstatix)
library(survival)
library(rms) # for Kapplan Meier
library(survminer)
library(geco.utils)
library(readxl)
library(data.table)

### 1. Load data ---------------------------------------------------------------
annot <- read_xlsx("D:/Dropbox/11p15.5 mosaicism/Capsule_temporelle_Integrated_table/Pediatric_table.xlsx")
annot. = as.data.frame(annot)

dim(annot.)
table(annot.$Type)

### 2. Add new annotations ------------------------------------------------------

`%nin%` = Negate(`%in%`)

annot. = annot. %>%
  mutate(Mosaic = ifelse(is.na(Mosaic), "wt", Mosaic))

# Add missing data
annot.[annot.$CHCID=="CHC3538S", "Age.at.surgery.recode.months"] <- 126.2
annot.[annot.$CHCID=="CHC3122S", "Age.at.surgery.recode.months"] <- 190.9

# Create Histological diagnosis variable with mosaic status
annot. = annot. %>%
  mutate(Histological.Diagnosis = case_when(Histological.Diagnosis=="HB.META" ~"HB", 
                                             Histological.Diagnosis=="HB.RELAPSE" ~"HB",
                                             Mosaic == "Mosaic" ~ "Mosaic",
                                             Histological.Diagnosis=="FF" ~"FF",
                                             Histological.Diagnosis=="FLC" ~"FLC",
                                             Histological.Diagnosis=="HCA" ~"HCA",
                                             Histological.Diagnosis=="HB" ~"HB",
                                             Histological.Diagnosis=="HCC" ~"HCC"),
         Age.surg.years = case_when(grepl("FF", CHCID) ~ (as.numeric(sub("SA,", "", sub("Foie fetal ", "", Miror.Bloc)))-40)/12.25,
                                    !is.na(Age.at.surgery.recode.months) ~ as.numeric(Age.at.surgery.recode.months)/12.25,
                                    T ~ NA_real_))

# Add missing annotation for two patients without T  

annot. = annot. %>%
  group_by(Patient.identification) %>% 
  mutate(M_by_Tumor = ifelse(Patient.identification %in% c("3185", "3577"), last(M_by_Tumor), M_by_Tumor),
         E_by_Tumor = ifelse(Patient.identification %in% c("3185", "3577"), last(E_by_Tumor), E_by_Tumor),
         F_by_Tumor = ifelse(Patient.identification %in% c("3185", "3577"), last(F_by_Tumor), F_by_Tumor),
         RNAseq_clustering_3 = ifelse(Patient.identification %in% c("3185", "3577"), last(RNAseq_clustering_3), RNAseq_clustering_3),
         NFE2L2 = ifelse(Patient.identification %in% c("3185", "3577"), last(NFE2L2), NFE2L2),
         TERT = ifelse(Patient.identification %in% c("3185", "3577"), last(TERT), TERT),
         PRETEXT_Stade = ifelse(Patient.identification %in% c("3185", "3577"), last(PRETEXT_Stade), PRETEXT_Stade),
         Tumor.Size..50mm.or...50mm  = ifelse(Patient.identification %in% c("3185", "3577"), last(Tumor.Size..50mm.or...50mm ), Tumor.Size..50mm.or...50mm),
         PRETEXT_M  = ifelse(Patient.identification %in% c("3185", "3577"), last(PRETEXT_M ), PRETEXT_M),
         DSS.status  = ifelse(Patient.identification %in% c("3185", "3577"), last(DSS.status), DSS.status),
         PFS.status  = ifelse(Patient.identification %in% c("3185", "3577"), last(PFS.status), PFS.status),
         TTD  = ifelse(Patient.identification %in% c("3185", "3577"), last(TTD), TTD)) %>%
  ungroup()

dim(annot.)

### 3. Recode to test significant associations --------------------------------

annot.asso <- annot. %>% 
  rowwise() %>% 
  mutate(More_than_1_nodule =  ifelse(Number.of.nodules > 1, "yes", "no"),
         RNAseq_LP_vs_others = ifelse(RNAseq_clustering_4 == "E", "LP", ifelse(!is.na(RNAseq_clustering_4), "other", NA_character_)),
         RNAseq_H_vs_others = ifelse(RNAseq_clustering_4 %in% c("F1", "F2"), "H", ifelse(!is.na(RNAseq_clustering_4), "other", NA_character_)),
         RNAseq_M_vs_others = ifelse(RNAseq_clustering_4 == "M", "M", ifelse(!is.na(RNAseq_clustering_4), "other", NA_character_)),
         F_by_Tumor_sup30 = ifelse(F_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         E_by_Tumor_sup30 = ifelse(E_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no"),
         E_by_Tumor_sup20 = ifelse(E_by_Tumor %in% paste0("sup", c(20, 30, 60, 90)), "yes", "no"),
         M_by_Tumor_sup30 = ifelse(M_by_Tumor %in% paste0("sup", c(30, 60, 90)), "yes", "no")) %>%
  as.data.frame()


# Manual correction of 11p15 status from mosaic patients and series
annot.asso[which(annot.asso$CHCID=="CHC05303T"),"To_keep_for_survival_analysis_worst"] ="yes"
annot.asso[which(annot.asso$CHCID=="CHC3126T"), "To_keep_for_survival_analysis_worst"] = "yes"
annot.asso[which(annot.asso$CHCID=="CHC4087T"), "To_keep_for_survival_analysis_worst"] = "yes"

# Create list of non-tumor samples analyzed in paper
annot.NT = annot. %>%
  filter(Papier_mosaiques_2021=="yes",
         grepl("N|S|FF", CHCID),
         Histological.Diagnosis %in% c("HB", "Mosaic")) %>%
  as.data.frame()  
list.NT = annot.NT$Patient.identification
length(list.NT)

# Filter the annot table 
annot.asso = annot.asso %>%
  filter(Patient.identification %in% list.NT,  # Restriction to non-tumor samples analyzed
         grepl("T", CHCID))  # keep only T

### 4. Create drivers annotation per patient------------------------------------
annot.asso = annot.asso %>%
  group_by(Patient.identification) %>%
  mutate(status_11p_by_patient = case_when(any(status_11p15_ext=="cn-LOH") ~ "cn-LOH",
                                           any(status_11p15_ext=="DEL-LOH") ~ "DEL-LOH",
                                           any(status_11p15_ext=="CDKN1C_mut") ~ "CDKN1C_mut",
                                           any(status_11p15_ext=="GOM_IC1") ~ "GOM_IC1",
                                           any(status_11p15_ext=="LOM_IC2") ~ "LOM_IC2",
                                           any(status_11p15_ext=="GAIN") ~ "GAIN",
                                           any(status_11p15_ext=="wt") ~ "wt",
                                           is.na(status_11p15_ext) ~ NA_character_),
         status_11p_alt_wt = case_when(status_11p_by_patient %in% c("cn-LOH", "DEL-LOH", "CDKN1C_mut", "GOM_IC1", "LOM_IC2", "GAIN") ~ "alt",
                                       status_11p_by_patient =="wt"~  "wt",
                                       is.na(status_11p_by_patient) ~ NA_character_),
         TERT_by_patient = case_when(any(TERT=="promoter_c.-124") ~ "alt",
                                     any(TERT=="promoter_c.-146") ~ "alt",
                                     any(TERT=="wt") ~ "wt",
                                     is.na(TERT) ~ NA_character_),
         NFE2L2_by_patient = case_when(any(NFE2L2=="missense") ~ "missense",
                                       any(NFE2L2=="wt") ~ "wt",
                                       is.na(NFE2L2) ~ NA_character_)) %>%
  ungroup()



# Add 11p15.5 timing info 

timing.11p.dir <- "D:/Dropbox/11p15.5 mosaicism/Timing_11p/"
timing.11p <- read_xlsx(file.path(timing.11p.dir, "Timing_11p_from_SuppFig5_Cancer_Discovery.xlsx")) %>% as.data.frame()
timing.11p$Patient.identification = as.character(timing.11p$Patient.identification)
annot.asso$Patient.identification = as.character(annot.asso$Patient.identification)
annot.asso = left_join(annot.asso, timing.11p, by="Patient.identification")

annot.asso = annot.asso %>%
  mutate(timing_11p = case_when(Mosaic=="Mosaic" ~ "non_private",
                                timing_11p=="trunk" ~  "non_private",
                                timing_11p=="private" ~  "private",
                                timing_11p == "wt" ~ "non_private",
                                is.na(timing_11p) ~ "non_private")) %>%
  
  as.data.frame()

dim(annot.asso)
annot.asso$timing_11p
annot.asso %>% select(c("CHCID", "timing_11p", "Timing_alteration"))


# Add 2 missing data 
annot.asso[which(annot.asso$CHCID=="CHC3082T"),"timing_11p"] ="private"  # Add 2 private beacause 11p alteration in one but not in another sample but no tree
annot.asso[which(annot.asso$CHCID=="CHC3575T"),"timing_11p"] ="private"  # Add 2 private beacause 11p alteration in one but not in another sample but no tree


# final filtering 
annot.asso = annot.asso %>%
  filter(Histological.Diagnosis %in% c("HB", "Mosaic"), # Analyze only HB
         To_keep_for_survival_analysis_worst %in% c("yes", "yes?"))  # Keep one sample per patient

dim(annot.asso)
table(annot.asso$Histological.Diagnosis)

table(annot.asso$timing_11p)
table(annot.asso$Mosaic)
table(annot.asso$status_11p_alt_wt)

# Create table for comparision late vs early alteration
annot.assoT = annot.asso %>% filter(Papier_mosaiques_2021=="yes") 
dim(annot.assoT)

#  create heterogeneity annot 
patient_het <- annot. %>% 
  filter(grepl("T", CHCID), grepl("HB", Histological.Diagnosis)) %>% 
  group_by(Patient.identification) %>% 
  summarise(anyNA = any(is.na(M_by_Tumor)),
            anyF = any(F_by_Tumor != "0"),
            anyE = any(E_by_Tumor != "0"),
            anyM = any(M_by_Tumor != "0"),
            anyC = any(C_by_Tumor != "0"),
            anySCUD = any(SCUD_by_Tumor != "0")) %>% 
  rowwise() %>% 
  mutate(Patient_plasticity = ifelse(sum(c(anyF, anyE, anyM, anyC, anySCUD)) > 1, "hete", "homo"),
         Patient_plasticity = ifelse(anyNA, "nana", Patient_plasticity),
         Patient_plasticity = ifelse(Patient.identification %in% c(4992, # na in table but hete in fiche anapath
                                                                   3122, # meta fetal pure but primary (no NGS) is hete
                                                                   3168),  # biopsy fetal pure but primary (no NGS) is hete
                                     "hete", Patient_plasticity))

tumor_het <- annot. %>% 
  filter(grepl("T", CHCID), 
         !is.na(To_keep_for_survival_analysis_worst) | Patient.identification %in% c(3511, 3970)) %>% # 3511 and 3970 yes are on the sample without RNAseq 
  rowwise() %>% 
  mutate(Tumor_plasticity = ifelse(sum(c(F_by_Tumor != "0", E_by_Tumor != "0", 
                                         M_by_Tumor != "0", C_by_Tumor != "0", 
                                         SCUD_by_Tumor != "0")) > 1, "hete", "homo"),
         Tumor_plasticity = ifelse(Patient.identification == "4992", "hete", Tumor_plasticity)) %>% # fiche anapath
  select(Patient.identification, Tumor_plasticity)

# remove 3959/3122 (only recurrent, no primary) and 3168/3966 (only biopsy, not possible to asses heterogeneity on the whole primary tumor)

tumor_het = tumor_het %>%
  mutate(Tumor_plasticity = ifelse(Patient.identification %in% c("3959", "3122", "3168", "3966"), NA_character_, Tumor_plasticity))
tumor_het$Patient.identification = as.character(tumor_het$Patient.identification)
patient_het$Patient.identification = as.character(patient_het$Patient.identification)

annot.het = left_join(annot.asso, tumor_het, by="Patient.identification")
annot.het = left_join(annot.het, patient_het, by="Patient.identification")
annot.het = annot.het[which(!duplicated(annot.het$Patient.identification)),]


# Export Source data files 
Source_supp_table_1 <- annot.asso %>%
  select(c("CHCID", 
           "Mosaic",
           "timing_11p",
           "status_11p_alt_wt",
           "Age.surg.years", 
           "Gender", 
           "E_by_Tumor_sup30", 
           "F_by_Tumor_sup30",
           "M_by_Tumor_sup30",
           "RNAseq_LP_vs_others",
           "RNAseq_H_vs_others",
           "RNAseq_M_vs_others",
           "NFE2L2_by_patient",
           "TERT_by_patient",
           "PRETEXT_Stade",
           "Tumor.Size..50mm.or...50mm",
           "More_than_1_nodule",
           "PRETEXT_M",
           "TTD",
           "TTP",
           "DSS.status",
           "PFS.status"))

write.table(Source_supp_table_1, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/Source_supp_table_1.txt", sep="\t")


### 5. Make statistical tests  -------------------------------------------------
  # 5.0 Age --------------------------------------------------------------------
t.test(annot.asso$Age.surg.years ~ annot.asso$Mosaic)
annot.test = annot.asso %>% filter(Mosaic=="Mosaic")
median(annot.test$Age.surg.years)
annot.test = annot.asso %>% filter(Mosaic=="wt")
median(annot.test$Age.surg.years)

annot.test = annot.assoT %>% filter(status_11p_alt_wt=="wt")
median(annot.test$Age.surg.years)
annot.test = annot.assoT %>% filter(status_11p_alt_wt=="alt")
median(annot.test$Age.surg.years)
table(annot.assoT$status_11p_alt_wt)
t.test(annot.assoT$Age.surg.years ~ annot.assoT$status_11p_alt_wt)

annot.test = annot.assoT %>% filter(timing_11p=="private")
median(annot.test$Age.surg.years)
annot.test = annot.assoT %>% filter(timing_11p=="non_private")
median(annot.test$Age.surg.years)
t.test(annot.assoT$Age.surg.years ~ annot.assoT$timing_11p)

wilcox.test(annot.asso$Age.surg.years ~ annot.asso$Gender)

test = annot.asso %>% filter(!is.na(E_by_Tumor_sup30))
wilcox.test(test$Age.surg.years ~ test$E_by_Tumor_sup30)

test = annot.asso %>% filter(!is.na(F_by_Tumor_sup30))
wilcox.test(test$Age.surg.years ~ test$F_by_Tumor_sup30)

test = annot.asso %>% filter(!is.na(M_by_Tumor_sup30))
wilcox.test(test$Age.surg.years ~ test$M_by_Tumor_sup30)

test = annot.asso %>% filter(!is.na(Tumor.Size..50mm.or...50mm))
wilcox.test(test$Age.surg.years ~ test$Tumor.Size..50mm.or...50mm)

test = annot.asso %>% filter(!is.na(More_than_1_nodule))
wilcox.test(test$Age.surg.years ~ test$More_than_1_nodule)

test = annot.asso %>% filter(!is.na(PRETEXT_M))
wilcox.test(test$Age.surg.years ~ test$PRETEXT_M)

test = annot.asso %>% filter(!is.na(NFE2L2_by_patient))
wilcox.test(test$Age.surg.years ~ test$NFE2L2_by_patient)

test = annot.asso %>% filter(!is.na(TERT_by_patient))
wilcox.test(test$Age.surg.years ~ test$TERT_by_patient)

test = annot.het %>% filter(!is.na(Tumor_plasticity))
wilcox.test(test$Age.surg.years ~ test$Tumor_plasticity)

test = annot.asso %>% filter(!is.na(RNAseq_LP_vs_others))
wilcox.test(test$Age.surg.years ~ test$RNAseq_LP_vs_others)

test = annot.asso %>% filter(!is.na(RNAseq_M_vs_others))
wilcox.test(test$Age.surg.years ~ test$RNAseq_M_vs_others)

test = annot.asso %>% filter(!is.na(RNAseq_H_vs_others))
wilcox.test(test$Age.surg.years ~ test$RNAseq_H_vs_others)

kruskal.test(annot.asso$Age.at.surgery.recode.months ~ annot.asso$PRETEXT_Stade)

  # 5.1 Sex --------------------------------------------------------------------
    # 5.1.1 Compare Mosaic and wt samples --------------------------------------
annot.test = annot.asso %>%
  group_by(Mosaic, Gender) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))
annot.test

# Age adjusted
annot.stat = annot.asso
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$Gender = factor(annot.stat$Gender, levels = c("F", "M"))
model2 <- glm(Gender ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(Gender))

annot.stat = as.data.frame(table(annot.stat$Gender, annot.stat$Mosaic))
colnames(annot.stat) = c("Gender", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Gender
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct = F)

    # 5.1.2 Compare late 11p ---------------------------------------------------

annot.test = annot.assoT %>%
  group_by(timing_11p, Gender) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT
annot.stat = as.data.frame(table(annot.stat$Gender, annot.stat$timing_11p))
colnames(annot.stat) = c("Gender", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Gender
annot.stat = annot.stat %>% 
  select(c("non_private", "private"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$Gender = factor(annot.stat$Gender, levels = c("F", "M"))
model2 <- glm(Gender ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.1.3 Compare 11palt vs wt -----------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, Gender) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat = as.data.frame(table(annot.stat$Gender, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("Gender", "status_11p_alt_wt", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Gender
annot.stat = annot.stat %>% 
  select(c("wt", "alt"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("wt", "alt"))
annot.stat$Gender = factor(annot.stat$Gender, levels = c("F", "M"))
model2 <- glm(Gender ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.2 Histology E tumor  -----------------------------------------------------
    # 5.2.1 Compare mosaic and wt ----------------------------------------------
annot.test = annot.asso %>% 
  filter(!is.na(E_by_Tumor_sup30)) %>% 
  group_by(Mosaic, E_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot.test

annot.stat = annot.asso %>%
  filter(!is.na(E_by_Tumor_sup30))
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$E_by_Tumor_sup30 = factor(annot.stat$E_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(E_by_Tumor_sup30 ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(E_by_Tumor_sup30))
annot.stat = as.data.frame(table(annot.stat$E_by_Tumor_sup30, annot.stat$Mosaic))
colnames(annot.stat) = c("E_by_Tumor_sup30", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$E_by_Tumor_sup30
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct = F)

    # 5.2.2 Compare timing 11p -------------------------------------------------

annot.test = annot.assoT %>%
  group_by(timing_11p, E_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.asso 
annot.stat = as.data.frame(table(annot.stat$E_by_Tumor_sup30, annot.stat$timing_11p))
colnames(annot.stat) = c("E_by_Tumor_sup30", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$E_by_Tumor_sup30
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$E_by_Tumor_sup30 = factor(annot.stat$E_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(E_by_Tumor_sup30 ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.2.3 Compare 11palt vs wt -----------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, E_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat = as.data.frame(table(annot.stat$E_by_Tumor_sup30, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("E_by_Tumor_sup30", "status_11p_alt_wt", "Freq")


annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 
rownames(annot.stat) = annot.stat$E_by_Tumor_sup30

annot.stat = annot.stat %>% 
  select(c("wt", "alt"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("wt", "alt"))
annot.stat$E_by_Tumor_sup30 = factor(annot.stat$E_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(E_by_Tumor_sup30 ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.3 Histology F tumor  -----------------------------------------------------
    # 5.3.1 Compare mosaic and wt ----------------------------------------------
annot.test = annot.asso %>% 
  filter(!is.na(F_by_Tumor_sup30)) %>% 
  group_by(Mosaic, F_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot.test

annot.stat = annot.asso %>%
  filter(!is.na(F_by_Tumor_sup30))
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$F_by_Tumor_sup30 = factor(annot.stat$F_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(F_by_Tumor_sup30 ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(F_by_Tumor_sup30))
annot.stat = as.data.frame(table(annot.stat$F_by_Tumor_sup30, annot.stat$Mosaic))
colnames(annot.stat) = c("F_by_Tumor_sup30", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$F_by_Tumor_sup30
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct = F)

    # 5.3.2 Compare timing 11p -------------------------------------------------

annot.test = annot.assoT %>%
  group_by(timing_11p, F_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.asso 
annot.stat = as.data.frame(table(annot.stat$F_by_Tumor_sup30, annot.stat$timing_11p))
colnames(annot.stat) = c("F_by_Tumor_sup30", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$F_by_Tumor_sup30
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$F_by_Tumor_sup30 = factor(annot.stat$F_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(F_by_Tumor_sup30 ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.3.3 Compare 11palt vs wt -----------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, F_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat = as.data.frame(table(annot.stat$F_by_Tumor_sup30, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("F_by_Tumor_sup30", "status_11p_alt_wt", "Freq")


annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 
rownames(annot.stat) = annot.stat$F_by_Tumor_sup30

annot.stat = annot.stat %>% 
  select(c("wt", "alt"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("wt", "alt"))
annot.stat$F_by_Tumor_sup30 = factor(annot.stat$F_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(F_by_Tumor_sup30 ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.4 Histology M tumor  -----------------------------------------------------
    # 5.4.1 Compare mosaic and wt ----------------------------------------------
annot.test = annot.asso %>% 
  filter(!is.na(M_by_Tumor_sup30)) %>% 
  group_by(Mosaic, M_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot.test

annot.stat = annot.asso %>%
  filter(!is.na(M_by_Tumor_sup30))
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$M_by_Tumor_sup30 = factor(annot.stat$M_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(M_by_Tumor_sup30 ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(M_by_Tumor_sup30))
annot.stat = as.data.frame(table(annot.stat$M_by_Tumor_sup30, annot.stat$Mosaic))
colnames(annot.stat) = c("M_by_Tumor_sup30", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$M_by_Tumor_sup30
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct = F)

    # 5.4.2 Compare timing 11p -------------------------------------------------

annot.test = annot.assoT %>%
  group_by(timing_11p, M_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.asso 
annot.stat = as.data.frame(table(annot.stat$M_by_Tumor_sup30, annot.stat$timing_11p))
colnames(annot.stat) = c("M_by_Tumor_sup30", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$M_by_Tumor_sup30
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$M_by_Tumor_sup30 = factor(annot.stat$M_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(M_by_Tumor_sup30 ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


    # 5.4.3 Compare 11palt vs wt -----------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, M_by_Tumor_sup30) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat = as.data.frame(table(annot.stat$M_by_Tumor_sup30, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("M_by_Tumor_sup30", "status_11p_alt_wt", "Freq")


annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 
rownames(annot.stat) = annot.stat$M_by_Tumor_sup30

annot.stat = annot.stat %>% 
  select(c("wt", "alt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("wt", "alt"))
annot.stat$M_by_Tumor_sup30 = factor(annot.stat$M_by_Tumor_sup30, levels = c("yes", "no"))
model2 <- glm(M_by_Tumor_sup30 ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.5 RNAseq M vs others -----------------------------------------------------
    # 5.5.1 mosaic vs non mosaic -----------------------------------------------
annot.test = annot.asso %>%
  filter(!is.na(RNAseq_M_vs_others)) %>%
  group_by(Mosaic, RNAseq_M_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))
annot.test

annot.stat = annot.asso %>%
  filter(!is.na(RNAseq_M_vs_others))

annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$RNAseq_M_vs_others = factor(annot.stat$RNAseq_M_vs_others, levels = c("M", "other"))
model2 <- glm(RNAseq_M_vs_others ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(RNAseq_M_vs_others))
annot.stat = as.data.frame(table(annot.stat$RNAseq_M_vs_others, annot.stat$Mosaic))
colnames(annot.stat) = c("RNAseq_M_vs_others", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$RNAseq_M_vs_others
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct =F)

    # 5.5.2 Compare timing 11p -------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(RNAseq_M_vs_others)) %>%
  group_by(timing_11p, RNAseq_M_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT  %>%
  filter(!is.na(RNAseq_M_vs_others))

annot.stat = as.data.frame(table(annot.stat$RNAseq_M_vs_others, annot.stat$timing_11p))
colnames(annot.stat) = c("RNAseq_M_vs_others", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$RNAseq_M_vs_others
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT  %>%
  filter(!is.na(RNAseq_M_vs_others))
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$RNAseq_M_vs_others = factor(annot.stat$RNAseq_M_vs_others, levels = c("M", "other"))
model2 <- glm(RNAseq_M_vs_others ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.5.3 11p alt vs wt  -----------------------------------------------------
annot.test = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt),
         !is.na(RNAseq_M_vs_others)) %>%
  group_by(status_11p_alt_wt, RNAseq_M_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt),
         !is.na(RNAseq_M_vs_others))
annot.stat = as.data.frame(table(annot.stat$RNAseq_M_vs_others, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("RNAseq_M_vs_others", "status_11p_alt_wt", "Freq")


annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 
rownames(annot.stat) = annot.stat$RNAseq_M_vs_others

annot.stat = annot.stat %>% 
  select(c("wt", "alt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("wt", "alt"))
annot.stat$RNAseq_M_vs_others = factor(annot.stat$RNAseq_M_vs_others, levels = c("M", "other"))
model2 <- glm(RNAseq_M_vs_others ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.6 RNAseq LP vs others -----------------------------------------------------------
    # 5.6.1 mosaic vs non mosaic -----------------------------------------------
annot.test = annot.asso %>%
  filter(!is.na(RNAseq_LP_vs_others)) %>%
  group_by(Mosaic, RNAseq_LP_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))
annot.test

annot.stat = annot.asso %>%
  filter(!is.na(RNAseq_LP_vs_others))

annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$RNAseq_LP_vs_others = factor(annot.stat$RNAseq_LP_vs_others, levels = c("LP", "other"))
model2 <- glm(RNAseq_LP_vs_others ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(RNAseq_LP_vs_others))
annot.stat = as.data.frame(table(annot.stat$RNAseq_LP_vs_others, annot.stat$Mosaic))
colnames(annot.stat) = c("RNAseq_LP_vs_others", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$RNAseq_LP_vs_others
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct =F)

    # 5.6.2 Compare timing 11p -------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(RNAseq_LP_vs_others)) %>%
  group_by(timing_11p, RNAseq_LP_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT  %>%
  filter(!is.na(RNAseq_LP_vs_others))

annot.stat = as.data.frame(table(annot.stat$RNAseq_LP_vs_others, annot.stat$timing_11p))
colnames(annot.stat) = c("RNAseq_LP_vs_others", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$RNAseq_LP_vs_others
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT  %>%
  filter(!is.na(RNAseq_LP_vs_others))
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$RNAseq_LP_vs_others = factor(annot.stat$RNAseq_LP_vs_others, levels = c("LP", "other"))
model2 <- glm(RNAseq_LP_vs_others ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.6.3 11p alt vs wt  -----------------------------------------------------
annot.test = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt),
         !is.na(RNAseq_LP_vs_others)) %>%
  group_by(status_11p_alt_wt, RNAseq_LP_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt),
         !is.na(RNAseq_LP_vs_others))
annot.stat = as.data.frame(table(annot.stat$RNAseq_LP_vs_others, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("RNAseq_LP_vs_others", "status_11p_alt_wt", "Freq")


annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 
rownames(annot.stat) = annot.stat$RNAseq_LP_vs_others

annot.stat = annot.stat %>% 
  select(c("wt", "alt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("wt", "alt"))
annot.stat$RNAseq_LP_vs_others = factor(annot.stat$RNAseq_LP_vs_others, levels = c("LP", "other"))
model2 <- glm(RNAseq_LP_vs_others ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.7 RNAseq H vs others -------------------------------------------------------
    # 5.7.1 mosaic vs non mosaic ---------------------------------------------
annot.test = annot.asso %>%
  filter(!is.na(RNAseq_H_vs_others)) %>%
  group_by(Mosaic, RNAseq_H_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))
annot.test

annot.stat = annot.asso %>%
  filter(!is.na(RNAseq_H_vs_others))

annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$RNAseq_H_vs_others = factor(annot.stat$RNAseq_H_vs_others, levels = c("H", "other"))
model2 <- glm(RNAseq_H_vs_others ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(RNAseq_H_vs_others))
annot.stat = as.data.frame(table(annot.stat$RNAseq_H_vs_others, annot.stat$Mosaic))
colnames(annot.stat) = c("RNAseq_H_vs_others", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$RNAseq_H_vs_others
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct =F)

    # 5.7.2 Compare timing 11p -----------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(RNAseq_H_vs_others)) %>%
  group_by(timing_11p, RNAseq_H_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT  %>%
  filter(!is.na(RNAseq_H_vs_others))

annot.stat = as.data.frame(table(annot.stat$RNAseq_H_vs_others, annot.stat$timing_11p))
colnames(annot.stat) = c("RNAseq_H_vs_others", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$RNAseq_H_vs_others
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct = F)

# Age adjusted
annot.stat = annot.assoT  %>%
  filter(!is.na(RNAseq_H_vs_others))
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$RNAseq_H_vs_others = factor(annot.stat$RNAseq_H_vs_others, levels = c("H", "other"))
model2 <- glm(RNAseq_H_vs_others ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.7.3 11p alt vs wt  -----------------------------------------------------
annot.test = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt),
         !is.na(RNAseq_H_vs_others)) %>%
  group_by(status_11p_alt_wt, RNAseq_H_vs_others) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt),
         !is.na(RNAseq_H_vs_others))
annot.stat = as.data.frame(table(annot.stat$RNAseq_H_vs_others, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("RNAseq_H_vs_others", "status_11p_alt_wt", "Freq")


annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 
rownames(annot.stat) = annot.stat$RNAseq_H_vs_others

annot.stat = annot.stat %>% 
  select(c("wt", "alt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT %>%
  filter(!is.na(status_11p_alt_wt))
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("wt", "alt"))
annot.stat$RNAseq_H_vs_others = factor(annot.stat$RNAseq_H_vs_others, levels = c("H", "other"))
model2 <- glm(RNAseq_H_vs_others ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.8 NFE2L2------------------------------------------------------------------
    # 5.8.1 Mosaic vs wt -------------------------------------------------------
annot.test = annot.asso %>%
  filter(!is.na(NFE2L2_by_patient)) %>%
  group_by(Mosaic, NFE2L2_by_patient) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot.test

annot.stat = annot.asso %>%
  filter(!is.na(NFE2L2_by_patient)) 

annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$NFE2L2_by_patient = factor(annot.stat$NFE2L2_by_patient, levels = c("missense", "wt"))
model2 <- glm(NFE2L2_by_patient ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(NFE2L2_by_patient)) 
annot.stat = as.data.frame(table(annot.stat$NFE2L2_by_patient, annot.stat$Mosaic))
colnames(annot.stat) = c("NFE2L2_by_patient", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$NFE2L2_by_patient
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct = F)


    # 5.8.2 Compare timing 11p -------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(NFE2L2_by_patient)) %>%
  group_by(timing_11p, NFE2L2_by_patient) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT 
annot.stat = as.data.frame(table(annot.stat$NFE2L2_by_patient, annot.stat$timing_11p))
colnames(annot.stat) = c("NFE2L2_by_patient", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$NFE2L2_by_patient
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$NFE2L2_by_patient = factor(annot.stat$NFE2L2_by_patient, levels = c("missense", "wt"))
model2 <- glm(NFE2L2_by_patient ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.8.3 11palt v wt -------------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(NFE2L2_by_patient),
         !is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, NFE2L2_by_patient) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT 
annot.stat = as.data.frame(table(annot.stat$NFE2L2_by_patient, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("NFE2L2_by_patient", "status_11p_alt_wt", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$NFE2L2_by_patient
annot.stat = annot.stat %>% 
  select(c("alt", "wt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("alt", "wt"))
annot.stat$NFE2L2_by_patient = factor(annot.stat$NFE2L2_by_patient, levels = c("missense", "wt"))
model2 <- glm(NFE2L2_by_patient ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.9 TERT -------------------------------------------------------------------
    # 5.9.1 Mosaic vs wt -------------------------------------------------------
annot.test = annot.asso %>%
  filter(!is.na(TERT_by_patient)) %>%
  group_by(Mosaic, TERT_by_patient) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))
annot.test

annot.stat = annot.asso %>%
  filter(!is.na(TERT_by_patient))
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$TERT_by_patient = factor(annot.stat$TERT_by_patient, levels = c("alt", "wt"))
model2 <- glm(TERT_by_patient ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(TERT_by_patient))
annot.stat = as.data.frame(table(annot.stat$TERT_by_patient, annot.stat$Mosaic))
colnames(annot.stat) = c("TERT_by_patient", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$TERT_by_patient
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct =F)

    # 5.9.2 Compare timing 11p -------------------------------------------------
annot.test = annot.assoT %>%
  filter(!is.na(TERT_by_patient)) %>%
  group_by(timing_11p, TERT_by_patient) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT 
annot.stat = as.data.frame(table(annot.stat$TERT_by_patient, annot.stat$timing_11p))
colnames(annot.stat) = c("TERT_by_patient", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$TERT_by_patient
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$TERT_by_patient = factor(annot.stat$TERT_by_patient, levels = c("alt", "wt"))
model2 <- glm(TERT_by_patient ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.9.3 11palt v wt -------------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(TERT_by_patient),
         !is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, TERT_by_patient) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT 
annot.stat = as.data.frame(table(annot.stat$TERT_by_patient, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("TERT_by_patient", "status_11p_alt_wt", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$TERT_by_patient
annot.stat = annot.stat %>% 
  select(c("alt", "wt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("alt", "wt"))
annot.stat$TERT_by_patient = factor(annot.stat$TERT_by_patient, levels = c("alt", "wt"))
model2 <- glm(TERT_by_patient ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


  # 5.10 PRETEXT stade ---------------------------------------------------------
    # 5.10.1 mos vs non mos ----------------------------------------------------
annot.test = annot.asso %>%
  filter(!is.na(PRETEXT_Stade)) %>%
  group_by(Mosaic, PRETEXT_Stade) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot.test
annot.test$PRETEXT_Stade = factor(annot.test$PRETEXT_Stade, levels = c("Pretext_I", "Pretext_II", "Pretext_III", "Pretext_IV"))

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(PRETEXT_Stade))

annot.stat = as.data.frame(table(annot.stat$Mosaic, annot.stat$PRETEXT_Stade))

colnames(annot.stat) = c("Mosaic", "Type", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mosaic,  
              values_from=Freq) 

annot.stat = as.data.frame(annot.stat)
row.names(annot.stat) = annot.stat$Type
annot.stat = annot.stat %>%
  select(c("wt", "Mosaic"))
annot.stat = t(annot.stat)
Test = prop_trend_test(annot.stat)
Test

    # 5.10.2 Compare timing 11p ------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(PRETEXT_Stade)) %>%
  group_by(timing_11p, PRETEXT_Stade) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(PRETEXT_Stade))
annot.stat = as.data.frame(table(annot.stat$timing_11p, annot.stat$PRETEXT_Stade))

colnames(annot.stat) = c("timing_11p", "Type", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

annot.stat = as.data.frame(annot.stat)
row.names(annot.stat) = annot.stat$Type
annot.stat = annot.stat %>%
  select(c("private", "non_private"))
annot.stat = t(annot.stat)
Test = prop_trend_test(annot.stat)
Test


    # 5.10.3 11p alt vs wt -----------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(PRETEXT_Stade),
         !is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, PRETEXT_Stade) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT %>%
  filter(!is.na(PRETEXT_Stade))
annot.stat = as.data.frame(table(annot.stat$status_11p_alt_wt, annot.stat$PRETEXT_Stade))

colnames(annot.stat) = c("status_11p_alt_wt", "Type", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

annot.stat = as.data.frame(annot.stat)
row.names(annot.stat) = annot.stat$Type
annot.stat = annot.stat %>%
  select(c("alt", "wt"))
annot.stat = t(annot.stat)
Test = prop_trend_test(annot.stat)
Test

  # 5.11 Tumor volume-----------------------------------------------------------
    # 5.11.1 Mosaic vs wt ------------------------------------------------------

annot.test = annot.asso %>%
  filter(!is.na(Tumor.Size..50mm.or...50mm)) %>%
  group_by(Mosaic, Tumor.Size..50mm.or...50mm) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))
annot.test

annot.stat = annot.asso  %>%
  filter(!is.na(Tumor.Size..50mm.or...50mm))
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$Tumor.Size..50mm.or...50mm = factor(annot.stat$Tumor.Size..50mm.or...50mm, levels = c(">50mm", "<=50mm"))
model2 <- glm(Tumor.Size..50mm.or...50mm ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(Tumor.Size..50mm.or...50mm))
annot.stat = as.data.frame(table(annot.stat$Tumor.Size..50mm.or...50mm, annot.stat$Mosaic))
colnames(annot.stat) = c("Size", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

row.names(annot.stat) = annot.stat$size
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct =F)

    # 5.11.2 Compare timing 11p ------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(Tumor.Size..50mm.or...50mm)) %>%
  group_by(timing_11p, Tumor.Size..50mm.or...50mm) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT 
annot.stat = as.data.frame(table(annot.stat$Tumor.Size..50mm.or...50mm, annot.stat$timing_11p))
colnames(annot.stat) = c("Tumor.Size..50mm.or...50mm", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Tumor.Size..50mm.or...50mm
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$Tumor.Size..50mm.or...50mm = factor(annot.stat$Tumor.Size..50mm.or...50mm, levels = c(">50mm", "<=50mm"))
model2 <- glm(Tumor.Size..50mm.or...50mm ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.11.3 11palt v wt -------------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(Tumor.Size..50mm.or...50mm),
         !is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, Tumor.Size..50mm.or...50mm) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.assoT 
annot.stat = as.data.frame(table(annot.stat$Tumor.Size..50mm.or...50mm, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("Tumor.Size..50mm.or...50mm", "status_11p_alt_wt", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Tumor.Size..50mm.or...50mm
annot.stat = annot.stat %>% 
  select(c("alt", "wt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("alt", "wt"))
annot.stat$Tumor.Size..50mm.or...50mm = factor(annot.stat$Tumor.Size..50mm.or...50mm, levels = c(">50mm", "<=50mm"))
model2 <- glm(Tumor.Size..50mm.or...50mm ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.12 >1 Nodule -------------------------------------------------------------
    # 5.12.1 mosaic vs wt ------------------------------------------------------
annot.test = annot.asso %>% 
  filter(!is.na(More_than_1_nodule)) %>% 
  group_by(Mosaic, More_than_1_nodule) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot.test

annot.stat = annot.asso %>%
  filter(!is.na(More_than_1_nodule))
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$More_than_1_nodule = factor(annot.stat$More_than_1_nodule, levels = c("yes", "no"))
model2 <- glm(More_than_1_nodule ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(More_than_1_nodule))
annot.stat = as.data.frame(table(annot.stat$More_than_1_nodule, annot.stat$Mosaic))
colnames(annot.stat) = c("More_than_1_nodule", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$More_than_1_nodule
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct =F)
  
    # 5.12.2 Compare timing 11p ------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(More_than_1_nodule)) %>%
  group_by(timing_11p, More_than_1_nodule) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.asso 
annot.stat = as.data.frame(table(annot.stat$More_than_1_nodule, annot.stat$timing_11p))
colnames(annot.stat) = c("More_than_1_nodule", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$More_than_1_nodule
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$More_than_1_nodule = factor(annot.stat$More_than_1_nodule, levels = c("yes", "no"))
model2 <- glm(More_than_1_nodule ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.12.3 HB alt v wt -------------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(status_11p_by_patient)) %>%
  filter(!is.na(More_than_1_nodule)) %>%
  group_by(status_11p_alt_wt, More_than_1_nodule) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.asso 
annot.stat = as.data.frame(table(annot.stat$More_than_1_nodule, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("More_than_1_nodule", "status_11p_alt_wt", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$More_than_1_nodule
annot.stat = annot.stat %>% 
  select(c("alt", "wt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("alt", "wt"))
annot.stat$More_than_1_nodule = factor(annot.stat$More_than_1_nodule, levels = c("yes", "no"))
model2 <- glm(More_than_1_nodule ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.13 PRETEXT_M ------------------------------------------------------------
    # 5.13.1 mosaic vs non mos -------------------------------------------------
annot.test = annot.asso %>%
  filter(!is.na(PRETEXT_M)) %>%
  group_by(Mosaic, PRETEXT_M) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))

annot.test

annot.stat = annot.asso  %>%
  filter(!is.na(PRETEXT_M)) 
annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$PRETEXT_M = factor(annot.stat$PRETEXT_M, levels = c("no", "yes"))
model2 <- glm(PRETEXT_M ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


#Do statistical test : Chi2
annot.stat = annot.asso %>%
  filter(!is.na(PRETEXT_M))
annot.stat = as.data.frame(table(annot.stat$PRETEXT_M, annot.stat$Mosaic))
colnames(annot.stat) = c("Meta", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Meta
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct = F)

    # 5.13.2 Compare timing 11p ------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(PRETEXT_M)) %>%
  group_by(timing_11p, PRETEXT_M) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.asso 
annot.stat = as.data.frame(table(annot.stat$PRETEXT_M, annot.stat$timing_11p))
colnames(annot.stat) = c("PRETEXT_M", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$PRETEXT_M
annot.stat = annot.stat %>% 
  select(c("private", "non_private"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$PRETEXT_M = factor(annot.stat$PRETEXT_M, levels = c("yes", "no"))
model2 <- glm(PRETEXT_M ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

    # 5.13.3 HB alt v wt -------------------------------------------------------

annot.test = annot.assoT %>%
  filter(!is.na(status_11p_by_patient)) %>%
  filter(!is.na(PRETEXT_M)) %>%
  group_by(status_11p_alt_wt, PRETEXT_M) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test

#Do statistical test : Chi2
annot.stat = annot.asso 
annot.stat = as.data.frame(table(annot.stat$PRETEXT_M, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("PRETEXT_M", "status_11p_alt_wt", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$PRETEXT_M
annot.stat = annot.stat %>% 
  select(c("alt", "wt"))
chisq.test(annot.stat, correct =F)

# Age adjusted
annot.stat = annot.assoT 
annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("alt", "wt"))
annot.stat$PRETEXT_M = factor(annot.stat$PRETEXT_M, levels = c("yes", "no"))
model2 <- glm(PRETEXT_M ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)

  # 5.14 Histological heterogeneity---------------------------------------------
    # 5.14.1 mosaic vs non mos  ------------------------------------------------
# TUMOR PLASTICITY
annot.test = annot.het %>%
  filter(!is.na(Tumor_plasticity)) %>%
  group_by(Mosaic, Tumor_plasticity) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(Mosaic) %>%
  summarise(nb = sum(count))
annot.test


annot.stat = annot.het  %>%
  filter(!is.na(Tumor_plasticity)) 

annot.stat$Mosaic = factor(annot.stat$Mosaic, levels = c("Mosaic", "wt"))
annot.stat$Tumor_plasticity = factor(annot.stat$Tumor_plasticity, levels = c("homo", "hete"))
model2 <- glm(Tumor_plasticity ~ Mosaic + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


#Do statistical test : Chi2
annot.stat = annot.het %>%
  filter(!is.na(Tumor_plasticity)) 

annot.stat = as.data.frame(table(annot.stat$Tumor_plasticity, annot.stat$Mosaic))
colnames(annot.stat) = c("Tumor_plasticity", "Mos", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = Mos,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Tumor_plasticity
annot.stat = annot.stat %>% 
  select(c("wt", "Mosaic"))
chisq.test(annot.stat, correct = F)

    # 5.14.2 timing 11p  -------------------------------------------------------
annot.test = annot.het %>%
  filter(!is.na(Tumor_plasticity)) %>%
  group_by(timing_11p, Tumor_plasticity) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(timing_11p) %>%
  summarise(nb = sum(count))
annot.test


annot.stat = annot.het  %>%
  filter(!is.na(Tumor_plasticity)) 

annot.stat$timing_11p = factor(annot.stat$timing_11p, levels = c("private", "non_private"))
annot.stat$Tumor_plasticity = factor(annot.stat$Tumor_plasticity, levels = c("homo", "hete"))
model2 <- glm(Tumor_plasticity ~ timing_11p + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


#Do statistical test : Chi2
annot.stat = annot.het %>%
  filter(!is.na(Tumor_plasticity)) 

annot.stat = as.data.frame(table(annot.stat$Tumor_plasticity, annot.stat$timing_11p))
colnames(annot.stat) = c("Tumor_plasticity", "timing_11p", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = timing_11p,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Tumor_plasticity
annot.stat = annot.stat %>% 
  select(c("non_private", "private"))
chisq.test(annot.stat, correct = F)


    # 5.14.3 11p alt vs wt------------------------------------------------------

annot.test = annot.het %>%
  filter(!is.na(Tumor_plasticity),
         !is.na(status_11p_alt_wt)) %>%
  group_by(status_11p_alt_wt, Tumor_plasticity) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) 

annot.test %>%
  group_by(status_11p_alt_wt) %>%
  summarise(nb = sum(count))
annot.test


annot.stat = annot.het  %>%
  filter(!is.na(Tumor_plasticity)) 

annot.stat$status_11p_alt_wt = factor(annot.stat$status_11p_alt_wt, levels = c("alt", "wt"))
annot.stat$Tumor_plasticity = factor(annot.stat$Tumor_plasticity, levels = c("homo", "hete"))
model2 <- glm(Tumor_plasticity ~ status_11p_alt_wt + Age.surg.years, data = annot.stat, family = "binomial")
summary(model2)


#Do statistical test : Chi2
annot.stat = annot.het %>%
  filter(!is.na(Tumor_plasticity),
         !is.na(status_11p_alt_wt))

annot.stat = as.data.frame(table(annot.stat$Tumor_plasticity, annot.stat$status_11p_alt_wt))
colnames(annot.stat) = c("Tumor_plasticity", "status_11p_alt_wt", "Freq")

annot.stat = annot.stat %>%
  pivot_wider(names_from = status_11p_alt_wt,  
              values_from=Freq) 

rownames(annot.stat) = annot.stat$Tumor_plasticity
annot.stat = annot.stat %>% 
  select(c("alt", "wt"))
chisq.test(annot.stat, correct = F)

### 6. Survival curves mosaic vs non mosaic ------------------------------------
  # 6.1 prepare dataframe-------------------------------------------------------
stop_survival_data <- 60 # choose how many months you want to represent

annot.asso$TTD = as.numeric(annot.asso$TTD)
annot.asso$TTP = as.numeric(annot.asso$TTP)
annot.asso$DSS.status = as.numeric(annot.asso$DSS.status)
annot.asso$PFS.status = as.numeric(annot.asso$PFS.status)

annot_surv = annot.asso %>% 
  filter(Age.surg.years<=3.3) %>%
  filter(!is.na(DSS.status),
         !is.na(PFS.status),
         !is.na(TTD)) %>%
  mutate(cnLOH_vs_no_cnLOH = ifelse(status_11p_by_patient=="cn-LOH", "cn-LOH", "no_cn-LOH"),
         OSS.status_truncated = ifelse(TTD > stop_survival_data,
                                       0,
                                       OSS.status),
         DSS.status_truncated = ifelse(TTD > stop_survival_data,
                                       0,
                                       DSS.status),
         TTD_truncated = ifelse(TTD > stop_survival_data,
                                stop_survival_data,
                                TTD),
         PFS.status_truncated = ifelse(TTP > stop_survival_data,
                                       0,
                                       PFS.status),
         TTP_truncated = ifelse(TTP > stop_survival_data,
                                stop_survival_data,
                                TTP)) %>% 
  as.data.frame()

# Export Source data files 
Source_figure_6d <- annot_surv %>%
  select(c("CHCID", 
           "Age.surg.years",
           "Mosaic",
           "cnLOH_vs_no_cnLOH",
           "DSS.status_truncated",
           "PFS.status_truncated",
           "TTP_truncated"))

write.table(Source_figure_6d, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Items/Source_figure_6d.txt", sep="\t")

table(annot_surv$DSS.status_truncated)
table(annot_surv$PFS.status_truncated)
dim(annot_surv)
GivePatientsSamples <- function(annot) {
  print(glue::glue("{length(annot_surv$CHCID)} samples from {length(unique(annot_surv$Patient.identification))} patients"))
}

  # 6.2 Create survival objects ------------------------------------------------

annot_surv$DSS_object = with(annot_surv, 
                             Surv(TTD_truncated, DSS.status_truncated == 1))

annot_surv$OS_object = with(annot_surv, 
                            Surv(TTD_truncated, OSS.status_truncated == 1))

annot_surv$PFS_object = with(annot_surv, 
                             Surv(TTP_truncated, PFS.status_truncated == 1))

  # 6.3 Spot significant association Mosaic-------------------------------------
surv_table <- data.frame(categ = "test", DSS = 1, PFS = 1, effectif_size = "100/100", big_enough = T) %>% factoall()

surv_test_DSS <- survdiff(DSS_object ~ annot_surv[, "Mosaic"], data = annot_surv, rho = 0) # rho = 0 for logrank test
log_rank_pval_DSS <- 1 - pchisq(surv_test_DSS$chisq, length(surv_test_DSS$n) - 1) 
print(glue::glue("DSS p-value = {log_rank_pval_DSS}"))

surv_test_PFS <- survdiff(PFS_object ~ annot_surv[, "Mosaic"], data = annot_surv, rho = 0) # rho = 0 for logrank test
log_rank_pval_PFS <- 1 - pchisq(surv_test_PFS$chisq, length(surv_test_PFS$n) - 1)
print(glue::glue("PFS p-value = {log_rank_pval_PFS}"))

effectif_s <- table(annot_surv[, "Mosaic"])

surv_table <- rbind(surv_table, 
                    c("Mosaic", log_rank_pval_DSS, log_rank_pval_PFS,
                      paste(effectif_s, collapse = "/"), sum(effectif_s > 3) > 1))

surv_table <- surv_table %>% 
  rowwise() %>% 
  mutate(DSS = as.numeric(DSS),
         PFS = as.numeric(PFS),
         big_enough = big_enough == "TRUE",
         min_p_val = min(DSS, PFS)) %>% 
  arrange(min_p_val) %>% 
  as.data.frame()



  # 6.4 Plot survival vs Mosaic status------------------------------------------

annot_var = annot_surv %T>% GivePatientsSamples() %>% 
  .[!is.na(annot_surv[, "Mosaic"]), ] %T>% GivePatientsSamples() # remove samples where the variable is unknown or not part of the groups you're interested in

DSS_my_categ = npsurv(DSS_object ~ annot_var[, "Mosaic"], data = annot_var)
PFS_my_categ = npsurv(PFS_object ~ annot_var[, "Mosaic"], data = annot_var)

gg_DSS <- survminer::ggsurvplot(DSS_my_categ, risk.table = T, pval = T, pval.method = F, conf.int =F, risk.table.height=0.3, 
                                 ylab = "Disease-specific survival", break.time.by = 12,
                                palette =  c("indianred", "black"),
                                risk.table.fontsize = 5,
                                legend.labs = c("Mosaic", "Non_Mosaic"), size = 1.5,  risk.table.y.text = FALSE) 

gg_PFS<- survminer::ggsurvplot(PFS_my_categ, risk.table = T, pval = T, pval.method = F, conf.int =F, risk.table.height=0.3,
                               risk.table.fontsize = 5,
                               ylab = "", break.time.by = 12, 
                               palette = c("indianred", "black"),
                               legend.labs = c("Mosaic", "Non_Mosaic"), size = 2, risk.table.y.text = FALSE)

print(gg_DSS)
print(gg_PFS)

### 7. Survival curves cn-LOH vs non cn-LOH ------------------------------------
  # 7.1 prepare dataframe-------------------------------------------------------
stop_survival_data <- 60 # choose how many months you want to represent

annot.asso$TTD = as.numeric(annot.asso$TTD)
annot.asso$TTP = as.numeric(annot.asso$TTP)
annot.asso$DSS.status = as.numeric(annot.asso$DSS.status)
annot.asso$PFS.status = as.numeric(annot.asso$PFS.status)

annot_surv = annot.asso %>% 
  filter(Age.surg.years<=3.3) %>%   # filter or not for all ages
  filter(!is.na(DSS.status),
         !is.na(PFS.status),
         !is.na(TTD),
         !is.na(status_11p_by_patient)) %>%
  mutate(cnLOH_vs_no_cnLOH = ifelse(status_11p_by_patient=="cn-LOH", "cn-LOH", "no_cn-LOH"),
          OSS.status_truncated = ifelse(TTD > stop_survival_data,
                                       0,
                                       OSS.status),
         DSS.status_truncated = ifelse(TTD > stop_survival_data,
                                       0,
                                       DSS.status),
         TTD_truncated = ifelse(TTD > stop_survival_data,
                                stop_survival_data,
                                TTD),
         PFS.status_truncated = ifelse(TTP > stop_survival_data,
                                       0,
                                       PFS.status),
         TTP_truncated = ifelse(TTP > stop_survival_data,
                                stop_survival_data,
                                TTP)) %>% 
  as.data.frame()



table(annot_surv$DSS.status_truncated)
table(annot_surv$PFS.status_truncated)
table(annot_surv$cnLOH_vs_no_cnLOH)
dim(annot_surv)
GivePatientsSamples <- function(annot) {
  print(glue::glue("{length(annot_surv$CHCID)} samples from {length(unique(annot_surv$Patient.identification))} patients"))
}

  # 7.2 Create survival objects ------------------------------------------------

annot_surv$DSS_object = with(annot_surv, 
                             Surv(TTD_truncated, DSS.status_truncated == 1))

annot_surv$OS_object = with(annot_surv, 
                            Surv(TTD_truncated, OSS.status_truncated == 1))

annot_surv$PFS_object = with(annot_surv, 
                             Surv(TTP_truncated, PFS.status_truncated == 1))

  # 7.3 Spot significant association cn-LOH -------------------------------------
surv_table <- data.frame(categ = "test", DSS = 1, PFS = 1, effectif_size = "100/100", big_enough = T) %>% factoall()

surv_test_DSS <- survdiff(DSS_object ~ annot_surv[, "cnLOH_vs_no_cnLOH"], data = annot_surv, rho = 0) # rho = 0 for logrank test
log_rank_pval_DSS <- 1 - pchisq(surv_test_DSS$chisq, length(surv_test_DSS$n) - 1) 
print(glue::glue("DSS p-value = {log_rank_pval_DSS}"))

surv_test_PFS <- survdiff(PFS_object ~ annot_surv[, "cnLOH_vs_no_cnLOH"], data = annot_surv, rho = 0) # rho = 0 for logrank test
log_rank_pval_PFS <- 1 - pchisq(surv_test_PFS$chisq, length(surv_test_PFS$n) - 1)
print(glue::glue("PFS p-value = {log_rank_pval_PFS}"))

effectif_s <- table(annot_surv[, "cnLOH_vs_no_cnLOH"])

surv_table <- rbind(surv_table, 
                    c("cnLOH_vs_no_cnLOH", log_rank_pval_DSS, log_rank_pval_PFS,
                      paste(effectif_s, collapse = "/"), sum(effectif_s > 3) > 1))

surv_table <- surv_table %>% 
  rowwise() %>% 
  mutate(DSS = as.numeric(DSS),
         PFS = as.numeric(PFS),
         big_enough = big_enough == "TRUE",
         min_p_val = min(DSS, PFS)) %>% 
  arrange(min_p_val) %>% 
  as.data.frame()



  # 7.4 Plot survival vs cn-LOH status------------------------------------------

annot_var = annot_surv %T>% GivePatientsSamples() %>% 
  .[!is.na(annot_surv[, "cnLOH_vs_no_cnLOH"]), ] %T>% GivePatientsSamples() # remove samples where the variable is unknown or not part of the groups you're interested in

DSS_my_categ = npsurv(DSS_object ~ annot_var[, "cnLOH_vs_no_cnLOH"], data = annot_var)
PFS_my_categ = npsurv(PFS_object ~ annot_var[, "cnLOH_vs_no_cnLOH"], data = annot_var)

gg_DSS <- survminer::ggsurvplot(DSS_my_categ, risk.table = T, pval = T, pval.method = F, conf.int =F, risk.table.height=0.3, 
                                ylab = "Disease-specific survival", break.time.by = 12,
                                palette =  c("indianred", "black"),
                                risk.table.fontsize = 5,
                                legend.labs = c("cn-LOH", "no cn-LOH"), size = 1.5,  risk.table.y.text = FALSE) 

gg_PFS<- survminer::ggsurvplot(PFS_my_categ, risk.table = T, pval = T, pval.method = F, conf.int =F, risk.table.height=0.3,
                               risk.table.fontsize = 5,
                               ylab = "", break.time.by = 12, 
                               palette = c("indianred", "black"),
                               legend.labs = c("cn-LOH", "no cn-LOH"), size = 2, risk.table.y.text = FALSE)

print(gg_DSS)
print(gg_PFS)

