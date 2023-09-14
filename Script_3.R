library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(readxl)

# Load table  ------------------------------------------------------------------
ana.table <- read.delim("D:/Dropbox/11p15.5 mosaicism/Anapath/for_R/IGF2_expression_relecture_for_R.txt", sep="\t", as.is=T, na.strings=c("",".","NA","na","#N/D","<NA>"))

# Create a new table with one line per patient for non tumor analysis
ana.unique = ana.table %>%
  filter(One_N_per_patient =="yes") 

# recode to max value
ana.unique[ana.unique == "0/x"] <- "x"
ana.unique[ana.unique == "x/xx"] <- "xx"
ana.unique[ana.unique == "xx/xxx"] <- "xxx"

ana.unique[ana.unique == "i0 ? i1"] <- "i1"
ana.unique[ana.unique == "i1 ? i2"] <- "i2"
ana.unique[ana.unique == "i2 ? i3"] <- "i3"

### 1. Non tumor analysis stats for mosaic vs non mosaic -----------------------
  # 1.1 Pattern ----------------------------------------------------------------
ana.test = as.data.frame(table(ana.unique$Histological.Diagnosis, ana.unique$Pattern_Zone_3_IGF2))

colnames(ana.test) = c("Histological.Diagnosis", "Type", "Freq")

# Do statistical test : Chi2
ana.stat = ana.test %>%
  pivot_wider(names_from = Histological.Diagnosis,  
  values_from=Freq) 

ana.stat = as.data.frame(ana.stat)
row.names(ana.stat) = ana.stat$Type
ana.stat = ana.stat %>%
  select(c("HB", "Mosaic"))
ana.stat = t(ana.stat)
Test = prop_trend_test(ana.stat)
Test

  # 1.2 Intensity ------------------------------------------------------------

ana.test = as.data.frame(table(ana.unique$Histological.Diagnosis, ana.unique$Intensity_Zone_3_IGF2))

colnames(ana.test) = c("Histological.Diagnosis", "Type", "Freq")

#Do statistical test : Chi2
ana.stat = ana.test %>%
  pivot_wider(names_from = Histological.Diagnosis,  
              values_from=Freq) 

ana.stat = as.data.frame(ana.stat)
row.names(ana.stat) = ana.stat$Type
ana.stat = ana.stat %>%
  select(c("HB", "Mosaic"))
ana.stat = t(ana.stat)
Test = prop_trend_test(ana.stat)
Test

### 2. Pattern graphs + stats --------------------------------------------------
  # 2.1 Non mosaic--------------------------------------------------------------
cols_pattern = c("0" = "white",  "x" = "#EF9191", "xx" = "#B71C1C", "xxx" = "#420000")

# Filter only non tumor non mosaic samples
ana.test = ana.unique %>%
  filter(Histological.Diagnosis=="HB") %>%
  pivot_longer(names_to = "Zone", cols = c("Pattern_Zone_1_IGF2", "Pattern_Zone_2_IGF2", "Pattern_Zone_3_IGF2"), values_to = "Pattern")

Source_figure_3a <- ana.unique %>%
  select(c("Histological.Diagnosis", 
           "Region", 
           "Patient.ID", 
           "Bloc", 
           "Pattern_Zone_1_IGF2", 
           "Pattern_Zone_2_IGF2", 
           "Pattern_Zone_3_IGF2",
           "Intensity_Zone_1_IGF2",
           "Intensity_Zone_2_IGF2",
           "Intensity_Zone_3_IGF2"))
write.table(Source_figure_3a, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_3a.txt", sep="\t")

ana.test = as.data.frame(table(ana.test$Zone, ana.test$Pattern))
colnames(ana.test) = c("Zone", "Type", "Freq")

#Do statistical test : Chi2
ana.stat = ana.test %>%
  pivot_wider(names_from = Zone,  
              values_from=Freq) 

ana.stat = as.data.frame(ana.stat)
row.names(ana.stat) = ana.stat$Type
ana.stat = ana.stat %>%
  select(c("Pattern_Zone_1_IGF2", "Pattern_Zone_2_IGF2", "Pattern_Zone_3_IGF2"))
ana.stat = t(ana.stat)
Test = prop_trend_test(ana.stat)
Test

#Prepare representation
ana.test = ana.test %>%
  group_by(Zone) %>%
  mutate(Freq = Freq/sum(Freq)*100)
ana.test$Zone = factor(ana.test$Zone, levels = c("Pattern_Zone_1_IGF2", "Pattern_Zone_2_IGF2", "Pattern_Zone_3_IGF2"))
ana.test$Type = factor(ana.test$Type, levels = c("xxx",  "xx",  "x", "0"))

#Graph
ggplot(ana.test, aes(x=Zone, y=Freq,fill = Type)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  scale_fill_manual(values = cols_pattern) +
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
        axis.text.x=element_blank()
        # axis.text.x = element_text(angle=90) 
  )

  # 2.2 Mosaic samples ---------------------------------------------------------
cols_pattern = c("0" = "white",  "x" = "#EF9191", "xx" = "#B71C1C", "xxx" = "#420000")

# Filter only non tumor non mosaic samples
ana.test = ana.unique %>%
  filter(Histological.Diagnosis=="Mosaic") %>%
  pivot_longer(names_to = "Zone", cols = c("Pattern_Zone_1_IGF2", "Pattern_Zone_2_IGF2", "Pattern_Zone_3_IGF2"), values_to = "Pattern")

ana.test = as.data.frame(table(ana.test$Zone, ana.test$Pattern))
colnames(ana.test) = c("Zone", "Type", "Freq")

#Do statistical test : Chi2
ana.stat = ana.test %>%
  pivot_wider(names_from = Zone,  
              values_from=Freq) 

ana.stat = as.data.frame(ana.stat)
row.names(ana.stat) = ana.stat$Type
ana.stat = ana.stat %>%
  select(c("Pattern_Zone_1_IGF2", "Pattern_Zone_2_IGF2", "Pattern_Zone_3_IGF2"))
ana.stat = t(ana.stat)
Test = prop_trend_test(ana.stat)
Test

#Prepare representation
ana.test = ana.test %>%
  group_by(Zone) %>%
  mutate(Freq = Freq/sum(Freq)*100)
ana.test$Zone = factor(ana.test$Zone, levels = c("Pattern_Zone_1_IGF2", "Pattern_Zone_2_IGF2", "Pattern_Zone_3_IGF2"))
ana.test$Type = factor(ana.test$Type, levels = c("xxx",  "xx",  "x", "0"))

#Graph
ggplot(ana.test, aes(x=Zone, y=Freq,fill = Type)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  scale_fill_manual(values = cols_pattern) +
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
        axis.text.x=element_blank()
        # axis.text.x = element_text(angle=90) 
  )

### 3. Intensity graphs + stats --------------------------------------------------
  # 2.1 Non mosaic--------------------------------------------------------------
cols_intensity = c( "i0" = "white", "i1" = "#C5CAE9", "i2"= "#303F9F", "i3" = "#00085c", "i0 ? i3" = "#7aba7a", "i1 ? i3" = "#228c22")

# Filter only non tumor non mosaic samples
ana.test = ana.unique %>%
  filter(Histological.Diagnosis=="HB") %>%
  pivot_longer(names_to = "Zone", cols = c("Intensity_Zone_1_IGF2", "Intensity_Zone_2_IGF2", "Intensity_Zone_3_IGF2"), values_to = "Intensity")


ana.test = as.data.frame(table(ana.test$Zone, ana.test$Intensity))
colnames(ana.test) = c("Zone", "Type", "Freq")

#Do statistical test : Chi2
ana.stat = ana.test %>%
  pivot_wider(names_from = Zone,  
              values_from=Freq) 

ana.stat = as.data.frame(ana.stat)
row.names(ana.stat) = ana.stat$Type
ana.stat = ana.stat %>%
  select(c("Intensity_Zone_1_IGF2", "Intensity_Zone_2_IGF2", "Intensity_Zone_3_IGF2"))
ana.stat = t(ana.stat)
Test = prop_trend_test(ana.stat)
Test

#Prepare representation
ana.test = ana.test %>%
  group_by(Zone) %>%
  mutate(Freq = Freq/sum(Freq)*100)
ana.test$Zone = factor(ana.test$Zone, levels = c("Intensity_Zone_1_IGF2", "Intensity_Zone_2_IGF2", "Intensity_Zone_3_IGF2"))
ana.test$Type = factor(ana.test$Type, levels = c("i3",  "i2",  "i1", "i0", "i1 ? i3", "i0 ? i3"))

#Graph
ggplot(ana.test, aes(x=Zone, y=Freq, fill = Type)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  scale_fill_manual(values = cols_intensity) +
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
        axis.text.x=element_blank()
        # axis.text.x = element_text(angle=90) 
  )

  # 2.2 Mosaic samples ---------------------------------------------------------
cols_intensity = c( "i0" = "white", "i1" = "#C5CAE9", "i2"= "#303F9F", "i3" = "#00085c", "i0 ? i3" = "#7aba7a", "i1 ? i3" = "#228c22")

# Filter only non tumor non mosaic samples
ana.test = ana.unique %>%
  filter(Histological.Diagnosis=="Mosaic") %>%
  pivot_longer(names_to = "Zone", cols = c("Intensity_Zone_1_IGF2", "Intensity_Zone_2_IGF2", "Intensity_Zone_3_IGF2"), values_to = "Intensity")


ana.test = as.data.frame(table(ana.test$Zone, ana.test$Intensity))
colnames(ana.test) = c("Zone", "Type", "Freq")

#Do statistical test : Chi2
ana.stat = ana.test %>%
  pivot_wider(names_from = Zone,  
              values_from=Freq) 

ana.stat = as.data.frame(ana.stat)
row.names(ana.stat) = ana.stat$Type
ana.stat = ana.stat %>%
  select(c("Intensity_Zone_1_IGF2", "Intensity_Zone_2_IGF2", "Intensity_Zone_3_IGF2"))
ana.stat = t(ana.stat)
Test = prop_trend_test(ana.stat)
Test

#Prepare representation
ana.test = ana.test %>%
  group_by(Zone) %>%
  mutate(Freq = Freq/sum(Freq)*100)
ana.test$Zone = factor(ana.test$Zone, levels = c("Intensity_Zone_1_IGF2", "Intensity_Zone_2_IGF2", "Intensity_Zone_3_IGF2"))
ana.test$Type = factor(ana.test$Type, levels = c("i3",  "i2",  "i1", "i0", "i1 ? i3", "i0 ? i3"))

#Graph
ggplot(ana.test, aes(x=Zone, y=Freq, fill = Type)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  scale_fill_manual(values = cols_intensity) +
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
        axis.text.x=element_blank()
        # axis.text.x = element_text(angle=90) 
  )

### 3. Veins and GS quantification Barkha --------------------------------------
  # 3.1 Load quantification data -----------------------------------------------

data.dir = "D:/Dropbox/11p15.5 mosaicism/Anapath/for_R"
veins_count = read_xlsx(file.path(data.dir, "veins_count.xlsx"))

veins_count = veins_count %>%
  mutate(group = case_when(grepl("Mos", Relecture_rank) &  Region == "normal" ~ "Mosaic_non_patch",
                           Region =="patch_IGF2" ~ "Mosaic_patch",
                           grepl("Normal", Relecture_rank) ~ "normal"),
         veins_by_area = Confirmed_veins / Area_analysed_QuPath,
         Doubtful_veins_by_area = Doubtful_veins / Area_analysed_QuPath)

Source_figure_3b <- veins_count %>%
  filter(group %in% c("Mosaic_patch", "normal")) %>%
  select(c("Region", 
           "Patient.ID", 
           "Bloc", 
           "group", 
           "veins_by_area",
           "Intensity_GS_staining_majority"))

write.table(Source_figure_3b, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_3b.txt", sep="\t")


  # 3.2 Intensity of GS staining in mosaic and normal --------------------------

table(veins_count$Intensity_GS_staining_majority)

veins_test = as.data.frame(table(veins_count$group, veins_count$Intensity_GS_staining_majority))
colnames(veins_test) = c("group", "Type", "Freq")
veins_test = veins_test %>%
  arrange(factor(Type, levels = c("absent", "weak", "medium", "strong")))

#Do statistical test : Chi2
veins_stat = veins_test %>%
  pivot_wider(names_from = group,  
              values_from=Freq) 

veins_stat = as.data.frame(veins_stat)
rownames(veins_stat) = veins_stat$Type
veins_stat = veins_stat %>%
  select(c("Mosaic_patch", "normal"))
veins_stat = t(veins_stat)
Test = prop_trend_test(veins_stat)
Test

#Prepare representation
veins_test = veins_test %>%
  group_by(group) %>%
  mutate(Freq = Freq/sum(Freq)*100) %>%
  filter(group %in% c("Mosaic_patch", "normal")) # Represent only mosaic area and normal liver

veins_test$group = factor(veins_test$group, levels = c("Mosaic_patch", "Mosaic_non_patch", "normal"))
veins_test$Type = factor(veins_test$Type, levels = c("strong","medium", "weak","absent"))

cols_intensity = c( "absent" = "white", "weak" = "#C5CAE9", "medium"= "#303F9F", "strong" = "#00085c")

#Graph
ggplot(veins_test, aes(x=group, y=Freq, fill = Type)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  scale_fill_manual(values = cols_intensity) +
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
       axis.text.x=element_blank())


  # 3.3 Number of veins  -------------------------------------------------------

veins_count$group = factor(veins_count$group, levels = c("Mosaic_patch", "Mosaic_non_patch", "normal"))
cols <- c("Mosaic_patch" = "indianred", "Mosaic_non_patch" = "darkblue", "normal" = "blue")

# All veins --------------------------------------------------------------------

veins_count = veins_count %>%
  mutate(all_veins = (Confirmed_veins + Doubtful_veins)/Area_analysed_QuPath) %>%
  filter(group %in% c("Mosaic_patch", "normal")) # Represent only mosaic area and normal liver


ggplot(veins_count, aes(x=group, y=all_veins, fill = group, color = group)) + 
  geom_jitter(width = 0.1, size=2)+
  scale_color_manual(values =cols) +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5) +
  theme_bw() +
  theme(legend.key = element_blank(),
        legend.background = element_rect(colour = 'black'),
        legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0.1, "cm"),
        panel.background = element_blank(),
            axis.text = element_blank()
  )

veins.stat = veins_count %>%
  filter(group %in% c("Mosaic_patch", "normal"))
wilcox.test(veins.stat$all_veins ~ veins.stat$group)

