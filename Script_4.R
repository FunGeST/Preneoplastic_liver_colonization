
# Load packages of interest ----------------------------------------------------

library(ggplot2)
library(dplyr)
library(viridis)
library(geco.utils)
library(tidyr)

# Load methylation IC1 x IC2 ---------------------------------------------------

IC1 <- geco.load("D:/Dropbox/11p15.5 mosaicism/Fetal_livers_Bonder/meth_fetal_IC1.RData")
IC2 <- geco.load("D:/Dropbox/11p15.5 mosaicism/Fetal_livers_Bonder/meth_fetal_IC2.RData")

Source_supp_figure_6 <- IC1
write.table(Source_supp_figure_6, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_6.txt", sep="\t")

Source_supp_figure_6bis <- IC2
write.table(Source_supp_figure_6bis, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_supp_figure_6bis.txt", sep="\t")


dim(IC1)
dim(IC2)

rownames(IC1) = IC1$ID_REF
rownames(IC2) = IC2$ID_REF

samp_of_interest = colnames(IC1)[grepl("SAMPLE", colnames(IC1))]

Mean_IC1 = data.frame(colMeans(IC1[,samp_of_interest]))
Mean_IC1$Sample = rownames(Mean_IC1)
colnames(Mean_IC1) = c("Mean_IC1", "Sample")

Mean_IC2 = data.frame(colMeans(IC2[,samp_of_interest]))
Mean_IC2$Sample = rownames(Mean_IC2)
colnames(Mean_IC2) = c("Mean_IC2", "Sample")

IC1 = IC1 %>%
  select(c("ID_REF", "pos",colnames(IC1)[grepl("SAMPLE", colnames(IC1))])) %>%
  pivot_longer(cols= samp_of_interest, names_to = "Sample", values_to = "Meth_value")

IC1 = left_join(IC1, Mean_IC1, by="Sample")

IC2 = IC2 %>%
  select(c("ID_REF","pos", colnames(IC2)[grepl("SAMPLE", colnames(IC2))])) %>%
  pivot_longer(cols= colnames(IC2)[grepl("SAMPLE", colnames(IC2))], names_to = "Sample", values_to = "Meth_value")

IC2 = left_join(IC2, Mean_IC2, by="Sample")

IC1 = IC1 %>%
  arrange(Mean_IC1)

IC1$Sample=factor(IC1$Sample, levels=unique(IC1$Sample))

cols = viridis(n=14, option="viridis", direction = -1, alpha=0.7)

# Plot IC1 and IC2 methylation -------------------------------------------------

IC1$Meth_value = as.numeric(IC1$Meth_value)
IC2$Meth_value = as.numeric(IC2$Meth_value)

ggplot(IC1, aes(x = Sample, y=Meth_value, fill=Sample)) + 
  geom_boxplot(aes(fill=Sample), outlier.color = "gray80", outlier.alpha = 0.5, color="black")+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(x = "", y="IC1 Methylaion value") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(size=8, angle=90,hjust = 1),
        axis.line = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(0.2, "cm")) 

samp.order <- unique(IC1$Sample)

IC2$Sample=factor(IC2$Sample, levels=samp.order)


ggplot(IC2, aes(x = Sample, y=Meth_value, fill=Sample)) + 
  geom_boxplot(aes(fill=Sample), outlier.color = "gray80", outlier.alpha = 0.5, color="black")+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(x = "", y="IC2 Methylaion value") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(size=8, angle=90,hjust = 1),
        axis.line = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(0.2, "cm")) 

# Calculate mosaic per sample percentage ---------------------------------------

merged.tab <- left_join(IC1, IC2, by="Sample")
merged.tab = merged.tab %>%
  filter(!duplicated(Sample)) %>%
  mutate(percent.mos = 100*(Mean_IC1-Mean_IC2)/(Mean_IC1 + Mean_IC2),
         MethRatio = Mean_IC1/Mean_IC2)

dim(merged.tab)

tab.to.print = merged.tab %>%
  select(c(Sample, Mean_IC1, Mean_IC2))
colnames(tab.to.print) = c("CHCID", "IC1_bval_merged", "IC2_bval_merged")

#write.table(tab.to.print, file="D:/Dropbox/11p15.5 mosaicism/Fetal_livers_Bonder/Fetal_livers_Bonder.txt", sep="\t")

merged.tab$Sample=factor(merged.tab$Sample, levels=samp.order)


ggplot(merged.tab, aes(x = Sample, y=percent.mos, fill=Sample, color = Sample)) + 
  geom_point(aes(fill=Sample))+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(-10,10)) +
  geom_hline(yintercept = 2.5, linetype=2, color="gray70")+
  geom_hline(yintercept = -2.5, linetype=2, color="gray70") +
  labs(x = "", y="Mosaic cells (%)") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(size=8, angle=90,hjust = 1),
        axis.line = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(0.2, "cm")) 

ggplot(merged.tab, aes(x = Sample, y=MethRatio, fill=Sample, color = Sample)) + 
  geom_point(aes(fill=Sample), size=3)+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(0,4)) +
  labs(x = "", y="IC1/IC2 methylation ratio") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(size=8, angle=90,hjust = 1),
        axis.line = element_line(colour = "black", linewidth = 1),
        axis.ticks.length = unit(0.2, "cm")) 

nb.mosaic <- data.frame(geco.load("D:/Dropbox/11p15.5 mosaicism/Fetal_livers_Bonder/nb_mosaic_cells.RData"))
