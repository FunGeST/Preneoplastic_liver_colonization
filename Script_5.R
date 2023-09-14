# Load installed packages-------------------------------------------------------
library(ggplot2)
library(readxl)
library(dplyr)
library(ggcorrplot)
library(geco.utils)

# Make corr plot multi samples -------------------------------------------------
#Directory-------------------------------------------------------------------
data.dir <- "D:/Dropbox/11p15.5 mosaicism/Summary_multi_techniques"
Multi_samp = read_xlsx(file.path(data.dir, "Summary_multi_techniques.xlsx"))

Source_figure_2b <- Multi_samp 
write.table(Source_figure_2b, "D:/Dropbox/11p15.5 mosaicism/MANUSCRIPT/Nat_com_revisions/Round_2/Source_files/Source_figure_2b.txt", sep="\t")


Multi_samp_plot = Multi_samp %>%
  select(c("Sample_ID", "WGS_WES", "IC1IC2_percentage", "RNAscope_IGF2H19"))

names = Multi_samp_plot$Sample_ID


Multi_samp_plot = Multi_samp_plot %>%
  select(-Sample_ID)

rownames(Multi_samp_plot) = names

Multi_samp_plot = as.data.frame(Multi_samp_plot)
Multi_samp_plot$WGS_WES = as.numeric(Multi_samp_plot$WGS_WES)
Multi_samp_plot$IC1IC2_percentage = as.numeric(Multi_samp_plot$IC1IC2_percentage)
Multi_samp_plot$RNAscope_IGF2H19 = as.numeric(Multi_samp_plot$RNAscope_IGF2H19)


Multi_samp_plot = Multi_samp_plot %>%
  mutate(WGS_WES = WGS_WES/100,
       IC1IC2_percentage = IC1IC2_percentage/100,
       RNAscope_IGF2H19 = RNAscope_IGF2H19/100)


corr.mos = ggcorrplot(Multi_samp_plot, 
           type = "full", 
           outline.color = "black",
           lab = FALSE,  
           show.legend = T,
           ggtheme=theme_void) +
  theme(axis.ticks.x = element_blank()) 

corr.mos + scale_fill_gradient2(limit=c(0,1), low="blue", high="red", mid = "yellow", midpoint = 0.5)

