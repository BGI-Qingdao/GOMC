rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 8g. Variance of biosynthetic diversity for genomes at different taxonomic rank from phylum to genus. 
# =========================================================================================================================

# Load libraries
library(ggplot2)
library(scales)
data <- read.table("data/ex-figure8g.Variance_Info.xls",header = T,sep = "\t")
ggplot(data,aes(Level,variance))+
  geom_boxplot(position = position_nudge(x=-0.15),width=0.3,outlier.size = 1.5,outlier.shape = 18)+
  stat_summary(fun="mean", geom="point", shape=20, size=2.5, color="red", fill="red",alpha=0.7,position = position_nudge(x=-0.15))+
  geom_point(position = position_nudge(x=0.15),size=0.5)+
  scale_y_log10(breaks = 10^(0:7),labels = trans_format("log10",math_format(10^.x)))+
  scale_x_discrete(limits = c("Phyla","Classes","Orders","Families","Genera"))+
  theme_bw()+
  labs(x=NULL,y="Variance of biosynthetic diversity")

if(!dir.exists("figs")) dir.create("figs")
ggsave("figs/Extended-Figure8g.pdf",width=90,height=150,units="mm")
