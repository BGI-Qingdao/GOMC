rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 6a. The predicted OGT of genomes with or without Cas operon. 
# =========================================================================================================================

# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(scales)

my_comparison <- list(c("Present","Absent"))

dat <- read.table("data/ex-figure6a.drep95.Cas_ARG_type.MGE.stat.add",header = T,sep = "\t",quote="")
dat$phylum <- factor(dat$phylum,levels=c("Cyanobacteria","Actinobacteriota","Armatimonadota","Nitrospirota","Zetaproteobacteria","Planctomycetota","Desulfobacterota","Bipolaricaulota","Acidobacteriota","Aquificota","Chloroflexota","Deinococcota","Thermotogota","Myxococcota","Campylobacterota","Desulfobacterota_I","Spirochaetota","Verrucomicrobiota","Poribacteria","Firmicutes","Firmicutes_A","Firmicutes_B","Marinisomatota","Gemmatimonadota","WOR-3","Patescibacteria","Thiomicrospirales","Burkholderiales","Pseudomonadales","Chromatiales","Enterobacterales","Methylococcales","Nevskiales","Rhodobacterales","Rhizobiales","Caulobacterales","Rhodospirillales","Flavobacteriales","Rhodothermales","Ignavibacteriales","Cytophagales","Chitinophagales","Bacteroidales","Asgardarchaeota","Aenigmatarchaeota","Methanobacteriota","Thermoproteota","Nanoarchaeota","Halobacteriota","Thermoplasmatota"))

if(!dir.exists("figs")) dir.create("figs")
pdf('figs/Extended-Figure6a.pdf',width=30,height=55.5)
ggplot(data=dat, aes(x=phylum,y=OGT,fill=Cas)) +
  geom_boxplot(outlier.size=0.5) +
  theme_bw() +
  labs(y="Optimal growth temperature (°C)") +
  theme(axis.text.x = element_text(size =20,angle=0,color = "black"),
        axis.title.y=element_blank(),
        axis.text.y= element_text(size = 50,color = "black"),
        axis.title.x=element_text(size=50,face="bold"),
        legend.title=element_blank(),
        legend.position="top",
        legend.text = element_text(size =20,color = "black")) +
  scale_fill_manual(values=c(Present = "#F8B2AA", Absent = "#D9D6E9")) +
  coord_flip() +
  stat_compare_means(aes(label = paste0(..p.signif..)),color="black",method = "wilcox.test",size=10,label.y = 95)
dev.off()

