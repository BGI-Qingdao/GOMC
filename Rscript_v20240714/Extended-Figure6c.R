rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 6c. The trend indicates a decrease in the upper limit number of MGEs with increased number of Cas operons.
# =========================================================================================================================

# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(scales)

dat <- read.table("data/ex-figure6c.cas_num.stat3",header = T,sep = "\t",quote="")
dat$Cas <- factor(dat$Cas,levels=c("1","2","3","4","5"))

if(!dir.exists("figs")) dir.create("figs")
pdf('figs/Extended-Figure6c.pdf',width=15,height=12)
ggplot(data=dat, aes(x=Cas,y=MGE,fill=Cas)) +
  geom_violin(position=position_dodge(0.5)) +
  geom_point(color = "black", size = 4) +
  theme_bw() +
  labs(y="MGE number") +
  theme(axis.text.x = element_text(size =20,angle=0,color = "black"),
        axis.title.y=element_text(size=25,face="bold"),
        axis.text.y= element_text(size = 20,color = "black"),
        axis.title.x=element_text(size=25,face="bold"),
        legend.title=element_blank(),
        legend.position="top",
        legend.text = element_text(size =20,color = "black")) +
  scale_fill_manual(values=c("1"="#FFEFCC", "2"="#FECCCC", "3"="#99CCCC", "4"="#CCCCFF", "5"="#3399CC")) +
  scale_y_continuous(trans=log2_trans(),breaks=trans_breaks("log2",function(x) 2^x),labels = trans_format("log2", math_format(2^.x))) 
dev.off()
#scale_y_continuous(expand = c(0, 0))

