rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 8h. cAMP prediction using deep-learning models. The bar chart shows the novel or known number of RiPPs subtypes of 133 cAMPs.
# =========================================================================================================================

# Load libraries
library(ggplot2)

# Read data
data <- read.table("data/ex-figure8g.amp.txt",header = T,sep = "\t")
data$Novetly <- factor(data$Novetly,levels = c("Novel","Known"))

if(!dir.exists("figs")) dir.create("figs")
ggplot(data,aes(x=Subtype,y=AMP_number,fill=Novetly))+
  geom_histogram(stat = "identity",width=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
        legend.position = c(0.9,0.8))+
  labs(y="c_AMP number")
  # scale_fill_brewer(type = "qual", palette = "Set2")
ggsave("figs/Extended-Figure8h.pdf", width = 4.96,height=4.96,units = "in")