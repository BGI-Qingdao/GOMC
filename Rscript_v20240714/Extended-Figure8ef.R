rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 8e and f, Rarefaction curves of the top 4 phyla (the other 16 phyla of the top 20 in embedded figure) and
# top 20 bacterial genera with most predicted biosynthetic potential, respectively. 
# =========================================================================================================================

# Load libraries
library(iNEXT)
library(ggplot2)
library(scales)

# Pass parameters
args<-commandArgs(T)

options(scipen=200)

# Read data
data <- list()
a <- list.files(paste0(args[1],"/iNEXT_Input"),pattern=".txt")
name <- gsub(".txt$","",a)
dir <- paste0(args[1],"/iNEXT_Input/",a)
for (i in 1:length(a)){
	data[[i]] <- read.table(dir[i], header = T, as.is = F)
}
names(data)<-name
# Run iNEXT
out <- iNEXT(data,datatype="incidence_raw",endpoint=as.numeric(args[2]),nboot=as.numeric(args[3]))
save(out, file="iNEXT_phylum.rds")

#load("iNEXT_genus.rds")
x <- 1:20
y <- dpois(x, lambda = 10)
data <- data.frame(X=x,y=y)
data$type <- as.factor(x)
shape_level <- nlevels(data[["type"]])
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}

p <- ggiNEXT(out, type=1, se = args[4])+
theme_bw()+
labs(x="Genomes",y="GCFs")+
scale_y_continuous(trans = log2_trans(),
    # breaks = trans_breaks("log2", function(x) 2^x),
    breaks = 2^(0:14),
    labels = trans_format("log2", math_format(2^.x))
)+
scale_shape_manual(values=shapes)+
coord_cartesian(ylim = c(2^5,2^13))
# theme(legend.position = "none")

if(!dir.exists("figs")) dir.create("figs")
ggsave("figs/Extended-Figure8ef.pdf", p, height=120, width=180,units="mm")

for (i in 1:length(out$iNextEst)){
  write.table(out$iNextEst[i],file=paste0("./temp/",names(out$iNextEst)[i],".txt"),quote = F,sep = "	",row.names = F)
}
