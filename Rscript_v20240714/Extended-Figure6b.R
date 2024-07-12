rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 6b. The uneven distribution patterns of defense systems.
# =========================================================================================================================

# Load libraries
library("gplots")
library("pheatmap")
library("RColorBrewer")

data.inter<-read.table("data/ex-figure6b.ecosystem.merge.ratio.sort.xls",header=T,sep = "\t",row.names=1)
data.inter <- as.matrix(data.inter)

#mycolors = colorRampPalette(colors = c("white","navy","firebrick3"))(100)
mycolors = colorRampPalette(colors = c("#F5F5F5","#1F695F"))(100)

#breaks = seq(0,160,length.out=200)
#gradient1 = colorpanel( sum( breaks[-1]<=5), "white", "navy")
#gradient2 = colorpanel( sum( breaks[-1]>5 ), "navy", "firebrick3")
#mycolors = c(gradient1,gradient2)
annotation_col=data.frame( ecosystem=c("Marine_GOMC","Marine_GOMC","Marine_GOMC","Marine_GOMC","Marine_GOMC","Aquatic_Thermal_springs","Aquatic_Thermal_springs","Aquatic_Thermal_springs","Aquatic_Non","Aquatic_Non","Aquatic_Non","Aquatic_Freshwater","Aquatic_Freshwater","Aquatic_Freshwater","Human_Digestive_system","Human_Digestive_system","Human_Digestive_system","Wastewater_Anaerobic_digestor","Wastewater_Anaerobic_digestor","Wastewater_Anaerobic_digestor","Mammals_Digestive_system","Mammals_Digestive_system","Mammals_Digestive_system","Terrestrial_Deep_subsurface","Terrestrial_Deep_subsurface","Terrestrial_Deep_subsurface","Terrestrial_Soil","Terrestrial_Soil","Terrestrial_Soil","Tibetan_Glacier","Tibetan_Glacier","Tibetan_Glacier"),type=c("ARG_MGE","ARG_MGE","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr","cas","crispr","Acr"))
rownames(annotation_col)=colnames(data.inter)
annotation_row=data.frame(phylum=c("Archaea","Archaea","Archaea","Archaea","Archaea","Bacteroidota","Bacteroidota","Bacteroidota","Bacteroidota","Bacteroidota","Bacteroidota","Bacteroidota","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Archaea","Archaea","Archaea","Bacteroidota","Bacteroidota","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Alphaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Gammaproteobacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"))
rownames(annotation_row)=rownames(data.inter)
ann_colors = list (phylum=c(Archaea="#3399FF",Bacteroidota="#FFC000",Alphaproteobacteria="#0BE33E",Gammaproteobacteria="#33A02C",Bacteria="#CCFF00"),type=c(ARG_MGE="white",cas="#66CCFF",crispr="#9999FF",Acr="#66FFFF"),ecosystem=c(Marine_GOMC="#E38303",Aquatic_Thermal_springs="#6699FF",Aquatic_Non="#FFFF00",Aquatic_Freshwater="#99CC00",Human_Digestive_system="#FF0000",Wastewater_Anaerobic_digestor="#00B050",Mammals_Digestive_system="#990099",Terrestrial_Deep_subsurface="#1F7884",Terrestrial_Soil="#00CCFF",Tibetan_Glacier="#FB8072"))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure6b.pdf", width=10,height=18)
obj <- pheatmap(
        data.inter,
        show_rownames=T,
        show_colnames=F,
		annotation_col = annotation_col,
		annotation_row = annotation_row,
		annotation_colors = ann_colors,
#       breaks=breaks,
#		annotation_row = annotation_row,
		col=mycolors,
        cluster_rows=F,
        cluster_cols=F,
		gaps_col=c(2,5,8,11,14,17,20,23,26,29),
		gaps_row=c(46),
#		cutree_rows=4,
        legend=T,
        fontsize=4,
		na_col="#FFFFFF",
#        main="\t\t\t\t\t  Hierarchical Clustering of DEGs(UV)",
#        display_numbers=T, number_format="%.2f",
		display_numbers = matrix(ifelse(data.inter == "NA","NA",""),nrow(data.inter)),
		cellwidth = 10,
		cellheight = 10,
        border_color="Grey60", #"WhiteSmoke",
)
dev.off()
