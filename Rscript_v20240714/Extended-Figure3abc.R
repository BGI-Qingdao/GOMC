rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 3abc. Microbial biogeography and metagenomic provinces (MPs).
# =================================================================================================================================

# Load libraries
suppressMessages(library(phyloseq)); packageVersion("phyloseq")
suppressMessages(library(microbiome)); packageVersion("microbiome")
suppressMessages(library(tidyverse)); packageVersion("tidyverse")
suppressMessages(library(plyr)); packageVersion("plyr")
suppressMessages(library(cluster)); packageVersion("cluster")
suppressMessages(library(ggplot2)); packageVersion("ggplot2")
suppressMessages(library(ggplotify)); packageVersion("ggplotify")
suppressMessages(library(ggpubr)); packageVersion("ggpubr")
suppressMessages(library(ggalluvial)); packageVersion("ggalluvial")
suppressMessages(library(gridExtra)); packageVersion("gridExtra")
suppressMessages(library(grid)); packageVersion("grid")
suppressMessages(library(RColorBrewer)); packageVersion("RColorBrewer")

# Read data
phylum2cols <- read.table("data/ex-figure3.umap_stat.cols4phyla.txt", header = TRUE, sep = "\t", comment.char = "") # colors for phyla

ps <- readRDS("data/ex-figure3.umap_stat.braken.RDS") # abundance table

if(file.exists("data/ex-figure3.umap_stat.jsd_dist.RDS")){
  ps.dist <- readRDS("data/ex-figure3.umap_stat.jsd_dist.RDS") %>% as.matrix %>% data.frame
}else{
  ps.dist <- phyloseq::distance(ps, method = "jsd")
  saveRDS(dist, "data/ex-figure3.umap_stat.jsd_dist.RDS")
}

umap.coordinates <- read.table("data/ex-figure3.umap_stat.coords.txt", sep = "\t", header = FALSE, skip = 1, row.names = 1) %>%
  select(V2,V3,V4) %>%
  `colnames<-`(c("group","UMAP1","UMAP2")) %>%
  mutate(group = paste("U", group, sep = ""))

# Visualize
if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure3abc.pdf", width = 18, height = 16)
define_region <- function(row, col){viewport(layout.pos.row = row, layout.pos.col = col)}
pushViewport(viewport(layout = grid.layout(nrow = 19, ncol = 18)))

# a. Distribution by depth
fig.by.depth <- list()

smp.types <- sample_data(ps) %>%
  data.frame() %>%
  pull(smp.type) %>%
  unique()

for(smp.type.slt in smp.types){

  print(paste("Start processing ", smp.type.slt, " ...", sep = ""))

  smp.type.slt.title <- case_when(smp.type.slt=="Host" ~ 'Host-associated', smp.type.slt!="Host" ~ smp.type.slt)

  ps.subset <- ps %>% subset_samples(smp.type == smp.type.slt)
  ps.subset <- prune_samples(sample_sums(ps.subset) > 0, ps.subset)
  ps.subset <- prune_taxa(taxa_sums(ps.subset) > 0, ps.subset)

  ps.taxa <- tax_table(ps.subset) %>% as.data.frame() %>% select(kingdom, phylum, class) %>%
             mutate(taxa = case_when(grepl("Crenarchaeota",phylum) ~ "Archaea|Crenarchaeota",
                                     grepl("Thermoplasmatota",phylum) ~ "Archaea|Thermoplasmatota",
                                     kingdom == "Archaea" & ! grepl("Thermoplasmatota|Crenarchaeota",phylum) ~ "Archaea|Other",
                                     grepl("Actinobacteriota",phylum) ~ "Bacteria|Actinobacteriota",
                                     grepl("Bacteroidota",phylum) ~ "Bacteria|Bacteroidota",
                                     grepl("Campylobacterota",phylum) ~ "Bacteria|Campylobacterota",
                                     grepl("Chloroflexota",phylum) ~ "Bacteria|Chloroflexota",
                                     grepl("Cyanobacteria",phylum) ~ "Bacteria|Cyanobacteria",
                                     grepl("Desulfobacterota",phylum) ~ "Bacteria|Desulfobacterota",
                                     grepl("Firmicutes",phylum) ~ "Bacteria|Firmicutes",
                                     grepl("Marinisomatota",phylum) ~ "Bacteria|Marinisomatota",
                                     grepl("Planctomycetota",phylum) ~ "Bacteria|Planctomycetota",
                                     class == "Alphaproteobacteria" ~ "Bacteria|Alphaproteobacteria",
                                     class == "Gammaproteobacteria" ~ "Bacteria|Gammaproteobacteria",
                                     phylum == "Proteobacteria" & ! class %in% c("Alphaproteobacteria","Gammaproteobacteria") ~ "Bacteria|Other Proteobacteria",
                                     grepl("SAR324",phylum) ~ "Bacteria|SAR324",
                                     grepl("Verrucomicrobiota",phylum) ~ "Bacteria|Verrucomicrobiota",
                                     kingdom == "Bacteria" & ! grepl("Actinobacteriota|Bacteroidota|Campylobacterota|Chloroflexota|Cyanobacteria|Desulfobacterota|Firmicutes|Marinisomatota|Planctomycetota|Proteobacteria|SAR324|Verrucomicrobiota",phylum) ~ "Bacteria|Other"))
  tax_table(ps.subset) <- as(ps.taxa, "matrix")

  ps.subset <- ps.subset %>% 
    aggregate_taxa(level = "taxa") %>% 
    transform(transform = "compositional")

  ps.meltd <- psmelt(ps.subset)[,c("Sample","Abundance","taxa")]
  ps.meltd <- ps.meltd %>% dplyr::group_by(Sample,taxa) %>% 
    dplyr::summarise(Abundance = sum(Abundance))

  smps.order <- ps.meltd %>% 
    filter(taxa == "Bacteria|Other") %>% 
    arrange(desc(Abundance)) %>% 
    pull(Sample)

  ps.meltd <- ps.meltd %>%
              mutate(Sample = factor(Sample, levels = smps.order)) %>%
              mutate(taxa = factor(taxa, levels = c("Archaea|Crenarchaeota","Archaea|Thermoplasmatota","Archaea|Other","Bacteria|Actinobacteriota",
                                                    "Bacteria|Bacteroidota","Bacteria|Campylobacterota","Bacteria|Chloroflexota","Bacteria|Cyanobacteria",
                                                    "Bacteria|Desulfobacterota","Bacteria|Firmicutes","Bacteria|Marinisomatota","Bacteria|Planctomycetota",
                                                    "Bacteria|Alphaproteobacteria","Bacteria|Gammaproteobacteria","Bacteria|Other Proteobacteria",
                                                    "Bacteria|SAR324","Bacteria|Verrucomicrobiota","Bacteria|Other")))

  fig.by.depth[[smp.type.slt]] <- ggplot() +
    geom_bar(data = ps.meltd, aes(fill = taxa, y = Abundance, x = Sample), position_stack(reverse = TRUE), stat = "identity", width = 0.8, linewidth = 0) + 
    scale_fill_manual(values = phylum2cols$col) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    theme_bw() +
    labs(x = "", y = "", fill = "") +
    theme(plot.margin = unit(c(0.1,0,0,0), "cm"),
          axis.line = element_blank(),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.01),
          axis.text.x = element_blank(),
          axis.text.y = element_text(margin = margin(t=0,r=0.5,b=0,l=0), size = 12),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth = 0.1),
          axis.ticks.length = unit(0.1,"cm"))

  if(smp.type.slt == "Sediment"){
    fig.by.depth[[smp.type.slt]] <- fig.by.depth[[smp.type.slt]] +
      theme(legend.position = c(1,-0.45),
            legend.title = element_blank(),
            legend.direction = "horizontal",
            legend.text = element_text(size = 12),
            legend.box.background = element_blank(),
            legend.background = element_blank()) +
      guides(fill = guide_legend(ncol = 3, byrow = TRUE))
  }else{
    fig.by.depth[[smp.type.slt]] <- fig.by.depth[[smp.type.slt]] + theme(legend.position = "none")
  }

  if(smp.type.slt == "SRF"){fig.by.depth[[smp.type.slt]] <- fig.by.depth[[smp.type.slt]] + labs(title = paste("a ",smp.type.slt.title,sep="")) + theme(plot.title = element_text(size = 24, face = "bold", hjust = -0.075))}
  if(smp.type.slt != "SRF"){fig.by.depth[[smp.type.slt]] <- fig.by.depth[[smp.type.slt]] + labs(title = smp.type.slt.title) + theme(plot.title = element_text(size = 24, face = "bold"))}
}

print(fig.by.depth[["SRF"]], vp = define_region(row = 1:4, col = 1:7))
print(fig.by.depth[["DCM"]], vp = define_region(row = 5:8, col = 1:7))
print(fig.by.depth[["MES"]], vp = define_region(row = 9:12, col = 1:5))
print(fig.by.depth[["BATH"]], vp = define_region(row = 9:12, col = 6:7))
print(fig.by.depth[["Sediment"]], vp = define_region(row = 13:16, col = 1:4))
print(fig.by.depth[["Host"]], vp = define_region(row = 13:16, col = 5:7))
# end of a.

# b. UMAP clusters
umaphylum2cols <- c("#3399FF","#0BE33E","#255AEF","#33A02C","#FFC000","#00CCCC","#DDABDD","#7BD37B","#B2DF8A","#CC33FF"
                    ,"#D09E00","#FF8001","#FF99CC","#93E9F7","#666699","#006600","#FF66FF","#FB9A99","#00CC99","#CC0099",
                    "#3366FF","#B15928","#6FAF2F","#9933FF","#CCFF00","#CCCCFF","#66CCFF","#CC3300","#9999FF","#66FFFF"
                    ,"#E38303","#6699FF","#FFFF00","#99CC00","#FF0000","#D60093","#00B050","#990099","#1F7884","#00CCFF"
                    ,"#67D5DB","#FB8072","#80B1D3","#A6D854","#DECBE4","#F2F2F2","#FCCDE5","#E7298A","#FC8D62","#7FC97F"
                    ,"#E6AB02","#8DD3C7","#E31A1C","#B3CDE3","#A6761D","#377EB8")
names(umaphylum2cols) <- c("U0","U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14","U15","U16",
                           "U17","U18","U19","U20","U21","U22","U23","U24","U25","U26","U27","U28","U29","U30","U31","U32",
                           "U33","U34","U35","U36","U37","U38","U39","U40","U41","U42","U43","U44","U45","U46","U47","U48",
                           "U49","U50","U51","U52","U53","U54","U55")

umap.labs <- data.frame()
for(umap.slt in unique(umap.coordinates$group)){
  umap.coordinates.slt <- umap.coordinates %>% filter(group == umap.slt)
  umap.medoid <- cluster::pam(umap.coordinates.slt[,c("UMAP1","UMAP2")], k = 1)$medoids
  umap.labs <- rbind(umap.labs, data.frame(group = umap.slt, UMAP1 = umap.medoid[1], UMAP2 = umap.medoid[2]))
}

fig.umap <- ggplot() +
  geom_point(data = umap.coordinates, aes(x = UMAP1, y = UMAP2, color = factor(group, levels = names(umaphylum2cols))), size = 1, alpha = 0.9) +
  scale_color_manual(values = umaphylum2cols) +
  geom_text(data = umap.labs, aes(x = UMAP1, y = UMAP2, label = group)) +
  theme_bw() +
  labs(title = "b", x = "UMAP1", y = "UMAP2") +
  scale_x_continuous(limits = c(min(umap.coordinates$UMAP1)-0.5, max(umap.coordinates$UMAP1)+0.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(min(umap.coordinates$UMAP2)-0.5, max(umap.coordinates$UMAP2)+0.5), expand = c(0,0)) +
  theme(plot.margin = unit(c(0,1,1,0.3), "cm"),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold", vjust = -3),
        plot.title = element_text(size = 24, face = "bold", hjust = -0.07, vjust = 0.1),
        axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.01),
        legend.position = "none")

print(fig.umap, vp = define_region(row = 1:11, col = 8:18))
# end of b.

# c. UMAP distribution
flow.df <- sample_data(ps) %>% 
  data.frame() %>%
  filter(smp.type %in% c("SRF","DCM","MES","BATH") & ! is.na(lat)) %>%
  select(c("umap","smp.type","lat")) %>%
  arrange(desc(lat)) %>%
  mutate(region = case_when(abs(lat) < 23.5 ~ 'T', between(lat, 23.5, 66.5) ~ 'NT', between(lat, -66.5, -23.5) ~ 'ST', lat > 66.5 ~ 'NP', lat < -66.5 ~ 'SP')) %>%
  mutate(smp.type = case_when(smp.type=="SRF" ~ 'S', smp.type=="DCM" ~ 'D', smp.type=="MES" ~ 'M', smp.type=="BATH" ~ "B"))

flow.df <- flow.df %>%
  dplyr::group_by(region, umap, smp.type) %>%
  dplyr::count(umap) %>%
  as.data.frame(.)

flow.df$region <- factor(flow.df$region, levels = c("NP","NT","T","ST","SP"))
flow.df$umap <- factor(flow.df$umap, levels = names(umaphylum2cols))
flow.df$smp.type <- factor(flow.df$smp.type, levels = c("S","D","M","B"))

is_alluvia_form(flow.df, axes = 1:3, silent = TRUE) # check format

fig.flow <- ggplot(flow.df, aes(y = n, axis1 = region, axis2 = umap, axis3 = smp.type)) +
  geom_alluvium(aes(fill = umap), alpha = 0.7) +
  geom_stratum(width = 0.08, reverse = TRUE, color = "grey", linewidth = 0) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = TRUE) +
  scale_x_continuous(breaks = 1:3, labels = c("Climate Zones", "UMAP", "Depth"), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = umaphylum2cols) +
  theme(plot.margin = unit(c(-2,1,1,2.5), "cm")) +
  labs(title = "c") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 24, face = "bold", hjust = -0.07, vjust = 0.1),
        axis.line = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.01),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

print(fig.flow, vp = define_region(row = 12:19, col = 8:18))
# end of c.

dev.off()
