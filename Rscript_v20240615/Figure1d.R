rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 1d. Contribution of this current study and extant published databases to each bacterial and archaeal phylum.
# =========================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# Read data
df <- read.table("data/figure1d.bin_taxonomy.txt", header = T, sep = "\t") %>%
  mutate(Type = factor(Type, levels = c("our", "OMD", "OceanDNA", "NCBI", "Interset"))) %>%
  mutate(Phylum = factor(Phylum, levels = c("Pelagibacterales", "Rhodobacterales", "Other Alphaproteobacteria", "Pseudomonadales",
                            "Enterobacterales", "Other Gammaproteobacteria", "Zetaproteobacteria", "Acidobacteriota",
                            "Actinobacteriota", "Bacteroidota", "Bdellovibrionota", "Campylobacterota", "Chloroflexota",
                            "Cyanobacteria", "Desulfobacterota", "Firmicutes", "Firmicutes_A", "Gemmatimonadota",
                            "Marinisomatota", "Myxococcota*", "Patescibacteria", "Planctomycetota", "Spirochaetota",
                            "Verrucomicrobiota", "Other Bacteria", "Halobacteriota", "Nanoarchaeota", "Thermoplasmatota",
                            "Thermoproteota", "Other Archaea")))

p <- ggplot(data = df, aes(x = Phylum, y = detail, fill = Type)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.8) +
  geom_segment(x = 24.8, y = 2090, xend = 24.8, yend = 3280) +
  annotate("text", label = "Specific", x = 24.3, y = 2700, size = 6, colour = "black", angle = 90) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#EB6060", "#41B263", "#468CCB", "#B4AFD5", "#D2D2D2"), labels = c("Newly recovered", "OMD", "OceanDNA", "NCBI", "Intersection")) +
  labs(title = "", x = "", y = "Genome number", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        axis.title = element_text(size = 20, face = "bold"),
        axis.ticks.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.3),
        axis.text.y = element_text(size = 15),
        legend.position.inside = c(0.9, 0.8),
        legend.text = element_text(size = 13),
        legend.background = element_blank(),
        legend.box.background = element_blank()) +
  guides(fill = guide_legend(position = "inside", override.aes = list(shape = 15, size = 7.5)))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure1d.pdf", width = 11, height = 6.5)
print(p)
dev.off()

