rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 1a. Geographic distribution of 43,191 newly recovered metagenome-assembled genomes (MAGs).
# =========================================================================================================================

# Load libraries
suppressMessages(library(scatterpie))
suppressMessages(library(maps))
suppressMessages(library(maptools))
suppressMessages(library(dplyr))

# Read data
world <- map_data("world")
names <- c("SRF","DCM","MES","BATH","Sediment","Host","Others")
colors <- c("forestgreen", "skyblue", "blue", "purple", "brown", "orange","steelblue")

df <- read.table("data/figure1a.global_distribution.txt", head = TRUE, sep = "\t") %>%
  mutate(radius_scaled_for_plot = 0.75*log(radius)+1) # Scale the radius (= number of MAGs) for visualization

# Visualize
plot.breaks <- 0.75*log(c(1, 20, 800, 12000))+1; names(plot.breaks) <- c(1, 20, 800, 12000) # Set breaks of radius (= number of MAGs) for legend

p <- ggplot(world, aes(x = long, y= lat)) +
  geom_map(map = world, aes(map_id = region), fill = "white", color = NA, linewidth = 0.001) +
  geom_scatterpie(data = df, aes(x = long, y = lat, group = region, r = radius_scaled_for_plot), color = NA, cols = names, alpha = 0.7, legend_name = "MAGs") +
  geom_scatterpie_legend(breaks = plot.breaks, labeller = function(x){names(plot.breaks[which(plot.breaks == x)])}, x = -160, y = -55, n = 4) +
  labs(title = "", x = "", y = "", fill = "Sample Type") +
  scale_fill_manual(values = colors) +
  coord_quickmap() +
  theme_bw() +
  theme(panel.background = element_rect(fill = '#EDEDEC', colour = 'NA'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1),
        legend.position = "right")

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure1a.pdf", width = 10, height = 8)
print(p)
dev.off()