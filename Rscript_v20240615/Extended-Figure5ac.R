rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 5. Bioprospecting of IsPETase candidates.
# a. The distribution of 1,598 IsPETase candidates with Ser-Asp-His catalytic triad across varying marine ecosystems.
# c. Ecosystem origins of different clades of the IsPETase candidates
# =========================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# Read data
df.cnt <- read.table("data/ex-figure5.pet.IsPETase_stat.txt", header = TRUE, sep = "\t") %>% # Number of IsPETase candidates retrieved from distinct ecosystems
  mutate(Type = factor(Type, levels = c("SRF", "DCM", "MES", "BATH", "Hydrothermal", "Sediment", "Host-associated", "Others")))

df.eco <- read.table("data/ex-figure5.pet.ecosystem_stats.txt", header = TRUE, sep = "\t") %>% # The ecological makeup within individual clades of PETase candidates
  mutate(Type = factor(Type, levels = c("SRF", "DCM", "MES", "BATH", "Hydrothermal", "Sediment", "Host-associated", "Others"))) %>%
  mutate(Clade = factor(Clade, levels = paste("Clade", 1:7, sep = "")))
  
# Visualize
p.cnt <- ggplot(data = df.cnt, aes(x = "", y = Number, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.8) +
  geom_text(aes(label = paste(Type, "(", Number, ")", sep = "")), position = position_stack(vjust = 0.5)) +
  coord_polar("y", start = 1.56) +
  scale_fill_manual(values = c("#BFE1DC", "#F8F8D7", "#D4D3E4", "#FABCB0", "#B7CCDD", "#F9D1AD", "#D0E6B8", "#F7DEEA")) +
  labs(title = "a", x = "", y = "", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 24, face = "bold", vjust = 0, hjust = 0.1), 
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

p.eco <- ggplot(data = df.eco, aes(x = Clade, y = detail, fill = Type)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.8, color = "black", linewidth = 0.5) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("#BFE1DC", "#F8F8D7", "#D4D3E4", "#FABCB0", "#B7CCDD", "#F9D1AD", "#D0E6B8", "#F7DEEA")) +
  labs(title = "c", x = "", y = "Percentage (%) ", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 24, face = "bold", vjust = 0, hjust = 0), 
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.3),
        axis.text.y = element_text(size = 14, ),
        legend.position = "right")

if(!dir.exists("figs")) dir.create("figs")
pdf('figs/Extended-Figure5ac.pdf', width = 10, height = 6, onefile = F)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p.cnt), x = -0.1, y = 0, width = 0.7, height = 1) +
  draw_plot(ggplotify::as.ggplot(p.eco), x = 0.5, y = 0, width = 0.5, height = 1)
dev.off()
