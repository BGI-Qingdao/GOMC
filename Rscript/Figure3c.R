rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 3c. The number (the middle heatmap) and fraction of MAGs encoding Cas operon (the upper part) or ARG (the lower part)
# in metagenomic provinces represent various marine ecosystems.
# =========================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(aplot))
suppressMessages(library(cowplot))
suppressMessages(library("RColorBrewer"))

# Read data
env.cols <- c("#3399FF90", "#FB9A9990", "#33A02C90", "#FFC00090", "#CCCCFF", "#E3830390", "#FF000090")
names(env.cols) <- c("G1", "G2", "G3", "G4", "G5", "G6", "core")

df.mag <- read.table("data/figure3c.cas_stat.mag_in_umap.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  mutate(umap = factor(row.names(.), levels = row.names(.)), group = 1)

df.arg <- read.table("data/figure3c.cas_stat.arg_in_umap.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  mutate(group = factor(row.names(.), levels = rev(names(env.cols))))

df.cas <- read.table("data/figure3c.cas_stat.cas_in_umap.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  mutate(group = factor(row.names(.), levels = rev(names(env.cols))))

# Format the data
df.arg.long <- tidyr::gather(df.arg, umap, percent, 1:(ncol(df.arg)-1), factor_key = TRUE) %>%
  mutate(umap = factor(umap, levels = df.mag$umap))

df.cas.long <- tidyr::gather(df.cas, umap, percent, 1:(ncol(df.arg)-1), factor_key = TRUE) %>%
  mutate(umap = factor(umap, levels = df.mag$umap))

df.group <- df.arg.long %>% filter(percent > 0 & group != "core") %>% select(group, umap) %>% unique

# Visualize
p1 <- ggplot(data = df.mag, aes(x = umap, y = group, fill = MAGs)) +
  geom_tile(color = "grey", linewidth = 1) +
  geom_text(aes(label = MAGs), angle = 90) +
  scale_fill_gradient(low = "white", high = "#1F695F") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "", x = "", y = "Genomes", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

p2 <- ggplot(data = df.cas.long) +
  geom_bar(aes(x = umap, y = -percent, fill = group), stat = "identity", color = NA, linewidth = 0) +
  annotate("rect", xmin = 30, xmax = 31.5, ymin = -55, ymax = -48, fill = "#FF000090", color = NA) +
  annotate("text", x = 39, y = -51.5, label = "Genomes with Cas & ARG (%)", size = 6) +
  scale_y_continuous(expand = c(0, 0), breaks = -seq(0, 60, by = 20), labels = seq(0, 60, by = 20)) +
  scale_fill_manual(values = env.cols, labels = c("Cas & ARG", "Oil degrading", "Marine/Seawater", "Sediment", "Tidal flats", "Host-associated", "Hydrothermal vent")) +
  labs(title = "", x = "", y = "Genomes with Cas (%)", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1))


guides(colour = guide_legend())

p3 <- ggplot(data = df.arg.long) +
  geom_bar(aes(x = umap, y = percent, fill = group), stat = "identity", color = NA, linewidth = 0) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 80, by = 20)) +
  scale_fill_manual(values = env.cols) +
  labs(title = "", x = "", y = "Genomes with ARG (%)") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 0.3),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.position = "none")

p4 <- ggplot(data = df.group %>% mutate(name = "group")) +
  geom_bar(aes(x = name, fill = group), stat = "count", position = "stack", color = "white", linewidth = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "", x = "", y = "", fill = "") +
  scale_fill_manual(values = env.cols) +
  coord_flip() +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

p <- p1 %>% insert_top(p2, height = 6)
p <- p %>% insert_bottom(p3, height = 8)

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure3c.pdf", width = 14, height = 9, onefile = F)
print(p)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p), x = 0, y = 0, width = 1, height = 0.95) +
  draw_plot(ggplotify::as.ggplot(p4), x = 0.0345, y = 0.905, width = 0.8065, height = 0.09)
dev.off()