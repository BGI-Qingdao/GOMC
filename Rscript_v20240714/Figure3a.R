rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# ==================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 3a. Barplot indicates the frequency of Cas operons in different lineages of all GOMC genomes.
# ==================================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(aplot))
suppressMessages(library(cowplot))

# Read data
cas.cols <- c("#5c91cc", "#5caf2f", "#309737", "#fac000", "#43b9c2", "#deb4d3", "#bdd584", "#9f5097", "#ed7900", "#f09cc3", "#addded",
              "#cd78ad", "#f19c9e", "#46b38a", "#c51a82", "#4d61aa", "#b4541b", "#6fa92c", "#804f9a", "#84caf0", "grey97")
names(cas.cols) <- c("I-A", "I-B", "I-C", "I-D", "I-E", "I-F", "I-G", "II-A", "II-B", "II-C", "III-A", 
                     "III-B", "III-C", "III-D", "III-E", "III-F", "IV", "V", "VI", "Hybrid", "None")

df.cas <- read.table("data/figure3a.cas_stat.ratio.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  mutate(Type = row.names(.)) %>%
  filter(Type != "MAGs")

df.mag <- read.table("data/figure3a.cas_stat.ratio.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  mutate(Type = row.names(.)) %>%
  filter(Type == "MAGs")

# Format the data
df.cas.long <- tidyr::gather(df.cas, taxon, percent, 1:(ncol(df.cas)-1), factor_key = TRUE) %>%
  mutate(taxon = factor(taxon, levels = colnames(df.cas))) %>%
  mutate(Type = factor(Type, levels = rev(names(cas.cols))))

df.mag.long <- tidyr::gather(df.mag, taxon, count, 1:(ncol(df.mag)-1), factor_key = TRUE) %>%
  mutate(taxon = factor(taxon, levels = colnames(df.mag)))

# Visualize
p1 <- ggplot(data = df.cas.long) +
  geom_bar(aes(x = taxon, y = percent, fill = Type), stat = 'identity', position = 'stack', width = 0.8) +
  scale_fill_manual(values = cas.cols) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "", x = "", y = "Cas operons fration (%)", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y.left = element_text(size = 24, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1),
        axis.ticks.length.y = unit(0.2, "cm"),
        axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.2),
        axis.text.y = element_text(size = 14),
        legend.position = "right",
        legend.box.background = element_blank(),
        legend.box.margin = unit(c(0, 0, 0, 2), "cm")) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1))

p2 <- ggplot(data = df.mag.long) +
  geom_segment(aes(x = taxon, xend = taxon, y = 0, yend = count), color = "#1d4ea1", linewidth = 1) +
  geom_point(aes(x = taxon, y = count), color = "#1d4ea1", size = 2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500), position = "right") +
  labs(title = "", x = "", y = "Genome count") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 1, color = "#1d4ea1"),
        axis.title.x = element_blank(),
        axis.title.y.right = element_text(size = 24, face = "bold", angle = 90, color = "#1d4ea1"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1, color = "#1d4ea1"),
        axis.ticks.length.y = unit(0.2, "cm"),
        axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.1, color = NA),
        axis.text.y = element_text(size = 14, color = "#1d4ea1"))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure3a.pdf", width = 24, height = 8, onefile = F)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p1), x = 0, y = 0, width = 1, height = 1) +
  draw_plot(ggplotify::as.ggplot(p2), x = 0.0325, y = 0, width = 0.931, height = 1)
dev.off()