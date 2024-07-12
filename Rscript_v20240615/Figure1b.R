rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 1b. A total of 43,191 MAGs had at least medium to higher quality
# =========================================================================================================================

# Load libraries
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(aplot))
suppressMessages(library(cowplot))

# Read data
df.quality <- read.table("data/figure1b.bin_quality.list", header = T, sep = "\t") %>%
  mutate(Type = factor(Type, levels = c("near", "high", "medium")))

df.unclassified <- read.table("data/figure1b.bin_unclassified.txt", header = T, sep = "\t") %>%
  mutate(phyla = factor(phyla, levels = c("Species","Genus","Family","Order","Class","Phylum")))

# Visualize
df.labels <- c()
df.labels <- c(df.labels, paste("Near complete (", nrow(df.quality %>% filter(Type == "near")), ", ", round(100 * nrow(df.quality %>% filter(Type == "near")) / nrow(df.quality), 2), "%)", sep = ""))
df.labels <- c(df.labels, paste("High quality (", nrow(df.quality %>% filter(Type == "high")), ", ", round(100 * nrow(df.quality %>% filter(Type == "high")) / nrow(df.quality), 2), "%)", sep = ""))
df.labels <- c(df.labels, paste("Medium quality (", nrow(df.quality %>% filter(Type == "medium")), ", ", round(100 * nrow(df.quality %>% filter(Type == "medium")) / nrow(df.quality), 2), "%)", sep = ""))

p1 <- ggplot(data = df.quality, aes(x = Completeness, y = Contamination, colour = Type)) +
  geom_point(size = 0.2) +
  labs(x = "Completeness (%)", y = "Contamination (%)", colour = "") +
  scale_color_manual(values = c("#619cff", "#f8766d", "#00ba38"), labels = df.labels) +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill = 'transparent', color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.ticks = element_line(linewidth = 1),
        legend.position = "right",
        legend.text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 10)))

p2 <- ggplot(data = df.quality) +
  geom_histogram(aes(x = Completeness, y = after_stat(count / sum(count)), fill = Type), color = "white", breaks = seq(50, 100, by = 2), width = 0.8, linewidth = 0) +
  scale_fill_manual(values = c("#619cff", "#f8766d", "#00ba38")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.105), breaks = seq(0, 0.1, by = 0.02), labels = c("", "", "", "", "", "10%")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(linewidth = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1),
        axis.ticks.length.y = unit(0.2, "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold"),
        legend.position = "none")

p3 <- ggplot(data = df.quality) +
  geom_histogram(aes(x = Contamination, y = after_stat(count / sum(count)), fill = Type), color = "white", breaks = seq(0, 10, by = 0.5), width = 0.8, linewidth = 0) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.25), breaks = seq(0, 0.25, by = 0.05), labels = c("", "", "", "", "", "25%")) +
  scale_fill_manual(values = c("#619cff", "#f8766d", "#00ba38")) +
  coord_flip() +
  labs(x = "", y = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color= NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 0.5),
        axis.ticks.x = element_line(linewidth = 1),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.2, "cm"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none")

p4 <- ggplot(data = df.unclassified) +
  geom_bar(aes(x = phyla, y = -log10(Number_of_genomes)), stat = "identity", color = "grey", fill = "grey") +
  geom_text(aes(x = phyla, y = -log10(Number_of_genomes), label = Number_of_genomes), size = 5, position = position_dodge(0.9), vjust = 1.5) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(limits = c(-log10(max(df.unclassified$Number_of_genomes)) - 0.3, 0), expand = c(0, 0), breaks = -4:0, labels = 10^(4:0)) +
  labs(title = "", x = "", y = "Unclassified genomes") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(linewidth = 1),
        axis.title = element_text(size = 16, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.position = "none")

p <- p1 %>% insert_top(p2, height = 0.25)
p <- p %>% insert_right(p3, width = 0.22)

if(!dir.exists("figs")) dir.create("figs")
pdf('figs/Figure1b.pdf', width = 15, height = 8, onefile = F)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p), x = 0, y = 0, width = 1, height = 1) +
  draw_plot(ggplotify::as.ggplot(p4), x = 0.075, y = 0.33, width = 0.35, height = 0.5)
dev.off()
