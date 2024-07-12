rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# ==================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 3d. Boxplots display the incidence rate of Cas operon in genomes from different ecosystems.
# ==================================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(aplot))
suppressMessages(library(cowplot))

# Read data
df <- read.table("data/figure3d.cas_stat.ecosystem.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  mutate(group = row.names(.))

# Format the data
df.long <- tidyr::gather(df, taxon, ratio, 1:(ncol(df)-1), factor_key = TRUE) %>%
  mutate(taxon = factor(taxon, levels = colnames(df)))

# Visualize
p <- ggplot() +
  geom_bar(data = df.long %>% filter(group == "Cas"), aes(x = taxon, y = ratio, fill = group), stat = 'identity', position = 'stack', width = 0.8) +
  geom_text(data = df.long %>% filter(group == "Cas"), aes(x = taxon, y = ratio + 1.5, label = round(ratio, 2)), size = 5) +
  geom_line(data = df.long %>% filter(group == "Crispr"), aes(x = taxon, y = ratio), color = "#dbd9da", linewidth = 1, group = 1) +
  geom_point(data = df.long %>% filter(group == "Crispr"), aes(x = taxon, y = ratio, color = group), size = 4) +
  geom_vline(xintercept = c(1.5, 4.5, 7.5), color = "#4169e0") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(df.long$ratio)+5)) +
  scale_fill_manual(values = c("#84a0d2"), labels = c("Cas operons")) +
  scale_color_manual(values = c("#dbd9da"), labels = c("CRISPR arrays")) +
  labs(title = "", x = "", y = "Incidence rate (%)", fill = "", color = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y.left = element_text(size = 20, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1),
        axis.ticks.length.y = unit(0.2, "cm"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.box.spacing = unit(0, "cm"),
        legend.spacing = unit(0, "cm")) +
  guides(fill = guide_legend(position = "inside", override.aes = list(shape = 15, size = 6), order = 1),
         color = guide_legend(position = "inside", override.aes = list(shape = 16, size = 7.5), order = 2))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure3d.pdf", width = 7, height = 7, onefile = F)
print(p)
dev.off()
