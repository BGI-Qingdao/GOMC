rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 3d. Line plots show the fractions of genomes coding ARG with (red line) or without (blue line) the presence of Cas operon.
# =================================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# Read data
df.ratio <- read.table("data/figure3d.cas_stat.cas_vs_arg.txt", header = TRUE, row.names = 1, sep = "\t") %>%
  mutate(diff = with_cas - without_cas) %>%
  mutate(group = case_when(diff > 0 ~ "with_higher", TRUE ~ "without_higher")) %>%
  mutate(group = factor(group, levels = c("with_higher", "without_higher"))) %>%
  mutate(upper_limit = case_when(with_cas > without_cas ~ with_cas, TRUE ~ without_cas)) %>%
  mutate(lower_limit = case_when(with_cas < without_cas ~ with_cas, TRUE ~ without_cas)) %>%
  mutate(taxon = factor(row.names(.), levels = row.names(.)))

# Visualize
p <- ggplot(data = df.ratio) +
  geom_bar(aes(x = taxon, y = upper_limit, fill = group), color = NA, width = 0.8, stat = 'identity', position = 'stack') +
  geom_bar(aes(x = taxon, y = lower_limit), fill = "white", color = NA, width = 0.8, stat = 'identity', position = 'stack') +
  geom_line(aes(x = taxon, y = without_cas), color = "#6495ED", group = 1) +
  geom_line(aes(x = taxon, y = with_cas), color = "#F08080", group = 1) +
  annotate("rect", xmin = 34.5, xmax = 35.5, ymin = 82, ymax = 87, fill = "white", color = "black") +
  annotate("text", x = 42.85, y = 84.7, label = "The difference between blue ratio and red ratio", size = 3) +
  annotate("rect", xmin = 34.4, xmax = 35.5, ymin = 91, ymax = 91.3, fill = "#F08080", color = "#F08080") +
  annotate("text", x = 42.9, y = 91, label = "ARG frequency for genomes with Cas operons", size = 3) +
  annotate("rect", xmin = 34.4, xmax = 35.5, ymin = 97.3, ymax = 97.6, fill = "#6495ED", color = "#6495ED") +
  annotate("text", x = 43.4, y = 97.3, label = "ARG frequency for genomes without Cas operons", size = 3) +
  scale_fill_manual(values = c("#F08080", "#6495ED")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(title = "", x = "", y = "Genomes with ARG (%)", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.3),
        axis.text.y = element_text(size = 14),
        legend.position = "none")

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure3d.pdf", width = 9, height = 3.8, onefile = F)
print(p)
dev.off()