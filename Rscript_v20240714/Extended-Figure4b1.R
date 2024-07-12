rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 4b. The workflow to identify Pfam domains potentially underpinning genome size expansion
# Part 1. The result of the phylogenetic regression analyses
# =========================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(ggpubr))
suppressMessages(library(ape))
suppressMessages(library(ggtree))
suppressMessages(library(cowplot))

# Read data
df <- read.table("data/ex-figure4b.counts.of.genes.with.pfams.in.each.genome.vs.gsize.pr.tsv", header = TRUE, sep = "\t", quote = "", comment.char = "")
df <- df %>% filter(!is.na(lambda)) %>% mutate(psignal = case_when(lambda <= 0.1 ~ '(0,0.1]', lambda > 0.1 & lambda <= 0.5 ~ '(0.1,0.5]',  lambda > 0.5 ~ '(0.5,1]'))

p <- ggplot(data = df, aes(x = fam.coefficient, y = -log10(fam.fdr), fill = factor(psignal, levels = c("(0,0.1]","(0.1,0.5]","(0.5,1]")))) +
  geom_point(aes(size = r2), alpha = 0.5, shape = 21, stroke = 0.1) +
  geom_vline(xintercept = 0, col = "black", linetype = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 2) +
  scale_size(range = c(0, 10)) +
  scale_fill_manual(name = "Lambda", values = c("darkred", "forestgreen", "blue")) +
  theme_bw() +
  labs(x = "Regression Coefficient (Gene Counts)", y = "-log10(FDR)", title = "") +
  guides(fill = guide_legend(title = "Lambda", override.aes = list(size = 8), nrow = 1, direction = "horizontal"),
         size = guide_legend(title = "R2", nrow = 1, direction = "horizontal")) +
  theme(plot.title = element_text(size = 32, hjust = -0.08, vjust = -4),
        strip.text.x = element_text(size = 32),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = 24),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.box.background = element_rect(fill = 'transparent'),
        legend.position = c(0.24, 0.92),
        legend.box = "vertical",
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure4b1.pdf", width = 13, height = 15, onefile = FALSE)
print(p)
dev.off()
