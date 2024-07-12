rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 4b. The workflow to identify Pfam domains potentially underpinning genome size expansion
# Part 2. The result of the ancestral proteome reconstruction analyses
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
df.votes <- read.table("data/ex-figure4b.votes.to.each.pfam.tsv", header = FALSE, sep = "\t", quote = "") %>% `colnames<-`(c("pfam","vote"))
df.votes <- df.votes %>% arrange(desc(vote))

p <- ggplot(data = df.votes, aes(x = factor(pfam, levels = pfam), y = vote)) +
  geom_bar(fill = "grey", width = 0.4, color = NA, stat = "identity") +
  geom_bar(data = df.votes %>% filter(vote >= 3), aes(x = factor(pfam, levels = pfam), y = vote), fill = "darkred", width = 0.4, color = NA, stat = "identity") +
  geom_hline(yintercept = 3, linetype = 2, linewidth = 0.5, color = "black") +
  scale_y_continuous(expand = c(0,0), breaks = c(0,3,10,20)) +
  theme_bw() +
  labs(x = "Pfam domains ranked by votes", y = "Vote", title = "") +
  theme(plot.title = element_text(size = 32, hjust = -0.085),
        axis.title.x = element_text(size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 32, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24)) +
  theme(panel.background = element_rect(fill = 'transparent', color = NA),
        panel.border = element_blank(),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.1, linetype = "solid", color = "black"),
        axis.line.y = element_line(linewidth = 0.1, linetype = "solid", color = "black"))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure4b2.pdf", width = 10, height = 9, onefile = FALSE)
print(p)
dev.off()