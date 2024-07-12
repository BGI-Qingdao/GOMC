rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =====================================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 5a. Statistics of the phylogenetic regression analyses between the selected 77 Pfam domains and bacterial genome sizes across multiple phyla
# =====================================================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))

# Read data
df.gsize <- read.table("data/ex-figure5.pr.subset.genome.size.of.phyla.with.large.genomes.xls", header = TRUE, row.names = 1, sep = "\t", quote = "", comment.char = "")
df.gsize <- aggregate(df.gsize$estimate_size, list(df.gsize$tgroup), FUN = mean) %>%
  `colnames<-`(c("tgroup", "mean")) %>%
  mutate(tgroup = case_when(grepl('d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria', tgroup) ~ "d__Bacteria;c__Gammaproteobacteria", TRUE ~ tgroup)) %>%
  mutate(tgroup = gsub(".*__", "", tgroup)) %>%
  arrange(desc(mean))

df <- read.table("data/ex-figure5.pr.stats.of.key.pfams.over.gomc.xls", header = TRUE, sep = "\t", quote = "", comment.char = "") %>%
  mutate(cont = case_when(r2 >= 0.5 ~ "R2 >= 0.5", r2 < 0.5 & r2 >= 0.3 ~ "0.3 <= R2 < 0.5", TRUE ~ "R2 < 0.3")) %>%
  mutate(cont = factor(cont, levels = c("R2 < 0.3", "0.3 <= R2 < 0.5", "R2 >= 0.5"))) %>%
  mutate(sig = case_when(fam.fdr < 0.01 ~ '**', fam.fdr >= 0.01 & fam.fdr < 0.05 ~ '*', TRUE ~ ' ')) %>%
  mutate(fam.fdr = format(fam.fdr, digits = 2)) %>%
  mutate(fam.fdr = case_when(grepl("NA", fam.fdr) ~ "", TRUE ~ fam.fdr)) %>%
  mutate(tgroup = case_when(grepl('d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria', tgroup) ~ "d__Bacteria;c__Gammaproteobacteria", TRUE ~ tgroup)) %>%
  mutate(tgroup = gsub(".*__", "", tgroup)) %>%
  filter(tgroup %in% df.gsize$tgroup) %>%
  mutate(tgroup = factor(tgroup, levels = df.gsize$tgroup)) %>%
  mutate(fam = gsub('\\.[0-9]+$', "", fam)) %>%
  mutate(coef.flag = case_when(fam.coefficient > 0 ~ "+", TRUE ~ "-")) %>%
  mutate(coef.flag = factor(coef.flag, levels = c("+", "-")))

p <- ggplot(data = df, aes(x = tgroup, y = fam))+
  geom_tile(aes(fill = cont), color = "black", alpha = 0.7) +
  geom_text(aes(label = sig, alpha = coef.flag), size = 2, vjust = 0.7) +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_alpha_manual(values = c(1, 0)) +
  scale_fill_manual(values = c("white", "#C5AB89", "#B95756")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0),position = "right")+
  theme(axis.text.x = element_text(color="black", angle = -45, hjust = 0.12, vjust = -0.1, size = 5),
        axis.text.y = element_text(color="black", size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1))


if(!dir.exists("figs")) dir.create("figs")
pdf('figs/Extended-Figure5a.pdf', width = 1.5, height = 8, onefile = F)
print(p)
dev.off()