rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =====================================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 3c and d. Quantification of editing efficiency for five selected editing sites of HBG gene (c) and BCL11a enhancer (d), respectively. 
# =====================================================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))

# Read data
df.hbg <- read.table("data/ex-figure3.crispr_cas.HBG.indel_freq.txt", header= TRUE, sep = "\t") %>%
  mutate(Group = factor(Group, levels = c("T", "M")))

df.bcl <- read.table("data/ex-figure3.crispr_cas.BCL11a.indel_freq.txt", header= TRUE, sep = "\t") %>%
  mutate(Group = factor(Group, levels = c("T", "M")))

# Visualize
my_comparison <- list(c("T", "M"))

p.hbg <- ggplot(data = df.hbg, aes(x = Group, y = Indel)) +
  geom_bar(aes(fill=Group), stat = "summary", fun = "mean", position = position_dodge()) +
  stat_summary(fun.data = "mean_sd", geom = "errorbar", colour = "black", width = 0.5, position = position_dodge(.9)) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = my_comparison) +
  # stat_compare_means(method = "t.test", label = "p.signif", comparisons = my_comparison, method.args = list(alternative = "greater")) +
  scale_fill_manual(values = c(T ="#CFB6D8", M = "#D2D2D2"), labels = c("Treatment", "Negative control")) +
  facet_wrap(~Site, ncol = 5) +
  labs(title = "", x = "", y = "HBG Indels (%)", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x  = element_text(size = 14, family = "serif"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 1),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "top",
        legend.text = element_text(size = 14))

p.bcl <- ggplot(data = df.bcl, aes(x = Group, y = Indel)) +
  geom_bar(aes(fill=Group), stat = "summary", fun = "mean", position = position_dodge()) +
  stat_summary(fun.data = "mean_sd", geom = "errorbar", colour = "black", width = 0.5, position = position_dodge(.9)) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = my_comparison) +
  # stat_compare_means(method = "t.test", label = "p.signif", comparisons = my_comparison, method.args = list(alternative = "greater")) +
  scale_fill_manual(values = c(T ="#CFB6D8", M = "#D2D2D2"), labels = c("Treatment", "Negative control")) +
  facet_wrap(~Site, ncol = 5) +
  labs(title = "", x = "", y = "BCL11a Indels (%)", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x  = element_text(size = 14, family = "serif"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 1),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "top",
        legend.text = element_text(size = 14))

if(!dir.exists("figs")) dir.create("figs")
pdf('figs/Extended-Figure3cd.pdf', width = 10, height = 12, onefile = F)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p.hbg), x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(ggplotify::as.ggplot(p.bcl), x = 0, y = 0, width = 1, height = 0.5)
dev.off()