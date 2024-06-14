rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 3f. The trend indicates a decrease in the upper limit number of ARGs with an increase in the number of Cas operons.
# =================================================================================================================================

# Load libraries
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

# Read data
df.count <- read.table("data/figure3f.cas_stat.cas_num.txt", header = TRUE, sep = "\t") %>%
  mutate(Cas = factor(Cas, levels = 1:5))

# Visualize
set.seed(12345)

p <- ggplot(data = df.count, aes(x = Cas, y = ARG_Strict, fill = Cas)) +
  geom_violin(position = position_dodge(0.5), linewidth = 0.5) +
  geom_boxplot(width = 0.1, position = position_dodge(0.5), linewidth = 0.5) +
  stat_compare_means(method = "kruskal.test", label.x = 2, label.y = 3.8, size = 6) +
  scale_x_discrete(breaks = 1:5, labels = c(1:4, ">=5")) +
  scale_y_continuous(trans = log2_trans(), breaks = c(1,2,4,8,16), labels = c(1,2,4,8,">=16")) +
  scale_fill_manual(values = c("1" = "#FFEFCC", "2" = "#FECCCC", "3" = "#99CCCC", "4" = "#CCCCFF", "5" = "#3399CC"))+
  labs(title = "", x = "Cas opreons count", y = "ARG count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold", vjust = -5, hjust = 0.5),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none")

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure3f.pdf", width = 6, height = 7, onefile = F)
print(p)
dev.off()