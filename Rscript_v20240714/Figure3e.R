rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 3e. The trend indicates a decrease in the upper limit number of ARGs with an increase in the number of Cas operons.
# =================================================================================================================================

# Load libraries
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(cowplot))

# Read data
df.count <- read.table("data/figure3e.cas_stat.cas_num.txt", header = TRUE, sep = "\t") %>%
  mutate(Cas = factor(Cas, levels = 1:5)) %>%
  mutate(group = case_when(Cas == 1 ~ "group1", TRUE ~ "group2"))

df.msd <- df.count %>%
  group_by(Cas) %>%
  summarise_at(vars(ARG_Strict), list(max = max, mean = mean, sd = sd, median = median, skew = moments::skewness)) %>% 
  as.data.frame() %>%
  mutate(Cas = as.numeric(Cas))
df.msd

# Visualize
set.seed(12345)

p1 <- ggplot() +
  geom_violin(data = df.count, aes(x = Cas, y = ARG_Strict, fill = Cas), position = position_dodge(0.5), linewidth = 0.5) +
  geom_jitter(data = df.count, aes(x = Cas, y = ARG_Strict), size = 0.01, color = "grey") +
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
        axis.title.x = element_text(size = 32, face = "bold"),
        axis.title.y = element_text(size = 32, face = "bold"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.position = "none")

p2 <- ggscatter(df.msd, x = "Cas", y = "max", add = "reg.line", add.params = list(color = "blue", fill = NA, linetype = "dashed")) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 1, label.y = 3) +
  stat_regline_equation(label.x = 1, label.y = 4) +
  labs(title = "", x = "Cas opreons count", y = "Max ARG count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 0.5),
        axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure3e.pdf", width = 7, height = 8, onefile = F)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p1), x = 0, y = 0, width = 1, height = 1) +
  draw_plot(ggplotify::as.ggplot(p2), x = 0.57, y = 0.605, width = 0.41, height = 0.38)
dev.off()