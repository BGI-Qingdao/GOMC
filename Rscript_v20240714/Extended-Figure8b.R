rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 8b. Comparison of BGCs in GOMC against BIG-FAM database. 
# =========================================================================================================================

# Load libraries
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# Read data
df.count <- read.table("data/ex-figure8b.novetly-phylum_input.txt", header = T, sep = "\t")

Phylum.order <- df.count %>%
  group_by(Phylum) %>%
  summarise(sum = sum(BGC_number)) %>%
  arrange(desc(sum)) %>%
  pull(Phylum) %>%
  unique()
Phylum.order <- Phylum.order[which(Phylum.order != "other")]
Phylum.order <- c(Phylum.order, "other")

df.count <- df.count %>%
  mutate(Phylum = factor(Phylum, levels = Phylum.order))

# Visulaize
p <- ggplot(data = df.count, aes(x = Cos_dist, y = BGC_number, fill = Phylum)) +
  geom_histogram(stat = "identity") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = c("0.0", "(0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]", "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1.0]")) +
  labs(x = "Minimum cosine distance of BGCs", y = "BGC number", fill = "") +
  scale_fill_manual(values = c(pals::brewer.paired(12), "grey")) +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.ticks.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure8b.pdf", width = 10, height = 4)
print(p)
dev.off()