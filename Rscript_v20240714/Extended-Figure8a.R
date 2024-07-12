rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =========================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 8a. BGCs predicted from GOMC genomes.
# =========================================================================================================================

# Load libraries
suppressMessages(library(ggplot2))
suppressMessages(library(aplot))

# Read data
df.count <- read.table("data/ex-figure8a.barplot_input.txt", header = T, sep = "\t")

Phylum.order <- df.count %>%
  group_by(Phylum) %>%
  summarise(sum = sum(Number)) %>%
  arrange(desc(sum)) %>%
  pull(Phylum) %>%
  unique()
Phylum.order <- Phylum.order[which(Phylum.order != "Other_Bacteria")]
Phylum.order <- c(Phylum.order, "Other_Bacteria")

Subtype.order <- df.count %>%
  group_by(Subtype) %>%
  summarise(sum = sum(Number)) %>%
  arrange(desc(sum)) %>%
  pull(Subtype) %>%
  unique()
Subtype.order <- Subtype.order[which(Subtype.order != "others")]
Subtype.order <- c(Subtype.order, "others")

df.count <- df.count %>%
  mutate(Phylum = factor(Phylum, levels = Phylum.order)) %>%
  mutate(Subtype = factor(Subtype, levels = Subtype.order))

df.gsize <- read.table("data/ex-figure8a.boxplot_mag_input.txt", header = T, sep = "\t") %>%
  mutate(Subtype = factor(Subtype, levels = Subtype.order))

df.length <- read.table("data/ex-figure8a.boxplot_bgc_input.txt",header = T, sep="\t") %>%
  mutate(Subtype = factor(Subtype, levels = Subtype.order))

# Visualize
p1 <- ggplot(data = df.count, aes(y = Subtype, x = Number, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  labs(title = "", x = "BGC Number", y = NULL) +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 20))

p2 <- ggplot(data = df.length, aes(y = Subtype, x = BGC_length)) +
  geom_boxplot(fill = "#376FB4", outlier.size = 0.5) +
  scale_x_log10(breaks = c(1, 10, 100)) +
  labs(title = "", x = "BGC Length (Kb)", y = NULL) +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 1),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank())

p3 <- ggplot(data = df.gsize, aes(y = Subtype, x = MAG_size)) +
  geom_boxplot(fill = "#376FB4", outlier.size = 0.5) +
  scale_x_log10(breaks = c(0.3, 1.0, 3.0, 10.0)) +
  labs(title = "", x = "Genome Size (Mb)", y = NULL) +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 1),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank())

p <- p1 %>% insert_right(p2, width = 0.25)
p <- p %>% insert_right(p3, width = 0.25)

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-figure8a.pdf", width = 14, height = 10)
print(p)
dev.off()