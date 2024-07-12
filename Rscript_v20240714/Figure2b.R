rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# ==================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Figure 2b. The heatmap illustrates the distribution of the top 33 functional domains across genomes in the Planctomycetota phylum.
# ==================================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(aplot))
suppressMessages(library(cowplot))

# Read data
# genome size; phyla ordered by mean
df.gsize <- read.table("data/figure2b.large_genome_data_genome_size.txt", header = TRUE, row.names = 1, sep = "\t", quote = "", comment.char = "")
df.gsize <- df.gsize %>% mutate(tgroup = case_when(grepl('d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria', tgroup) ~ "d__Bacteria;c__Gammaproteobacteria", TRUE ~ tgroup))
df.gsize$tgroup <- gsub(".*__", "", df.gsize$tgroup)
df.mean <- aggregate(df.gsize$estimate_size, list(df.gsize$tgroup), FUN = mean) %>%
  `colnames<-`(c("tgroup", "mean")) %>%
  arrange(desc(mean))
df.gsize <- df.gsize %>% mutate(tgroup = factor(tgroup, levels = df.mean$tgroup))

# count of pfams
df.count <- read.table("data/figure2b.large_genome_data_gene_count.txt", header = TRUE, row.names = 1, sep = "\t", quote = "", comment.char = "")
df.count <- df.count %>% filter(rownames(df.count) %in% rownames(df.gsize))
df.count <- df.count[match(rownames(df.gsize), rownames(df.count)),]
colnames(df.count) <- gsub('\\.[0-9]+$', "", colnames(df.count))

# phylogenetic regression stats
df.preg <- read.table("data/figure2b.large_genome_data_key_pfams.txt", header = TRUE, sep = "\t", quote = "", comment.char = "")
df.preg$tgroup <- gsub(".*__", "", df.preg$tgroup)
df.preg$fam <- gsub('\\.[0-9]+$', "", df.preg$fam)

# Extract data of Planctomycetota
df.gsize <- df.gsize %>% filter(tgroup == "Planctomycetota") %>%
  arrange(size) %>%
  mutate(genome = factor(rownames(.), levels = rownames(.)))

df.count <- df.count[rownames(df.gsize),]

df.preg <- df.preg %>% filter(tgroup == "Planctomycetota") %>%
  arrange(desc(r2)) %>%
  filter(r2 >= 0.3) %>%
  mutate(sig = case_when(fam.fdr < 0.01 ~ '**', fam.fdr >= 0.01 & fam.fdr < 0.05 ~ '*', TRUE ~ ' ')) 

# Visualize
genome.order <- rownames(df.gsize)
fam.order <- rev(df.preg$fam)

df.var <- log2(df.count + 1)
df.var <- df.var %>% select(all_of(fam.order)) %>% mutate(genome = factor(rownames(.), levels = genome.order))
df.var.long <- tidyr::gather(df.var, pfam, count, 1:length(fam.order), factor_key = TRUE)
df.var.long$pfam <- factor(df.var.long$pfam, levels = fam.order)

p1 <- ggplot(df.var.long, aes(genome, pfam, fill = count)) + 
  geom_tile(alpha = 1) +
  scale_fill_gradient2(low = "#5590CC", mid = "white", high = "#F6999A", midpoint = median(df.var.long$count), breaks = c(0, 2.5, 5, 7.5)) +
  labs(title = "", x = "", y = "", fill = "log2(Gene Counts)") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.1),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 21),
        legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24, face = "bold", angle = 90),
        legend.title.position = "left",
        legend.key.height = unit(1.6, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.ticks = element_line(color = "black"),
        legend.background = element_blank(),
        legend.box.background = element_blank())

p2 <- ggplot(df.preg, aes(fam, r2)) +
  geom_bar(stat = "identity", color = "white", fill = "#c9c0d3", alpha = 1, width = 0.8, linewidth = 0) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.2,0.4,0.6)) +
  coord_flip() +
  labs(title = "", x = "", y = bquote(""~bolditalic(R)^2), fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(linewidth = 1),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 1),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.2, "cm"),
        axis.text.x = element_text(size = 20, angle = -90, vjust = 0, hjust = 0),
        axis.text.y = element_blank(),
        legend.position = "none")

p <- p1 %>% insert_right(p2, width = 0.06)

p3 <- ggplot(df.gsize, aes(genome, size)) +
  geom_bar(stat = "identity", color = "#d4b8b4", fill = "#d4b8b4", alpha = 1, width = 1, linewidth = 0) +
  scale_y_continuous(expand = c(0,0), breaks = c(5,10,15,20), limits = c(0,20), position = "right") +
  labs(title = "", x = "", y = "Genome Size\n(Mb)", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 1),
        axis.title.x = element_blank(),
        axis.title.y.right = element_text(size = 24, face = "bold", angle = 90),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 1),
        axis.ticks.length.y = unit(0.2, "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        legend.position = "none")

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Figure2b.pdf", width = 25, height = 20, onefile = F)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p), x = 0, y = 0.02, width = 1, height = 0.85) +
  draw_plot(ggplotify::as.ggplot(p3), x = 0.053, y = 0.855, width = 0.879, height = 0.13)
dev.off()