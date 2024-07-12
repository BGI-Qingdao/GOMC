rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 4b. The workflow to identify Pfam domains potentially underpinning genome size expansion
# Part 3. The result of the exploration of associations between genome size and gene copies
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
df.info <- read.table("data/ex-figure4b.tableS1A.update.xls", header = TRUE, sep = "\t", quote = "", comment.char = "")

df.count <- read.table("data/ex-figure4b.counts.of.genes.with.pfams.in.each.genome.xls", header = TRUE, row.names = 1, sep = "\t", quote = "", comment.char = "")
fam <- "PF13360.9"
df.count <- df.count %>% select(all_of(fam)) %>% `colnames<-`(c("count"))
rownames(df.info) <- df.info$Genome

df.count <- cbind(df.info[match(rownames(df.count), rownames(df.info)),], df.count) %>% select(Genome, Genome_size.bp., count) %>% `colnames<-`(c("genome","size","number_of_genes"))

df.rank <- df.count %>% filter(number_of_genes == 0) %>% mutate(group = "= 0") %>% select(genome, size, number_of_genes, group)
df.rank <- rbind(df.rank, df.count %>% filter(number_of_genes == 1) %>% mutate(group = "= 1") %>% select(genome, size, number_of_genes, group))
df.rank <- rbind(df.rank, df.count %>% filter(number_of_genes == 2) %>% mutate(group = "= 2") %>% select(genome, size, number_of_genes, group))
df.rank <- rbind(df.rank, df.count %>% filter(number_of_genes >= 3) %>% mutate(group = ">= 3") %>% select(genome, size, number_of_genes, group))

group.levels <- c("= 0", "= 1", "= 2", ">= 3")

df.mean <- data.frame()

for(group.slt in group.levels){
  df.mean <- rbind(df.mean, data.frame(group = group.slt, value = mean(df.rank %>% filter(group == group.slt) %>% pull(size))))
}

p <- ggplot(data = df.rank, aes(x = factor(group, levels = group.levels), y = size/1000000)) +
  geom_violin(color = "blue", linetype = 2, fill = NA, linewidth = 0.1)+
  geom_point(data = df.mean, aes(x = factor(group, levels = group.levels), y = value/1000000), color = "blue", size = 5) +
  geom_line(data = df.mean, aes(x = factor(group, levels = group.levels), y = value/1000000), color = "blue", linewidth = 1, linetype = 2, group = 1) +
  stat_compare_means(comparisons = list(c("= 1", "= 0")), method = "wilcox.test", method.args = list(alternative = "greater", exact = FALSE), label.y = 16.6, size = 5, color = "black") +
  stat_compare_means(comparisons = list(c("= 2", "= 1")), method = "wilcox.test", method.args = list(alternative = "greater", exact = FALSE), label.y = 17.8, size = 5, color = "black") +
  stat_compare_means(comparisons = list(c(">= 3", "= 2")), method = "wilcox.test", method.args = list(alternative = "greater", exact = FALSE), label.y = 19.2, size = 5, color = "black") +
  ylim(0,21.5) +
  theme_bw() +
  labs(title = "", x = "Number of genes with PF13360.9\ndomain(s) in a genome", y = "Genome size (Mbp)", color = "") +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure4b3.pdf", width = 4, height = 8, onefile = FALSE)
print(p)
dev.off()

