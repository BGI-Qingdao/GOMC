rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =====================================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 5b. Distribution of Pfam domains within each phylum as genome size increases.
# Only Pfams with R² ≥ 0.5 from the regression analyses are shown for each phylum.
# =====================================================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(vegan))
suppressMessages(library(ape))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(aplot))
suppressMessages(library(cowplot))

# Read data
df.gsize <- read.table("data/ex-figure5.pr.subset.genome.size.of.phyla.with.large.genomes.xls", header = TRUE, row.names = 1, sep = "\t", quote = "", comment.char = "") %>%
  mutate(tgroup = case_when(grepl('d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria', tgroup) ~ "d__Bacteria;c__Gammaproteobacteria", TRUE ~ tgroup)) %>%
  mutate(tgroup = gsub(".*__", "", tgroup))
df.mean <- aggregate(df.gsize$estimate_size, list(df.gsize$tgroup), FUN = mean) %>%
  `colnames<-`(c("tgroup", "mean")) %>%
  arrange(desc(mean))
df.gsize <- df.gsize %>% mutate(tgroup = factor(tgroup, levels = df.mean$tgroup))
df.count <- read.table("data/ex-figure5.pr.subset.counts.of.genes.with.pfams.in.each.genome.xls", header = TRUE, row.names = 1, sep = "\t", quote = "", comment.char = "") %>%
  filter(rownames(.) %in% rownames(df.gsize))
df.count <- df.count[match(rownames(df.gsize), rownames(df.count)),]
colnames(df.count) <- gsub('\\.[0-9]+$', "", colnames(df.count))

df.preg <- read.table("data/ex-figure5.pr.stats.of.key.pfams.over.gomc.xls", header = TRUE, sep = "\t", quote = "", comment.char = "") %>%
  filter(r2 >= 0.5 & fam.coefficient > 0 & fam.fdr < 0.05)
df.preg$tgroup <- gsub(".*__", "", df.preg$tgroup)
df.preg$fam <- gsub('\\.[0-9]+$', "", df.preg$fam)

# Visualize
p.bar <- list()
p.hmap <- list()

for(i in 1:9){
  tgroup.slt <- levels(df.gsize$tgroup)[i]
  
  print(paste(i, tgroup.slt, sep = "|"))
  
  fam.keep <- df.preg %>% filter(tgroup == tgroup.slt) %>% arrange(r2) %>% select(fam) %>% unlist %>% unique
  
  df.var <- df.count[rownames(df.gsize %>% filter(tgroup == tgroup.slt)),] %>% select(all_of(fam.keep))
  
  genome.order <- df.gsize %>% filter(tgroup == tgroup.slt) %>% arrange(size) %>% rownames(.)
  
  df.var <- log2(df.var + 1)
  df.var <- df.var %>% mutate(genome = factor(rownames(.), levels = genome.order))
  df.var.long <- tidyr::gather(df.var, pfam, count, 1:length(fam.keep), factor_key = TRUE)
  
  df.res <- df.gsize %>%
    filter(tgroup == tgroup.slt) %>%
    mutate(genome = factor(rownames(.), levels = genome.order))
    
  p.bar[[i]] <- ggplot(df.res, aes(genome, size)) +
    geom_bar(stat = "identity", color = "#d4b8b4", fill = "#d4b8b4", alpha = 1, width = 1, linewidth = 0) +
    scale_y_continuous(expand = c(0,NA), breaks = c(5,10,15,20), limits = c(0,20)) +
    labs(title = tgroup.slt, x = "", y = "Size\n(Mb)", fill = "") +
    theme_bw() +
    theme() +
    theme(plot.title = element_text(size = 30, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 24),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth = 0.01),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1),
          legend.position = "none")
  
  p.hmap[[i]] <- ggplot(df.var.long, aes(genome, pfam, fill = count)) + 
    geom_tile(alpha = 1) +
    scale_fill_gradient2(low = "#5590CC", mid = "white", high = "#F6999A", midpoint = median(df.var.long$count)) +
    labs(title = "", x = "", y = "", fill = "log2(Gene Counts)") +
    theme_bw() +
    theme() +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1),
          legend.position = "none")
  
  if(i == 9) p.hmap[[i]] <- p.hmap[[i]] + theme(legend.position = "bottom")
}

if(!dir.exists("figs")) dir.create("figs")
pdf('figs/Extended-Figure5b.pdf', width = 20, height = 35, onefile = F)
ggdraw() +
  draw_plot(ggplotify::as.ggplot(p.hmap[[9]]), x = 0, y = 0, width = 1, height = 0.1) +
  draw_plot(ggplotify::as.ggplot(p.bar[[9]]), x = 0.014, y = 0.096, width = 0.986, height = 0.04) +
  draw_plot(ggplotify::as.ggplot(p.hmap[[8]]), x = 0, y = 0.135, width = 1, height = 0.11) +
  draw_plot(ggplotify::as.ggplot(p.bar[[8]]), x = 0.014, y = 0.24, width = 0.986, height = 0.04) +
  draw_plot(ggplotify::as.ggplot(p.hmap[[7]]), x = 0, y = 0.28, width = 1, height = 0.035) +
  draw_plot(ggplotify::as.ggplot(p.bar[[7]]), x = 0.014, y = 0.31, width = 0.986, height = 0.04) +
  draw_plot(ggplotify::as.ggplot(p.hmap[[6]]), x = 0, y = 0.35, width = 1, height = 0.03) +
  draw_plot(ggplotify::as.ggplot(p.bar[[6]]), x = 0.014, y = 0.375, width = 0.986, height = 0.04) +
  draw_plot(ggplotify::as.ggplot(p.hmap[[5]]), x = 0, y = 0.415, width = 1, height = 0.095) +
  draw_plot(ggplotify::as.ggplot(p.bar[[5]]), x = 0.014, y = 0.505, width = 0.986, height = 0.04) +
  draw_plot(ggplotify::as.ggplot(p.hmap[[4]]), x = 0, y = 0.545, width = 1, height = 0.07) +
  draw_plot(ggplotify::as.ggplot(p.bar[[4]]), x = 0.014, y = 0.61, width = 0.986, height = 0.04) +
  draw_plot(ggplotify::as.ggplot(p.hmap[[3]]), x = 0, y = 0.65, width = 1, height = 0.11) +
  draw_plot(ggplotify::as.ggplot(p.bar[[3]]), x = 0.014, y = 0.755, width = 0.986, height = 0.04) +
  draw_plot(ggplotify::as.ggplot(p.hmap[[1]]), x = 0, y = 0.8, width = 1, height = 0.16) +
  draw_plot(ggplotify::as.ggplot(p.bar[[1]]), x = 0.014, y = 0.955, width = 0.986, height = 0.04)
dev.off()