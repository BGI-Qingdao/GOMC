rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 4a. Correlations between genome features and genome size in the phyla with large genomes
# =========================================================================================================================

# Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# Read data
df.gsize <- read.table("data/ex-figure4a.tableS1A.update.final.xls", header = TRUE, sep = "\t", quote = "", comment.char = "") %>%
  select(Genome, Genome_size.bp., GC..., Completeness..., Contamination..., GTDB_taxonomy.v207., Estimated_genome_size.bp., Gene_number, Gene_length.bp., Average_gene_length.bp., Gene_GC..., Coding_density...) %>%
  `colnames<-`(c("genome", "size", "gc", "completeness", "contamination", "taxonomy", "estimate_size", "gene_number", "gene_length", "average_gene_length", "gene_gc", "coding_density")) %>%
  mutate(tgroup = case_when(! grepl("p__Proteobacteria", taxonomy) ~ gsub(";c__.*", "", taxonomy), 
                            grepl("c__Alphaproteobacteria|c__Gammaproteobacteria", taxonomy) ~ gsub(";o__.*", "", taxonomy), 
                            grepl("p__Proteobacteria", taxonomy) & ! grepl("c__Alphaproteobacteria|c__Gammaproteobacteria", taxonomy) ~ "d__Bacteria;p__Proteobacteria;c__Other")) %>%
  mutate(tgroup = case_when(grepl("d__Bacteria;p__Myxococcota_A", tgroup) ~ "d__Bacteria;p__Myxococcota", TRUE ~ tgroup)) %>%
  select(-taxonomy) %>%
  mutate(size = size / 1000000) %>%
  mutate(gene_length = gene_length / 1000000) %>%
  mutate(intergenic = size - gene_length)

df.large <- table(df.gsize %>% filter(size > 8) %>% select(tgroup))
df.large <- df.large[which(df.large >= 5)]

df.gsize <- df.gsize %>% filter(tgroup %in% names(df.large) & completeness >= 70 & contamination <= 10) %>% mutate(tgroup = gsub("d__Bacteria;.*__","",tgroup))

df.mean <- aggregate(df.gsize$estimate_size, list(df.gsize$tgroup), FUN = mean) %>%
  `colnames<-`(c("tgroup", "mean")) %>%
  arrange(desc(mean))

df.gsize <- df.gsize %>% mutate(tgroup = factor(tgroup, levels = df.mean$tgroup)) 

df.gsize <- df.gsize %>%
  mutate(intergenic = size - gene_length) %>%
  select(tgroup, size, gc, intergenic, average_gene_length, coding_density)

df.long <- tidyr::gather(df.gsize, group, value, 3:6, factor_key = TRUE) %>%
  mutate(group = factor(group, levels = c("gc", "intergenic", "average_gene_length", "coding_density")))

# Visualize
p <- ggplot(data = df.long, aes(x = size, y = value)) +
  geom_point(color = "grey", alpha = 0.8, size = 0.1) +
  geom_smooth(color = "forestgreen", fill = scales::muted("green")) +
  facet_grid(group ~ tgroup, scales = "free_y", labeller = labeller(group = c(gc = "GC (%)", intergenic = "Intergenic\nRegion (Mb)", average_gene_length = "Average Gene\nLength (base)", coding_density = "Coding Density"))) +
  theme_bw() +
  labs(x = "Genome Size (Mb)", y = "", title = "") +
  theme(strip.text.x = element_text(size= 18), 
        strip.text.y = element_text(size= 20), 
        axis.title.x = element_text(size = 20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_line(linewidth = 0.1),
        axis.text = element_text(size = 16),
        panel.background = element_rect(fill = 'transparent', color = NA),
        panel.border = element_rect(linewidth = 0.1),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.1, linetype = "solid", color = "black"),
        axis.line.y = element_line(linewidth = 0.1, linetype = "solid", color = "black"))

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure4a.pdf", width = 26, height = 10, onefile = FALSE)
print(p)
dev.off()