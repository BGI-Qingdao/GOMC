rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set working directory to the path of currently opened R script.

# =================================================================================================================================
# Code used to produce the figures and analysis of "Global marine microbial diversity and its potential in bioprospecting"
# Extended Figure 4b. The workflow to identify Pfam domains potentially underpinning genome size expansion
# Part 4. The VennDiagram
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
suppressMessages(library(VennDiagram))
suppressMessages(library(ggplotify))

# Read data
df.lg.angst <- read.table("data/ex-figure4b.votes.to.each.pfam.tsv", header = FALSE, row.names = 1, sep = "\t", quote = "") %>% `colnames<-`(c("vote"))
df.lg.pr <- read.table("data/ex-figure4b.counts.of.genes.with.pfams.in.each.genome.vs.gsize.pr.tsv", header = TRUE, row.names = 1, sep = "\t", quote = "")
df.gg.trend <- read.table("data/ex-figure4b.pfam.slt.sig.v2.xls", header = TRUE, row.names = 1, sep = "\t", quote = "")

pfam.angst <- df.lg.angst %>% filter(vote >= 3) %>% rownames(.)
pfam.pr <- df.lg.pr %>% filter(adjusted.r2 >= 0.5 & fam.coefficient > 0 & fam.fdr < 0.05) %>% rownames(.)
pfam.trend <- df.gg.trend %>% filter(fdr10 < 0.05 & fdr21 < 0.05 & fdr32 < 0.05) %>% rownames(.)

myCol <- c("#43B5C3", "#82CCBB", "#C8E3B5")

p <- venn.diagram(
  x = list(pfam.pr, pfam.angst, pfam.trend),
  category.names = c("Set 1", "Set 2", "Set 3"),
  filename = NULL,
  output = TRUE ,
  imagetype="png" ,
  height = 100, 
  width = 100, 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col = "black",
  fill = myCol,
  cex = 3,
  fontfamily = "sans",
  cat.cex = 3,
  cat.default.pos = "outer",
  cat.pos = c(-28, 28, 136),
  cat.dist = c(0.055, 0.06, 0.07),
  cat.fontfamily = "sans",
  cat.col = myCol,
  rotation = 1
)

p <- grobTree(p)
p <- as.ggplot(p)

if(!dir.exists("figs")) dir.create("figs")
pdf("figs/Extended-Figure4b4.pdf", width = 4, height = 4, onefile = FALSE)
print(p)
dev.off()