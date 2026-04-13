library(circlize)
args<-commandArgs(T)
pdf(args[3])
chr <- read.table(args[1],header = T,sep="\t",stringsAsFactors = F)
circos.genomicInitialize(chr,plotType = NULL)

circos.track(
  ylim = c(0, 1), 
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 0.8, col = "#0000FF30")
    circos.text(
      mean(xlim), mean(ylim), chr, cex = 0.7,
      col = "black", facing = "outside",
      niceFacing = TRUE
    )
  }, 
  track.height = 0.05, bg.border = NA
)
link <- read.table(args[2],header = F, sep = "\t")
from <- link[,c(1:3)]
to <- link[,c(4:6)]
# circos.genomicLink(from, to, col = rand_color(1,transparency = 0.7), lwd = 1)
circos.genomicLink(from, to, col = "#f4cccc50", lwd = 1)
dev.off()