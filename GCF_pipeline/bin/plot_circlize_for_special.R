library(circlize)

args<-commandArgs(T)

outfile <- paste0(args[1],"/",args[2])
pdf(outfile)
circos.clear()
chr <- read.table("Circlize_Input/Chr_len.txt",header = T,sep="\t",stringsAsFactors = F)
chr_num <- dim(chr)[1]
color_r <- rand_color(chr_num, transparency=0.7)
co <- paste(chr[,1],color_r,sep="\t")
write.table(file="Color_Info.txt", co, quote=F, row.names=F, col.names=F)
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

link <- read.table("Circlize_Input/Links_input.txt",header = F, sep = "\t")
chr_u <- unique(link[,1])
for (i in 1:length(chr_u)){
  zj <- link[grep(pattern = chr_u[i],link[,1]),]
  f <- zj[,c(1:3)]
  t <- zj[,c(4:6)]
  circos.genomicLink(f, t, col = color_r[i], lwd = 1)
}

dev.off()
