# 加载包
library(iNEXT)
library(ggplot2)

#传参
args<-commandArgs(T)

options(scipen=200)

# 创建输入矩阵
data <- list()
a <- list.files(paste0(args[1],"/iNEXT_Input"),pattern=".txt")
name <- gsub(".txt$","",a)
dir <- paste0(args[1],"/iNEXT_Input/",a)
for (i in 1:length(a)){
	data[[i]] <- read.table(dir[i], header = T, as.is = F)
}
names(data)<-name
# 运行iNEXT
out <- iNEXT(data,datatype="incidence_raw",endpoint=as.numeric(args[2]),nboot=as.numeric(args[3]))
save(out, file=paste0(args[1],"/iNEXT.rds"))
# 绘图
p <- ggiNEXT(out, type=1, se = args[4])+
theme_bw()+
labs(x="Genomes",y="GCFs")+
theme(legend.position = "none")
ggsave(paste0(args[1],"/Rarefaction_result.pdf"), p)

# 输出数据文件
for (i in 1:length(out$iNextEst)){
  write.table(out$iNextEst[i],file=paste0(args[1],"/temp/",names(out$iNextEst)[i],".txt"),quote = F,sep = "	",row.names = F)
}