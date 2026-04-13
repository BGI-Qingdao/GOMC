#!/dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/envs/gcf/bin/Rscript
library(ggplot2)
library(scales)
data <- read.table("./Variance_Info.txt",header = T,sep="\t")
ggplot(data,aes(variance,Level))+
  geom_boxplot(position = position_nudge(y=0.05),width=0.1,outlier.size = 1.5,outlier.shape = 18)+
  stat_summary(fun="mean", geom="point", shape=20, size=2.5, color="red", fill="red",alpha=0.7,position = position_nudge(y=0.05))+
  geom_point(position = position_nudge(y=-0.05),size=0.5)+
  scale_x_log10(breaks = 10^(0:7),labels = trans_format("log10",math_format(10^.x)))+
  scale_y_discrete(limits = c("Genera","Families","Orders","Classes","Phyla"))+
  theme_bw()+
  labs(x="Variance of Biosynthetic diversity",y=NULL)
ggsave("./Variance_Boxplot.pdf")
ggsave("./Variance_Boxplot.png")
