rm(list=ls())

##普通制图
setwd("E:/Xiong/data/20170706/final/distance")
table <- read.csv(file ='identification_rate.csv')

table2 <- as.matrix(table) ##原始表格需转化为矩阵，要求height为向量或矩阵
rownames(table2) <- c("tree","distance","BLAST","character")
colnames(table2) <- c("A","M","P","Y","AM","AP","AY","MP","MY","PY","AMP","AMY","APY","MPY","AMPY")

##ggplot物种鉴定率
library(ggplot2)
table <- read.csv(file ='resolution.csv')
table$Marker=factor(table$Marker, levels=c("A","M","T","Y","AM","AT","AY","MT","MY","TY","AMT","AMY","ATY","MTY","AMTY"))
table$method=factor(table$method, levels=c("BLOG","Tree","Distance","BLAST"))
discrimination_rate_plot <-ggplot(table) + 
  geom_bar(aes(x=Marker,y=value,fill=method),stat="identity",position="dodge")+ 
  scale_fill_brewer(palette="gray") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(y="Discrimination rate (%)")
ggsave("discrimination_rate_plot.pdf",plot=discrimination_rate_plot)
  
rm(list=ls())
##ggplot barcoding gap
##数据准备
##提取种间距离，去除NA  
  dis_mat_int <- read.csv(file ='matK_inter_distance2.csv')
  dis_mat_int<-na.omit(dis_mat_int)
  type<- rep("inter",69557)
  dis_mat_int<-cbind(dis_mat_int,type) ### 2列合并，形成新的数据框
  colnames(dis_mat_int) <-c("distance","type") 
##提取种内距离，去除NA  
  dis_mat_itr <- read.csv(file ='matK_intra_distance.csv')
  dis_mat_itr<-na.omit(dis_mat_itr)
  type<- rep("intra",1318)
  dis_mat_itr<-cbind(dis_mat_itr,type)
  colnames(dis_mat_itr) <-c("distance","type")
##合并种内、种间数据框
  dis_mat <- rbind(dis_mat_int, dis_mat_itr)
  write.csv(dis_mat,"matK_distance2.csv")
##统计范围内的频数，写入终表格
 table_mat <- table(dis_mat)
 write.csv(table_mat,"table_mat.csv") ##计算概率
table_mat2<-read.csv("table_mat2.csv")

dis_psb_int <- read.csv(file ='psbA_inter_distance2.csv')
dis_psb_int<-na.omit(dis_psb_int)
type<- rep("inter",62424)
dis_psb_int<-cbind(dis_psb_int,type)
colnames(dis_psb_int) <-c("distance","type")


dis_psb_itr <- read.csv(file ='psbA_intra_distance.csv')
dis_psb_itr<-na.omit(dis_psb_itr)
type<- rep("intra",1139)
dis_psb_itr<-cbind(dis_psb_itr,type)
colnames(dis_psb_itr) <-c("distance","type")

table_psb <- table(dis_psb)
write.csv(table_psb,"table_psb.csv") ##计算概率
table_psb2<-read.csv("table_psb2.csv")

dis_ycf_int <- read.csv(file ='ycf_inter_distance2.csv')
dis_ycf_int<-na.omit(dis_ycf_int)
type<- rep("inter",64473)
dis_ycf_int<-cbind(dis_ycf_int,type)
colnames(dis_ycf_int) <-c("distance","type")


dis_ycf_itr <- read.csv(file ='ycf_intra_distance.csv')
dis_ycf_itr<-na.omit(dis_ycf_itr)
type<- rep("intra",1229)
dis_ycf_itr<-cbind(dis_ycf_itr,type)
colnames(dis_ycf_itr) <-c("distance","type")

dis_ycf <- rbind(dis_ycf_int, dis_ycf_itr)
write.csv(dis_ycf,"ycf_distance2.csv")


table_ycf <- table(dis_ycf)
write.csv(table_ycf,"table_ycf.csv") ##计算概率
table_ycf2<-read.csv("table_ycf2.csv")
rm(list=ls())


##制图
##彩色拼合图
library(ggplot2)
table_atp2<-read.csv("table_atp2.csv")
p1 <- ggplot(table_atp2, aes(x=distance,frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="atpI-atpH") +
  xlim(-0.001,0.03)+
  ylim(0,0.6)

table_mat2<-read.csv("table_mat2.csv")
p2 <- ggplot(table_mat2, aes(x=distance,frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="matK") +
  xlim(-0.001,0.03)+
  ylim(0,0.6)

table_psb3<-read.csv("table_psb3.csv")
p3 <- ggplot(table_psb3, aes(x=distance,y=frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="trnH-psbA") +
  xlim(-0.001,0.06)+
  ylim(0,0.6)


table_ycf2<-read.csv("table_ycf2.csv")
p4 <- ggplot(table_ycf2, aes(x=distance,y=frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="ycf1") +
  xlim(-0.001,0.03)+
  ylim(0,0.6)

library(gridExtra)
bar.gap.colour<-grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
ggsave("barcode_gap_colour.pdf",bar.gap.colour)

##黑白拼合图
table_atp2<-read.csv("table_atp2.csv")
p1 <- ggplot(table_atp2, aes(x=distance,frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  scale_fill_manual(values=c("grey","black")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="atpI-atpH") +
  xlim(-0.001,0.03)+
  ylim(0,0.6)

table_mat2<-read.csv("table_mat2.csv")
p2 <- ggplot(table_mat2, aes(x=distance,frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  scale_fill_manual(values=c("grey","black")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="matK") +
  xlim(-0.001,0.03)+
  ylim(0,0.6)

table_psb3<-read.csv("table_psb3.csv")
p3 <- ggplot(table_psb3, aes(x=distance,y=frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  scale_fill_manual(values=c("grey","black")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="trnH-psbA") +
  xlim(-0.001,0.06)+
  ylim(0,0.6)


table_ycf2<-read.csv("table_ycf2.csv")
p4 <- ggplot(table_ycf2, aes(x=distance,y=frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  scale_fill_manual(values=c("grey","black")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="ycf1") +
  xlim(-0.001,0.03)+
  ylim(0,0.6)

library(gridExtra)
bargap.black<-grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
ggsave("barcoding_gap_black.pdf",bargap.black)

##自动组合图
##黑白
table_cob <- read.csv("table_combination.csv")
ggplot(table_cob, aes(x=distance,y=frequency)) + 
  geom_bar(aes(color=type),alpha=.6,position="identity",stat="identity") +
  scale_color_manual(values=c("white","black","grey"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  facet_wrap(~marker,nrow=2)

##黑白选用
ggplot(table_cob, aes(x=distance,y=frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  scale_fill_manual(values=c("grey","black")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  facet_wrap(~marker,nrow=2)

##彩色
ggplot(table_cob, aes(x=distance,y=frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  facet_wrap(~marker,nrow=2)

ggplot(table_cob, aes(x=distance,y=frequency)) + 
  geom_bar(aes(fill=type),alpha=.6,position="identity",stat="identity") +
  scale_fill_manual(values=c("lightgreen","seagreen")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  facet_wrap(~marker,nrow=2)



colors()

##查看颜色集
library(RColorBrewer)
display.brewer.all()

##黑白模式
+scale_colour_grey()
+scale_color_manual(values=c("white","black","grey"))

#### 其他命令
ggplot(dis_atp3, aes(x=distance,..ncount..)) + 
    geom_histogram(aes(fill=type),binwidth=0.0005,alpha=.6,position="identity") +
    stat_density(geom = 'line', position = 'identity',aes(colour = type)) +
    theme_bw()
    geom_vline(aes(xintercept=0.01), linetype="dashed", size=1) 

  
   ggplot(dis_atp3, aes(x=distance,y=..count../sum(..count..))) + 
    geom_histogram(aes(fill=type),bins=200,alpha=.8,position="identity") +
     scale_y_continuous(breaks = "pp", labels = "percent")
     theme_bw()
     install.packages('dplyr')
  library()
  ggplot(mpg, aes(class)) +
    geom_bar(aes(fill=drv), stat="identity",position="dodge")
  mpg

par(mfrow=c(1,1))
barplot(table2,beside=T)
legend(x=-1,y=35,legend=rownames(table2),ylab="discrimination rate")
barplot(table2,beside=T,legend.text=T)


legend("topleft",legend=rownames(table2))

#####
setwd("E:/Xiong/data/20170706/final/distance")
dis_cyc <- read.csv(file ='4ampy_cyc.csv', header = FALSE, sep = ",") ##设置文件名
dis_ilex <- read.csv(file ='4ampy_ilex.csv', header = FALSE, sep = ",")
dis_cer <- read.csv(file ='4ampy_cer.csv', header = FALSE, sep = ",")
dis_que <- read.csv(file ='4ampy_que.csv', header = FALSE, sep = ",")
dis_lob <- read.csv(file ='4ampy_lob.csv', header = FALSE, sep = ",")

dis_cyc2 <- unlist(dis_cyc[,2:147])
dis_ilex2 <- unlist(dis_ilex[,2:125])

boxplot(unlist(dis_cyc[,2:147]),unlist(dis_ilex[,2:125]),unlist(dis_cer[,2:18]),unlist(dis_que[,2:11]),unlist(dis_lob[,2:7]),width=rep(0.3,5),ylim=c(0,0.020),space = rep (0.5,5),na.action = TRUE)
boxplot(unlist(dis_ilex[,2:125]),na.rm = TRUE)
boxplot(unlist(dis_cer[,2:18]),na.rm = TRUE)
boxplot(unlist(dis_que[,2:11]),na.rm = TRUE)
boxplot(unlist(dis_lob[,2:7]),na.rm = TRUE)
rm(mylist)
mylist <- list()
mylist[[1]] <- unlist(dis_ilex[,2:125])
name <- c("cyc","ilex")
names(mylist[[1]]) <- name[1]
mylist[[2]] <- unlist(dis_cyc[,2:147])
names(mylist[[2]]) <- ilex
mylist <- unlist(dis_ilex[,2:125])
boxplot(mylist[[1]],mylist[[2]],ylim=c(0,0.020),space = rep (0.5,5),na.action = TRUE)
ggplot(data = data2) + geom_violin() + geom_boxplot(width = 0.3, outlier.colour = NA, fill = 'blue') + stat_summary(fun.y = 'median', geom = 'point', shape = 18, colour = 'orange')


num <- 7
summary(unlist(distance[,2:num]),na.rm = TRUE)
sd(unlist(distance[,2:num]),na.rm = TRUE)
write.csv(table, file = "species_list.csv")


####################
ggplot各组距离箱线图
dis_que <- read.csv(file ='4ampy_que2.csv', header = FALSE, sep = ",")
dis_que<-na.omit(dis_que)
write.csv(dis_que,"4ampy_que3.csv")

dis_cyc <- read.csv(file ='4ampy_cyc2.csv', header = FALSE, sep = ",")
dis_cyc<-na.omit(dis_cyc)
write.csv(dis_cyc,"4ampy_cyc3.csv")

dis_ilex <- read.csv(file ='4ampy_ilex.csv', header = FALSE, sep = ",")
dis_ilex<-na.omit(dis_ilex)
write.csv(dis_ilex,"4ampy_ilex3.csv")

dis_lob <- read.csv(file ='4ampy_lob2.csv', header = FALSE, sep = ",")
dis_lob<-na.omit(dis_lob)
write.csv(dis_lob,"4ampy_lob3.csv")

distance <- read.csv("4ampy_distance.csv",sep=",")
colnames(distance) <- c("distance","groups")

ggplot(distance, aes(groups, distance)) + 
  geom_violin(fill="lightblue") +
  geom_boxplot(width=.05,fill="white",alpha=.8) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(x="sections", y="p-distance")
ggsave("p_dis_sections.pdf")
