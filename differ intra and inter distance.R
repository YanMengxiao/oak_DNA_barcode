library(stringr)

rm(list=ls())
##文件读入
setwd("E:/Xiong/data/20170706/final/distance") ##设置文件夹
distance <- read.csv(file ='4ampy_lob.csv', header = FALSE, sep = ",") ##设置文件名
##定义变量
species <- c()
row_start <- 1 ##设定第一个种起始为第一行
row_end <- c()
name <- "[A-Z].*[a-z]"  ##截取种名部分

##查找每个种的起始行号
num <- 7
for (i in 1:num)
  if (str_match(distance[i,1], name) == str_match(distance[i+1,1], name)) {
    m <- i
  } else {
    species <- append(species, str_match(distance[i,1], name))
    row_end <- append(row_end,i)
    row_start <-append(row_start,i+1)
  }
row_end <- append(row_end,num) ##设定最后一个种为end最后一行
species <- append(species, str_match(distance[num,1], name)) 


species
row_start
row_end

##写入数据框
table <- data.frame(species,row_start,row_end)
table
write.csv(table, file = "total_specieslist.csv")

##
num <- 7
summary(unlist(distance[,2:num]),na.rm = TRUE)
sd(unlist(distance[,2:num]),na.rm = TRUE)
write.csv(table, file = "species_list.csv")

##种内距离
##读取各物种种内距离，提取数值，存成总表
spe_num <- 7
intra_dist <- list()
intra_theta_full <- data.frame()
coales_depth_full <- data.frame()
for (i in 1:spe_num) ## 设置物种个数
  {
  SN <- table[i,1]  ##物种名
  RS <- table[i,2]  ##行起始数字
  RE <- table[i,3]  ##行结束数字
  CS <- RS+1  ##列起始数字，多名字列故加一
  if (RE<num)   {CE <- RE+1} ##列结束数字，因距离表格最后一列空行未能读入，最后一列不加一
  else {CE <- RE}
  intra_dist[[i]] <- as.numeric(unlist(str_match_all(unlist(distance[RS:RE,CS:CE]),"0.*"))) ##提取距离矩阵，去列表；提取含0数字，去列表；转化为数值型
  intra_theta_full[i,1] <- mean(intra_dist[[i]],na.rm = TRUE)
  if (max(intra_dist[[i]],na.rm = TRUE) == -Inf) {intra_dist[[i]] <- 0}
  else {coales_depth_full[i,1] <- max(intra_dist[[i]],na.rm = TRUE)}
   ##将精简后的各物种向量加入数据框
   # fliename <- paste(table[i,1],".xlsx",sep="") ##粘结文件名和扩展名，合成完整文件名
  #  namepath[[i]] <- paste("E:/Xiong/data/20170706/final/distance/",fliename,sep="")  ##粘结存储路径和完整文件名
  #  write.xlsx(intra_dist[[i]],namepath[[i]]) ## write.xlsx(x, file），列表实现批量操作
}

summary(unlist(intra_dist),na.rm = TRUE)
mean(unlist(intra_dist),na.rm = TRUE)
avg_intra_dist <- mean(unlist(intra_dist), na.rm = TRUE)
sd(unlist(intra_dist), na.rm = TRUE)
avg_intra_dist


rownames(intra_theta_full) <- table[,1]
summary(intra_theta_full[,1],na.rm = TRUE)
sd(intra_theta_full[,1],na.rm = TRUE)

rownames(coales_depth_full) <- table[,1]
summary(coales_depth_full[,1],na.rm = TRUE)
sd(coales_depth_full[,1],na.rm = TRUE)

write.csv(unlist(intra_dist), file = "atpIH_intra_distance.csv")


## 将每个种的种内距离存成xlsx子表格
intra_dist <- list()
namepath <- list()
for (i in 1:spe_num) {
  SN <- table[i,1]  ##物种名
  RS <- table[i,2]  ##行起始数字
  RE <- table[i,3]  ##行结束数字
  CS <- RS+1  ##列起始数字，多名字列故加一
  if (RE<num)   {CE <- RE+1} ##列结束数字，因距离表格最后一列空行未能读入，最后一列不加一
  else {CE <- RE}
  intra_dist[[i]]  <- distance[RS:RE,CS:CE] ##将子表格加入数据框
  # fliename <- paste(table[i,1],".xlsx",sep="") ##粘结文件名和扩展名，合成完整文件名
  # namepath[[i]] <- paste("E:/Xiong/data/20170706/final/distance/",fliename,sep="")  ##粘结存储路径和完整文件名
  # write.xlsx(intra_dist[[i]],namepath[[i]]) ## write.xlsx(x, file），列表实现批量操作
}


##种间距离
##读取各物种种内距离，提取数值，存成总表
inter_dist <- distance
spe_num <-7
for (i in 1:spe_num) ## 设置物种个数
{
  SN <- table[i,1]  ##物种名
  RS <- table[i,2]  ##行起始数字
  RE <- table[i,3]  ##行结束数字
  CS <- RS+1  ##列起始数字，多名字列故加一
  if (RE < num)  {CE <- RE+1} ##列结束数字，因距离表格最后一列空行未能读入，最后一列不加一
  else {CE <- RE}
  inter_dist[RS:RE,CS:CE] = NA
}
write.csv(inter_dist, file = "distance2.csv") ##手动删除行列

##scan读入为字符型向量，转化为数值型
inter_dist2 <- as.numeric(scan("distance2.csv", what="", sep = ",", na.strings = "NA")) 
summary(inter_dist2, na.rm = TRUE)
min(inter_dist2, na.rm = TRUE)
max(inter_dist2, na.rm = TRUE)
mean(inter_dist2, na.rm = TRUE)

inter_dist2 <- na.omit(inter_dist2)
write.csv(inter_dist2, file = "matK_inter_distance2.csv")
intra_dist2 <- na.omit(unlist(intra_dist))
write.csv(intra_dist2, file = "matK_intra_distance2.csv")

##
summary(unlist(distance3[,2:166]),na.rm = TRUE)
sd(unlist(distance[,2:141]),na.rm = TRUE)
write.csv(table, file = "species_list.csv")


##秩和检验
wilcox.test(intra_dist2,inter_dist2,paired=FALSE, conf.level = 0.95)
t.test(intra_dist2,inter_dist2,paired=FALSE, conf.level = 0.95)


##频数图
library(ggplot2)
ggplot(diamonds, aes(carat)) +  geom_histogram()
par(new=TRUE)
hist(inter_dist2, freq = FALSE, breaks= seq(0,0.06,0.001),col="grey")
par(new=TRUE)
hist(intra_dist2, freq = FALSE, breaks= seq(0,0.06,0.001),col="white")

par(new=TRUE)
intra_dist_cut<-cut(intra_dist2,breaks= seq(0,0.06,0.001)) ##binning data
intra_dist_cut_tab <- table(intra_dist_cut) ## binned data table
intra_dist_cut
barplot(intra_dist_cut_tab)

intra_dist3 <- read.csv(file ="matK_intra_distance2.csv", header = FALSE, sep = ",")
try <- intra_dist3[,1]
hist(inter_dist2, freq = FALSE, breaks= seq(0,0.06,0.001))
lines(density(intra_dist2))
str(inter_dist2)
barplot(inter_dist3)

inter_dist3 <- c()
inter_dist3 <- append(try, intra_dist2)

names(intra_dist[[1]]) <- table[1,1]
names(intra_dist[[1]])
names(intra_dist)
intra_dist
rm(try)

inter_dist3<-as.numeric(na.omit(inter_dist2))
mean(inter_dist3)
op <- par(mfrow = c(2, 2))
hist(islands,freq=FALSE)
