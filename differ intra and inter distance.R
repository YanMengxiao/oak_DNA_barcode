library(stringr)

rm(list=ls())
##�ļ�����
setwd("E:/Xiong/data/20170706/final/distance") ##�����ļ���
distance <- read.csv(file ='4ampy_lob.csv', header = FALSE, sep = ",") ##�����ļ���
##�������
species <- c()
row_start <- 1 ##�趨��һ������ʼΪ��һ��
row_end <- c()
name <- "[A-Z].*[a-z]"  ##��ȡ��������

##����ÿ���ֵ���ʼ�к�
num <- 7
for (i in 1:num)
  if (str_match(distance[i,1], name) == str_match(distance[i+1,1], name)) {
    m <- i
  } else {
    species <- append(species, str_match(distance[i,1], name))
    row_end <- append(row_end,i)
    row_start <-append(row_start,i+1)
  }
row_end <- append(row_end,num) ##�趨���һ����Ϊend���һ��
species <- append(species, str_match(distance[num,1], name)) 


species
row_start
row_end

##д�����ݿ�
table <- data.frame(species,row_start,row_end)
table
write.csv(table, file = "total_specieslist.csv")

##
num <- 7
summary(unlist(distance[,2:num]),na.rm = TRUE)
sd(unlist(distance[,2:num]),na.rm = TRUE)
write.csv(table, file = "species_list.csv")

##���ھ���
##��ȡ���������ھ��룬��ȡ��ֵ������ܱ�
spe_num <- 7
intra_dist <- list()
intra_theta_full <- data.frame()
coales_depth_full <- data.frame()
for (i in 1:spe_num) ## �������ָ���
  {
  SN <- table[i,1]  ##������
  RS <- table[i,2]  ##����ʼ����
  RE <- table[i,3]  ##�н�������
  CS <- RS+1  ##����ʼ���֣��������йʼ�һ
  if (RE<num)   {CE <- RE+1} ##�н������֣������������һ�п���δ�ܶ��룬���һ�в���һ
  else {CE <- RE}
  intra_dist[[i]] <- as.numeric(unlist(str_match_all(unlist(distance[RS:RE,CS:CE]),"0.*"))) ##��ȡ�������ȥ�б�����ȡ��0���֣�ȥ�б���ת��Ϊ��ֵ��
  intra_theta_full[i,1] <- mean(intra_dist[[i]],na.rm = TRUE)
  if (max(intra_dist[[i]],na.rm = TRUE) == -Inf) {intra_dist[[i]] <- 0}
  else {coales_depth_full[i,1] <- max(intra_dist[[i]],na.rm = TRUE)}
   ##�������ĸ����������������ݿ�
   # fliename <- paste(table[i,1],".xlsx",sep="") ##ճ���ļ�������չ�����ϳ������ļ���
  #  namepath[[i]] <- paste("E:/Xiong/data/20170706/final/distance/",fliename,sep="")  ##ճ��洢·���������ļ���
  #  write.xlsx(intra_dist[[i]],namepath[[i]]) ## write.xlsx(x, file�����б�ʵ����������
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


## ��ÿ���ֵ����ھ�����xlsx�ӱ���
intra_dist <- list()
namepath <- list()
for (i in 1:spe_num) {
  SN <- table[i,1]  ##������
  RS <- table[i,2]  ##����ʼ����
  RE <- table[i,3]  ##�н�������
  CS <- RS+1  ##����ʼ���֣��������йʼ�һ
  if (RE<num)   {CE <- RE+1} ##�н������֣������������һ�п���δ�ܶ��룬���һ�в���һ
  else {CE <- RE}
  intra_dist[[i]]  <- distance[RS:RE,CS:CE] ##���ӱ���������ݿ�
  # fliename <- paste(table[i,1],".xlsx",sep="") ##ճ���ļ�������չ�����ϳ������ļ���
  # namepath[[i]] <- paste("E:/Xiong/data/20170706/final/distance/",fliename,sep="")  ##ճ��洢·���������ļ���
  # write.xlsx(intra_dist[[i]],namepath[[i]]) ## write.xlsx(x, file�����б�ʵ����������
}


##�ּ����
##��ȡ���������ھ��룬��ȡ��ֵ������ܱ�
inter_dist <- distance
spe_num <-7
for (i in 1:spe_num) ## �������ָ���
{
  SN <- table[i,1]  ##������
  RS <- table[i,2]  ##����ʼ����
  RE <- table[i,3]  ##�н�������
  CS <- RS+1  ##����ʼ���֣��������йʼ�һ
  if (RE < num)  {CE <- RE+1} ##�н������֣������������һ�п���δ�ܶ��룬���һ�в���һ
  else {CE <- RE}
  inter_dist[RS:RE,CS:CE] = NA
}
write.csv(inter_dist, file = "distance2.csv") ##�ֶ�ɾ������

##scan����Ϊ�ַ���������ת��Ϊ��ֵ��
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


##�Ⱥͼ���
wilcox.test(intra_dist2,inter_dist2,paired=FALSE, conf.level = 0.95)
t.test(intra_dist2,inter_dist2,paired=FALSE, conf.level = 0.95)


##Ƶ��ͼ
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