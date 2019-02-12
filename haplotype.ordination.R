
rm(list=ls())
setwd("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype")

##按照名字自动读取haplotype
haplo.info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/hap.info.abrev.csv",as.is=T)
seq.name<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/haplotype.seq.name.csv",as.is=T)
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
n=401
for (i in 1: n) {
  seq.name$hap[i]<-haplo.info$hap[grep(seq.name[i,2],haplo.info$name)]
}
#本步结束会报错，但是执行正常
write.csv(seq.name,"haplotype.info.csv")

#将haplotype写入总表
haplo.info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/haplotype.info.csv")
#haplo.info$seq.name预先重命名，_|_替换为___，-替换为_
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
info$hap<-haplo.info$hap[match(info$sequence,haplo.info$seq.name)]
write.csv(info,"E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv")

#读入序列，并将序列名和序列组建表格
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings") 
library("Biostrings")
fastaFile <- readDNAStringSet("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/4ampy.add.outgroup.allsec.total.fasta")
seq<- data.frame(name = names(fastaFile), sequence = paste(fastaFile))
seq$name<-gsub("_[|]_","___",seq$name)
seq$name<-gsub("-","_",seq$name)

##根据info中use列中TRUE标注，抽取序列
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
simple.seq<-data.frame(name=info$tip.label[which(info$rasp)],
                       sequence=seq$sequence[match(info$sequence[which(info$rasp)],seq$name)])
simple.seq<-data.frame(name=info$sequence[grep("Quercus|Protobalanus",info$section)],
                       sequence=seq$sequence[match(info$sequence[grep("Quercus|Protobalanus",info$section)],seq$name)])
pq.gps<-info[grep("Quercus|Protobalanus",info$section),c("sequence","longitude","latitude")]
write.csv(pq.gps,"E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/pq.gps.csv")
#simple.seq$name<-info$tip.label[match(simple.seq$name,info$sequence)]
##转为fasta
hap.seq<-paste(">",simple.seq$name,sep="")
hap.seq<-paste(hap.seq,simple.seq$sequence,sep="###")
write.csv(hap.seq,"E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/pq.csv") ##注意检查序列的错误，手动将###替换为换行符\n
#存储文件位置错乱library(seqinr) write.fasta(sequence,name,'simple.fasta')

#调取rasp所用信息
rasp.info<-info[which(info$rasp),]
write.csv(rasp.info,"haplotype.info.csv")
setwd("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/rasp/")


##共享单倍型信息提取
haploinfo<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/hap.info.shared.csv")
hap_share<-haploinfo$hap[which(haploinfo$shared)]
hap_share
cyc_share<-hap_share[1:7]
que_share<-hap_share[9:19]
ilex_share<-hap_share[22:29]

info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
cyc_hap<-info[which(info$hap %in% cyc_share),]
write.csv(cyc_hap,"cyc_hap.csv")
que_hap<-info[which(info$hap %in% que_share),]
write.csv(que_hap,"que_hap.csv")
ilex_hap<-info[which(info$hap %in% ilex_share),]
write.csv(ilex_hap,"ilex_hap.csv")
