
rm(list=ls())
setwd("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype")

##���������Զ���ȡhaplotype
haplo.info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/hap.info.abrev.csv",as.is=T)
seq.name<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/haplotype.seq.name.csv",as.is=T)
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
n=401
for (i in 1: n) {
  seq.name$hap[i]<-haplo.info$hap[grep(seq.name[i,2],haplo.info$name)]
}
#���������ᱨ��������ִ������
write.csv(seq.name,"haplotype.info.csv")

#��haplotypeд���ܱ�
haplo.info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/haplotype.info.csv")
#haplo.info$seq.nameԤ����������_|_�滻Ϊ___��-�滻Ϊ_
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
info$hap<-haplo.info$hap[match(info$sequence,haplo.info$seq.name)]
write.csv(info,"E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv")

#�������У������������������齨����
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings") 
library("Biostrings")
fastaFile <- readDNAStringSet("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/4ampy.add.outgroup.allsec.total.fasta")
seq<- data.frame(name = names(fastaFile), sequence = paste(fastaFile))
seq$name<-gsub("_[|]_","___",seq$name)
seq$name<-gsub("-","_",seq$name)

##����info��use����TRUE��ע����ȡ����
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
simple.seq<-data.frame(name=info$tip.label[which(info$rasp)],
                       sequence=seq$sequence[match(info$sequence[which(info$rasp)],seq$name)])
simple.seq<-data.frame(name=info$sequence[grep("Quercus|Protobalanus",info$section)],
                       sequence=seq$sequence[match(info$sequence[grep("Quercus|Protobalanus",info$section)],seq$name)])
pq.gps<-info[grep("Quercus|Protobalanus",info$section),c("sequence","longitude","latitude")]
write.csv(pq.gps,"E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/pq.gps.csv")
#simple.seq$name<-info$tip.label[match(simple.seq$name,info$sequence)]
##תΪfasta
hap.seq<-paste(">",simple.seq$name,sep="")
hap.seq<-paste(hap.seq,simple.seq$sequence,sep="###")
write.csv(hap.seq,"E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/pq.csv") ##ע�������еĴ����ֶ���###�滻Ϊ���з�\n
#�洢�ļ�λ�ô���library(seqinr) write.fasta(sequence,name,'simple.fasta')

#��ȡrasp������Ϣ
rasp.info<-info[which(info$rasp),]
write.csv(rasp.info,"haplotype.info.csv")
setwd("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/rasp/")


##������������Ϣ��ȡ
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