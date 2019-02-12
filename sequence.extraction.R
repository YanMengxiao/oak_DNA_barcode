## extract sequence
# 
# source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings") 
library("Biostrings")
fastaFile <- readDNAStringSet("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/4ampy.add.outgroup.allsec.total.fasta")
seq<- data.frame(name = names(fastaFile), sequence = paste(fastaFile))
# seq$name<-gsub("_[|]_","___",seq$name)
# seq$name<-gsub("-","_",seq$name)

##����info��use����TRUE��ע����ȡ����
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/climate/mydata/ilex with climate data2.5.csv",as.is=T)
# simple.seq<-data.frame(name=info$tip.label[which(info$rasp)],
#                        sequence=seq$sequence[match(info$sequence[which(info$rasp)],seq$name)])
simple.seq<-data.frame(name=seq$name[match(info$sequence,seq$name)],
                       sequence=seq$sequence[match(info$sequence,seq$name)])
#simple.seq$name<-info$tip.label[match(simple.seq$name,info$sequence)]
##תΪfasta
hap.seq<-paste(">",simple.seq$name,sep="")
hap.seq<-paste(hap.seq,simple.seq$sequence,sep="###")
write.csv(hap.seq,"E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/ilex.seq.csv") ##ע�������еĴ����ֶ���###�滻Ϊ���з�\n
#�洢�ļ�λ�ô���library(seqinr) write.fasta(sequence,name,'simple.fasta')