## extract sequence
# 
# source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings") 
library("Biostrings")
fastaFile <- readDNAStringSet("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/4ampy.add.outgroup.allsec.total.fasta")
seq<- data.frame(name = names(fastaFile), sequence = paste(fastaFile))
# seq$name<-gsub("_[|]_","___",seq$name)
# seq$name<-gsub("-","_",seq$name)

##根据info中use列中TRUE标注，抽取序列
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/climate/mydata/ilex with climate data2.5.csv",as.is=T)
# simple.seq<-data.frame(name=info$tip.label[which(info$rasp)],
#                        sequence=seq$sequence[match(info$sequence[which(info$rasp)],seq$name)])
simple.seq<-data.frame(name=seq$name[match(info$sequence,seq$name)],
                       sequence=seq$sequence[match(info$sequence,seq$name)])
#simple.seq$name<-info$tip.label[match(simple.seq$name,info$sequence)]
##转为fasta
hap.seq<-paste(">",simple.seq$name,sep="")
hap.seq<-paste(hap.seq,simple.seq$sequence,sep="###")
write.csv(hap.seq,"E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/ilex.seq.csv") ##注意检查序列的错误，手动将###替换为换行符\n
#存储文件位置错乱library(seqinr) write.fasta(sequence,name,'simple.fasta')