#!/usr/bin/env Rscript
library(seqinr)
library(data.table)
options(warn=2)

kofam.path = "output/BinAnnotation/kofam/Larrelt/Larrelt.1.tsv"
cds.path = "output/BinAnnotation/metaerg/Larrelt/Larrelt.1/data/cds.faa"
out.path = "test.tsv"
out.ko.sure.path = "test2.tsv"
out.ko.path = "test3.tsv"
unknownthrshldeval = 1.0e-10


if(exists("snakemake")){
  kofam.path = snakemake@input[["kofam"]]
  cds.path = snakemake@input[["query"]]
  unknownthrshldeval = snakemake@params[["ifunsure_eval"]]
  out.path = snakemake@output[["gffbase"]]
  out.ko.path = snakemake@output[["kos"]]
  out.ko.sure.path = snakemake@output[["kos_onlysure"]]
}



cols = c("GeneID","KO","thrshld","score","eval","descr")
kofam = system(paste0("grep -v '#' ",shQuote(kofam.path)),intern=T)
kofam = sapply(kofam,trimws)
kofam = sub("^(?:\\*\\s+)?(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(.+)$","\\1\t\\2\t\\3\t\\4\t\\5\t\\6",kofam,perl=T)
kofam = paste0(kofam,collapse="\n")
kofam = fread(kofam,header=F,col.names=cols,na.strings=c(getOption("datatable.na.strings","NA"),"-"))
kofam = kofam[(score>=thrshld)|(is.na(thrshld)&eval<=unknownthrshldeval),]

fwrite(kofam[!is.na(thrshld),.(GeneID,KO)],out.ko.sure.path,sep="\t",col.names=F)
fwrite(kofam[,.(GeneID,KO)],out.ko.path,sep="\t",col.names=F)

cds.lengths = rbindlist(lapply(read.fasta(cds.path),function(entry)data.table(GeneID=attr(entry,"name"),length=length(entry))))
kofam = cds.lengths[kofam,on="GeneID"]
kofam[grepl(".+\\[EC:([ -\\.\\d]+)\\]",descr),EC:=sub(".+\\[EC:([ -\\.\\d]+)\\]","\\1",descr)]
kofam[,txt:=sub("(.+)\\[EC:([ -\\.\\d]+)\\]","\\1",descr)]
kofam[!is.na(EC),EC:=gsub(" ",",",EC,fixed=T)]

kofam[,attributes:=paste0("Name=",lapply(txt,trimws))]
kofam[,attributes:=paste0(attributes,";Ontology_term=",KO)]
kofam[,attributes:=paste0(attributes,";KO_thrshld=",thrshld)]
kofam[,attributes:=paste0(attributes,";KO_score=",score)]
kofam[!is.na(EC),attributes:=paste0(attributes,";Dbxref=",gsub("(\\d{1,2}(\\.(\\-|\\d{1,2})){3})","EC:\\1",EC))]

kofam = kofam[,.(seqid=GeneID,source=ifelse(!is.na(thrshld),"kofamkoala-sure","kofamkoala-unsure-eval-sel"),type="protein match",start=1,end=length,score=eval,strand="+",phase=".",attributes)]

fwrite(kofam,out.path,sep="\t")


