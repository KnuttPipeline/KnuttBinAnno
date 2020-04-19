#!/usr/bin/env Rscript
library(seqinr)
library(data.table)
options(warn=2)
source("scripts/commonGFF3Tools.R")

#save.image()

kofam.path <- snakemake@input[["kofam"]]
cds.path <- snakemake@input[["query"]]
unknownthrshldeval <- snakemake@params[["ifunsure_eval"]]
out.path <- snakemake@output[["gffbase"]]
out.ko.path <- snakemake@output[["kos"]]
out.ko.sure.path <- snakemake@output[["kos_onlysure"]]
sequenceontology_path <- snakemake@input[["so"]]


kofam <- fread(kofam.path,na.strings=c(getOption("datatable.na.strings","NA"),"-"))
kofam <- kofam[(score>=thrshld)|(is.na(thrshld)&eval<=unknownthrshldeval),]

fwrite(kofam[!is.na(thrshld),.(seqid,KO)],out.ko.sure.path,sep="\t",col.names=F)
fwrite(kofam[,.(seqid,KO)],out.ko.path,sep="\t",col.names=F)

cds.lengths <- rbindlist(lapply(read.fasta(cds.path),function(entry)data.table(seqid=attr(entry,"name"),length=length(entry))))
kofam <- cds.lengths[kofam, on="seqid"]
kofam[, KO:=paste0("ko:", KO)]
kofam[grepl(".+\\[EC:([ -\\.\\d]+)\\]",descr),EC:=sub(".+\\[EC:([ -\\.\\d]+)\\]","\\1",descr)]
kofam[,descr:=trimws(sub("(.+)\\[EC:([ -\\.\\d]+)\\]","\\1",descr))]
kofam[,EC:=gsub("(\\d{1,2}(\\.(\\-|\\d{1,2})){3})","EC:\\1",EC)]
kofam[,EC:=strsplit(EC," ",EC,fixed=T)]
combineCols <- function(dat){
   dat <- do.call(function(...)mapply(function(...)unlist(unname(c(...))),...),dat)
   dat <- lapply(dat,function(subels)subels[!is.na(subels)])
   dat[sapply(dat,length)==0] <- "NA"
   dat
}
kofam[,EC:=combineCols(list(EC))]
setnames(kofam,c("descr","KO","thrshld","EC","length","score","eval"),c("Note","Ontology_term","KO_thrshld","Dbxref", "end","KO_score","score"))
kofam$hit <- NULL
kofam[,source:=ifelse(!is.na(KO_thrshld),"kofamkoala","kofamkoala-eval-sel")]
kofam[,type:="protein_match"]
kofam[,strand:="+"]
kofam[,phase:="."]
kofam[,start:=1]

kofam <- formatDataFrame(kofam, sequenceontology_path)

fwrite(kofam,out.path,sep="\t")


