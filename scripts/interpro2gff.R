#!/usr/bin/env Rscript
library(data.table)
options(warn=2)
source("scripts/commonGFF3Tools.R")

#save.image()

interprofile <- snakemake@input[["dat"]]
out.path <- snakemake@output[[1]]
sequenceontology_path <- snakemake@input[["so"]]


interpro <- fread(interprofile,fill=T)
for(col in names(interpro)) set(interpro, i=which(interpro[[col]]==""), j=col, value=NA)
interpro[,len:=NULL]
interpro[,md5:=NULL]
interpro[,date:=NULL]
interpro[,status:=NULL]
interpro[,target_strand:=NA]
interpro[,target_start:=1]
interpro[,target_end:=1]
interpro[,GO:=strsplit(GO, "|", fixed=T)]
interpro[sapply(GO,length)==0,GO:=NA]
interpro[,pathways:=gsub(" ", "", pathways, fixed=T)]
interpro[,pathways:=strsplit(pathways, "|", fixed=T)]
interpro[interpro!="",interpro:=paste0("InterPro:", interpro)]
interpro[,Dbxref:=combineCols(.SD),.SD=c("interpro", "pathways")]
interpro[,interpro:=NULL]
interpro[,pathways:=NULL]
interpro[,type:="protein_match"]
interpro[,phase:="."]
interpro[,strand:="+"]
setnames(interpro,c("analysis","sigacc","interprodescr","stop","GO"),c("source","target_id","Note","end","Ontology_term"))
interpro <- formatDataFrame(interpro, sequenceontology_path)
fwrite(interpro,out.path,sep="\t")


