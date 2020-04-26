#!/usr/bin/env Rscript
library(data.table)
options(warn=2)
source("scripts/commonGFF3Tools.R")

hyddbfile <- snakemake@input[["dat"]]
out.path <- snakemake@output[[1]]
sequenceontology_path <- snakemake@input[["so"]]


hyddb <- fread(hyddbfile,fill=T)
if(nrow(hyddb)!=0){
   hyddb <- hyddb[,.(seqid=query, start=querystart, end=queryend, score=eval, source="HydDBClient", 
       type="protein_match", strand="+", phase=".", target_id=paste0(center,"|",subject),
       target_start=subjectstart, target_end=subjectend, target_strand=NA, Dbxref=paste0("HydDB:",Classification))]
   hyddb <- formatDataFrame(hyddb, sequenceontology_path)
}else{
#seqid	source	type	start	end	score	strand	phase	attributes
   hyddb <- data.table(seqid=character(),source=character(),type=character(),start=numeric(),end=numeric(),score=numeric(),strand=character(),phase=character(),attributes=character())
}
fwrite(hyddb,out.path,sep="\t",quote=F)


