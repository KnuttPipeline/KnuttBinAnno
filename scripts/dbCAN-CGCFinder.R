#!/usr/bin/env Rscript
options(warn=2)

library(data.table)
library(seqinr)


dbcan.path = "output/BinAnnotation/dbCAN/HotVent/HotVent_dbcan_all.tsv"
cds.path = "output/BinAnnotation/metaerg/HotVent/data/cds_BinIds.faa"
maxnotsiggenes = 1
out.path = "test.tsv"

if(exists("snakemake")){
  dbcan.path = snakemake@input[["dbcan"]]
  cds.path = snakemake@input[["query"]]
  maxnotsiggenes = snakemake@params[["spacersig"]]
  out.path = snakemake@output[[1]]
}

dbcan = fread(dbcan.path)
cds=read.fasta(cds.path)
contig2cds=sub(".+cid=([^ ]+).*","\\1",sapply(cds[1:5],attr,which="Annot"))
contig2cds=as.data.table(cbind(contig=contig2cds,seqid=names(contig2cds)))
dbcan.markers = dbcan[contig2cds,on="seqid"]
dbcan.markers[is.na(source),source:=""]

dbcan.markers = dbcan.markers[,.(source=paste0(source,collapse=",")),keyby=c("contig","seqid")]
dbcan.markers[,source:=strsplit(source,",",fixed=T)]
contigs.withhits = dbcan.markers[,sum(sapply(source,length)),by=contig][V1>0,contig]
dbcan.markers[,index:=1:.N,by=contig]
dbcan.markers[sapply(source,length)>0,spaces:=c(0,diff(index))-1,by=contig]
if(length(contigs.withhits)==0){
  file.create(out.path)
  quit()
}
dbcan.markers = dbcan.markers[contig %in% contigs.withhits,.SD[first(which(!is.na(spaces))):.N,],by=contig]

assignCluster <- function(vector,maxdist){
  curcluster = 0
  clusterassignment = rep(0,length(vector))
  for(i in seq_along(vector)){
     if(!is.na(vector[[i]]) && vector[[i]]>maxdist){
       curcluster = curcluster + 1
     }
     clusterassignment[[i]]=curcluster
  }
  clusterassignment
}

dbcan.markers[,cluster:=assignCluster(spaces,maxnotsiggenes),by=contig]
dbcan.markers[,clusterid:=.GRP,by=c("contig","cluster")]
dbcan.markers[,clustersize:=.N,by=clusterid]
dbcan.markers[,importantgenecount:=sum(sapply(source,length)>0),by=clusterid]
dbcan.markers = dbcan.markers[importantgenecount>1,.(clustersize=first(clustersize),dbCANgenecount=first(importantgenecount),
   members=paste0(seqid,collapse="|;|"),sources=paste0(dbcan[seqid,source,on="seqid"],collapse="|;|"),
   attributes=paste0(dbcan[seqid,attributes,on="seqid"],collapse="|;|")
  ),
  keyby=c("contig","clusterid")]


fwrite(dbcan.markers,out.path,sep="\t")

