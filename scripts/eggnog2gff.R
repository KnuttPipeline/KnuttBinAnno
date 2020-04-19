#!/usr/bin/env Rscript

options(warn=2)
library(seqinr)
library(data.table)
source("scripts/commonGFF3Tools.R")

eggnog.anno.path <- snakemake@input[["anno"]]
cds.path <- snakemake@input[["query"]]
out.path <- snakemake@output[[1]]
sequenceontology_path <- snakemake@input[["so"]]



annocols <- c("GeneID","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","best_tax_level","Preferred_name","GOs","EC",
"KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","tax_scope","eggnogOGS","bestorsmallestog","COG","descr")

eggnog.anno <- fread(cmd=paste0("grep -v '#' ",shQuote(eggnog.anno.path)), header=F, col.names=annocols)
cds.lengths <- rbindlist(lapply(read.fasta(cds.path), function(entry)data.table(GeneID=attr(entry,"name"), length=length(entry))))
eggnog.anno <- cds.lengths[eggnog.anno, on="GeneID"]
setnames(eggnog.anno, c("best_tax_level", "Preferred_name", "descr"), c("tax","Name","Note"))
eggnog.anno[Name=="", Name:=NA]
eggnog.anno[Note=="",Note:=NA]
eggnog.anno[tax=="",tax:=NA]
eggnog.anno[,GOs:=strsplit(GOs,",",fixed=T)]
eggnog.anno[,KEGG_ko:=strsplit(KEGG_ko,",",fixed=T)]
eggnog.anno[,EC:=strsplit(gsub("(\\d{1,2}(\\.(\\-|\\d{1,2})){3})","EC:\\1",EC),",",fixed=T)]
eggnog.anno[,KEGG_Pathway:=strsplit(gsub("((ko|map)\\d+)","KEGG_PATHWAY:\\1",KEGG_Pathway),",",fixed=T)]
eggnog.anno[,KEGG_Module:=strsplit(gsub("(M\\d+)","KEGG_MODULE:\\1",KEGG_Module),",",fixed=T)]
eggnog.anno[,KEGG_Reaction:=strsplit(gsub("(R\\d+)","KEGG_REACTION:\\1",KEGG_Reaction),",",fixed=T)]
eggnog.anno[,KEGG_rclass:=strsplit(gsub("(RC\\d+)","KEGG_RCLASS:\\1",KEGG_rclass),",",fixed=T)]
eggnog.anno[,BRITE:=strsplit(gsub("((ko|br)\\d+)","BRITE:\\1",BRITE),",",fixed=T)]
eggnog.anno[,KEGG_TC:=strsplit(gsub("([0-9]+\\.[A-Z](\\.[0-9]+)+)","TC:\\1",KEGG_TC),",",fixed=T)]
eggnog.anno[,CAZy:=strsplit(gsub("(([A-Z]{2,3}\\d{1,4})((?:_\\d{1,3})?))","CAZY:\\1",CAZy),",",fixed=T)]
eggnog.anno[,BiGG_Reaction:=strsplit(gsub("([^,]+)","BIGG_REACTION:\\1",BiGG_Reaction),",",fixed=T)]
eggnog.anno[,eggnogOGS:=strsplit(gsub("([^,]+)","EGGNOG:\\1",eggnogOGS),",",fixed=T)]
eggnog.anno[,COG:=strsplit(gsub("([^,]+)","EGGNOG_FUNCTION:\\1",COG),",",fixed=T)]
combineCols <- function(dat){
   dat <- do.call(function(...)mapply(function(...)unlist(unname(c(...))),...),dat)
   dat <- lapply(dat,function(subels)subels[!is.na(subels)])
   dat[sapply(dat,length)==0] <- "NA"
   dat
}
eggnog.anno[, Ontology_term:=combineCols(.SD),.SD=c("GOs", "KEGG_ko")]
eggnog.anno[, Dbxref:=combineCols(.SD),.SD=c("EC","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","eggnogOGS","COG")]
eggnog.anno <- eggnog.anno[,.(seqid=GeneID,source="eggnog-mapper-v2",type="protein_match",start=1,end=length,score=seed_ortholog_evalue,strand="+",phase=".",
   target_id=paste0("eggNOG|",tax_scope,"|",seed_eggNOG_ortholog), target_start=1, target_end=1, target_strand=NA, tax, Name, Note, Ontology_term, Dbxref)]
eggnog.anno[sapply(Ontology_term,length)==0,Ontology_term:="NA"]
eggnog.anno[sapply(Dbxref,length)==0,Dbxref:="NA"]
eggnog.anno <- formatDataFrame(eggnog.anno, sequenceontology_path)

fwrite(eggnog.anno,out.path,sep="\t")
