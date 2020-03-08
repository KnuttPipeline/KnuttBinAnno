#!/usr/bin/env Rscript
library(seqinr)
library(data.table)
options(warn=2)

eggnog.anno.path = "output/BinAnnotation/eggnog/Larrelt/Larrelt.1.emapper.annotations"
cds.path = "output/BinAnnotation/metaerg/Larrelt/Larrelt.1/data/cds.faa"
out.path = "test.tsv"

if(exists("snakemake")){
  eggnog.anno.path = snakemake@input[["anno"]]
  cds.path = snakemake@input[["query"]]
  out.path = snakemake@output[[1]]
}



gfftextescaper <- function(textvector){
  textvector = gsub("\t","%09",textvector,fixed=T)
  textvector = gsub("\n","%0A",textvector,fixed=T)
  textvector = gsub("\r","%0D",textvector,fixed=T)
  textvector = gsub("%","%25",textvector,fixed=T)
  for(char in c(1:31,127)){
    textvector = gsub(intToUtf8(char),paste0("%",format(as.hexmode(char),width=2,upper.case=T)),textvector,fixed=T,useBytes=T)
  }
}

gffcol9escaper <- function(textvector){
 # Todo fixme
  #textvector = gfftextescaper(textvector)
  textvector = gsub(";","%3B",textvector,fixed=T)
  textvector = gsub("=","%3D",textvector,fixed=T)
  textvector = gsub("&","%26",textvector,fixed=T)
  gsub(",","%2C",textvector,fixed=T)
}

gffseqidsingleEncoder <- function(text){
    OK <- paste0("[^","ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz0123456789",".:^*$@!+_?-|","]")
    x <- strsplit(text, "")[[1L]]
    z <- grep(OK, x)
    if (length(z)) {
        y <- sapply(x[z], function(x) paste0("%", toupper(as.character(charToRaw(x))),
            collapse = ""))
        x[z] <- y
    }
    paste(x, collapse = "")
}

gffseqidEncoder <- function(textvector){
   unlist(lapply(textvector,gffseqidsingleEncoder))
}


annocols=c("GeneID","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","best_tax_level","Preferred_name","GOs","EC",
"KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","tax_scope","eggnogOGS","bestorsmallestog","COG","descr")

eggnog.anno = fread(cmd=paste0("grep -v '#' ",shQuote(eggnog.anno.path)),header=F,col.names=annocols)
cds.lengths = rbindlist(lapply(read.fasta(cds.path),function(entry)data.table(GeneID=attr(entry,"name"),length=length(entry))))
eggnog.anno = cds.lengths[eggnog.anno,on="GeneID"]

eggnog.anno[,attributes:=paste0("target=eggNOG|",tax_scope,"|",seed_eggNOG_ortholog)]
eggnog.anno[,attributes:=paste0(attributes,";tax=",best_tax_level)]
eggnog.anno[Preferred_name!="",attributes:=paste0(attributes,";Name=",Preferred_name)]
eggnog.anno[descr!="",attributes:=paste0(attributes,";descr=",descr)]



eggnog.anno[,GOs:=strsplit(GOs,",",fixed=T)]
eggnog.anno[,KEGG_ko:=strsplit(KEGG_ko,",",fixed=T)]
eggnog.anno[,onts:=rowSums(as.data.frame(lapply(.SD,function(col)sapply(col,length)))),.SD=c("GOs","KEGG_ko")]
eggnog.anno[onts>0,attributes:=paste0(attributes,";Ontology_term=",apply(.SD,1,function(stringlists)paste0(unlist(stringlists),collapse=","))),.SD=c("GOs","KEGG_ko")]

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
eggnog.anno[,xrefs:=rowSums(as.data.frame(lapply(.SD,function(col)sapply(col,length)))),.SD=c("EC","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","eggnogOGS","COG")]
eggnog.anno[xrefs>0,attributes:=paste0(attributes,";Dbxref=",apply(.SD,1,function(stringlists)paste0(unlist(stringlists),collapse=","))),.SD=c("EC","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","eggnogOGS","COG")]




eggnog.anno = eggnog.anno[,.(seqid=gffseqidEncoder(GeneID),source="eggnog-mapper-v2",type="protein match",start=1,end=length,score=seed_ortholog_evalue,strand="+",phase=".",attributes)]

fwrite(eggnog.anno,out.path,sep="\t")
