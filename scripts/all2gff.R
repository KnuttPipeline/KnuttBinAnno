options(warn=2)
source("scripts/commonGFF3Tools.R")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(seqinr))


infiles <- snakemake@input$gffs
sofile <- snakemake@input$so
cdsfile <- snakemake@input$cds
gfffile <- snakemake@output$gff
tsvfile <- snakemake@output$tsv


gffbase <- rbindlist(lapply(infiles, function(f)as.data.table(readBaseGFF(fread(f)))),fill=T)
setorder(gffbase, seqid,start,end,score)
tsvcols <- c("source","target","Note", "TC", "CAZY", "EC", "InterPro", "Reactome", "MetaCyc",
     "GO", "tax", "EGGNOG", "EGGNOG_FUNCTION", "BIGG_REACTION", "ko", "HydDB", "Name")
tsvcols <- tsvcols[tsvcols %in% colnames(gffbase)]
collapseels <- function(els){
    els <- unique(unlist(els))
    if(is.character(els))
        els <- gsub(",", ";", els)
    els <- els[!is.na(els) & els!=""]
    paste0(els,collapse=",")
}
tsvdat <- gffbase[, lapply(.SD,collapseels), .SD=tsvcols, by="seqid"]
setkey(tsvdat, seqid)
fwrite(tsvdat,tsvfile,sep="\t")
writeGFF(gffbase, cdsfile, gfffile, sofile)

