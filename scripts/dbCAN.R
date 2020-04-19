#!/usr/bin/env Rscript
options(warn=2)

#options(echo=T)
#save.image()

library(data.table)
library(tidyr)
library(seqinr)

if(exists("snakemake")){
  hotpep.path <- snakemake@input[["hotpep"]]
  hotpep.ref.dirpath <- snakemake@params[["hotpepref"]]
  diamond.path <- snakemake@input[["diamond"]]
  diamond.ref.data.path <- snakemake@input[["diamond_ref"]]
  hmmscan.path <- snakemake@input[["hmmscan"]]
  cds.path <- snakemake@input[["query"]]
  gff.tsv.path <- snakemake@output[[1]]
  hmmscan.eval <- snakemake@params[["hmmscaneval"]]
  hmmscan.coverage <- snakemake@params[["hmmscancov"]]
  tf1.hmmscan.path <- snakemake@input[["tf1hmmscan"]]
  tf2.hmmscan.path <- snakemake@input[["tf2hmmscan"]]
  stp.hmmscan.path <- snakemake@input[["stphmmscan"]]
  tf.hmmscan.eval <- snakemake@params[["tf_hmmscaneval"]]
  tf.hmmscan.coverage <- snakemake@params[["tf_hmmscancov"]]
  stp.hmmscan.eval <- snakemake@params[["stp_hmmscaneval"]]
  stp.hmmscan.coverage <- snakemake@params[["stp_hmmscancov"]]
  tcdb.path <- snakemake@input[["tcdb"]]
  tcdb.ref.path <- snakemake@input[["tcdb_ref"]]
  sequenceontology_path <- snakemake@input[["so"]]
}


setDTthreads(1)

hmmscan.regex <- "^([A-Z]{2,3}\\d{1,4}).*?((?:_\\d{1,3})?)$"
hotpepcols <- c("CAZyFamily","PPRSubfamily","GeneID","Frequency","Hits","SignaturePeptides")
if(file.info(hotpep.path)$size > 0){
  hotpep <-  fread(hotpep.path, col.names=hotpepcols, header=F, sep="\t")
}else{
  hotpep <- data.table(CAZyFamily = character(), PPRSubfamily = character(), GeneID = character(), Frequency = character(), Hits = character(), SignaturePeptides = character())
}
hotpep.ref.files <- list.files(paste0(hotpep.ref.dirpath,"/CAZY_PPR_patterns"), pattern="*_group_ec.txt", recursive=T, full.names=T)
hotpep.ref <- rbindlist(lapply(hotpep.ref.files, function(file)cbind(file=basename(file),fread(file,header=F,sep="\t"))))
setnames(hotpep.ref,c("V1","V2"),c("PPRSubfamily","EC"))
#hotpep.ref = hotpep.ref[!is.na(EC)&EC!="",]
hotpep.ref[,EC:=gsub("x", "-", EC, fixed=T)]
hotpep.ref[,EC:=strsplit(EC,",",fixed=T)]
hotpep.ref[,EC:=lapply(EC,function(ECs)unlist(lapply(ECs,function(ec)gsub("([^:]+):.+","\\1",ec))))]
hotpep.ref[,CAZyFamily:=gsub("(.+)_group_ec\\.txt","\\1",file)]
hotpep <- hotpep.ref[hotpep,on=c("CAZyFamily","PPRSubfamily")]
hotpep[,CAZyFamily:=strsplit(CAZyFamily,",",fixed=T)]

if(file.info(diamond.path)$size>0){
  diamond <- fread(diamond.path,col.names=c("GeneID","CAZyID","Identity","Length","Mismatches","Gap open","Query.Start","Query.End","Subject.Start","Subject.End","eval","score"), header=F)
}else{
  diamond <- data.table(GeneID=character(), CAZyID=character(), Identity=numeric(), Length=numeric(),`Gap open`=numeric(),Query.Start=numeric(),Query.End=numeric(),Subject.Start=numeric(),Subject.End=numeric(),eval=numeric(),score=numeric())
}
diamond.ref.data <- fread(diamond.ref.data.path,header=T)
setnames(diamond.ref.data,c("sseqid","CAZyECs"),c("CAZyID","data"))

#sseqid	taxid	CAZyECs	superkingdom	phylum	class	order	family	genus	species

diamond.ref.data[,data:=strsplit(data," ",fixed=T)]
ECpattern <- "[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+"
diamond.ref.data[grepl(".",data,fixed=T),EC:=lapply(data,grep,pattern=ECpattern,perl=T,value=T)]
diamond.ref.data[,CAZyFamily:=lapply(data,grep,pattern=ECpattern,invert=T,perl=T,value=T)]
diamond.ref.data[,data:=NULL]

setkey(diamond.ref.data,"CAZyID")
setkey(diamond,"CAZyID")
diamond <- diamond.ref.data[diamond]

read.hmm <- function(file,feval,fcov){
  lines <- as.integer(system2("grep",args = c("-o","^[^#]*", file," | wc -l | awk '{print $1}'"),stdout = TRUE))
  if(lines>0){
    hmmscan.cmd <- paste0("grep -v '^#' ",shQuote(file)," |  awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_.\"\n\";}}'")
    hmmscan <- fread(cmd=hmmscan.cmd,header=F,col.names=c("HMMProfile","ProfileLength","GeneID","GeneLength","eval","Profile.Start","Profile.End","Gene.Start","Gene.End"))
    hmmscan[,Coverage:=(Profile.End-Profile.Start)/ProfileLength]
    hmmscan <- hmmscan[eval<feval & Coverage>=fcov]
  }else{
    hmmscan <- data.table(HMMProfile=character(),ProfileLength=numeric(),GeneID=character(),GeneLength=numeric(),eval=numeric(),Profile.Start=numeric(),Profile.End=numeric(),Gene.Start=numeric(),Gene.End=numeric())
  }
  hmmscan
}
hmmscan <- read.hmm(hmmscan.path,hmmscan.eval,hmmscan.coverage)
hmmscan[,basename:=gsub("(.+)\\.p?hmm","\\1",HMMProfile)]
hmmscan[grepl(hmmscan.regex,basename),CAZyFamily:=strsplit(sub(hmmscan.regex,"\\1\\2",basename),",",fixed=T)]
hmmscan[!grepl(hmmscan.regex,basename),CAZyFamily:=strsplit("",",",fixed=T)]


tf1.hmmscan <- read.hmm(tf1.hmmscan.path,tf.hmmscan.eval,tf.hmmscan.coverage)
tf2.hmmscan <- read.hmm(tf2.hmmscan.path,tf.hmmscan.eval,tf.hmmscan.coverage)
stp.hmmscan <- read.hmm(stp.hmmscan.path,stp.hmmscan.eval,stp.hmmscan.coverage)

if(file.info(tcdb.path)$size>0){
  tcdb <- fread(tcdb.path,col.names=c("GeneID","TCDBID","Identity","Length","Mismatches","Gap open","Query.Start","Query.End","Subject.Start","Subject.End","eval","score"), header=F)
}else{
  tcdb <-  data.table(GeneID=character(),TCDBID=character(),Identity=numeric(),Length=numeric(),Mismatches=numeric(),`Gap open`=numeric(),Query.Start=numeric(),Query.End=numeric(),Subject.Start=numeric(),Subject.End=numeric(),eval=numeric(),score=numeric())
}
tcdb.ref <- data.table(txt=system(paste0("grep '>' ",shQuote(tcdb.ref.path)),intern=T))
tcdb.ref[,TCDBID:=sub("^>([^ ]+)( ?.*)$","\\1",txt)]
tcdb.ref[,descr:=sub("^>([^ ]+)( ?.*)$","\\2",txt)]
tcdb.ref[,c("asc","db","target","TCFamily"):=lapply(tstrsplit(TCDBID,"|",fixed=T), trimws)]
tcdb <- tcdb.ref[tcdb,on="TCDBID"]



diamond.hits <- diamond[,.(GeneID,CAZyFamily),] #Only one result per query allowed
hotpep.hits <- hotpep[,.(CAZyFamily=paste(CAZyFamily,collapse=",")),by=GeneID]
hmmscan.hits <- hmmscan[,.(CAZyFamily=paste(CAZyFamily,collapse=",")),by=GeneID]

hits <- rbindlist(list(diamond.hits,hotpep.hits,hmmscan.hits))
hits <- hits[,.(CAZyFamily=paste(unique(unlist(CAZyFamily)),collapse=",")),by=GeneID]
setkey(hits,"GeneID")

source("scripts/commonGFF3Tools.R")

addECandCAZyDbxRef <- function(CAZyFamily, EC){
  cazydb <- lapply(CAZyFamily,function(cazylist)unname(sapply(cazylist,function(cazy)paste0("CAZY:",trimws(cazy)))))
  ecdb <- lapply(EC,function(cazylist)unname(sapply(cazylist,function(cazy)paste0("EC:",trimws(cazy)))))
  res <- mapply(function(ca,ec)unlist(c(ca,ec)), cazydb,ecdb, SIMPLIFY=FALSE)
  res
}

if(nrow(diamond)>0){
  # Default dbCAN always reports the strand being +, which makes sense for BLASTP
  diamond.gff <- diamond[,.(seqid=GeneID, source="Knutt-dbCAN-diamond",type="protein_match", start=Query.Start,end=Query.End,
        score=eval,strand="+",phase=NA, target_id=paste0("CAZyDB.fa|",CAZyID), target_start=Subject.Start, target_end=Subject.End,
        target_strand=NA, Dbxref=as.list(addECandCAZyDbxRef(CAZyFamily, EC)))]
  diamond.gff <- formatDataFrame(diamond.gff, sequenceontology_path)
}else{
  diamond.gff <- NULL
}
if(nrow(hotpep)>0){
  cds.lengths <- rbindlist(lapply(read.fasta(cds.path),function(entry)data.table(GeneID=attr(entry,"name"),length=length(entry))))
  hotpep.gff <- cds.lengths[hotpep, on="GeneID"]
  hotpep.gff <- hotpep.gff[,.(seqid=GeneID, Note=ifelse(sapply(EC,length)==0,NA,"ECs derived from the CAZy (Sub)Family+PPR as a whole!"),
  source="Knutt-dbCAN-HotPep",type="protein_match",start=1,end=length,score=".",strand="+",phase=".",Dbxref=as.list(addECandCAZyDbxRef(CAZyFamily, EC)))]
  hotpep.gff <- formatDataFrame(hotpep.gff, sequenceontology_path)
}else{
  hotpep.gff <- NULL
}
if(nrow(hmmscan)>0){
  hmmscan.gff <- hmmscan[,.(seqid=GeneID,source="Knutt-dbCAN-hmmscan",type="protein_match",start=Gene.Start,end=Gene.End,
  score=eval, strand="+", phase=".",target_id=paste0("CAZyDB.hmm|",HMMProfile), target_start=Profile.Start, target_end=Profile.End, 
  target_strand=NA, Note=ifelse(sapply(CAZyFamily,length)!=0,NA,"Model name was not parseable"), Dbxref=lapply(CAZyFamily,function(cazylist)sapply(cazylist,function(cazy)paste0("CAZY:",trimws(cazy)))))]
  hmmscan.gff <- formatDataFrame(hmmscan.gff, sequenceontology_path)
}else{
  hmmscan.gff <- NULL
}
if(nrow(tf1.hmmscan)>0){
  tf1.hmmscan.gff <- tf1.hmmscan[,.(seqid=GeneID, source="Knutt-dbCAN-tf1-hmmscan",type="protein_match",start=Gene.Start,end=Gene.End,
    score=eval, strand="+", phase=".",target_id=paste0("DBD-Pfam|", HMMProfile), target_start=Profile.Start, target_end=Profile.End,target_strand=NA)]
  tf1.hmmscan.gff <- formatDataFrame(tf1.hmmscan.gff, sequenceontology_path)
}else{
  tf1.hmmscan.gff <- NULL
}
if(nrow(tf2.hmmscan)>0){
  tf2.hmmscan.gff <- tf2.hmmscan[,.(seqid=GeneID, source="Knutt-dbCAN-tf2-hmmscan",type="protein_match",start=Gene.Start,end=Gene.End,
    score=eval, strand="+", phase=".",target_id=paste0("DBD-SUPERFAMILY|", HMMProfile), target_start=Profile.Start, target_end=Profile.End,target_strand=NA)]
  tf2.hmmscan.gff <- formatDataFrame(tf2.hmmscan.gff, sequenceontology_path)
}else{
  tf2.hmmscan.gff <- NULL
}
if(nrow(stp.hmmscan)>0){
  stp.hmmscan.gff <- stp.hmmscan[,.(seqid=GeneID, source="Knutt-dbCAN-stp-hmmscan",type="protein_match",start=Gene.Start,end=Gene.End,
    score=eval, strand="+", phase=".",target_id=paste0("STP|", HMMProfile), target_start=Profile.Start, target_end=Profile.End,target_strand=NA)]
  stp.hmmscan.gff <- formatDataFrame(stp.hmmscan.gff, sequenceontology_path)
}else{
  stp.hmmscan.gff <- NULL
}
if(nrow(tcdb)>0){
  tcdb.gff <- tcdb[,.(seqid=GeneID, source="Knutt-tcdb-diamond",type="protein_match",start=Query.Start,end=Query.End,
    score=eval, strand="+", phase=".",target_id=paste0(db,"|",target), target_start=Subject.Start, target_end=Subject.End,
    target_strand=NA, Note=trimws(descr), Dbxref=paste0("TC:", TCFamily))]
  tcdb.gff <- formatDataFrame(tcdb.gff, sequenceontology_path)
}else{
  tcdb.gff <- NULL
}

gff <- as.data.table(rbind(diamond.gff,hotpep.gff,hmmscan.gff,tf1.hmmscan.gff,tf2.hmmscan.gff,stp.hmmscan.gff,tcdb.gff))
setkey(gff,seqid)
fwrite(gff,gff.tsv.path,sep="\t")



















#diamond.gff[,index:=.I]
#dbtargetcols = c("db", "eval", "score","Identity")
#diamond.gff = diamond.gff[,c("index","GeneID","Query.Start","Query.End","CAZyFamily",dbtargetcols),with=F]
#TODO qcov  !qcov!, identity
#diamond.gff.targets = diamond.gff[,c("index",dbtargetcols),with=F]
#diamond.gff.targets[,(cols):=lapply(.SD,as.character),.SD=dbtargetcols]
#setcolorder(diamond.gff.targets, c("index",dbtargetcols))
#diamond.gff.targets = melt(diamond.gff.targets,id.vars=c("index"))
#diamond.gff.targets[,element:=paste0(variable,":",value)]
#diamond.gff.targets[,element:=gffcol9escaper(element)]
#diamond.gff.targets=diamond.gff.targets[,.(CAZyDBdiamondtarget=paste(element,collapse=" ")),by=index]


#diamond.gff.score = diamond.gff[,c("index","eval"),with=F]

#diamond.gff = diamond.gff.targets[diamond.gff[,.SD,.SD=setdiff(colnames(diamond.gff),dbtargetcols)],on="index"]
#diamond.gff = melt(diamond.gff,id.vars=c("index","GeneID","Query.Start","Query.End"))
#diamond.gff = diamond.gff[,.(value=paste0(value,collapse=",")),by=c("GeneID","Query.Start","Query.End","variable")]
#diamond.gff[,element:=paste0(variable,"=",value)]
#diamond.gff = diamond.gff[,.(attributes=paste0(element,collapse=";")),by=c("index","GeneID","Query.Start","Query.End")]
#diamond.gff = diamond.gff.score[diamond.gff,on="index"]

