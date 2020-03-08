#!/usr/bin/env Rscript
options(warn=2)

#options(echo=T)

library(data.table)
library(tidyr)
library(seqinr)

hotpep.path = "output/BinAnnotation/dbCAN/Larrelt/Larrelt.1_hotpep.txt"
hotpep.ref.dirpath = "output/BinAnnotation/dbCAN/HotPep/CAZY_PPR_patterns"
diamond.path = "output/BinAnnotation/dbCAN/Larrelt/Larrelt.1_diamond.tsv"
diamond.ref.data.path = "output/refDBs/readanno/prot/CAZyDB.tsv"
hmmscan.path = "output/BinAnnotation/dbCAN/Larrelt/Larrelt.1_hmmscan.out"
cds.path = "output/BinAnnotation/metaerg/Larrelt/Larrelt.1/data/cds.faa"
gff.tsv.path = "test.tsv" 

tf1.hmmscan.path = "output/BinAnnotation/dbCAN/Larrelt/Larrelt.1_tf1_hmmscan.out"
tf2.hmmscan.path = "output/BinAnnotation/dbCAN/Larrelt/Larrelt.1_tf2_hmmscan.out"
stp.hmmscan.path = "output/BinAnnotation/dbCAN/Larrelt/Larrelt.1_stp_hmmscan.out"

tcdb.path = "output/BinAnnotation/dbCAN/Larrelt/Larrelt.1_tcdb_diamond.tsv"
tcdb.ref.path = "output/refDBs/dbCAN/tcdb.fa"

hmmscan.eval = 1e-15
hmmscan.coverage = 0.35
tf.hmmscan.eval = 1e-4
tf.hmmscan.coverage = 0.35
stp.hmmscan.eval = 1e-4
stp.hmmscan.coverage = 0.3




if(exists("snakemake")){
  hotpep.path = snakemake@input[["hotpep"]]
  hotpep.ref.dirpath = snakemake@params[["hotpepref"]]
  diamond.path = snakemake@input[["diamond"]]
  diamond.ref.data.path = snakemake@input[["diamond_ref"]]
  hmmscan.path = snakemake@input[["hmmscan"]]
  cds.path = snakemake@input[["query"]]
  gff.tsv.path = snakemake@output[[1]]
  hmmscan.eval = snakemake@params[["hmmscaneval"]]
  hmmscan.coverage = snakemake@params[["hmmscancov"]]
  tf1.hmmscan.path = snakemake@input[["tf1hmmscan"]]
  tf2.hmmscan.path = snakemake@input[["tf2hmmscan"]]
  stp.hmmscan.path = snakemake@input[["stphmmscan"]]
  tf.hmmscan.eval = snakemake@params[["tf_hmmscaneval"]]
  tf.hmmscan.coverage = snakemake@params[["tf_hmmscancov"]]
  stp.hmmscan.eval = snakemake@params[["stp_hmmscaneval"]]
  stp.hmmscan.coverage = snakemake@params[["stp_hmmscancov"]]
  tcdb.path = snakemake@input[["tcdb"]]
  tcdb.ref.path = snakemake@input[["tcdb_ref"]]
}



hmmscan.regex = "^([A-Z]{2,3}\\d{1,4}).*?((?:_\\d{1,3})?)$"
hotpepcols = c("CAZyFamily","PPRSubfamily","GeneID","Frequency","Hits","SignaturePeptides")
if(file.info(hotpep.path)$size>0){
  hotpep =  fread(hotpep.path,col.names=hotpepcols,header=F,sep="\t")
}else{
  hotpep = data.table(CAZyFamily = character(), PPRSubfamily = character(), GeneID = character(), Frequency = character(), Hits = character(), SignaturePeptides = character())
}
hotpep.ref.files = list.files(paste0(hotpep.ref.dirpath,"/CAZY_PPR_patterns"),pattern="*_group_ec.txt",recursive=T,full.names=T)
hotpep.ref = rbindlist(lapply(hotpep.ref.files, function(file)cbind(file=basename(file),fread(file,header=F,sep="\t"))))
setnames(hotpep.ref,c("V1","V2"),c("PPRSubfamily","EC"))
#hotpep.ref = hotpep.ref[!is.na(EC)&EC!="",]
hotpep.ref[,EC:=gsub("x","-",EC,fixed=T)]
hotpep.ref[,EC:=strsplit(EC,",",fixed=T)]
hotpep.ref[,EC:=lapply(EC,function(ECs)unlist(lapply(ECs,function(ec)gsub("([^:]+):.+","\\1",ec))))]
hotpep.ref[,CAZyFamily:=gsub("(.+)_group_ec\\.txt","\\1",file)]
hotpep = hotpep.ref[hotpep,on=c("CAZyFamily","PPRSubfamily")]
hotpep[,CAZyFamily:=strsplit(CAZyFamily,",",fixed=T)]

if(file.info(diamond.path)$size>0){
  diamond = fread(diamond.path,col.names=c("GeneID","CAZyID","Identity","Length","Mismatches","Gap open","Query.Start","Query.End","Subject.Start","Subject.End","eval","score"), header=F)
}else{
  diamond = data.table(GeneID=character(), CAZyID=character(), Identity=numeric(), Length=numeric(),`Gap open`=numeric(),Query.Start=numeric(),Query.End=numeric(),Subject.Start=numeric(),Subject.End=numeric(),eval=numeric(),score=numeric())
}
diamond.ref.data = fread(diamond.ref.data.path,header=T)
setnames(diamond.ref.data,c("sseqid","CAZyECs"),c("CAZyID","data"))

#sseqid	taxid	CAZyECs	superkingdom	phylum	class	order	family	genus	species

diamond.ref.data[,data:=strsplit(data," ",fixed=T)]
ECpattern="[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+"
diamond.ref.data[grepl(".",data,fixed=T),EC:=lapply(data,grep,pattern=ECpattern,perl=T,value=T)]
diamond.ref.data[,CAZyFamily:=lapply(data,grep,pattern=ECpattern,invert=T,perl=T,value=T)]
diamond.ref.data[,data:=NULL]

setkey(diamond.ref.data,"CAZyID")
setkey(diamond,"CAZyID")
diamond = diamond.ref.data[diamond]

read.hmm <- function(file,feval,fcov){
  lines = as.integer(system2("grep",args = c("-o","^[^#]*", file," | wc -l | awk '{print $1}'"),stdout = TRUE))
  if(lines>0){
    hmmscan.cmd = paste0("grep -v '^#' ",shQuote(file)," |  awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_.\"\n\";}}'")
    hmmscan = fread(cmd=hmmscan.cmd,header=F,col.names=c("HMMProfile","ProfileLength","GeneID","GeneLength","eval","Profile.Start","Profile.End","Gene.Start","Gene.End"))
    hmmscan[,Coverage:=(Profile.End-Profile.Start)/ProfileLength]
    hmmscan = hmmscan[eval<feval & Coverage>=fcov]
  }else{
    hmmscan = data.table(HMMProfile=character(),ProfileLength=numeric(),GeneID=character(),GeneLength=numeric(),eval=numeric(),Profile.Start=numeric(),Profile.End=numeric(),Gene.Start=numeric(),Gene.End=numeric())
  }
  hmmscan
}
hmmscan = read.hmm(hmmscan.path,hmmscan.eval,hmmscan.coverage)
hmmscan[,basename:=gsub("(.+)\\.p?hmm","\\1",HMMProfile)]
hmmscan[grepl(hmmscan.regex,basename),CAZyFamily:=strsplit(sub(hmmscan.regex,"\\1\\2",basename),",",fixed=T)]
hmmscan[!grepl(hmmscan.regex,basename),CAZyFamily:=strsplit("",",",fixed=T)]


tf1.hmmscan = read.hmm(tf1.hmmscan.path,tf.hmmscan.eval,tf.hmmscan.coverage)
tf2.hmmscan = read.hmm(tf2.hmmscan.path,tf.hmmscan.eval,tf.hmmscan.coverage)
stp.hmmscan = read.hmm(stp.hmmscan.path,stp.hmmscan.eval,stp.hmmscan.coverage)

if(file.info(tcdb.path)$size>0){
  tcdb = fread(tcdb.path,col.names=c("GeneID","TCDBID","Identity","Length","Mismatches","Gap open","Query.Start","Query.End","Subject.Start","Subject.End","eval","score"), header=F)
}else{
  tcdb =  data.table(GeneID=character(),TCDBID=character(),Identity=numeric(),Length=numeric(),Mismatches=numeric(),`Gap open`=numeric(),Query.Start=numeric(),Query.End=numeric(),Subject.Start=numeric(),Subject.End=numeric(),eval=numeric(),score=numeric())
}
tcdb.ref = data.table(txt=system(paste0("grep '>' ",shQuote(tcdb.ref.path)),intern=T))
tcdb.ref[,TCDBID:=sub("^>([^ ]+)( ?.*)$","\\1",txt)]
tcdb.ref[,descr:=sub("^>([^ ]+)( ?.*)$","\\2",txt)]
tcdb.ref[,c("asc","db","target","TCFamily"):=tstrsplit(TCDBID,"|",fixed=T)]
tcdb=tcdb.ref[tcdb,on="TCDBID"]



diamond.hits = diamond[,.(GeneID,CAZyFamily),] #Only one result per query allowed
hotpep.hits = hotpep[,.(CAZyFamily=paste(CAZyFamily,collapse=",")),by=GeneID]
hmmscan.hits = hmmscan[,.(CAZyFamily=paste(CAZyFamily,collapse=",")),by=GeneID]

hits = rbindlist(list(diamond.hits,hotpep.hits,hmmscan.hits))
hits = hits[,.(CAZyFamily=paste(unique(unlist(CAZyFamily)),collapse=",")),by=GeneID]
setkey(hits,"GeneID")

# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

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


addECandCAZyDbxRef=function(dt){
   dt[,cazydb:=lapply(CAZyFamily,function(cazylist)paste0(lapply(cazylist,function(cazy)paste0("CAZY:",trimws(cazy))),collapse=","))] 
   dt[,ecdb:=lapply(EC,function(cazylist)paste0(lapply(cazylist,function(cazy)paste0("EC:",trimws(cazy))),collapse=","))]
   dt[sapply(CAZyFamily,length)>0|sapply(EC,length)>0,attributes:=paste0(attributes,";Dbxref=")]
   dt[sapply(CAZyFamily,length)>0,attributes:=paste0(attributes,cazydb)]
   dt[sapply(CAZyFamily,length)>0&sapply(EC,length)>0,attributes:=paste0(attributes,",")]
   dt[sapply(EC,length)>0,attributes:=paste0(attributes,ecdb)]
}

diamond.gff = copy(diamond)
diamond.gff[,db:=paste0("CAZyDB.fa|",gsub(" ","%20",CAZyID,fixed=T)," ",as.character(Subject.Start)," ",as.character(Subject.End))]
diamond.gff[,attributes:=paste0("target=",db)]
addECandCAZyDbxRef(diamond.gff)
if(nrow(diamond.gff)>0){
  diamond.gff = diamond.gff[,.(seqid=gffseqidEncoder(GeneID),source="MetaKnutt-dbCAN-diamond",type="protein match",start=Query.Start,end=Query.End,score=eval,strand="+",phase=".",attributes)]
}else{
  diamond.gff = NULL
}

cds.lengths = rbindlist(lapply(read.fasta(cds.path),function(entry)data.table(GeneID=attr(entry,"name"),length=length(entry))))
hotpep.gff = cds.lengths[hotpep,on="GeneID"]
hotpep.gff[,attributes:=paste0("HotPepHits=",as.character(Hits),";HotPepFrequency=",as.character(Frequency),";HotPepCAZyPPRSubfamily=",PPRSubfamily,";SignaturePeptides=",SignaturePeptides)]
addECandCAZyDbxRef(hotpep.gff)
hotpep.gff[sapply(EC,length)>0,attributes:=paste0(attributes,";Note=ECs derived from the CAZy (Sub)Family+PPR as a whole")]
if(nrow(hotpep.gff)>0){
  hotpep.gff = hotpep.gff[,.(seqid=gffseqidEncoder(GeneID),source="MetaKnutt-dbCAN-HotPep",type="protein match",start=1,end=length,score=".",strand="+",phase=".",attributes)]
}else{
  hotpep.gff = NULL
}

hmmscan.gff = copy(hmmscan)
hmmscan.gff[,db:=paste0("CAZyDB.hmm|",gsub(" ","%20",HMMProfile,fixed=T)," ",as.character(Profile.Start)," ",as.character(Profile.End))]
hmmscan.gff[,attributes:=paste0("target=",db)]
hmmscan.gff[,cazydb:=lapply(CAZyFamily,function(cazylist)paste0(lapply(cazylist,function(cazy)paste0("CAZY:",trimws(cazy))),collapse=","))]
hmmscan.gff[sapply(CAZyFamily,length)>0,attributes:=paste0(attributes,";Dbxref=",cazydb)]
hmmscan.gff[sapply(CAZyFamily,length)==0,attributes:=paste0(attributes,";Note=Model name was not parseable")]
if(nrow(hmmscan.gff)>0){
  hmmscan.gff = hmmscan.gff[,.(seqid=gffseqidEncoder(GeneID),source="MetaKnutt-dbCAN-hmmscan",type="protein match",start=Gene.Start,end=Gene.End,score=eval,strand="+",phase=".",attributes)]
}else{
  hmmscan.gff = NULL
}

tf1.hmmscan.gff = copy(tf1.hmmscan)
tf1.hmmscan.gff[,db:=paste0("DBD-Pfam|",gsub(" ","%20",HMMProfile,fixed=T)," ",as.character(Profile.Start)," ",as.character(Profile.End))]
tf1.hmmscan.gff[,attributes:=paste0("target=",db)]
if(nrow(tf1.hmmscan.gff)>0){
  tf1.hmmscan.gff = tf1.hmmscan.gff[,.(seqid=gffseqidEncoder(GeneID),source="MetaKnutt-dbCAN-tf1-hmmscan",type="protein match",start=Gene.Start,end=Gene.End,score=eval,strand="+",phase=".",attributes)]
}else{
  tf1.hmmscan.gff = NULL
}

tf2.hmmscan.gff = copy(tf2.hmmscan)
tf2.hmmscan.gff[,db:=paste0("DBD-SUPERFAMILY|",gsub(" ","%20",HMMProfile,fixed=T)," ",as.character(Profile.Start)," ",as.character(Profile.End))]
tf2.hmmscan.gff[,attributes:=paste0("target=",db)]
if(nrow(tf2.hmmscan.gff)>0){
  tf2.hmmscan.gff = tf2.hmmscan.gff[,.(seqid=gffseqidEncoder(GeneID),source="MetaKnutt-dbCAN-tf2-hmmscan",type="protein match",start=Gene.Start,end=Gene.End,score=eval,strand="+",phase=".",attributes)]
}else{
  tf2.hmmscan.gff = NULL
}

stp.hmmscan.gff = copy(stp.hmmscan)
stp.hmmscan.gff[,db:=paste0("STP|",gsub(" ","%20",HMMProfile,fixed=T)," ",as.character(Profile.Start)," ",as.character(Profile.End))]
stp.hmmscan.gff[,attributes:=paste0("target=",db)]
if(nrow(stp.hmmscan.gff)>0){
  stp.hmmscan.gff = stp.hmmscan.gff[,.(seqid=gffseqidEncoder(GeneID),source="MetaKnutt-dbCAN-stp-hmmscan",type="protein match",start=Gene.Start,end=Gene.End,score=eval,strand="+",phase=".",attributes)]
}else{
  stp.hmmscan.gff = NULL
}

tcdb.gff = copy(tcdb)
tcdb.gff[,db:=paste0(paste0(db,target,sep="|")," ",as.character(Subject.Start)," ",as.character(Subject.End))]
tcdb.gff[,attributes:=paste0("target=",db)]
tcdb.gff[,attributes:=paste0(attributes,";Name=",gffcol9escaper(descr))]
tcdb.gff[,attributes:=paste0(attributes,";Dbxref=TC:",TCFamily)]
if(nrow(tcdb.gff)>0){
  tcdb.gff = tcdb.gff[,.(seqid=gffseqidEncoder(GeneID),source="MetaKnutt-tcdb-diamond",type="protein match",start=Query.Start,end=Query.End,score=eval,strand="+",phase=".",attributes)]
}else{
  tcdb.gff = NULL
} 

gff = rbind(diamond.gff,hotpep.gff,hmmscan.gff,tf1.hmmscan.gff,tf2.hmmscan.gff,stp.hmmscan.gff,tcdb.gff)
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

