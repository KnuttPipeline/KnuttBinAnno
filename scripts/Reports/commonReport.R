#!/usr/bin/env Rscript

##
## commonReport.R - Utility functions used by the reports
##
## Knutt.org/KnuttReads2Bins


# Load packages every Report requires and set thread counts
commonOptions  <- function(threads){
    options(warn=2)
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(rlang))
    suppressPackageStartupMessages(require(data.table))
    suppressPackageStartupMessages(require(plotly))
    suppressPackageStartupMessages(require(parallel))
    suppressPackageStartupMessages(require(stringr))
    suppressPackageStartupMessages(require(DT))
    suppressPackageStartupMessages(require(tidyr))
    suppressPackageStartupMessages(require(tools))
    suppressPackageStartupMessages(require(matrixStats))
    if(is.na(threads))
      threads=detectCores()
    setDTthreads(threads)
    options(mc.cores = threads)
    theme_set(theme_linedraw())
    knitr::opts_chunk$set(fig.height=7) 
}

# From a data.table with elements like this "varX=valA;varY=valB"
# create a new one, recreating the original rows (like a join operation)
extractNamedValues <- function(data,colname="attributes",sep1=";",keyvalsep="=",valuepaste=","){
  require(data.table)
  dt=data.table(col=as.character(data[[colname]]))
  dt[,row__:=.I]
  dt=dt[!is.na(col),]
  dt=dt[,.(col=unlist(strsplit(col,split=sep1,fixed=T))),by=row__]
  dt[,c("key", "val"):=tstrsplit(col, split=keyvalsep,fixed=T)]
  dt[,col:=NULL]
  dt[,val:=trimws(val)]
  dt=dcast(dt,row__~key,value.var = "val",fun.aggregate = function(x)paste(x,collapse=valuepaste))
  data[,row__:=.I]
  res=dt[data,on="row__"]
  res[,row__:=NULL]
  return(res)
}

# Read the CDS file from metaerg
readCDS <- function(cds.file){
  suppressPackageStartupMessages(require(seqinr))
  cdsdat <- sapply(read.fasta(cds.file), attr, which="Annot")
  cdsdat <- sub(">([^ ]+) (.*)", "cds=\\1 metaergname=\\2", cdsdat)
  cdsdat <- strsplit(cdsdat, " (?=[^ ]+=)", perl=T)
  extractdat <- function(textvector){
    kvpairs <- strsplit(textvector, "=", fixed=T)
    v <- sapply(kvpairs, "[[", 2)
    names(v) <- sapply(kvpairs, "[[", 1)
    as.data.table(t(v))
  }
  cdsdat <- rbindlist(lapply(cdsdat, extractdat))
  cdsdat
}

# For DT tables create a link from a path
formatasLink <- function(x)gsub("^(.+/(.+)$)",paste0("<a  target=_blank href='\\1'>\\2</a>"),x,perl = T)

# Create a data.table from a list like this list(sample=list(R1=datatable,R2=datatable))
# (binds the data.tables and adds a sample and readdirection column)
deflateSampleReadData <- function(readlist,name){
  res <- rbindlist(unlist(lapply(names(readlist),unlist(function(sample)lapply(names(readlist[[sample]]),function(read)list(sample=sample,readdirection=read,readlist[[sample]][[read]])),recursive = F)),recursive = F))
  res
}

# Read data which is formatted as a list of file names and add a argument combination to each table 
# colvals is a data.frame like list (list of named lists)
readData <- function(filelist,colvals,readfun=fread){
  data <- lapply(filelist,readfun)
  for (i in seq_along(data)) {
     for (colname in names(colvals)){
       set(data[[i]], NULL ,colname,colvals[[colname]][[i]])
     }
  }
  data <- rbindlist(data)
  setcolorder(data, c(names(colvals)))
  data
}


# Split a string every "every" characters and insert sep
insertEvery <- function(strings,every=5,sep="\n")sapply(gsub(paste0("(.{",as.character(every),"})"), paste0("\\1",sep),strings),trimws)

# Calculate the boxplot stats
calculateBoxplotstats <- function(x){
  require(ggplot2)
  dat=boxplot.stats(x,do.out = F)
  list(ymin=dat$stats[[1]],lower=dat$stats[[2]],middle=dat$stats[[3]],upper=dat$stats[[4]],ymax=dat$stats[[5]],N=dat$n)
}

# Draw a htmltools compatible widget in a rmarkdown "asis" block
drawplotasis = function(title,subtitle,plot,height=NA,width=NA){
  require(rlang)
  if(!is.na(height))
    title=paste0(title," {data-height=",as.character(height),"}")
  if(!is.na(width))
    title=paste0(title," {data-width=",as.character(height),"}")
  cat(paste0("\n\n### ",title,"\n\n"))
  if("htmlwidget"%in% class(plot))
    plot=htmltools::tagList(list(a=plot))
  if(is.character(plot))
    cat(plot)
  else
    print(plot)
  cat(paste0("\n\n> ",subtitle,"\n\n"))
}
drawrow = function()cat("\n\nRow\n-------------------------------------\n\n")

# Generate the info for the table with the file info
genInfoBlock <- function(inputfiles){
  require(tools)
  finfos=file.info(inputfiles)
  infodf=data.table(`File Name`=sapply(inputfiles,basename),`File Path`=inputfiles,`MD5 Checksum`=md5sum(inputfiles),`Size (in MB)`=round(finfos$size/10^6,1),`Last modified`=finfos$mtime,`Created`=finfos$ctime,check.names = F)
  infodf
}

# KEGG Data
# Convert the given Ks to their KO counterpart (seperated by ;)
keggKOECConverter <- function(KOs){
  require(KEGGREST)
  require(data.table)
  tolookup=keggLink("enzyme",unique(KOs))
  tolookup=data.table(KO=sub("ko:(.+)","\\1",names(tolookup)),ECsin=sub("ec:(.+)","\\1",tolookup))
  tolookup=tolookup[,.(EC=paste0(ECsin,collapse=";")),by=KO]
  setkey(tolookup,KO)
  tolookup[KOs]$EC
}

# KEGG Data
# If an EC was moved to one other EC, rename it
cleanMovedECs <- function(ECs){
  require(KEGGREST)
  require(data.table)
  toreplace=sapply(unique(ECs),function(EC)tryCatch(keggGet(paste0("ec:",EC))[[1]]$NAME,error=function(e)NA,silent=T))
  toreplace=lapply(toreplace[!is.na(toreplace)],"[[",1)
  toreplace=toreplace[grep("Transferred to \\d+\\.\\d+\\.\\d+\\.\\d+$",toreplace)]
  ornames = names(toreplace)
  toreplace=gsub("Transferred to (\\d+\\.\\d+\\.\\d+\\.\\d+)$","\\1",toreplace)
  lookup=data.table(oldecs=ornames,newec=toreplace,key = c("oldecs"))
  res=lookup[ECs]
  res[is.na(newec),newec:=oldecs]
  res$newec
}

# KEGG Data
# Lookup the names of the given EC numbers
keggECNamelookup <- function(ECs,maxnamelen=60){
  require(KEGGREST)
  require(data.table)
  ecset=unique(ECs)
  entries=lapply(ecset,function(EC)tryCatch(keggGet(paste0("ec:",EC))[[1]],error=function(e)NULL,silent=T))
  sysnames=sapply(entries,function(en)en$SYSNAME[[1]])
  sysnames[sapply(sysnames,is.null)]=NA
  fullnames=sapply(entries,function(en)en$NAME[[1]])
  fullnames[sapply(fullnames,is.null)]=NA
  tolookup=data.table(EC=ecset,sysname=unlist(sysnames),name=unlist(fullnames))
  tolookup[is.na(sysnames),sysname:=name]
  tolookup[is.na(sysname)|sapply((sysname),nchar)>maxnamelen,sysname:=EC]
  setkey(tolookup,EC)
  tolookup[ECs]$sysname
}

# KEGG Data
# Lookup Ks for the given ECs (multiple results seperated by ;)
keggECKOConverter <- function(ECs){
  require(KEGGREST)
  require(data.table)
  tolookup=keggLink("ko",unique(ECs))
  tolookup=data.table(EC=sub("ec:(.+)","\\1",names(tolookup)),KOsin=sub("ko:(.+)","\\1",tolookup))
  tolookup=tolookup[,.(KO=paste0(KOsin,collapse=";")),by=EC]
  setkey(tolookup,EC)
  tolookup[ECs]$KO
}

# Extract the CAZy groups and ECs from the given vector
extractCAZy <- function(datacol){
  require(data.table)
  res=data.table(data=unique(datacol))
  res[,dataspl:=strsplit(data," ",fixed=T)]
  ECpattern="[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+"
  res[,ECs:=rep(list(),.N)]
  res[grepl(".",dataspl,fixed=T),ECs:=mclapply(dataspl,grep,pattern=ECpattern,perl=T,value=T)]
  res[,CazySubClasses:=mclapply(dataspl,grep,pattern=ECpattern,invert=T,perl=T,value=T)]
  res[,CazyClasses:=mclapply(CazySubClasses,sub,pattern="(\\D{2,3}).+",replacement="\\1")]
  res=res[datacol,on="data"]
  res[,dataspl:=NULL]
  res[,data:=NULL]
  res
}


unnestSingleColUnique<- function(dat,col){
  require(tidyr)
  require(dplyr)
  lookup = dat[[col]] %>% unique %>% tibble(datacol=sapply(.,paste0,collapse=";"),tounnest=.) %>% unnest("tounnest")
  result = dat %>% mutate(datacol=sapply(get(col),paste0,collapse=";")) %>% select(-all_of(col)) %>% full_join(lookup,by="datacol") %>% select(-c(datacol))%>% rename(!!col:="tounnest")
  result
}

# Order the "values" column as a factor using the values in "field".
# The "values" are also seperated by "grouping" (second axis in the heatmap)
orderHeatMap <- function(dat,values="CazySubClasses",field="rel",grouping="sample"){
  if(length(unique(unlist(dat[, c(values)])))==1){
    return(invisible(dat))
  }
  valmatrix = dat[,.(C__OL=sum(get(field))),by=c(grouping,values)]
  setnames(valmatrix,"C__OL",field)
  form=paste0(values,"~",grouping)
  valmatrix=as.data.frame(dcast(valmatrix,form,value.var=field,fill=0))
  rownames(valmatrix)=valmatrix[,1]
  valmatrix[,1]=NULL
  distmat=dist(valmatrix)
  clustering=as.dendrogram(hclust(distmat,method = "ward.D2"))
  #plot(clustering)
  dat[,(values):=factor(get(values),levels=c(labels(clustering)),ordered = T)]
  invisible(dat)
}


filtermax <- function(dt,groupcol=c("phylum"),filtercol="class",maxn=10,additionalsubcols=c("family"),additionalparentcols=c()){
  filtered = dt[,.(N=sum(N)),by=c(groupcol,filtercol)]
  setorder(filtered,-N)
  if(length(groupcol)>0){
    filtered = filtered[,head(get(filtercol),maxn-1),by=groupcol]
    setkey(filtered)
  }
  else{
    filtered = filtered[,head(get(filtercol),maxn-1)]
  }
  res = copy(dt)
  setkeyv(res,c(groupcol,filtercol))
  others = res[!filtered]
  res = res[filtered]
  for (toreplace in c(additionalsubcols,filtercol))  {
    others[,(toreplace):="Other"]
  }
  res = rbind(res,others)
  res = res[,.(N=sum(N)),by=c(additionalparentcols,groupcol,filtercol,additionalsubcols)]
}

drawSampleTaxplots <- function(taxdata_sample,baseheight=576){
  domaindata = taxdata_sample[,.(N=sum(N)),by=c("sample","superkingdom")]
  domaindata[,rel:=N/sum(N)*100,by="sample"]
  #domaindata=as.data.table(complete(domaindata,sample,superkingdom,fill=list(N=0,rel=0)))
  orderHeatMap(domaindata,values="superkingdom",field="rel",grouping="sample")
  orderHeatMap(domaindata,values="sample",field="rel",grouping="superkingdom")
  domainplot = ggplot(domaindata)+aes_string(x="sample",y="rel",fill="superkingdom")+ geom_bar(stat="identity")+ xlab("Sample") + scale_y_continuous(labels = function(x) paste0(x, "%"))+ ylab("Percentage of classified Reads")
  
  drawplotasis("Application Note",NA,"In the following plots the entries (samples and taxonomic levels) have been sorted using Ward linkage clustering with euclidian distance. This was always done after agglomeration.")
  drawrow()
  drawplotasis("Domains","The domain composition of each sample",ggplotly(domainplot),baseheight)
  
  phylumdata = taxdata_sample[,.(N=sum(N)),by=c("sample","superkingdom","phylum")]
  phylumdata = filtermax(phylumdata,groupcol = c("superkingdom"),filtercol = "phylum",maxn=20,additionalsubcols = c(),additionalparentcols=c("sample"))
  phylumdata[,rel:=N/sum(N)*100,by="sample"]
  phylumdata=as.data.table(complete(phylumdata,sample,nesting(superkingdom,phylum),fill=list(N=0,rel=0)))
  phylumdata = split(phylumdata,phylumdata$superkingdom)
  bla = lapply(phylumdata,function(pdat)if(length(unique(pdat$phylum))>1)orderHeatMap(pdat,values="phylum",field="rel",grouping="sample"))
  bla = lapply(phylumdata,function(pdat)if(length(unique(pdat$sample))>1)orderHeatMap(pdat,values="sample",field="rel",grouping="phylum"))
  phylumplots = lapply(phylumdata,function(phylumdata)ggplot(phylumdata)+aes_string(x="sample",y="phylum",fill="rel") +labs(fill="% of filtered\n Reads")+ geom_raster() + xlab("Sample") + ylab("Phylum"))
  
  for (i in seq_along(phylumplots)) {
    plotname = names(phylumplots)[[i]]
    if(plotname != "Unclassified"&&plotname != "(Taxonomy Gap)"&&plotname != "Viruses"){
      if(i%%2==0){
        drawrow()
        drawn=T
      }else
        drawn=F
      drawplotasis(plotname,"The 20 (or less) most abundant phyla in this domain (over all samples). Smaller phyla have been aggregated into \"Other\".",ggplotly(phylumplots[[i]]),baseheight)
    }
  }
  
  classdatat = taxdata_sample[,.(N=sum(N)),by=c("sample","superkingdom","phylum","class","family")] 
  classdatat = filtermax(classdatat,groupcol = c("sample","superkingdom"),filtercol = "phylum",maxn=10,additionalsubcols = c("class","family"),additionalparentcols = c())
  classdatat = filtermax(classdatat,groupcol = c("sample","superkingdom","phylum"),filtercol = "family",maxn=5,additionalsubcols = c(),additionalparentcols = c("phylum","class"))
  classdatat[,rel:=N/sum(N)*100,by="sample"]
  #classdatat=as.data.table(complete(classdatat,sample,nesting(superkingdom,phylum,class,family),fill=list(N=0,rel=0)))
  classdatat[,phylclass:=paste0(phylum,".",class)]
  classdatat = split(classdatat,list(classdatat$superkingdom))
  bla = lapply(classdatat,function(pdat)if(length(unique(pdat$phylclass))>1)orderHeatMap(pdat,values="phylclass",field="rel",grouping="sample"))
  bla = lapply(classdatat,function(pdat)if(length(unique(pdat$sample))>1)orderHeatMap(pdat,values="sample",field="rel",grouping="phylclass"))
  classplots = lapply(classdatat,function(dt)ggplot(dt)+aes_string(x="phylclass",y="rel",fill="family") + geom_bar(stat="identity",position = "stack") + xlab("Phylum.Class") + ylab("% of filtered Reads") + facet_grid(cols = vars(sample))+coord_flip())
  
  if(!drawn)
    drawrow()
  classplots=classplots[!names(classplots) %in% c("Unclassified","(Taxonomy Gap)")]
  for (i in seq_along(classplots)) {
    plotname = names(classplots)[[i]]
    drawplotasis(paste0("Families in ",plotname),"The 10 largest phyla in every sample and at most 5 familes in those phyla each. ",ggplotly(classplots[[i]]),baseheight*1.75)
    if(i!=length(classplots))
      drawrow()
  }
  
}




drawBlastPlots <- function(hitdata,readanno_sampled_overview,baseheight=576){
  hitcount.data = hitdata[,.N,by="sample"]
  hitcount.data = hitcount.data[readanno_sampled_overview,on="sample"]
  hitcount.data=hitcount.data[,.(rel=sum(N)/reads*100),by=sample]
  hitcount = ggplot(hitcount.data) + aes(x=sample,y=rel,fill=sample) +theme(legend.position = "none")+  geom_bar(stat="identity")+ xlab("Sample") + scale_y_continuous(labels = function(x) paste0(x, "%"))+ ylab("Percentage of classified Reads")
  
  pidentlength = ggplot(hitdata) + aes(x=length,y=pident) + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) + facet_grid(cols = vars(sample)) + xlab("Alignment Length") + ylab("Percent identity") + labs(fill="(Estimated)\nDensity")
  logevallength = ggplot(hitdata) + aes(x=length,y=log(evalue)) + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) + facet_grid(cols = vars(sample)) + xlab("Alignment Length") + ylab("ln(Evalue)") + labs(fill="(Estimated)\nDensity")
  hitdata[,scov:=length/slen*100]
  covlength = ggplot(hitdata) + aes(x=length,y=scov) + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) + facet_grid(cols = vars(sample)) + xlab("Alignment Length") + ylab("Subject coverage") + labs(fill="(Estimated)\nDensity")  + scale_y_continuous(labels = function(x) paste0(x, "%"))
  subjecthitcounts = hitdata[,.N,by=c("sname","sample")]
  subjecthitcounts[,rel:=N/sum(N)*100,by="sample"]
  persubjectrelcount =  ggplot(subjecthitcounts) + aes(x=sample,y=rel,color=sample) + geom_boxplot() + theme(legend.position = "none") + xlab("Sample") + ylab("Hits on a single subject of all") + scale_y_continuous(labels = function(x) paste0(x, "%"))
  
  
  drawplotasis("Number of hits","The percentage of the reads which mapped to this database.",ggplotly(hitcount))
  drawrow()
  drawplotasis("Percent identity to alignment length","The distribution of the read lengths strongly influences this and the following plots.",pidentlength,width=baseheight*2)
  drawplotasis("Log Evalue to alignment length","The log of the evalue corresponds to the bit score.",logevallength)
  drawrow()
  drawplotasis("Subject Coverage to alignment length","How much of the subject sequence length does the alignment cover?",covlength)
  drawplotasis("Hit count per subject","The number of reads in a sample that mapped to the same subject",persubjectrelcount)
}


# For Sourmash
fixTaxCounts <- function(taxdat){
  samplen <- unique(taxdat$sample)
  taxdat[, sample:=NULL]
  for(curdepth in 2:length(taxcols)){
    subtract <- taxdat[depth==curdepth, ]
    if(nrow(subtract)==0)
      next
    toreplace <- taxcols[curdepth:length(taxcols)]
    subtract[, c(toreplace):=""]
    subtract <- subtract[, .(minus=sum(count)), by=taxcols]
    taxdat <- merge(taxdat, subtract, on=taxcols, all.x=TRUE)
    taxdat[is.na(minus), minus:=0]
    taxdat[, count:=count-minus]
    taxdat[, minus:=NULL]
    taxdat <- taxdat[count!=0,]
  }
  taxdat[, depth:=NULL]
  setcolorder(taxdat, c("count"))
  taxdat[,(taxcols):=lapply(.SD, function(col)sub(".__(.+)","\\1",col)),.SD=taxcols]
  taxdat[, sample:=samplen]
  taxdat
}