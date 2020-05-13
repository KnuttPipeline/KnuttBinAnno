suppressPackageStartupMessages(library(ontologyIndex))

escapeGFF3Text <- function(textvector){
  textvector <- gsub("\t", "%09", textvector, fixed=T)
  textvector <- gsub("\n", "%0A", textvector, fixed=T)
  textvector <- gsub("\r" ,"%0D", textvector, fixed=T)
  textvector <- gsub("%", "%25", textvector, fixed=T)
  for(char in c(1:31,127)){
    textvector <- gsub(intToUtf8(char), paste0("%", format(as.hexmode(char), width=2, upper.case=T)), textvector, fixed=T, useBytes=T)
  }
  textvector
}

unescapeGFF3Text <- function(textvector){
  textvector <- gsub("%09", "\t", textvector, fixed=T)
  textvector <- gsub("%0A", "\n", textvector, fixed=T)
  textvector <- gsub("%0D", "\r" , textvector, fixed=T)
  textvector <- gsub("%25", "%", textvector, fixed=T)
  for(char in c(1:31,127)){
    textvector <- gsub(paste0("%", format(as.hexmode(char), width=2, upper.case=T)), intToUtf8(char), textvector, fixed=T, useBytes=T)
  }
  textvector
}

escapeGFF3Col9 <- function(textvector){
  textvector <- escapeGFF3Text(textvector)
  textvector <- gsub(";", "%3B", textvector, fixed=T)
  textvector <- gsub("=", "%3D", textvector, fixed=T)
  textvector <- gsub("&", "%26", textvector, fixed=T)
  textvector <- gsub(",", "%2C", textvector, fixed=T)
  textvector
}

unescapeGFF3Col9 <- function(textvector){
  textvector <- unescapeGFF3Text(textvector)
  textvector <- gsub("%3B", ";", textvector, fixed=T)
  textvector <- gsub("%3D", "=", textvector, fixed=T)
  textvector <- gsub("%26", "&", textvector, fixed=T)
  textvector <- gsub("%2C", ",", textvector, fixed=T)
  textvector
}

escapeGFF3Col1 <- function(textvector){
  charstoescape <- gregexpr("[^a-zA-Z0-9.:^*$@!+_?-| ]", textvector)
  charstoescape <- unique(unlist(regmatches(textvector, charstoescape)))
  for(char in charstoescape){
      textvector <- gsub(char, paste0("%", format(as.hexmode(utf8ToInt(char)), width=2, upper.case=T)), textvector, fixed=T, useBytes=T)
  }
  textvector
}

unescapeGFF3Col1 <- function(textvector){
  sapply(textvector, URLdecode)
}

fixTypeNames <- function(typevector){
    oldvals <- c("bac_23SrRNA","protein match")
    newvals <- c("rRNA_23S","protein_match")
    for(i in seq_along(oldvals)){
        typevector[typevector==oldvals[[i]]] <- newvals[[i]]
    }
    typevector
}

# Description from https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
formatDataFrame <- function(inputframe, sequenceontologyfile="so.obo"){
    result <- data.frame(seqid=escapeGFF3Col1(inputframe$seqid), stringsAsFactors=F)
    # The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature.
    result$source <- escapeGFF3Text(inputframe$source)
    # The type of the feature (previously called the "method"). This is constrained to be either a term from the Sequence Ontology or an SO accession number.
    sequenceontology <- get_ontology(sequenceontologyfile)
    valid_ids <- get_descendants(sequenceontology, "SO:0000110")
    valid_names <- unname(sequenceontology$name[valid_ids])
    inputframe$type <- fixTypeNames(inputframe$type)
    for (val in inputframe$type) {
       if(! val %in% c(valid_names, valid_ids)){
           stop(paste0(val, " is not a valid sequence ontology term or name!"))
       }
    }
    result$type <- inputframe$type
    # positive 1-based integer coordinate, Start is always less than or equal to end
    if(! all(inputframe$start <= inputframe$end)){
        stop("All start positions need to be smaller or equal to the ends!")
    }
    result$start <- inputframe$start
    result$end <- inputframe$end
    # Floating point number,  E-values be used for sequence similarity features, and that P-values be used for ab initio gene prediction features
    # Missing values: .
    result$score <- as.character(inputframe$score)
    # +, -, .(not stranded), ?(unknown)
    if(! all(inputframe$strand %in% c("+", "-", ".", "?"))){
        stop("Not allowed strand character found!")
    }
    result$strand <- inputframe$strand
    # the phase indicates where the next codon begins relative to the 5' end, CDS mandatory
    cdsterms <- c(sequenceontology$name[get_descendants(sequenceontology, "SO:0000316")], get_descendants(sequenceontology, "SO:0000316"))
    if(! all(as.character(inputframe$phase[inputframe$type %in% cdsterms]) %in% c("0", "1", "2"))){
        stop("All CDS entries require a valid phase!")
    }
    result$phase <- inputframe$phase
    if("Parent" %in% colnames(inputframe)){
        if(!all(unlist(attributes$Parent) %in% unlist(attributes$ID))){
            exit("All Parent values need to be IDs")
        }
    }
    if(! "attributes" %in% colnames(inputframe)){
        attributes <- setdiff(colnames(inputframe), c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "target_id", "target_start", "target_end", "target_strand"))
        attributes_vals <- inputframe[, attributes, with=F]
        if(any(c("target_id", "target_start", "target_end", "target_strand") %in% colnames(inputframe))){
            formatTarget <- function(rowi, df){
                if(is.na(df[rowi, "target_id"]))
                    return(list())
                target_id <- gsub(" ", "%20", as.character(unlist(df[rowi, "target_id"])), fixed=T)
                target_start <- as.character(unlist(df[rowi, "target_start"]))
                target_end <- as.character(unlist(df[rowi, "target_end"]))
                target_strand <- ifelse(is.na(unlist(df[rowi, "target_strand"])), "", paste0(" ", as.character(unlist(df[rowi, "target_strand"]))))
                paste0(target_id, " ", target_start, " ", target_end, target_strand)
            }
            attributes_vals <- cbind(attributes_vals,target=lapply(1:nrow(inputframe), formatTarget , df=inputframe))
        }
        attributes_vals <- lapply(attributes_vals, function(col)unname(sapply(col, function(row)paste0(escapeGFF3Col9(as.character(unlist(row))), collapse=","))))
        attributes_vals <- mapply(function(n, col)ifelse(col=="NA"|is.na(col)|col=="", NA, paste0(n, "=", col)), names(attributes_vals), attributes_vals, SIMPLIFY=FALSE)
        attributes_vals <- as.data.frame(x=attributes_vals)
        result$attributes <- apply(attributes_vals,1,function(vals)paste0(vals[!is.na(vals)], collapse=";"))
    }else{
        result$attributes <- inputframe$attributes
    }
    result[is.na(result)] <- "."
    result <- result[order(result$seqid, result$start, result$end),]
    result
}

parseDbxref <- function(textvector){
    parseDf <- function(el){
        if(length(el)==0)
            return(list())
        sp <- strsplit(el, ":", fixed=T)
        lens <- sapply(sp,length)
        key <- sapply(sp, function(keyval)trimws(keyval[[1]]))
        val <- sapply(sp, function(keyval)trimws(keyval[[2]]))
        res <- data.frame(key, val, stringsAsFactors=F)
        res <- split(res, res$key)
        res <- lapply(res, "[[", "val")
        res
    }
    rowlists <- lapply(textvector,parseDf)
    allcols <- unique(unlist(lapply(rowlists, names)))
    fixDF <- function(rowlist){
        missingcols <- setdiff(allcols, names(rowlist))
        if(length(missingcols)>0)
            rowlist[missingcols] <- list(character())
        rowlist <- rowlist[allcols]
        rowlist
    }
    result <- as.data.frame(do.call(rbind, lapply(rowlists, fixDF)))
}

parseCol9 <- function(textvector){
    getKeyVal <- function(el){
        sp <- strsplit(el, "=", fixed=T)
        res <- sapply(sp, function(keyval)trimws(keyval[[2]]))
        res <- lapply(res, function(val)unescapeGFF3Col9(trimws(strsplit(val, ",", fixed=T)[[1]])))
        names(res) <-sapply(sp, function(keyval)trimws(keyval[[1]]))
        res
    }
    rowlists <- lapply(strsplit(textvector, ";", fixed=T), getKeyVal)
    allcols <- unique(unlist(lapply(rowlists, names)))
    fixDF <- function(rowlist){
        missingcols <- setdiff(allcols, names(rowlist))
        if(length(missingcols)>0)
            rowlist[missingcols] <- list(character())
        rowlist <- rowlist[allcols]
        rowlist
    }
    result <- as.data.frame(do.call(rbind, lapply(rowlists, fixDF)))
    #dbCANres[is.na(Dbxref), Dbxref:=list(list(NULL))]
    #listcols <- colnames(result)[sapply(result, class)=="list"]
    #for(col in listcols){
    #    result[[col]][is.na(result[[col]])] <- list(NULL)
    #}
    if("Dbxref" %in% colnames(result)){
        result <- cbind(result, parseDbxref(result$Dbxref))
        result$Dbxref <- NULL
    }
    if("Ontology_term" %in% colnames(result)){
        result <- cbind(result, parseDbxref(result$Ontology_term))
        result$Ontology_term <- NULL
    }
    for(col in colnames(result)){
        lens <- sapply(result[[col]], length)
        if(all(lens == 0 | lens == 1) & is.list(result[[col]])){
            result[[col]][lens == 0] <- NA
            result[[col]] <- unlist(result[[col]])
        }
    }
    result
}

readBaseGFF <- function(df){
    if(nrow(df)==0){
      df$attributes <- NULL
      return(df)
    }
    df$seqid <- unescapeGFF3Col1(df$seqid)
    df$source <- unescapeGFF3Text(df$source)
    att <- df$attributes
    df$attributes <- NULL
    df <- cbind(df, parseCol9(att))
    df
}


readGFF <- function(file){
    tempfile <- tempfile()
    on.exit(unlink(tempfile))
    cmd <- paste0("sed -n '/##FASTA/q;p' ", file, " | sed /^#/d > ", tempfile)
    system(cmd)
    result <- read.table(tempfile, header=F, sep="\t", quote="", stringsAsFactors=F)
    unlink(tempfile)
    colnames(result) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    readBaseGFF(result)
}

#tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14])
#d <- readGFF("output/MetaErg/bin.12/data/master.gff.txt")
#sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt

# regiondf: seqid, start, end
writeGFF <- function(df, fastafile, file, sequenceontologyfile="so.obo"){
    df <- formatDataFrame(df, sequenceontologyfile=sequenceontologyfile)
    dflist <- split(df, gffbase$seqid)
    regiondf <- do.call(rbind, lapply(dflist,function(df)data.frame(seqid=df$seqid[[1]],start=min(df$end),end=max(df$end))))
    write("##gff-version 3.2.1", file=file)
    regionlines <- paste0("##sequence-region", " ", regiondf$seqid, " ", regiondf$start, " ", regiondf$end)
    for(region in seq_along(regionlines)){
        write(regionlines[[region]], file=file, append=TRUE)
        write.table(dflist[[region]], file, append=T, quote=F, sep="\t", row.names=F, col.names=F)
    }
    write("##FASTA", file=file, append=TRUE)
    system(paste0("cat ",fastafile," >> ",file))
}

combineCols <- function(dat){
   dat <- do.call(function(...)mapply(function(...)unlist(unname(c(...))),...),dat)
   dat <- lapply(dat,function(subels)subels[!is.na(subels)])
   dat[sapply(dat,length)==0] <- "NA"
   dat
}
