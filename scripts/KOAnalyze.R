#!/usr/bin/env Rscript

library(KEGGREST)

library(pbapply)
#options(echo=T)
options(warn=2)
#pboptions(type="timer")

ko.path = "Desktop/old.tsv"
out.path = "Desktop/old_Weurt_kos_all_modules.tsv"


if(exists("snakemake")){
  ko.path = snakemake@input[[1]]
  out.path = snakemake@output[[1]]
}

ko.data <- read.delim(ko.path, header = F)
colnames(ko.data) <- c("CDS", "KO")



linkAndWait <- function(chunk) {
  result <- keggLink("module", chunk)
  Sys.sleep(0.5)
  return(result)
}
downloadAndWait <- function(chunk) {
  result <- keggGet(chunk)
  Sys.sleep(0.5)
  return(result)
}

downloadSubModules <- function(modules) {
  #print("Checking for modules in module definitions")
  defs = unlist(sapply(modules, function(module)module$DEFINITION))
  submodules = unname(unlist(regmatches(defs,gregexpr("M\\d{5}",defs))))
  modules.names = unname(sapply(modules, function(module)
    module$ENTRY))
  submodules = setdiff(submodules, modules.names)
  if (length(submodules)) {
    submodules <-
      split(submodules,
            ceiling(seq_along(submodules) / 5))
    #print("Downloading additional modules...")
    submodules <-
      pblapply(submodules, downloadAndWait)
    submodules <- unlist(submodules, recursive = F)
    modules <- c(submodules, modules)
    modules <- downloadSubModules(modules)
  }
  modules
}

downloadModules <-
  function(kovector) {
    kovector <- unique(as.character(kovector))
    kovector <-
      split(kovector, ceiling(seq_along(kovector) / 3))
    #print("Linking modules...")
    modules.linked <-
      unique(unlist(unname((
        pbsapply(kovector, linkAndWait)
      ))))
    modules <-
      split(modules.linked, ceiling(seq_along(modules.linked) / 5))
    #print("Downloading modules...")
    modules <- pblapply(modules, downloadAndWait)
    modules <- unlist(modules, recursive = F)
    modules <- downloadSubModules(modules)
    modules
  }


replaceSubmodules <- function(modules.defs){
  defs = unlist(modules.defs)
  submodules = unname(unlist(regmatches(defs,gregexpr("M\\d{5}",defs))))
  toreplace = modules.defs[submodules]
  if(length(toreplace)){
    for (mod in names(toreplace))
      modules.defs <- sapply(modules.defs,function(modules.defs)gsub(mod,replacement = paste0("(",toreplace[[mod]],")"),modules.defs,fixed = T))
    return(modules.defs)
  }
  modules.defs
}


splitIntoBlocks <- function(defstring) {
  characters = strsplit(defstring, "", fixed = T)[[1]]
  openings = ifelse(characters == "(", 1, 0)
  closings = ifelse(characters == ")",-1, 0)
  balance = cumsum(openings + closings) == 0
  blockspaces = which(characters == " " & balance)
  blocks = substring(defstring, c(1, blockspaces + 1), c(blockspaces - 1, nchar(defstring)))
  blocks[blocks != "--"]
}

findNextClosingParenthesis <- function(string, startpos) {
  followingstr = substr(string, startpos, nchar(string))
  characters = strsplit(followingstr, "", fixed = T)[[1]]
  openings = ifelse(characters == "(", 1, 0)
  closings = ifelse(characters == ")",-1, 0)
  balance = cumsum(openings + closings) == 0
  endingpos = startpos + which(balance)[[1]]
  endingpos - 1
}

extractOptionalFromBlock <- function(block) {
  opt = "-K\\d{5}"
  woanyopt = gsub(opt, "", block)
  singleopt = unlist(regmatches(block, gregexpr(opt, block)))
  singleopt = sapply(singleopt, function(opt)
    substr(opt, start = 2, stop = nchar(opt)), USE.NAMES = F)
  groupopts.start =  as.numeric(gregexpr("-(", block, fixed = T))
  if (groupopts.start[[1]] == -1) {
    groupopts = list()
  } else{
    groupopts.start = groupopts.start + 1
    groupopts.stop = sapply(groupopts.start, findNextClosingParenthesis, string =
                              block)
    groupopts = sapply(seq_along(groupopts.start), function(i)
      substr(block, groupopts.start[[i]], groupopts.stop[[i]]))
    for (groupopt in groupopts) {
      woanyopt = sub(paste0("-", groupopt), "", woanyopt, fixed = T)
    }
  }
  opts = unique(c(singleopt, groupopts))
  list(wo = woanyopt, opts = opts)
}

replaceOperators <- function(block) {
  res = gsub(",", "|", block, fixed = T)
  res = gsub(" ", "&", res, fixed = T)
  gsub("+", "&", res, fixed = T)
}


sandbox <- new.env(parent = emptyenv())
for (fun in c("&", "|", "(")) {
  sandbox[[fun]] <- get(fun, "package:base")
}

checkBlock <- function(block, foundkterms) {
  neededkos = unique(unlist(regmatches(block, gregexpr("K\\d{5}", block))))
  neededkos.hits = as.character(neededkos %in% foundkterms)
  names(neededkos.hits) = neededkos
  for (ko in neededkos) {
    block = gsub(ko, neededkos.hits[[ko]], block, fixed = T)
  }
  res = eval(parse(text = block), env = sandbox)
  list(hits = neededkos.hits, complete = res)
}

checkModuleBlocks <- function(blocks, foundkterms) {
  modulehits = lapply(blocks, checkBlock, foundkterms = foundkterms)
  completion = unlist(lapply(modulehits, function(blockhit)
    blockhit$complete))
  terms = unlist(lapply(modulehits, function(blockhit)
    blockhit$hits))
  list(
    modulehits = modulehits,
    completion = completion,
    terms = paste0(unique(names(terms)), collapse = ","),
    hits = paste0(unique(names(
      which(terms == "TRUE")
    )), collapse = ",")
  )
}

checkModule <- function(blocks, optionals, foundkterms) {
  main = checkModuleBlocks(blocks, foundkterms)
  opt = checkModuleBlocks(optionals, foundkterms)
  data.frame(
    allterms = main$terms,
    hitterms = main$hits,
    alloptionals = opt$terms,
    hitoptionals = opt$hits,
    optionalblockcount = length(optionals),
    optionalblockhits = sum(opt$completion),
    blockcount = length(blocks),
    blockhits = sum(main$completion),
    completion = sum(main$completion) / length(blocks)
  )
}

checkModules <-
  function(modules.defs.clean,
           modules.defs.optionals,
           foundkterms) {
    module.checks = lapply(seq_along(modules.defs.clean), function(i)
      checkModule(modules.defs.clean[[i]], modules.defs.optionals[[i]], foundkterms))
    module.checks = cbind(module = names(modules.defs.clean), do.call(rbind, module.checks))
    
  }



# Get all unique modules present with the K terms
modules <- downloadModules(ko.data$KO)

if(length(modules)==0){
  result <- data.frame(module=character(),name=character(),Module.Type=character(),Upper.Class=character(),Lower.Class=character(),allterms=character(),hitterms=character(),alloptionals=character(),hitoptionals=character(),optionalblockcount=integer(),optionalblockhits=integer(),blockcount=integer(),blockhits=integer(),completion=double(),hits=character(),optionals=character()) 
}else{

modules.classes <-
  t(as.data.frame(strsplit(
    sapply(modules, "[[", "CLASS", USE.NAMES = F), "; ", fixed = T
  )))
colnames(modules.classes) <-
  c("Module.Type", "Upper.Class", "Lower.Class")
modules.data = cbind(data.frame(
  module = sapply(modules, function(module)
    module$ENTRY),
  name = sapply(modules, function(module)
    module$NAME)
), modules.classes)
rownames(modules.data) <- NULL

modules.defs = sapply(modules, function(module)
  module$DEFINITION)
names(modules.defs) <- sapply(modules, function(module)
  module$ENTRY)
# Fill submodules
modules.defs = replaceSubmodules(modules.defs)

# Split into blocks (including grouped entries)
modules.defs.blocks = lapply(modules.defs, function(defgroup)
  unlist(lapply(defgroup, splitIntoBlocks)))

# Replace operators
modules.defs.blocks = lapply(modules.defs.blocks, function(blocks)
  unlist(lapply(blocks, replaceOperators)))

# Remove optionals from every block
modules.defs.clean = lapply(modules.defs.blocks, function(blocks)
  lapply(blocks, extractOptionalFromBlock))
modules.defs.optionals = lapply(modules.defs.clean, function(blocks)
  unlist(lapply(blocks, function(block)
    block$opts)))
modules.defs.clean = lapply(modules.defs.clean, function(blocks)
  sapply(blocks, function(block)
    block$wo))



result <- checkModules(modules.defs.clean, modules.defs.optionals, ko.data$KO)
extractHitCds <-
  function(field)
    sapply(strsplit(as.character(result[[field]]), ",", fixed = T), function(kos)
      paste0(as.character(ko.data[ko.data$KO %in% kos, "CDS"]), collapse =
               ","))

result <-  cbind(result,
      hits = extractHitCds("hitterms"),
      optionals = extractHitCds("hitoptionals")
    )

result <- result[result$blockhits != 0,]
result <- merge(modules.data, result, by = "module")
}
write.table(result, out.path, row.names = F,sep = "\t")

