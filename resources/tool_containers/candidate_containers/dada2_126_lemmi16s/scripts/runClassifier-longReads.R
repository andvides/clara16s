#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
queryReads=args[1]
reference=args[2]
cpus=as.integer(args[3])
path_files=args[4]
aux_parameters=args[5]

library(dada2)
packageVersion("dada2")

#Default tool parameters
parameters <- list(
  maxN=0,
  truncQ=10,
  truncLen=0,
  phix=FALSE,
  compress=FALSE,
  removeBimeraDenovo_method="consensus",
  verbose=TRUE
)

#updating parameters from config file
for (values in unlist(strsplit(aux_parameters,","))) {
  key=unlist(strsplit(values,"="))[[1]]
  value=unlist(strsplit(values,"="))[[2]]
  parameters[[key]] <- value
}
#print(parameters)

path <- path_files
fnFs <- sort(list.files(path, pattern="queryReads.fq", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "."), `[`, 1)


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

names(filtFs) <- sample.names

print("filterAndTrim")
out <- filterAndTrim(fnFs, filtFs, multithread=TRUE, 
                     maxN=as.integer(parameters$maxN), rm.phix=as.logical(parameters$phix), compress=as.logical(parameters$compress),
                     truncLen = as.integer(parameters$truncLen), truncQ=as.integer(parameters$truncQ))

errF <- learnErrors(filtFs, multithread=TRUE)

print("derep")
derep_forward <- derepFastq(filtFs, verbose=as.logical(parameters$verbose))
names(derep_forward) <- sample.names

print("dada")
dada_forward <- dada(derep_forward, err=errF, multithread=TRUE)

seq_table <- makeSequenceTable(dada_forward)
seqtab.nochim <- removeBimeraDenovo(seq_table, method=parameters$removeBimeraDenovo_method, multithread=TRUE, verbose=as.logical(parameters$verbose))

print("assignTaxonomy")
taxa <- assignTaxonomy(seqtab.nochim,reference,taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),multithread=cpus)

reads  <-  t(seqtab.nochim)
results <- merge(reads,taxa,by=0)
results <- results[,-1]
write.csv(results,file="results.txt",row.names = FALSE)

