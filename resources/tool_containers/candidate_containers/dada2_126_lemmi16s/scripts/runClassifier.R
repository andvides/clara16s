#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
queryReads=args[1]
reference=args[2]
cpus=as.integer(args[3])
path_files=args[4]
aux_parameters=args[5]

library(dada2); 
packageVersion("dada2")

#Default tool parameters
parameters <- list(
  maxN=0,
  truncQ=10,
  trunc_len_forward=200,
  trunc_len_reverse=200,
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
fnFs <- sort(list.files(path, pattern="queryReads1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="queryReads2.fq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "."), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, truncQ=20, rm.phix=FALSE, compress=FALSE, multithread=TRUE)
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, truncLen = c(200, 200), rm.phix=FALSE, compress=FALSE, multithread=TRUE)
print("filterAndTrim")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, multithread=TRUE, 
                     maxN=as.integer(parameters$maxN), rm.phix=as.logical(parameters$phix), compress=as.logical(parameters$compress),
                     truncLen = c(as.integer(parameters$trunc_len_forward),as.integer(parameters$trunc_len_reverse)), truncQ=as.integer(parameters$truncQ))
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

print("derep")
derep_forward <- derepFastq(filtFs, verbose=as.logical(parameters$verbose))
derep_reverse <- derepFastq(filtRs, verbose=as.logical(parameters$verbose))
names(derep_forward) <- sample.names
names(derep_reverse) <- sample.names

print("dada")
dada_forward <- dada(derep_forward, err=errF, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=errR, multithread=TRUE)

print("merge")
merged_reads <- mergePairs(dada_forward, derep_forward, dada_reverse,derep_reverse, verbose=as.logical(parameters$verbose))
seq_table <- makeSequenceTable(merged_reads)
seqtab.nochim <- removeBimeraDenovo(seq_table, method=parameters$removeBimeraDenovo_method, multithread=TRUE, verbose=as.logical(parameters$verbose))

print("assignTaxonomy")
taxa <- assignTaxonomy(seqtab.nochim,reference,taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),multithread=cpus)
#write.table(taxa,file='results.txt',row.names = FALSE, col.names = FALSE, sep=';')
reads  <-  t(seqtab.nochim)
results <- merge(reads,taxa,by=0)
results <- results[,-1]
write.csv(results,file="results.txt",row.names = FALSE)
