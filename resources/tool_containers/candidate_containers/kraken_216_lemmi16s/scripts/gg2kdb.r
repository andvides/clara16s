#!/usr/bin/env Rscript
rm(list=ls()) ; gc() 

      
##   get sequences   =========================================================
      
      # attributes, names, length, hist etc. for DNA work
      library(optparse)        
      library("Biostrings")
      library("parallel")

      option_list <- list(
      make_option("--path_in",
              type = "character",
              default = "dataset",
              help = "Path to input directory [default = %default]"),
  
      make_option("--n_cores",
              type = "integer",
              default = 4,
              help = "Number of CPU cores to use [default = %default]")
      )

      opt <- parse_args(OptionParser(option_list = option_list))

      path_in  <- opt$path_in
      n_cores  <- opt$n_cores

      cat("INPUT DIRECTORY:", path_in, "\n")
      cat("USING CORES:", n_cores, "\n\n")
      
      path_out <- "tmp"
      length(seqs <- readDNAStringSet( filepath = paste0(path_in, "/reference.fasta" )))  #   23,467,470
      head(tax <- read.table(paste0( path_in, "/reference.tsv"), header = TRUE, sep = "\t"))   # 23M lines


  ##   make taxonomy   =========================================================
      
      rank_is <- c( "d__" = "domain", "p__" = "phylum", 
                    "c__" = "class", "o__" = "order", 
                    "f__" = "family", "g__" = "genus",
                    "s__" = "species")      
      
      # 3 taxonomic dummy levels added here for fun
      tax_dummies <- 3  
      taxs <- as.data.frame(do.call("rbind", mclapply( tax$Taxon, function(aa){  
        bb <- unlist(strsplit( aa, split = "\\; "))[1:(7+tax_dummies)]
        bb
      }, mc.cores = n_cores, mc.cleanup = TRUE)), stringsAsFactors = FALSE)
      colnames(taxs) <- rank_is
      taxs$header <- tax$Feature.ID
      rownames(taxs) <- taxs$header
    
    ## crucial - accession - taxon map
      taxs$effective_taxon <- unlist( mclapply( 1:nrow(taxs), function(aa){  
        bb <- taxs[ aa, 1:(7+tax_dummies) ]
        bb[ max(which( !( grepl("__$", bb, perl = TRUE) | is.na(bb) ) ))  ]
      },  mc.cores = n_cores, mc.cleanup = TRUE))

      uniq_tax <- apply( taxs[, 1:(7+tax_dummies)], 2, function(aa){
        bb <- aa[ !( is.na(aa) ) ]
        sort(unique(unlist(bb)))
        })


  ## subset to seqs in 1.4-1.65kbp range    ====================================
      
      s_widths <- width(seqs)
      # 332,588 in full version
      length(seqs_w <- seqs[ s_widths > 1400 & s_widths < 1650 ])
      # match with those in taxs -  312,905 in full version
      length( usable_seqs <- seqs_w[ names(seqs_w) %in% taxs$header ] )
      
      
  ## taxid proto   -----------------------------------------------------------

    ## root is 1, so floop/domains start from 2...
      floop <- 1    
    ## counter with universal assignment to var
      taxid_proto <- lapply( names(uniq_tax)[ 
        sapply( uniq_tax, function(aa){ !identical( aa, character(0)) })
        ], function(aa){   
          if( aa == "domain"){ floop <<- 1 }
          bb <- uniq_tax[ aa ][[1]]
          cc <- (floop+1):(floop+length(bb))    #  stringr::str_pad( (floop+1):(floop+length(bb)), width = 7, pad = "0" )
          floop <<- floop + length(cc)
          list( "name_txt" = bb, "taxid" = cc)
      })
      names(taxid_proto) <- rank_is

    ## no mult lookups - decide rank ahead of time, use rank_is & grep through that alone (e.g. based on prefix)
      identify <- function(aa){
        bb_rank <- rank_is[ gsub("__.*", "__", aa )]                           # aa <- "g__GWA2-37-10"
        if( aa %in% taxid_proto[[bb_rank]][[1]] ){
          taxid_proto[[bb_rank]][[2]][ match(aa, taxid_proto[[bb_rank]][[1]] ) ]
        }else{
          print(paste0("taxon ", aa, " wasn't in the taxid map"))
        }
      }

      identify_parent <- function(aa){   # aa <- "g__UBA6532"         aa <- "s__"
        if( grepl("d__", aa)){
          "1"                      # snoonos get unos
        }else{
          bb_rank <- rank_is[ gsub("__.*", "__", aa) ]
          cc <- taxs[ match(aa, taxs[ , bb_rank ]) , ( which( names(taxs) == bb_rank) -1) ] 
          identify( cc )
          }   
        }


  ## names proto   -----------------------------------------------------------
  
      names_proto <- data.frame(
        "tax_id" = unlist( mclapply( unlist( uniq_tax), identify, mc.cores = n_cores, mc.cleanup = TRUE )),
        "name_txt" = unlist( uniq_tax),
        "unique name" = "",                   # safety net ignored.
        "name class" = "scientific name",     # synonym, common name, ... lower case!
        stringsAsFactors = FALSE
      )

      
  ## proto nodes   -----------------------------------------------------------
      
      nodes_proto <- data.frame(
        "tax_id" = unlist(mclapply( names_proto$name_txt, identify, mc.cores = n_cores, mc.cleanup = TRUE )),                                           ## parallel
        "parent tax_id" = unlist( mclapply( names_proto$name_txt, identify_parent, mc.cores = n_cores, mc.cleanup = TRUE )),                             ## parallel
        "rank" = rank_is[ unlist( mclapply( names_proto$name_txt, function(aa){ gsub("(^\\w__).*", "\\1", aa) }, mc.cores = n_cores, mc.cleanup = TRUE )) ] ,     ## parallel
        "embl code" = "",                  # no one seems to know about this one
        "division id" = 0,                  # 0 = Bacteria (but also apparently Archaea...)
        "inherited div flag" = 1,                            
        "genetic code id" = 11,                   # 11 = Bact/Arch/plastid
        "inherited GC" = 1,                     # 1 if node inherits genetic code from parent
        "mitochondrial genetic code" = 0,       # see gencode.dmp file
        "inherited MGC" = 1,                # 1 if node inherits mitochondrial gencode from parent
        "GenBank hidden" = 0,               # 1 if name is suppressed in GenBank entry lineage
        "hidden subtree root" = 0,          # 1 if this subtree has no sequence data yet
        "comments" = "")

            
    ## rooting the tree - must be appended to :
      # - nodes.dmp (note 8,0,1, instead of 11,1,0)
      nodes_proto <- rbind( 
        c( "1","1","no rank","","8","0","1","0","0","0","0","0",""),
        nodes_proto
      )
      # - and names.dmp
      names_proto <- rbind( 
        c( "1", "all", "all", "synonym"),
        c( "1", "root", "root", "scientific name"),
        names_proto
      )

      
  ## modify the seq headers    =================================================
      
      names(usable_seqs) <- paste0(
        names(usable_seqs), 
        "|kraken:taxid|", 
        unlist( mclapply( taxs[ names(usable_seqs) , "effective_taxon"], identify, mc.cores = n_cores, mc.cleanup = TRUE )
                )
        )
            

##   out that!   ---------------------------------------------------------------

      # head( nodes_proto )
      # head( names_proto )
      # head( usable_seqs )

      ## note sep for nodes & names
      writeXStringSet( usable_seqs, filepath = paste0( path_out, "/gg2_2024.09__FullLength313k__seqs.fna"), 
                       compress = FALSE, format = "fasta" )
      
      write.table(names_proto, file = paste0( path_out, "/gg2_2024.09__FullLength313k__names.dmp"), 
                  sep = '\t|\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
      write.table(nodes_proto, file = paste0( path_out, "/gg2_2024.09__FullLength313k__nodes.dmp"), 
                  sep = '\t|\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
