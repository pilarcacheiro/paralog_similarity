################################################################################
### Project: Paralogs 
### Script: GO_Similarity_2.R 
### Purpose: Compute similarity between paralogs based on GO annotations
################################################################################


# load packages -----------------------------------------------------------
library(dplyr);library(tidyr);library(stringr);library(readr)
library(GOSemSim)

if (!require("dplyr")) install.packages("dplyr")
library("dplyr")

if (!require("tidyr")) install.packages("tidyr")
library("tidyr")

if (!require("string")) install.packages("stringr")
library("stringr")

if (!require("readr")) install.packages("readr")
library("readr")

if (!require("GOSemSim")) install.packages("GOSemSim")
library("GOSemSim")


# input files -------------------------------------------------------------

filename <- "gene_with_protein_product.txt"

if (!file.exists(filename)) {
  filename <- paste("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types",
                    filename,
                    sep = "/")
}

protein.coding.genes <- readr::read_delim(filename,
                                          delim = "\t",
                                          col_names = TRUE)



paralogs <- read_delim("./data_output/paralogs_ensembl_proteincoding.txt.gz",delim="\t") %>%
  inner_join(protein.coding.genes,by=c("hgnc_gene_id" = "hgnc_id")) %>%
  dplyr::rename(entrez_gene = entrez_id) %>%
  inner_join(protein.coding.genes,by=c("hsapiens_paralog_hgnc_gene_id" = "hgnc_id")) %>%
  dplyr::rename(entrez_paralog = entrez_id)

paralog.pairs <- paralogs %>%
  dplyr::select(entrez_gene,entrez_paralog) %>%
  distinct()



# import GO annotations ---------------------------------------------------

## Biological process, molecular function and cellular component

hsbp <- godata('org.Hs.eg.db', ont="BP")
hsmf <- godata('org.Hs.eg.db', ont="MF")
hscc <- godata('org.Hs.eg.db', ont="CC")
hsbp.go <- hsbp@geneAnno

## compute pairwise similarity scores

pair.go.similarity.df <- list()

#for(i in 1:dim(paralog.pairs)[1]){

## this takes a long time to run
## only the 10 first rows

for(i in 1:10){
  
  gene = paralog.pairs[i,1]
  paralog = paralog.pairs[i,2]
  
  sim.bp.reisnik = mgeneSim(genes = c(gene,paralog), semData =hsbp, 
                            measure = "Resnik", drop = "IEA", combine = "BMA",verbose = TRUE) 
  sim.mf.reisnik = mgeneSim(genes = c(gene,paralog), semData =hsmf,
                            measure = "Resnik", drop = "IEA", combine = "BMA",verbose = TRUE) 
  sim.cc.reisnik = mgeneSim(genes = c(gene,paralog), semData =hscc, 
                            measure = "Resnik", drop = "IEA", combine = "BMA",verbose = TRUE) 
  
  if(length(sim.bp.reisnik)>1) {sim.bp.reisnik.score = sim.bp.reisnik[1,2] 
  } else {
    sim.bp.reisnik.score = NA
  }
  
  if(length(sim.mf.reisnik)>1) {sim.mf.reisnik.score = sim.mf.reisnik[1,2] 
  } else {
    sim.mf.reisnik.score = NA
  }
  
  if(length(sim.cc.reisnik)>1) {sim.cc.reisnik.score = sim.cc.reisnik[1,2] 
  } else {
    sim.cc.reisnik.score = NA
  }
  
  
  pair.go.similarity.df[[i]] = 
    data.frame(gene,paralog,sim.bp.reisnik.score,sim.mf.reisnik.score,sim.cc.reisnik.score)
  
  cat("iteration:", i, " \n") 
  
  flush.console()
  
}

#}

## merge the three similarity scores

pair.go.similarity.df.all <- do.call(rbind,pair.go.similarity.df)