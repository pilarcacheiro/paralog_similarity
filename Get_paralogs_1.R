################################################################################
### Project: Paralogs 
### Script: Get_paralogs_1.R 
### Purpose: Get paralogs from Ensembl95 (biomaRt)
### Notes: 1) Only protein coding genes are considered; 
### 2) mapping paralog subtype to  time of duplication event (1:oldest);
### 3) Bidirectional % aa similarity
################################################################################


# load packages -----------------------------------------------------------


if (!require("dplyr")) install.packages("dplyr")
library("dplyr")

if (!require("plyr")) install.packages("plyr")
library("plyr")

if (!require("tidyr")) install.packages("tidyr")
library("tidyr")

if (!require("string")) install.packages("stringr")
library("stringr")

if (!require("readr")) install.packages("readr")
library("readr")

if (!require("Hmisc")) install.packages("Hmisc")
library("Hmisc")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library("biomaRt")


# input files -------------------------------------------------------------

## hgnc protein coding genes
## filename <- "gene_with_protein_product.txt"
## this piece of code is from previous code review feedback (thanks @iain)

filename <- "gene_with_protein_product.txt"

if (!file.exists(filename)) {
  filename <- paste("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types",
                    filename,
                    sep = "/")
}

protein.coding.genes <- readr::read_delim(filename,
                                          delim = "\t",
                                          col_names = TRUE)




# retrieve paralogs from ensembl ------------------------------------------

## retrieve ensembl paralogs 
## define Mart: homo sapiens
mart.human.ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

## define attributes to get

ensembl.paralogs <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                         "hsapiens_paralog_ensembl_gene",
                                         "hsapiens_paralog_associated_gene_name",
                                         "hsapiens_paralog_chromosome",
                                         "hsapiens_paralog_subtype",
                                         "hsapiens_paralog_orthology_type",
                                         "hsapiens_paralog_perc_id",
                                         "hsapiens_paralog_perc_id_r1"),
                          filters = "ensembl_gene_id",
                          values = protein.coding.genes$ensembl_gene_id,
                          mart = mart.human.ensembl)


## merge with hgnc file to get hgnc id 

ensembl.paralogs.hgnc <-  ensembl.paralogs %>%
  left_join(protein.coding.genes %>% dplyr::select(hgnc_id,ensembl_gene_id),
            by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  rename(hgnc_gene_id = hgnc_id) %>%
  left_join(protein.coding.genes %>% dplyr::select(hgnc_id,ensembl_gene_id),
            by = c("hsapiens_paralog_ensembl_gene" = "ensembl_gene_id")) %>%
  rename(hsapiens_paralog_hgnc_gene_id = hgnc_id) %>%
  dplyr::select(ensembl_gene_id,
                hgnc_gene_id,external_gene_name,
                hsapiens_paralog_ensembl_gene,
                hsapiens_paralog_hgnc_gene_id,
                hsapiens_paralog_associated_gene_name,
                hsapiens_paralog_subtype,
                hsapiens_paralog_orthology_type,
                hsapiens_paralog_perc_id,
                hsapiens_paralog_perc_id_r1)

## add ordinal phylogeny

subtype <- Cs(Opisthokonta,Bilateria,Chordata,Vertebrata,
              Euteleostomi,Sarcopterygii,Tetrapoda,Amniota,Mammalia,
              Theria,Eutheria,Boreoeutheria,Euarchontoglires,Primates,
              Haplorrhini,Simiiformes,Catarrhini,Hominoidea,Hominidae,
              Homininae,Homo.sapiens)


## keep only those paralogues with hgnc id according to hgnc file
## (restrict the associated paralogues to protein coding genes)

ensembl.paralogs.hgnc.proteincoding <- ensembl.paralogs.hgnc %>%
  mutate(hsapiens_paralog_subtype = ifelse(hsapiens_paralog_subtype =="Homo sapiens",
                                           "Homo.sapiens",hsapiens_paralog_subtype)) %>%
  mutate(hsapiens_paralog_subtype_ordinal = plyr::mapvalues(hsapiens_paralog_subtype,subtype,c(1:21))) %>%
  dplyr::select(hgnc_gene_id,hsapiens_paralog_hgnc_gene_id,
                hsapiens_paralog_subtype,hsapiens_paralog_subtype_ordinal,
                hsapiens_paralog_orthology_type:hsapiens_paralog_perc_id_r1) %>%
  dplyr::filter(!is.na(hsapiens_paralog_hgnc_gene_id))


## function to compute number of paralogues for different % of sequence identity and 
## different times of duplication event

paralog.identity.count.unidir <-  function(paralogues,identity) {
  
  count <- paralogues %>% 
    filter(hsapiens_paralog_perc_id >= identity) %>%
    group_by(hgnc_gene_id,hsapiens_paralog_subtype_ordinal) %>% 
    tally() 
  
  return(count)
  
}

paralog.identity.count.bidir <-  function(paralogues,identity) {
  
  count <- paralogues %>% 
    filter(hsapiens_paralog_perc_id >= identity) %>%
    filter(hsapiens_paralog_perc_id_r1  >= identity) %>%
    group_by(hgnc_gene_id,hsapiens_paralog_subtype_ordinal) %>% 
    tally()
  
  return(count)
  
}

# (messy) bunch of code to get to the final dataset

ensembl.paralogs.hgnc.proteincoding.count.unidir <- protein.coding.genes %>% 
  left_join(paralog.identity.count.unidir(ensembl.paralogs.hgnc.proteincoding ,0) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)),
            by=c("hgnc_id"="hgnc_gene_id")) %>%
  dplyr::rename(n.paralogs.all = n) %>%
  left_join(paralog.identity.count.unidir(ensembl.paralogs.hgnc.proteincoding ,10) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.10 = n) %>%
  left_join(paralog.identity.count.unidir(ensembl.paralogs.hgnc.proteincoding ,20) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.20 = n) %>%
  left_join(paralog.identity.count.unidir(ensembl.paralogs.hgnc.proteincoding ,30) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.30 = n) %>%
  left_join(paralog.identity.count.unidir(ensembl.paralogs.hgnc.proteincoding ,40) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.40 = n) %>%
  left_join(paralog.identity.count.unidir(ensembl.paralogs.hgnc.proteincoding ,50) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.50 = n) %>%
  dplyr::select(-gene_paralog_subtype) %>%
  as.data.frame() %>%
  replace(is.na(.),0)



ensembl.paralogs.hgnc.proteincoding.count.bidir <- protein.coding.genes %>%
  dplyr::select(hgnc_id) %>%
  distinct() %>%
  left_join(paralog.identity.count.bidir(ensembl.paralogs.hgnc.proteincoding ,0) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)),
            by=c("hgnc_id"="hgnc_gene_id")) %>%
  dplyr::rename(n.paralogs.all = n) %>%
  left_join(paralog.identity.count.bidir(ensembl.paralogs.hgnc.proteincoding ,10) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.10 = n) %>%
  left_join(paralog.identity.count.bidir(ensembl.paralogs.hgnc.proteincoding ,20) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.20 = n) %>%
  left_join(paralog.identity.count.bidir(ensembl.paralogs.hgnc.proteincoding ,30) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.30 = n) %>%
  left_join(paralog.identity.count.bidir(ensembl.paralogs.hgnc.proteincoding ,40) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.40 = n) %>%
  left_join(paralog.identity.count.bidir(ensembl.paralogs.hgnc.proteincoding ,50) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",
                                                   hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::rename(n.paralogs.50 = n) %>%
  dplyr::select(-gene_paralog_subtype) %>%
  as.data.frame() %>%
  replace(is.na(.),0)


# export files ------------------------------------------------------------

write.table(ensembl.paralogs.hgnc.proteincoding,
            gzfile("./data_output/paralogs_ensembl_proteincoding.txt.gz"),
            quote = F, sep = "\t", row.names = F)

write.table(ensembl.paralogs.hgnc.proteincoding.count.unidir,
            gzfile("./data_output/paralogs_count_ensembl_proteincoding.unidirectional.txt.gz"),
            quote = F, sep = "\t", row.names = F)

write.table(ensembl.paralogs.hgnc.proteincoding.count.bidir,
            gzfile("./data_output/paralogs_count_ensembl_proteincoding.bidirectional.txt.gz"),
            quote = F, sep = "\t", row.names = F)




