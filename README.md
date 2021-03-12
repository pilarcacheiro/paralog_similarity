## Paralog similarity

Code to retrieve paralogs from Ensembl BioMart using the Bioconductor BiomaRt package
and compute functional similarity between pairs of paralogues based on Gene Ontology
annotations.

It does not require any local input file - all the data can be retrieve from
ftp repositories or trough biomaRt.

### Get_paralog_1.R

It requires to connect to the biomaRt server - sometimes unavailable and time consuming-  
Feedback I'm most interested in:

* functions: `paralog.identity.count.unidir` & `paralog.identity.count.bidir` and 
also `ensembl.paralogs.hgnc.proteincoding.count.unidir` and 
`ensembl.paralogs.hgnc.proteincoding.count.bidir` objects. I realise this piece
of code is messy and I guess it needs to be split / improved somehow.

### GO_Similariy_2.R

It relies on the output from  Get_paralog_1.R. It commputes pairwise GO similarity scores
for different types of GO annotations.  
Feedback I'm most interested in:
* the infinite loop! How to optimise this!!!

Any other feedback, comments more than welcome. Thanks.

