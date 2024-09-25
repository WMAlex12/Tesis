############CARGAR LIBRERIAS###################
library(mia)
library(dplyr)
library(knitr)
library(readr)
library(phyloseq)
############ Manipulación de datos ###################
## Usando las funciones para acceder a la info. taxonomica ## 


library(mia)

mae_MGYS00001153

mae_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[1]]) 
mae_phyloseq
 
mae2_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[2]])
mae2_phyloseq

mae3_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[3]]) 
mae3_phyloseq

mae_phyloseq_MGYS00006110 <- save(c(mae_phyloseq, mae2_phyloseq, mae3_phyloseq), 
                                  file = "data/mae_phyloseq_proccess/mae_phyloseq_MGYS00006110.RData")


mae_TSE_MGYS00006110 <- makePhyloseqFromTreeSE(mae_phyloseq, tree_name = "phylo")

###### Visualización de datos ######
rank_names(mae_phyloseq)

## ATRIBUTOS ## 
sample_variables(mae_phyloseq)

checkTaxonomy(mae_phyloseq)
##  [1] TRUE
all(!taxonomyRankEmpty(mae_MGYS00006110, rank = "Kingdom"))
##  [1] TRUE
table(taxonomyRankEmpty(mae_phyloseq, rank = "Genus"))
##  
##  FALSE  TRUE 
##   8008 11208
table(taxonomyRankEmpty(mae_phyloseq, rank = "Species"))
##  
##  FALSE  TRUE 
##   1413 17803


## 
data("GlobalPatterns", package ="mia")
data 
class(GlobalPatterns)
class(mae_phyloseq)
