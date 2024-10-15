###################################################### 
library(mia)
library(dplyr)
library(knitr)
library(kableExtra)
library(readr)
library(phyloseq)
library(scater)
library(SpiecEasi)
library(SPRING)
library(NetCoMi)
library(igraph)
library(ComplexHeatmap)
library(circlize)

################Cargar Paquetes#######################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

BiocManager::install("limma")

############ Manipulación de datos ###################
## Usando las funciones para acceder a la info. taxonomica ## 
## Crea y almacena los datos del tipo MultiAssayExperiment a archivos tipo Phyloseq para su almacenamiento 
## El almacenamiento para cada Phyloseq lo hace por la lista que se genera en el paso anterior
## EStos valores siempre pueden cambiar

mae_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00005537[[1]]) 
mae_phyloseq

mae2_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00005537[[2]])
mae2_phyloseq

mae3_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[3]]) 
mae3_phyloseq

mae_phyloseq_MGYS0000 <- save(c(mae_phyloseq, mae2_phyloseq), 
                                  file = "data/mae_phyloseq_proccess/mae_phyloseq_MGYS00006110.RData")

##### Prepocesamiento de los datos para libreria y atmosfera mia ##### 
## Convertir los archivos phyloseq a TreeSummarizedExperiment

mae_TSE_MGYS00005537_1_1 <- makeTreeSEFromPhyloseq(mae_phyloseq)
mae_TSE_MGYS00005537_1_1 

mae_TSE_MGYS00005537_1_2 <- makeTreeSEFromPhyloseq(mae2_phyloseq)
mae_TSE_MGYS00005537_1_2

mae_TSE_MGYS00006110_3 <- makeTreeSEFromPhyloseq(mae3_phyloseq)
mae_TSE_MGYS00006110_3
