
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

mae_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00005139[[1]]) 
mae_phyloseq

mae2_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00005139[[2]])
mae2_phyloseq

mae3_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00005139[[3]]) 
mae3_phyloseq

save(mae_phyloseq, mae2_phyloseq, mae3_phyloseq, file = "data/mae_phyloseq_proccess/mae_phyloseq_MGYS00005139.RData")

rm(mae_phyloseq, mae2_phyloseq)
##### Prepocesamiento de los datos para libreria y atmosfera mia ##### 
## Convertir los archivos phyloseq a TreeSummarizedExperiment

TSE_MGYS00005537_1 <- makeTreeSEFromPhyloseq(mae_phyloseq)
TSE_MGYS00005537_1 

TSE_MGYS00005537_2 <- makeTreeSEFromPhyloseq(mae2_phyloseq)
TSE_MGYS00005537_2

TSE_MGYS00005139_3 <- makeTreeSEFromPhyloseq(mae3_phyloseq)
TSE_MGYS00005139_3

save(TSE_MGYS00005139_1,TSE_MGYS00005537_2, file = "data/mae_TSE/TSE_MGYS00005537.RData")

################ Data_Wrangling #############################
## Comprobar si la información es utilizable 
checkTaxonomy(TSE_MGYS00005139_1)

## Devolver las columnas que contienen  la información taxonómica
taxonomyRanks(TSE_MGYS00005139_1)

## Esto nos permite identificar si s posible realizar subconjuntos de los datos
rowData(TSE_MGYS00005139_1)[, taxonomyRanks(TSE_MGYS00005139_1)]

## Mostrar si hay Valores NA en la columna y contabilizar los NA
## En este caso se hace para las columnas de Family y Genus
all(!taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Family"))

table(taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Family"))

all(!taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Genus"))

table(taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Genus"))

## Crear sub conjuntos por caracateristica 
# Selección de una flia completa de datos con sus respectivo identificador 
# 
FamilySubset <- rowData(TSE_MGYS00005139_1)$Family |> table() |> head() |> kable()
FamilySubset










######################################################
#### Tabla de abundancias sin paqueteria mia ####
counts_data <- assay(TSE_MGYS00005537_2)

abundancia_resumen <- rowData(TSE_MGYS00005537_2)$Family %>%
  as.data.frame() %>%
  bind_cols(as.data.frame(counts_data)) %>%
  mutate(total_abundancia = rowSums(across(starts_with("MGYA")), na.rm = TRUE)) %>%
  ungroup()
abundancia_resumen


########### Aglomeración ###############
tse_family <- agglomerateByRank(TSE_MGYS00005139_1, rank = "Family", update.tree = TRUE)
tse_family

rowData(tse_family)

rowTree(tse_family)


TSE_MGYS00005139_1_prev_family <- agglomerateByPrevalence(
  TSE_MGYS00005139_1,
  rank = "Family",
  prevalence = 20 / 100,
  detection =  1
)
TSE_MGYS00005139_1_prev_family


# Transform "counts" assay to relative abundances ("relabundance"), with
# pseudocount 1
TSE_MGYS00005139_1 <- transformAssay(
  TSE_MGYS00005139_1, assay.type = "counts", method = "relabundance", pseudocount = 1)
TSE_MGYS00005139_1
# Transform relative abundance assay ("relabundance") to "clr", using
# pseudocount if necessary; name the resulting assay to "clr"
TSE_MGYS00005139_1 <- transformAssay(
  x = TSE_MGYS00005139_1, assay.type = "relabundance", method = "clr", pseudocount = TRUE,
  name = "clr")


assay(TSE_MGYS00005139_1, "relabundance") |> head()

TSE_MGYS00005139_1 <- transformAssay(TSE_MGYS00005139_1, method = "pa")
assay(TSE_MGYS00005139_1, "pa") |> head()



