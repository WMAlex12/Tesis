
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
library(RDS)
library(metagMisc)
################Cargar Paquetes#######################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

BiocManager::install("limma")
remotes::install_github("vmikk/metagMisc")
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
rowData(TSE_MGYS00005139_1)[, taxonomyRanks(TSE_MGYS00005139_1)]$Phylum

## Mostrar si hay Valores NA en la columna y contabilizar los NA
## En este caso se hace para las columnas de Family y Genus
all(!taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Phylum"))

table(taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Phylum"))

all(!taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Genus"))

table(taxonomyRankEmpty(TSE_MGYS00005139_1, rank = "Genus"))

## Crear sub conjuntos por caracateristica 
# Selección de una flia completa de datos con sus respectivo identificador 
# 
FamilySubset <- rowData(TSE_MGYS00005139_1)$Family |> table() |> head() |> kable()
FamilySubset


convertToPhyloseq
phyloseq_to_df()





######################################################
#### Tabla de abundancias por phylum####
###
counts_data <- assay(TSE_MGYS00005537_2)

abundancia_resumen <- rowData(TSE_MGYS00005537_2)$Family %>%
  as.data.frame() %>%
  bind_cols(as.data.frame(counts_data)) %>%
  mutate(total_abundancia = rowSums(across(starts_with("MGYA")), na.rm = TRUE)) %>%
  ungroup()
abundancia_resumen
####### ESta tabla vamos a usar #########

TSE_MGYS00005139_1 <- transformAssay(TSE_MGYS00005139_1, method = "relabundance")
TSE_MGYS00005139_1 <- agglomerateByPrevalence(
  TSE_MGYS00005139_1,
  rank = "Phylum",
  assay.type = "relabundance")
saveRDS(TSE_MGYS00005139_1, "data/TabalaAbundancia/TSE_Abundancia_MGYS00005139_1.rds")

TSE_MGYS00005139_2 <- transformAssay(TSE_MGYS00005139_2, method = "relabundance")
TSE_MGYS00005139_2 <- agglomerateByPrevalence(
  TSE_MGYS00005139_2,
  rank = "Phylum",
  assay.type = "relabundance")
saveRDS(TSE_MGYS00005139_2, "data/TabalaAbundancia/TSE_Abundancia_MGYS00005139_2.rds")

TSE_MGYS00005139_3 <- transformAssay(TSE_MGYS00005139_3, method = "relabundance")
TSE_MGYS00005139_3 <- agglomerateByPrevalence(
  TSE_MGYS00005139_3,
  rank = "Phylum",
  assay.type = "relabundance")
saveRDS(TSE_MGYS00005139_3, "data/TabalaAbundancia/TSE_Abundancia_MGYS00005139_3.rds")


TSE_MGYS00005139_1

########### Reorganización de datos para trabajar con SpiecEasi ###############
##############################################################
library(devtools)
install_github("zdk123/SpiecEasi")

# Other packages you need to install are
install.packages('igraph')
BiocManager::install("microbiome")

install.packages("intergraph")
install.packages("GGally")
devtools::install_github("briatte/ggnet")

install.packages("network")
install.packages("ggnetwork")
install.packages("ggpubr")

##############################################################
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(SpiecEasi) # Network analysis for sparse compositional data  
library(network)
library(intergraph)
#devtools::install_github("briatte/ggnet")
library(ggnet)
library(igraph)


mae_phyloseq

tax_table(mae_phyloseq)[, colnames(tax_table(mae_phyloseq))] <- gsub(tax_table(mae_phyloseq)[, colnames(tax_table(mae_phyloseq))], pattern = "[a-z]__", replacement = "")

tax_table(mae_phyloseq)[tax_table(mae_phyloseq)[, "Phylum"] == "", "Phylum"] <- "Unidentified"

# save the pseq object

saveRDS(mae_phyloseq, "./data/mae_phyloseq_proccess/mae_phyloseq_MGYS00005139.rds")

# check for features of data  
summarize_phyloseq(mae_phyloseq)
summarize_phyloseq(ps.otu)




###############Preparando  datos para spiec Easi##########################
### Convert a TSE en un phyloseq 
Phyloseq_MGYS00005139_1 <- makePhyloseqFromTreeSE(TSE_Abundancia_MGYS00005139_1)

### Convertir la tabla taxonomica a un dataframe 
tax_MGYS00005139_1 <- as.data.frame(tax_table(Phyloseq_MGYS00005139_1)@.Data)


### Convertir la tabla OTU a un dataframe 
OTU_MGYS00005139_1 <- t(otu_table(Phyloseq_MGYS00005139_1)@.Data)

# In practice, use more repetitions
set.seed(1204)

net.c_1 <- spiec.easi(OTU_MGYS00005139_1, method='mb', icov.select.params=list(rep.num=1)) # reps have to increases for real data

saveRDS(net.c_1, "Output/net.c_1.rds")

#please use more numebr of rep.num (99 or 999) the paraemters 

## Create graph object and get edge values  