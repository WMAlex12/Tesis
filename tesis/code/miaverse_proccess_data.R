############CARGAR LIBRERIAS###################
library(mia)
library(dplyr)
library(knitr)
library(kableExtra)
library(readr)
library(phyloseq)
############ Manipulación de datos ###################
## Usando las funciones para acceder a la info. taxonomica ## 
## Crea y almacena los datos del tipo MultiAssayExperiment a archivos tipo Phyloseq para su almacenamiento 

## El almacenamiento para cada Phyloseq lo hace por la lista que se genera en el paso anterior
## EStos valores siempre pueden cambiar

mae_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[1]]) 
mae_phyloseq
 
mae2_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[2]])
mae2_phyloseq

mae3_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[3]]) 
mae3_phyloseq

mae_phyloseq_MGYS00006110 <- save(c(mae_phyloseq, mae2_phyloseq, mae3_phyloseq), 
                                  file = "data/mae_phyloseq_proccess/mae_phyloseq_MGYS00006110.RData")

##### Prepocesamiento de los datos para libreria y atmosfera mia ##### 
## Convertir los archivos phyloseq a TreeSummarizedExperiment

mae_TSE_MGYS00006110_1 <- makeTreeSEFromPhyloseq(mae_phyloseq)
mae_TSE_MGYS00006110_1 

mae_TSE_MGYS00006110_2 <- makeTreeSEFromPhyloseq(mae2_phyloseq)
mae_TSE_MGYS00006110_2

mae_TSE_MGYS00006110_3 <- makeTreeSEFromPhyloseq(mae3_phyloseq)
mae_TSE_MGYS00006110_3
################ Data_Wrangling #############################

##### mae_TSE_MGYS00006110_1 ##### 
## Dimension de la tabla 
dim(mae_TSE_MGYS00006110_1)

## Unique te permite identificar los datos de la seccion colData > listData
## Recuerda que deben ser datos principalmente de la ubicación o del bioma 
## donde se recolectaron 
unique(mae_TSE_MGYS00006110_1$sample_environment..biome.)

## Crea una tabla que te muestra cuantos analisis se realizaron en ese bioma 
mae_TSE_MGYS00006110_1$sample_environment..biome. %>% table()


## Crear un Subconjunto con respecto a el phylum 
## phylum[1,n] donde n son todos los valores que necesites con respecto a tu base de datos
## Puedes utilziar cualquier valor diferente a phylum que se encuentre en la sección de: 
## rowRanges/listData
cat(paste(unique(rowData(mae_TSE_MGYS00006110_1)$Genus)[1:30], collapse = "\n"))

## Crea una tabla de frecuencia con la libreria kable_extra 
## de los 14 phylum existentes 
rowData(mae_TSE_MGYS00006110_1)$Genus %>% table() %>%
  kable() %>%
  kableExtra::kable_styling("striped", latex_options="scale_down") %>%
  kableExtra::scroll_box(width = "80%")

### Subconjunto por muestra y caracteristica ## 
mae_TSE_sub_MGYS00006110_1 <- mae_TSE_MGYS00006110_1[
  rowData(mae_TSE_MGYS00006110_1)$Phylum %in% c("Actinobacteria", "Firmicutes", "Proteobacteria") &
    !is.na(rowData(mae_TSE_MGYS00006110_1)$Phylum),
  mae_TSE_MGYS00006110_1$sample_environment..biome. %in% c("cheese")]
## Dimensiones del subset
dim(mae_TSE_sub_MGYS00006110_1)


### Limpieza de datos seleccionando por caracteristica 
## Rlimina columnas y filas vacias 
mae_TSE_sub_genus_MGYS00006110_1 <- agglomerateByRank(mae_TSE_MGYS00006110_1, rank = "Genus", na.rm = FALSE)
mae_TSE_sub_genus_MGYS00006110_1
# List bacteria that we want to include
genera <- c(
  "Class:bacilli",
  "Genus:Lactobacillus")
# Subset data
mae_TSE_sub1_genus_MGYS00006110_1 <- mae_TSE_sub_genus_MGYS00006110_1[genera, ]

mae_TSE_sub1_genus_MGYS00006110_1

## Lista de cada sample 
colSums(assay(mae_TSE_MGYS00006110_1, "counts"))[1:91]

## Hace la limpieza de estudios que tengan 0 o na 
mae_TSE_sub_genus_MGYS00006110_1 <- mae_TSE_sub_MGYS00006110_1[
  , colSums(assay(mae_TSE_sub_MGYS00006110_1, "counts")) != 0]
mae_TSE_sub_MGYS00006110_1

## Aqui se muestran sin na 
colSums(assay(mae_TSE_sub_genus_MGYS00006110_1, "counts"))[1:91]


##### mae_TSE_MGYS00006110_2 #####


##### mae_TSE_MGYS00006110_3 #####



###### Visualización de datos ######
rank_names(mae_TSE_MGYS00006110)

## ATRIBUTOS ## 
sample_variables(mae_TSE_MGYS00006110)

checkTaxonomy(mae_TSE_MGYS00006110)
  ##  [1] TRUE
all(!taxonomyRankEmpty(mae_TSE_MGYS00006110, rank = "Kingdom"))
##  [1] TRUE
table(taxonomyRankEmpty(mae_TSE_MGYS00006110, rank = "Genus"))
##  
##  FALSE  TRUE 
##   8008 11208
table(taxonomyRankEmpty(mae_TSE_MGYS00006110, rank = "Species"))
##  
##  FALSE  TRUE 
##   1413 17803


## 
mae_TSE_MGYS00006110 <- makeTreeSEFromPhyloseq(mae_phyloseq)



class(GlobalPatterns)
class(mae_TSE_MGYS00006110)
head(mae_TSE_MGYS00006110)






