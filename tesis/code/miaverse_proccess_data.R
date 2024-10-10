############CARGAR LIBRERIAS###################
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

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

BiocManager::install("limma")

############ Manipulación de datos ###################
## Usando las funciones para acceder a la info. taxonomica ## 
## Crea y almacena los datos del tipo MultiAssayExperiment a archivos tipo Phyloseq para su almacenamiento 

## El almacenamiento para cada Phyloseq lo hace por la lista que se genera en el paso anterior
## EStos valores siempre pueden cambiar
mae_MGYS00005537_1 <- mae
mae
mae_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00005537_1[[1]]) 
mae_phyloseq
 
mae2_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00005537_1[[2]])
mae2_phyloseq

mae3_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(mae_MGYS00006110[[3]]) 
mae3_phyloseq

mae_phyloseq_MGYS00006110 <- save(c(mae_phyloseq, mae2_phyloseq, mae3_phyloseq), 
                                  file = "data/mae_phyloseq_proccess/mae_phyloseq_MGYS00006110.RData")

##### Prepocesamiento de los datos para libreria y atmosfera mia ##### 
## Convertir los archivos phyloseq a TreeSummarizedExperiment

mae_TSE_MGYS00005537_1_1 <- makeTreeSEFromPhyloseq(mae_phyloseq)
mae_TSE_MGYS00005537_1_1 

mae_TSE_MGYS00005537_1_2 <- makeTreeSEFromPhyloseq(mae2_phyloseq)
mae_TSE_MGYS00005537_1_2

mae_TSE_MGYS00006110_3 <- makeTreeSEFromPhyloseq(mae3_phyloseq)
mae_TSE_MGYS00006110_3

## Almacenamiento de datos en formato RData
save(mae_TSE_MGYS00005537_1_2, file = "data/mae_TSE/mae_TSE_MGYS00005537_1_2.RData")

################ Data_Wrangling #############################

##### mae_TSE_MGYS00005537_1 ##### 
## Dimension de la tabla 
dim(mae_TSE_MGYS00005537_1_1)

## Unique te permite identificar los datos de la seccion colData > listData
## Recuerda que deben ser datos principalmente de la ubicación o del bioma 
## donde se recolectaron 
unique(mae_TSE_MGYS00005537_1$sample_environment..biome.)

## Crea una tabla que te muestra cuantos analisis se realizaron en ese bioma 
mae_TSE_MGYS00005537_1$sample_environment..biome. %>% table()


## Crear un Subconjunto con respecto a el phylum 
## phylum[1,n] donde n son todos los valores que necesites con respecto a tu base de datos
## Puedes utilziar cualquier valor diferente a phylum que se encuentre en la sección de: 
## rowRanges/listData
cat(paste(unique(rowData(mae_TSE_MGYS00005537_1)$Genus)[1:60], collapse = "\n"))

## Crea una tabla de frecuencia con la libreria kable_extra 
## de los 14 phylum existentes 
rowData(mae_TSE_MGYS00005537_1)$Genus %>% table() %>%
  kable() %>%
  kableExtra::kable_styling("striped", latex_options="scale_down") %>%
  kableExtra::scroll_box(width = "80%")

### Subconjunto por muestra y caracteristica ## 
mae_TSE_sub_MGYS00005537_1 <- mae_TSE_MGYS00005537_1[
  rowData(mae_TSE_MGYS00005537_1)$Phylum %in% c("Actinobacteria", "Firmicutes", "Proteobacteria") &
    !is.na(rowData(mae_TSE_MGYS00005537_1)$Phylum),
  mae_TSE_MGYS00005537_1$sample_environment..biome. %in% c("cheese")]
## Dimensiones del subset
dim(mae_TSE_sub_MGYS00005537_1)


### Limpieza de datos seleccionando por caracteristica 
## Rlimina columnas y filas vacias 
mae_TSE_sub_genus_MGYS00005537_1 <- agglomerateByRank(mae_TSE_MGYS00005537_1, rank = "Genus", na.rm = FALSE)
mae_TSE_sub_genus_MGYS00005537_1
# List bacteria that we want to include
genera <- c(
  "Class:bacilli",
  "Genus:Lactobacillus")
# Subset data
mae_TSE_sub1_genus_MGYS00005537_1 <- mae_TSE_sub_genus_MGYS00005537_1[genera, ]

mae_TSE_sub1_genus_MGYS00005537_1

## Lista de cada sample 
colSums(assay(mae_TSE_MGYS00005537_1, "counts"))[1:91]

## Hace la limpieza de estudios que tengan 0 o na 
mae_TSE_sub_genus_MGYS00005537_1 <- mae_TSE_sub_MGYS00005537_1[
  , colSums(assay(mae_TSE_sub_MGYS00005537_1, "counts")) != 0]
mae_TSE_sub_MGYS00005537_1

## Aqui se muestran sin na 
colSums(assay(mae_TSE_sub_genus_MGYS00005537_1, "counts"))[1:91]



##### Información Taxonomica ##### 

## Comprobar que la información taxonomica sea compatible con el paquete mia 
checkTaxonomy(mae_TSE_MGYS00005537_1)

## Devuelve las columnas que mia puede utilizar 
taxonomyRanks(mae_TSE_MGYS00005537_1)

# Muestra una posible opción para crear un subconjunto 
rowData(mae_TSE_MGYS00005537_1)[, taxonomyRanks(mae_TSE_MGYS00005537_1)]

## Muestra por medio de valores logicos cuantos datos estan vacios 
all(!taxonomyRankEmpty(mae_TSE_MGYS00005537_1, rank = "Kingdom"))

table(taxonomyRankEmpty(mae_TSE_MGYS00005537_1, rank = "Genus"))

table(taxonomyRankEmpty(mae_TSE_MGYS00005537_1, rank = "Species"))

phylum <- !is.na(rowData(mae_TSE_MGYS00005537_1)$Phylum) &
  vapply(data.frame(apply(
    rowData(mae_TSE_MGYS00005537_1)[, taxonomyRanks(mae_TSE_MGYS00005537_1)[3:7]], 1L, is.na)), all, logical(1))

head(getTaxonomyLabels(mae_TSE_MGYS00005537_1[phylum,]))

head(getTaxonomyLabels(mae_TSE_MGYS00005537_1[phylum,], with.rank = TRUE))

head(getTaxonomyLabels(mae_TSE_MGYS00005537_1[phylum,], with.rank = TRUE, make.unique = FALSE))


### Recolectar información sobre determinados taxones ### 

#Lista de taxones unicos paa el rango especificado
head(getUniqueTaxa(mae_TSE_MGYS00005537_1, rank = "Phylum"))

#Identifica todos los taxones que coinciden con "Acidobacteria" existentes en la tabla
mapTaxonomy(mae_TSE_MGYS00005537_1, taxa = "Acidobacteria")


## Arbol de gerarquias 

getHierarchyTree(mae_TSE_MGYS00005537_1)
mae_TSE_MGYS00005537_1 <- addHierarchyTree(mae_TSE_MGYS00005537_1)
mae_TSE_MGYS00005537_1
  

getTaxonomyRanks()

##### mae_TSE_MGYS00006110_2 #####


##### mae_TSE_MGYS00006110_3 #####



###### Visualización de datos ######
rowData(mae_TSE_MGYS00005537_1)

## ATRIBUTOS ## 
sample_variables(mae_TSE_MGYS00005537_1)

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





############ Creación de redes neuronales ############## 

## Cargar librerías ## 

if(!require(SpiecEasi)){
  devtools::install_github("zdk123/SpiecEasi")
}

if(!require(SPRING)){
  devtools::install_github("GraceYoon/SPRING")
}

if(!require(NetCoMi)){
  devtools::install_github(
    "stefpeschel/NetCoMi", force = TRUE, ref = "TSE",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    repos = c("https://cloud.r-project.org/", BiocManager::repositories()))
}


#######################################################

## Subir datos en formato TSE  


dim(mae_TSE_MGYS00005537_1)

checkTaxonomy(mae_TSE_MGYS00005537_1)

taxonomyRanks((mae_TSE_MGYS00005537_1))
# Agregación a nivel de género

mae_TSE_MGYS00005537_1 <- agglomerateByRank(mae_TSE_MGYS00005537_1, rank = "Genus")
mae_TSE_MGYS00005537_1

# Añadir ensayo de abundancia relativa

mae_TSE_MGYS00005537_1 <- transformAssay(
  mae_TSE_MGYS00005537_1,
  assay.type = "counts",
  method = "relabundance",
  MARGIN = "samples")
mae_TSE_MGYS00005537_1

# Filtrado de prevalencia (mantener los géneros con prevalencia > 20%)
mae_TSE_MGYS00005537_1 <- subsetByPrevalentFeatures(
  mae_TSE_MGYS00005537_1,
  prevalence = 0.2,
  detection = 0,
  assay.type = "relabundance")

# Añadir ensayo con abundancias transformadas en log10
mae_TSE_MGYS00005537_1 <- transformAssay(mae_TSE_MGYS00005537_1, method = "log10", pseudocount = 1)
mae_TSE_MGYS00005537_1

# Añadir ensayo con abundancias transformadas de clr
mae_TSE_MGYS00005537_1 <- transformAssay(mae_TSE_MGYS00005537_1, method = "clr", pseudocount = 1)
mae_TSE_MGYS00005537_1

####################### SPRING ############################################
dim(mae_TSE_MGYS00005537_1)

set.seed(13075)

spring_est <- SPRING(
  t(assay(mae_TSE_MGYS00005537_1, "counts")),
  Rmethod = "approx",
  thresh = 0.05,
  lambdaseq = "data-specific")

# Get index of the optimal lambda selected by StARS
opt.K <- spring_est$output$stars$opt.index

# Store partial correlation matrix belonging to the optimal lambda as matrix
spring_cor <- SpiecEasi::symBeta(as.matrix(spring_est$output$est$beta[[opt.K]]))
spring_cor <- as.matrix(spring_cor)
rownames(spring_cor) <- colnames(spring_cor) <- rownames(mae_TSE_MGYS00005537_1)
diag(spring_cor) <- 1

# Arguments:
# - assoMat: association matrix
# - threshold: associations below the threshold are set to zero
# - dissTrans: dissimilarity transformation ("signed" or "unsigned")

transform_asso <- function(assoMat, thresh = NULL, dissTrans = "signed") {
  # Sparsification
  if (!is.null(thresh)) {
    assoMat[abs(assoMat) < thresh] <- 0
  }
  
  # Compute dissimilarity matrix
  if (dissTrans == "signed") {
    dissMat <- sqrt(0.5 * (1 - assoMat))
  } else {
    dissMat <- sqrt(1 - assoMat^2)
  }
  
  # Dissimilarity between nodes with zero correlation is set to 1
  # (these nodes are unconnected and thus should have maximum dissimilarity)
  dissMat[assoMat == 0] <- 1
  
  # Compute similarity matrix
  simMat <- 1 - dissMat
  
  # Turn into igraph object
  graphObj <- SpiecEasi::adj2igraph(simMat)
  
  return(list(graph = graphObj, adja = simMat, asso = assoMat, diss = dissMat))
}

# Create graph object
spring_graph <- transform_asso(spring_cor)$graph
spring_graph

############################## NetCoMi #######################################################

## Creacion de modelo NetCoMi 
## Utilizando el objeto mae_TSE_MGYS00005537_1 sin filtrar y se eliminan
## los taxones que aparecen menos de un 20% en la muestra 
mae_TSE_MGYS00005537_2 <- renameTaxa(mae_TSE_MGYS00005537_2)

netcomi_net <- netConstruct(
  mae_TSE_MGYS00005537_2,
  taxRank = "Genus",
  filtTax = "numbSamp",
  filtTaxPar = list(numbSamp = 0.2),
  measure = "spring",
  measurePar = list(thresh = 0.05, Rmethod = "approx"),
  sparsMethod = "none",
  dissFunc = "signed",
  seed = 13075
)

## Ahora eñ objeto genera una lista de aristas. 
## EStas aristas pueden mostrar asociación, disimilitud y la adyacencia entre cada arista

edgelist <- netcomi_net$edgelist1[
  order(netcomi_net$edgelist1$adja, decreasing = TRUE), ]

head(edgelist)

tail(edgelist)


netcomi_graph <- SpiecEasi::adj2igraph(abs(netcomi_net$adjaMat1))




################ Visualizacion ############################################
# Node sizes

vsize <- (colMeans(t(assay(mae_TSE_MGYS00005537_1, "log10"))) + 1) * 3
vsize
# Fruchterman-Reingold layout from igraph package
set.seed(13075)
lay_fr <- layout_with_fr(spring_graph)
lay_fr
#par(mfrow = c(1,2))
plot(spring_graph, layout = lay_fr, vertex.size = vsize,
     vertex.label = NA, main = "SPRING network")
plot(netcomi_graph, layout = lay_fr, vertex.size = vsize,
     vertex.label = NA, main = "NetCoMi network\n(with SPRING associations)")

################## Medidas de centralidad ########################################

get_centr <- function(graph_obj) {
  # We access igraph directly with "::" because there are more packages loaded
  # in this chapter that contain a degree() function.
  df <- data.frame(Degree = igraph::degree(graph_obj))
  df$Betweenness <- betweenness(graph_obj)
  df$Closeness <- closeness(graph_obj, normalized = TRUE)
  df$Eigenvector <- eigen_centrality(graph_obj)$vector
  return(df)
}

centr_df <- get_centr(spring_graph)
rownames(centr_df) <- rownames(spring_cor)
head(centr_df, 15)
############### Escalar tamaños de nodos por grado ################################ 
get_vsizes <- function(centr_df) {
  df <- as.matrix(centr_df)
  df[, "Betweenness"] <- log(df[, "Betweenness"])
  df[, "Closeness"] <- df[, "Closeness"] * 10
  df[, "Eigenvector"] <- df[, "Eigenvector"] * 10
  df[is.infinite(df) | is.na(df)] <- 0
  return(df)
}

vsize_df <- get_vsizes(centr_df )
head(vsize_df)


par(mfrow = c(2,2))
for (i in seq_along(centr_df)) {
  plot(spring_graph, layout = lay_fr, vertex.size = vsize_df[, i],
       vertex.label = NA, main = colnames(centr_df)[i])
}


### mapas de calor agrupados ##

# Color vector
col <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("royalblue4", "lightblue", "white", "orange", "firebrick3"))

Heatmap(
  adja_sel,
  col = col,
  rect_gp = gpar(col = "gray", lwd = 1),
  show_row_names = FALSE,
  show_column_names = FALSE,
  name = "Association")

## Analisis de red con netcomi

netcomi_netprops <- netAnalyze(
  netcomi_net,
  clustMethod = "cluster_fast_greedy",
  hubPar = "eigenvector",
  normDeg = T)

summary(netcomi_netprops, numbNodes = 5)


plot(netcomi_netprops,
     repulsion = 0.98,
     rmSingles = TRUE,
     shortenLabels = "intelligent",
     labelScale = FALSE,
     nodeSize = "eigenvector",
     nodeSizeSpread = 3,
     hubBorderCol = "gray40",
     nodeColor = "cluster",
     cexNodes = 1.8,
     edgeTranspHigh = 20,
     title1 = "Network properties highlighted",
     showTitle = TRUE,
     cexTitle = 2,
     mar = c(1, 3, 4, 8) 
     )

legend(0.75, 1, cex = 1.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)


library(RColorBrewer)

# Generate vector with phylum names for node coloring
phyla <- as.factor(rowData(mae_TSE_MGYS00005537_1)$Phylum)
names(phyla) <- rowData(mae_TSE_MGYS00005537_1)$Genus

# Create color vector
colvec <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(levels(phyla)))


plot(netcomi_netprops,
     repulsion = 0.98,
     rmSingles = TRUE,
     shortenLabels = "intelligent",
     labelScale = FALSE,
     nodeSize = "mclr",
     nodeColor = "feature",
     featVecCol = phyla,
     colorVec =  colvec,
     nodeTransp = 20,
     highlightHubs = FALSE,
     cexNodes = 0,
     edgeTranspHigh = 20,
     title1 = "Data features highlighted",
     showTitle = TRUE,
     cexTitle = 2.3,
     mar = c(1, 10, 4, 6))

# Add legends
legend(0.75, 1, cex = 0.8, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)

# Colors used in the legend should be equally transparent as in the plot
col_transp <- colToTransp(colvec, 20)

legend(-1.8, 1, cex = 0.8, pt.cex = 2.5, title = "Phylum:",
       legend=levels(phyla), col = col_transp, bty = "n", pch = 16)
