##############################################################

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
###################################################### 
library(mia)
library(dplyr)
library(knitr)
library(kableExtra)
library(readr)
library(phyloseq)
library(scater)
library(SPRING)
library(NetCoMi)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(RDS)
library(metagMisc)
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
library(devtools)

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

study <- readline(prompt = "Introduce el estudio: ")
PATH <- paste("~/Documents/WMAlex/Tesis/tesis/data/MgnifyDownloads/mae_getresult/mae_", 
              study,".rds", sep =  "")

assign("physeq",readRDS(PATH) )


## Convertir los datos MAE a un phyloseq en caso de ser necesario 

physeq_1 <- makePhyloseqFromTreeSummarizedExperiment(physeq[[1]]) 
physeq_1

physeq_2 <- makePhyloseqFromTreeSummarizedExperiment(physeq[[2]])
physeq_2

physeq_3 <- makePhyloseqFromTreeSummarizedExperiment(physeq[[3]]) 
physeq_3

## Guarda los archivos en la carpeta data y especificamente a la de 
save(physeq_1, physeq_2, physeq_3, file = paste("data/mae_phyloseq_proccess/phyloseq_", study,".RData", sep = ""))

##### Prepocesamiento de los datos para libreria y atmosfera mia ##### 
## Convertir los archivos phyloseq a TreeSummarizedExperiment

TSE1 <- assign(paste("TSE_", study, "_1", sep = ""), makeTreeSEFromPhyloseq(physeq_1))

TSE2 <- assign(paste("TSE_", study, "_2", sep = ""), makeTreeSEFromPhyloseq(physeq_2))

TSE3 <- assign(paste("TSE_", study, "_3", sep = ""), makeTreeSEFromPhyloseq(physeq_2))


save(TSE1 ,TSE2 ,TSE3 , file = paste("data/TSE/TSE_", study,".RData", sep = "") )

################ Data_Wrangling #############################
## Comprobar si la información es utilizable 
checkTaxonomy(TSE1)

checkTaxonomy(TSE2)

checkTaxonomy(TSE3)

######################################################
####### Creación de tablas de abundancias #########
#' 
#' 
#' 

TSE_Abundance_1 <- transformAssay(TSE1, method = "relabundance")
TSE_Abundance_1 <- agglomerateByPrevalence(
  TSE_Abundance_1,
  rank = "Class",
  assay.type = "relabundance")

TSE_Abundance_2 <- transformAssay(TSE2, method = "relabundance")
TSE_Abundance_2 <- agglomerateByPrevalence(
  TSE_Abundance_2,
  rank = "Class",
  assay.type = "relabundance")

TSE_Abundance_3 <- transformAssay(TSE3, method = "relabundance")
TSE_Abundance_3 <- agglomerateByPrevalence(
  TSE_Abundance_3,
  rank = "Class",
  assay.type = "relabundance")

save(TSE_Abundance_1, TSE_Abundance_2, TSE_Abundance_3, file = paste("data/TablaAbundancia/TSE_Abundance_", study,".RData", sep = "") )

head(TSE_Abundance_1)

########### Reorganización de datos para trabajar con SpiecEasi ###############

tax_table(physeq_1)[, colnames(tax_table(physeq_1))] <- gsub(tax_table(physeq_1)[, colnames(tax_table(physeq_1))], pattern = "[a-z]__", replacement = "")

tax_table(physeq_1)[tax_table(physeq_1)[, "Class"] == "", "Class"] <- "Unidentified"

# save the pseq object

save(physeq_1, file = paste("data/phyobjects/taxtable_phy_", study,".RData", sep = "") )

# check for features of data  
summarize_phyloseq(physeq_1)

###############Preparando  datos para spiec Easi##########################
### Convert a TSE en un phyloseq 
Phyloseq_Abundance_1 <- makePhyloseqFromTreeSE(TSE1)

Phyloseq_Abundance_2 <- makePhyloseqFromTreeSE(TSE2)

Phyloseq_Abundance_3 <- makePhyloseqFromTreeSE(TSE3)

### Convertir la tabla taxonomica a un dataframe 
tax_1 <- as.data.frame(tax_table(Phyloseq_Abundance_1)@.Data)

### Convertir la tabla OTU a un dataframe 
OTU_1 <- t(otu_table(Phyloseq_Abundance_1)@.Data)

dim(OTU_1)
# In practice, use more repetitions
set.seed(1204)

netSE_1 <- spiec.easi(OTU_1, method='mb', icov.select.params=list(rep.num=999)) # reps have to increases for real data


save(netSE_1, file = paste("data/network_correlations/netSE_", study, ".RData" ,sep = ""))
######################################################################################
#please use more numebr of rep.num (99 or 999) the paraemters 

## Create graph object and get edge values  
class(netSE_1)

n.c <- symBeta(getOptBeta(netSE_1))

colnames(n.c) <- rownames(n.c) <- colnames(t(OTU_1))

vsize <- log2(apply(OTU_1, 2, mean)) # add log abundance as properties of vertex/nodes.
vsize


net_graph_1 <- graph.adjacency(n.c, mode='undirected', add.rownames = TRUE, weighted = TRUE)
net_graph_1 # we can see all the attributes and weights


## Pone un 0.1 a todos los valores menores a 0 
E(net_graph_1)$weight[E(net_graph_1)$weight <= 0] <- 0.1

coords.fdr = layout.fruchterman.reingold(net_graph_1)

E(net_graph_1)[weight > 0]$color<-"#321321AA" #now color the edges based on their values positive is steelblue
E(net_graph_1)[weight < 0]$color<-"red"  #now color the edges based on their values

plot(net_graph_1, layout=coords.fdr, vertex.size = 3, vertex.label.cex = 0.1)

#
net_graph_1 <- asNetwork(net_graph_1)
network::set.edge.attribute(net_graph_1, "color", ifelse(net_graph_1 %e% "weight" > 0, "steelblue", "orange"))

colnames(tax_table(Phyloseq_Abundance_1))

summary(net_graph_1)

phyla <- tax_table(Phyloseq_Abundance_1)[, "Family"]
phyla <- as.character(phyla)

net_graph_1 %v% "Family" <- phyla
net_graph_1 %v% "nodesize" <- vsize

mycolors <- scale_color_manual(values = rainbow(168))  # Usando una paleta rainbow para 62 colores

# Alternativamente, si prefieres colores con buen contraste, puedes usar colores en el espacio HCL (hue, chroma, luminance)
mycolors <- scale_color_manual(values = rainbow_hcl(168))

p <- ggnet2(net_graph_1, node.color = "Family", 
            label = F, node.size = "nodesize", 
            label.size = 0.5, edge.color = "color") + 
  guides(color = guide_legend(title = "Family"), size = "none") + mycolors

p

View(net_graph_1)
# If net_graph_1 is an adjacency matrix
net_graph_1 <- graph_from_adjacency_matrix(net_graph_1)

# If net_graph_1 is an edge list (data frame or matrix)
net_graph_1 <- graph_from_edgelist(net_graph_1)
head(net_graph_1)
fam_mb <- degree.distribution(net_graph_1)
plot(0:(length(fam_mb)-1), fam_mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
## Plot mejorado 
library(ggraph)
library(tidygraph)

set_graph_style(plot_margin = margin(1,1,1,2))
graph <- as_tbl_graph(net_graph_1)

# Not specifying the layout - defaults to "auto"
ggraph(graph, layout = coords.fdr, circular = T) + 
  geom_edge_link(aes(colour = "sample")) + 
  geom_node_point()


