##############################################################

install_github("zdk123/SpiecEasi")
devtools::install_github("vmikk/metagMisc")
# Other packages you need to install are
install.packages('igraph')
BiocManager::install("metagMisc")

install.packages("intergraph")
install.packages("GGally")
devtools::install_github("GraceYoon/SPRING")

install.packages("network")
install.packages("ggnetwork")
install.packages("ggpubr")


devtools::install_github("stefpeschel/NetCoMi", 
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
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
BiocManager::install("SPRING")

BiocManager::install("limma")
remotes::install_github("vmikk/metagMisc")