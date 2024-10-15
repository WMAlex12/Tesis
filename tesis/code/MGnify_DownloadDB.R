
###Descarga e instalacion de librerias##
#install.packages(c("devtools","mia","plyr","dplyr","reshape2","httr","urltools","tidyverse"))
  ## Instalación de Biomanagerpara descargar librerias externas al CRAN

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

## Versión actual funcional para r versión 4.4.1

#BiocManager::install(version = "3.19")

## Descarga de MGnifyR con Biomanager

#BiocManager::install("mia") 
#BiocManager::install("phyloseq")
BiocManager::install("MGnifyR")

# Librerias necesarias para utilizar MGnifyR

library(devtools) 
library(tidyverse) 
library(mia) 
library(plyr)
library(dplyr) 
library(reshape2) 
library(httr) 
library(urltools)
library(remotes) 
library(MGnifyR) 
library(phyloseq) 
############################Solucion error utf-8###############################

# En caso de generar un Error: no encoding suplied default utf-8
# Realizar los siguientes cambios
## Recolectar la ubicacion del archivo .Rprofile 
file_path <- file.path(Sys.getenv("HOME"), ".Rprofile")
file_path

# Esta linea de codigo te va a abrir el archivo .Rprofile y te permitira editarlo
#file.edit(file.path(file_path))

# Solo se debe de poner la siguiente linea,guardar los cambios y reiniciar RStudio 
# options(encoding = "utf-8")

###############################################################################


study <- "MGYS00000747"
study

# Set up the MGnify client instance

mgclnt <- MgnifyClient(usecache = TRUE)

# Retrieve the list of analyses associated with a study

accession_list <- searchAnalysis(mgclnt, "studies", study)

# Download all associated study/sample and analysis metadata

meta_dataframe <- getMetadata(mgclnt, accession_list)

savePathMetadata <- paste("data/MgnifyDownloads/meta_data/meta_data_",study, ".RData", sep = "")
savePathMetadata

save(meta_dataframe, file = savePathMetadata)

# Convert analyses outputs to a single `MultiAssayExperiment` object

mae <- getResult(mgclnt, meta_dataframe$analysis_accession, usecache = TRUE)
savePathMae <- paste("data/MgnifyDownloads/mae_getresult/mae_",study, ".RData", sep = "")
savePathMae
save(mae, file = savePathMae)

 


