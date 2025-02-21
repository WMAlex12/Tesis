
############ Manipulación de datos ###################
## Usando las funciones para acceder a la info. taxonomica ## 
## Crea y almacena los datos del tipo MultiAssayExperiment a archivos tipo Phyloseq para su almacenamiento 
## El almacenamiento para cada Phyloseq lo hace por la lista que se genera en el paso anterior
## EStos valores siempre pueden cambiar

study <- readline(prompt = "Introduce el estudio: ")
PATH <- paste("~/Documentos/Tesis/tesis/data/MgnifyDownloads/mae_getresult/mae_", 
              study,".rds", sep =  "")

assign("physeq",readRDS(PATH))
physeq <- mae
## Convertir los datos MAE a un phyloseq en caso de ser necesario 

physeq_1 <- convertToPhyloseq(physeq[[1]]) 
physeq_1

physeq_2 <- convertToPhyloseq(physeq[[2]])
physeq_2

physeq_3 <- convertToPhyloseq(physeq[[3]]) 
physeq_3

## Guarda los archivos en la carpeta data y especificamente a la de 
save(physeq_1, physeq_2, physeq_3, file = paste("data/mae_phyloseq_proccess/phyloseq_", study,".RData", sep = ""))

##### Prepocesamiento de los datos para libreria y atmosfera mia ##### 
## Convertir los archivos phyloseq a TreeSummarizedExperiment

TSE1 <- assign(paste("TSE_", study, "_1", sep = ""), convertFromPhyloseq(physeq_1))

TSE2 <- assign(paste("TSE_", study, "_2", sep = ""), convertFromPhyloseq(physeq_2))

TSE3 <- assign(paste("TSE_", study, "_3", sep = ""), convertFromPhyloseq(physeq_2))


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

