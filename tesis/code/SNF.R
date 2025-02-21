
########### Reorganización de datos para trabajar con Redes de similitud ###############

tax_table(physeq_1)[, colnames(tax_table(physeq_1))] <- gsub(tax_table(physeq_1)[, colnames(tax_table(physeq_1))], pattern = "[a-z]__", replacement = "")

tax_table(physeq_1)[tax_table(physeq_1)[, "Class"] == "", "Class"] <- "Unidentified"

# save the pseq object

save(physeq_1, file = paste("data/phyobjects/taxtable_phy_", study,".RData", sep = "") )

# check for features of data  
summarize_phyloseq(physeq_1)

###############Preparando  datos para SNF##########################
### Convert a TSE en un phyloseq 
Phyloseq_Abundance_1 <- makePhyloseqFromTreeSummarizedExperiment(TSE_Abundance_1)

Phyloseq_Abundance_2 <- makePhyloseqFromTreeSummarizedExperiment(TSE_Abundance_2)

Phyloseq_Abundance_3 <- makePhyloseqFromTreeSummarizedExperiment(TSE_Abundance_3)

### Convertir la tabla taxonomica a un dataframe 
tax_1 <- as.data.frame(tax_table(Phyloseq_Abundance_1)@.Data)

### Convertir la tabla OTU a un dataframe 
OTU_1 <- t(otu_table(Phyloseq_Abundance_1)@.Data)

dim(OTU_1)
# In practice, use more repetitions
set.seed(1204)

########################## SNF ####################################

#' 
#' Calculamos la distancia entre los datos
#' 

dist1 <- as.matrix(dist(t(OTU_1)))
dist2 <- as.matrix(dist(OTU_2))

#' 
#' Primero se hacen las matrices de distancia 
#' Los valores de K y sigma se colocan empiricamente 
#'
  
W1 <- affinityMatrix(t(dist1), K = 6, sigma = 0.75)
W2 <- affinityMatrix(dist2, K = 25, sigma = 0.75)

#'
#' Generamos un  heatmap
#'


displayClustersWithHeatmap(W1, spectralClustering(W1, K = 2))
displayClustersWithHeatmap(W2, spectralClustering(W2, K = 2))


