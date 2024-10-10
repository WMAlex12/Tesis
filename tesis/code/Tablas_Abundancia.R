
Tabla_MGYS00006110_1 <- as.data.frame(otu_table(mae_phyloseq))
Tabla_MGYS00006110_1


Tabla_MGYS00006110_2 <- as.data.frame(otu_table(mae2_phyloseq))
Tabla_MGYS00006110_2

Tabla_MGYS00006110_3 <- as.data.frame(otu_table(mae3_phyloseq))
Tabla_MGYS00006110_3


seqs <- colnames(Tabla_MGYS00006110_1)
colnames(Tabla_MGYS00006110_1) <- paste0('ASV_',seq(1:ncol(Tabla_MGYS00006110_1)))
nams <- colnames(Tabla_MGYS00006110_1)
