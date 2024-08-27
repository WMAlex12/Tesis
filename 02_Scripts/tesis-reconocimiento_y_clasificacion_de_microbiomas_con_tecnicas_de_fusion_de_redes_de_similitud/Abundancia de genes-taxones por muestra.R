library(MultiAssayExperiment)

# Mostrar las clases de datos en el objeto MAE
assayNames <- names(assays(mae))
print(assayNames)

# Si el objeto es similar a un SummarizedExperiment
counts_data <- assays(mae)[[1]]  # Usa el índice apropiado si hay varios assays
head(counts_data)

# Muestra los nombres de los ensayos (e.g., "counts", "logcounts")

# Extraer datos de un assay específico, por ejemplo, "counts"
counts_data <- assays(mae)[["microbiota"]]

# Ver las primeras filas del assay
head(counts_data)
# Extraer los datos en formato largo (long format) para ggplot
library(tidyr)
# Verificar si el experimento "microbiota" es un SummarizedExperiment
if (inherits(microbiota_experiment, "SummarizedExperiment")) {
  # Extraer los datos de conteo
  counts_data <- assays(microbiota_experiment)[[1]]  # Usa el primer assay
  
  # Extraer los nombres de los genes o taxones
  gene_names <- rownames(microbiota_experiment)
  
  # Convertir los datos de conteo en un formato largo
  counts_long <- as.data.frame(counts_data)
  counts_long <- pivot_longer(counts_long, cols = everything(), names_to = "Sample", values_to = "Count")
  
  # Agregar los nombres de los genes/taxones
  counts_long$Gene <- rep(gene_names, each = ncol(counts_data))
  
  # Mostrar la tabla organizada
  head(counts_long)
}

library(ggplot2)

# Visualizar la abundancia de genes/taxones por muestra
ggplot(counts_long, aes(x = Sample, y = Count, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Abundancia de genes/taxones por muestra", x = "Muestra", y = "Abundancia") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 4),  # Tamaño del texto de las etiquetas del eje X
    axis.text.y = element_text(size = 10),  # Tamaño del texto de las etiquetas del eje Y
    axis.title = element_text(size = 12),   # Tamaño del texto de los títulos de los ejes
    plot.title = element_text(size = 14),   # Tamaño del texto del título de la gráfica
    legend.text = element_text(size = 8),   # Tamaño del texto de la leyenda
    legend.title = element_text(size = 10)  # Tamaño del título de la leyenda
  )

