
########### Reorganización de datos para trabajar con SpiecEasi ###############

tax_table(physeq_1)[, colnames(tax_table(physeq_1))] <- gsub(tax_table(physeq_1)[, colnames(tax_table(physeq_1))], pattern = "[a-z]__", replacement = "")

tax_table(physeq_1)[tax_table(physeq_1)[, "Class"] == "", "Class"] <- "Unidentified"

# save the pseq object

save(physeq_1, file = paste("data/phyobjects/taxtable_phy_", study,".RData", sep = "") )

# check for features of data  
summarize_phyloseq(physeq_1)

###############Preparando  datos para spiec Easi##########################
### Convert a TSE en un phyloseq 
Phyloseq_Abundance_1 <- makePhyloseqFromTreeSummarizedExperiment(TSE_Abundance_1)

Phyloseq_Abundance_2 <- makePhyloseqFromTreeSummarizedExperiment(TSE_Abundance_2)

Phyloseq_Abundance_3 <- makePhyloseqFromTreeSummarizedExperiment(TSE_Abundance_3)

### Convertir la tabla taxonomica a un dataframe 
tax_2 <- as.data.frame(tax_table(Phyloseq_Abundance_2)@.Data)

### Convertir la tabla OTU a un dataframe 
OTU_2 <- t(otu_table(Phyloseq_Abundance_2)@.Data)

dim(OTU_2)
# In practice, use more repetitions
set.seed(1204)

netSE_1 <- spiec.easi(OTU_1, method='mb', icov.select.params=list(rep.num=999)) # reps have to increases for real data


save(netSE_1, file = paste("data/network_correlations/netSE_", study, ".RData" ,sep = ""))
######################################################################################
#please use more numebr of rep.num (99 or 999) the paraemters 

## Create graph object and get edge values  
class(netSE_2)

n.c <- symBeta(getOptBeta(netSE_2))

colnames(n.c) <- rownames(n.c) <- colnames(t(OTU_2))

vsize <- log2(apply(OTU_2, 2, mean)) # add log abundance as properties of vertex/nodes.
vsize


net_graph_2 <- graph.adjacency(n.c, mode='undirected', add.rownames = TRUE, weighted = TRUE)
net_graph_2 # we can see all the attributes and weights


## Pone un 0.1 a todos los valores menores a 0 
E(net_graph_2)$weight[E(net_graph_2)$weight <= 0] <- 0.1

coords.fdr = layout.fruchterman.reingold(net_graph_2)

E(net_graph_2)[weight > 0]$color<-"#321321AA" #now color the edges based on their values positive is steelblue
E(net_graph_2)[weight < 0]$color<-"red"  #now color the edges based on their values

plot(net_graph_2, layout=coords.fdr, vertex.size = 3, vertex.label.cex = 0.1)

#
net_graph_2 <- asNetwork(net_graph_2)
network::set.edge.attribute(net_graph_2, "color", ifelse(net_graph_2 %e% "weight" > 0, "steelblue", "orange"))

colnames(tax_table(Phyloseq_Abundance_2))

summary(net_graph_2)

phyla <- tax_table(Phyloseq_Abundance_2)[, "Family"]
phyla <- as.character(phyla)

net_graph_2 %v% "Family" <- phyla
net_graph_2 %v% "nodesize" <- vsize

mycolors <- scale_color_manual(values = rainbow(168))  # Usando una paleta rainbow para 62 colores

# Alternativamente, si prefieres colores con buen contraste, puedes usar colores en el espacio HCL (hue, chroma, luminance)
mycolors <- scale_color_manual(values = rainbow_hcl(168))

p <- ggnet2(net_graph_2, node.color = "Family", 
            label = F, node.size = "nodesize", 
            label.size = 0.5, edge.color = "color") + 
  guides(color = guide_legend(title = "Family"), size = "none") + mycolors

p

View(net_graph_2)
# If net_graph_2 is an adjacency matrix
net_graph_2 <- graph_from_adjacency_matrix(net_graph_2)

# If net_graph_2 is an edge list (data frame or matrix)
net_graph_2 <- graph_from_edgelist(net_graph_2)
head(net_graph_2)
fam_mb <- degree.distribution(net_graph_2)
plot(0:(length(fam_mb)-1), fam_mb, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
## Plot mejorado 
library(ggraph)
library(tidygraph)

set_graph_style(plot_margin = margin(1,1,1,2))
graph <- as_tbl_graph(net_graph_2)

# Not specifying the layout - defaults to "auto"
ggraph(graph, layout = coords.fdr, circular = T) + 
  geom_edge_link(aes(colour = "sample")) + 
  geom_node_point()

