### Cargar libreria ###

library(RCy3)
library(igraph)

cytoscapePing()  # Verifica conexión

# Transicion de nodos con igraph
node_transitivity <- transitivity(W1, type = "local", isolates = "zero")

# Extracion de vertices 
names(W1) <- names(V(W1))

write.table(W1, file = "Output/Cy_nets/Nodes_transitivity/nodeTransitivity_MGYS0006110")

createNetworkFromIgraph(net_graph_1, "net_MGYS00005139_2")

