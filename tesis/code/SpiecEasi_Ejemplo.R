##############################################################
library(devtools)
install_github("zdk123/SpiecEasi")

# Other packages you need to install are
install.packages('igraph')
BiocManager::install("microbiomeutilities")

install.packages("intergraph")
install.packages("GGally")
devtools::install_github("briatte/ggnet")

install.packages("network")
install.packages("ggnetwork")
install.packages("ggpubr")

##############################################################
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

# install.packages("vctrs")
# devtools::install_github("microsud/microbiomeutilities")

ps1 <- readRDS("data/phyobjects/ps.ng.tax")
ps1
head(tax_table(ps1))

ps1.stool <- subset_samples(ps1, bodysite == "Stool")
ps1.stool


ps1.stool.otu <- prune_taxa(taxa_sums(ps1.stool) > 100, ps1.stool)
ps1.stool.otu
# Add taxonomic classification to OTU ID
ps1.stool.otu.f <- microbiomeutilities::format_to_besthit(ps1.stool.otu)

head(tax_table(ps1.stool.otu))

head(tax_table(ps1.stool.otu.f))
  

otu.c <- t(otu_table(ps1.stool.otu.f)@.Data) #extract the otu table from phyloseq object

tax.c <- as.data.frame(tax_table(ps1.stool.otu.f)@.Data)#extract the taxonomy information

head(tax.c)

phyla <- map_levels(colnames(otu.c), from = "best_hit", to = "Phylum", tax_table(ps1.stool.otu.f))
stool.net %v% "Phylum" <- phyla
stool.net %v% "nodesize" <- vsize
# use this only for first attempt to run it on server to save time
saveRDS(otu.c, "data/TabalaAbundancia/stool.otu.c.rds")

saveRDS(tax.c, "data/TabalaAbundancia/stool.tax.c.rds")


###########################################################################

install.packages("torch")
library (torch)
# In practice, use more repetitions
set.seed(1244)

net.c <- spiec.easi(otu.c, method='mb', icov.select.params=list(rep.num=50)) # reps have to increases for real data



saveRDS(net.c, "~input_data/net.c.rds")

#please use more numebr of rep.num (99 or 999) the paraemters 

## Create graph object and get edge values  



