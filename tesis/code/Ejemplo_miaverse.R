# Store data into se and check dimensions
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
# Show dimensions (features x samples)
dim(tse)
##  [1] 19216    26
unique(tse$SampleType)
##  [1] Soil               Feces              Skin              
##  [4] Tongue             Freshwater         Freshwater (creek)
##  [7] Ocean              Sediment (estuary) Mock              
##  9 Levels: Feces Freshwater Freshwater (creek) Mock ... Tongue
# Show the frequency of each value
tse$SampleType %>% table()
# Subset by sample
tse_sub <- tse[ , tse$SampleType %in% c("Feces", "Skin", "Tongue")]

# Show dimensions
dim(tse_sub)
##  [1] 19216     9
# Subset by feature
tse_sub <- tse[
  rowData(tse)$Phylum %in% c("Actinobacteria","Chlamydiae"), ]

# Show dimensions
dim(tse_sub)
##  [1] 1652   26
# Inspect possible values for phylum
# and show the first five values
tse_phylum <- cat(paste(unique(rowData(tse)$Phylum)[1:5], collapse = "\n"))
##  Crenarchaeota
##  Euryarchaeota
##  Actinobacteria
##  Spirochaetes
##  MVP-15
# Show the frequency of each value
rowData(tse)$Phylum %>% table() %>% head() %>%
  kable() %>%
  kableExtra::kable_styling("striped", latex_options="scale_down") %>%
  kableExtra::scroll_box(width = "100%")
# Subset by sample and feature and remove NAs

tse_sub <- tse[
  rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") &
    !is.na(rowData(tse)$Phylum),
  tse$SampleType %in% c("Feces", "Skin", "Tongue")]
# Show dimensions
dim(tse_sub)
##  [1] 1652    9
# Agglomerate data at Genus level
tse_genus <- agglomerateByRank(tse, rank = "Genus", na.rm = FALSE)
# List bacteria that we want to include
genera <- c(
  "Class:Thermoprotei",
  "Genus:Sulfolobus",
  "Genus:Sediminicola")
# Subset data
tse_sub <- tse_genus[genera, ]

tse_sub
##  class: TreeSummarizedExperiment 
##  dim: 3 26 
##  metadata(1): agglomerated_by_rank
##  assays(1): counts
##  rownames(3): Class:Thermoprotei Genus:Sulfolobus Genus:Sediminicola
##  rowData names(7): Kingdom Phylum ... Genus Species
##  colnames(26): CL3 CC1 ... Even2 Even3
##  colData names(7): X.SampleID Primer ... SampleType Description
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  rowLinks: a LinkDataFrame (3 rows)
##  rowTree: 1 phylo tree(s) (19216 leaves)
##  colLinks: NULL
##  colTree: NULL
# List total counts of each sample
colSums(assay(tse_genus_sub, "counts"))[1:5]
##      CL3     CC1     SV1 M31Fcsw M11Fcsw 
##        1       0       0       1       1

# Remove samples that do not contain any bacteria
tse_genus_sub <- tse_sub[
  , colSums(assay(tse_sub, "counts")) != 0]
tse_sub
##  class: TreeSummarizedExperiment 
##  dim: 3 26 
##  metadata(1): agglomerated_by_rank
##  assays(1): counts
##  rownames(3): Class:Thermoprotei Genus:Sulfolobus Genus:Sediminicola
##  rowData names(7): Kingdom Phylum ... Genus Species
##  colnames(26): CL3 CC1 ... Even2 Even3
##  colData names(7): X.SampleID Primer ... SampleType Description
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  rowLinks: a LinkDataFrame (3 rows)
##  rowTree: 1 phylo tree(s) (19216 leaves)
##  colLinks: NULL
##  colTree: NULL
