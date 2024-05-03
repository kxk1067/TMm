# TMm
DADA2 and Phyloseq analysis
install.packages("tidyverse")
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
install.packages("vegan")
install.packages("dendextend")
install.packages("viridis")
##install.packages("pairwiseAdonis")

install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("ggpubr")
##install.packages("ANCOMBC")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ANCOMBC", force = TRUE)
##install.packages("microbiome")
library(BiocManager)
BiocManager::install("microbiome", force = TRUE)
install.packages("eulerr")
##install.packages("phyloseq")
devtools::install_github("joey711/phyloseq")

install.packages("venn")
install.packages("VennDiagram")
library(VennDiagram)
install.packages("venn(list_core)")
install.packages("meta")
##install.packages("meta(phy)")
# Load necessary library
library(phyloseq)

# Assuming 'phy' is your phyloseq object

# Access metadata
metadata <- sample_data(phy)

# Print metadata to see its structure and contents
print(metadata)
##install.packages("core_members")
# Install BiocManager package if you haven't already
install.packages("BiocManager")

# Load BiocManager package
library(BiocManager)

# Install core_members from Bioconductor
BiocManager::install("core_members")

library(grid)
library(meta)
library(tidyverse) ; packageVersion("tidyverse") # 1.3.2
library(phyloseq) ; packageVersion("phyloseq") # 1.42.0
library(vegan) ; packageVersion("vegan") # 2.6.4
library(dendextend) ; packageVersion("dendextend") # 1.16.0
library(viridis) ; packageVersion("viridis") # 0.6.2
library("pairwiseAdonis"); packageVersion("pairwiseAdonis") # 0.4
library(ggpubr)
library(ANCOMBC)
library(microbiome)
library(eulerr)
library(phyloseq)
library(venn)

setwd ("/Users/kirtanakrishnakumar/Desktop/TMM Project/TMM Project/")

#uploadi/Users/kirtanakrishnakumar/Desktop/TMM Project/TMM Project/#uploading your taxonomy and count data
asv_tab = read.delim("ASVs_counts.tsv", sep = "\t")
asv_tax = read.delim("ASVs_taxonomy.tsv", sep="\t")

#need to construct around your data using your metadata table from Run Selector

sample_info_tab = read.csv("CleanedPRJNA51101_KrishnaKumar.csv")

#need to make columns the rownames all of them
#need to figure out what the row names are for all your data, likely if you followed the tutorial the tax and tab are "X")
#for your sample data file the first row needs to be the SRR that is listed in your ASV files
#load the dplyr package
library("dplyr")
asv_tab = asv_tab %>% tibble::column_to_rownames("X")
asv_tax = asv_tax %>% tibble::column_to_rownames("X")
sample_df = sample_info_tab %>% tibble::column_to_rownames("Run")



#transforming your taxonomy data into a matrix for entry into phyloseq
asv_mat = as.matrix(asv_tab)
tax_mat = as.matrix(asv_tax)

#transform into phyloseq objects
asv = otu_table(asv_mat, taxa_are_rows=TRUE)
tax=tax_table(tax_mat)
samples=sample_data(sample_df)


#make a phyloseq object
phy =phyloseq(asv, tax, samples)

#remove chloroplast and mitochondria from your phyloseq object
phy <- phy %>% subset_taxa( family!= "mitochondria" | is.na(family) & class!="Chloroplast" | is.na(class) ) 

#number of taxa
ntaxa(phy)

#number of samples

nsamples(phy)

#number of variables
sample_variables(phy)

#load required packages
install.packages("ggplot2")
library("ggplot2")
library("phyloseq")
#plotting richness, be sure to change x and color by what you're trying to evaluate
plot_richness(phy, x="Host",color = "env_medium", measure=c("Chao1", "Shannon","Fisher")) +
  geom_boxplot(alpha=0.2)
plot_richness 

#with/without Fisher
#plotting richness, be sure to change x and color by what you're trying to evaluate
plot_richness(phy, x="lat_lon",color = "Host", measure=c("Chao1", "Shannon", "Fisher")) +
  geom_boxplot(alpha=0.2)
plot_richness 

#if you want to show what is significant and then what you want to compare
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
comparisons_material = list( c("Thalassia testudinum", "Syringodium filiforme"), c("Syringodium filiforme", "not applicable"), c("Thalassia testudinum", "not applicable"))

install.packages("ggpubr")
library("ggpubr")
#then rerun richness
plot_richness(phy, x="Host", measure=c("Chao1", "Shannon")) +
  geom_boxplot(alpha=0.2)+
  sta2t_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)


plot_richness(phy, x="Host",color = "env_medium", measure=c("Chao1", "Shannon")) +
  geom_boxplot(alpha=0.2)+
  facet_grid(~env_medium)+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)

#make a stacked barplot of the data at the phylum level
phy_phylum = tax_glom(phy, taxrank = "phylum",NArm=FALSE)
rel_abund_phylum = phyloseq::transform_sample_counts(phy_phylum,
                                                           function(x){x / sum(x)})
phyloseq::plot_bar(rel_abund_phylum, fill = "phylum")+
  geom_bar(aes(color = phylum, fill = phylum), stat = "identity", position = "stack", color = "black") +
  labs(x = "env_medium", y = "Relative Abundance\n") +
  facet_wrap(~env_medium, scales = "free") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))

#find the core members of the microbiome
#convert to relative abundance
pseq_rel = microbiome::transform(phy, "compositional")

#make a variable for all of the conditions you want to compare
stuff = unique(as.character(meta(phy)$Species))
print(stuff)
# Load necessary libraries
library(phyloseq)

# Assuming 'phy' is your phyloseq object

# Access metadata
metadata <- sample_data(phy)
print(metadata)
# Extract unique values of the "Species" column
#stuff <- unique(as.character(metadata$Species))
env_medium <- metadata$env_medium
print(paste0("Stuff Variable:",env_medium))

core_members <- function(physeq_obj, detection = 0.001, prevalence = 0.1) {
  # Transform to relative abundance
  physeq_rel <- transform_sample_counts(physeq_obj, function(x) x / sum(x))
  
  # Identify core taxa
  core_taxa <- taxa_names(physeq_rel)[apply(otu_table(physeq_rel), 1, function(taxon_abundance) {
    # Calculate the proportion of samples where the taxon is detected above the threshold
    detected_in_samples <- mean(taxon_abundance > detection)
    # Return TRUE if taxon prevalence meets the threshold
    return(detected_in_samples >= prevalence)
  })]
  
  return(core_taxa)
}

#make a for loop to go through each bit of "stuff" one by one and combine identified core taxa into a list
list_core <- c() # an empty object to store information
for (n in env_medium){ # for each variable n in stuff
  print(paste0("Identifying Core Taxa for ", n))
   print(n)
  ps.sub <- subset_samples(pseq_rel, env_medium == n) # Choose sample from DiseaseState by n
  print(paste("Number of samples in subset:",nsamples(ps.sub)))
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.10) #prevelence really matters--how much of the sample can it make upt 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
print(list_core)
}
install.packages("devtools")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram", force = TRUE)
print(list_core[["rhizosphere"]])

# Collect these into a list
sets_list <- list(
  Rhizosphere = list_core[["rhizosphere"]],
  Phyllosphere = list_core[["phyllosphere"]],
  Rhizome = list_core[["rhizome"]],
  MarineWater = list_core[["marine water"]],
  RhizosphereSediment = list_core[["rhizosphere sediment"]],
  CoastalMarineSediment = list_core[["coastal marine sediment"]]
)
print(sets_list)
#plot(venn(list_core))

#venn.plot <- venn.diagram(x=sets_list[c("coastal marine sediment","rhizosphere sediment","marine water")], category.names = c("coastal marine sediment","rhizosphere sediment","marine water"),filename= NULL,output = TRUE)
#print(venn.plot)
#$grid.draw(venn.plot)
mycols = c("Rhizosphere"="#d6e2e9", "Phyllosphere"="#cbf3f0", "Rhizome"="#fcf5c7","MarineWater"="#FFA07A","RhizosphereSediment"="#D4EFDF","CoastalMarineSediment"="#F4D03F") 
# Assuming you are only comparing three sets for simplicity

venn.plot <- venn.diagram(
  x = sets_list[c("Rhizosphere", "Phyllosphere", "Rhizome","Marinewater","RhizosphereSediment","CoastaMarineSediment")],
  category.names = c("Rhizosphere", "Phyllosphere", "Rhizome","Marinewater","RhizosphereSediment","CoastalMarineSediment"),
  fill = mycols, #Fills colors
  filename = NULL,
  output = TRUE)
#set 2 
mycols = c("Rhizosphere"="#d6e2e9", "Phyllosphere"="#cbf3f0", "Rhizome"="#fcf5c7","MarineWater"="#FFA07A","RhizosphereSediment"="#D4EFDF","CoastalMarineSediment"="#F4D03F") 
venn.plot <- venn.diagram(
  x = sets_list[c("Marinewater","RhizosphereSediment","CoastaMarineSediment")],
  category.names = c("Marinewater","RhizosphereSediment","CoastalMarineSediment"),
  fill = mycols, #Fills colors
  filename = NULL,
  output = TRUE)
mycols = c("Rhizosphere"="#d6e2e9", "Phyllosphere"="#cbf3f0", "Rhizome"="#fcf5c7","marine water"="#FFA07A","rhizosphere sediment"="#D4EFDF","coastal marine sediment"="#F4D03F") 
# Assuming you are only comparing three sets for simplicity
venn.plot <- venn.diagram(
  x = sets_list[c("marine water","rhizosphere sediment","coastal marine sediment")],
  category.names = c("marine water", "rhizosphere sediment","coastal marine sediment"),
  fill = mycols, #Fills colors
  filename = NULL,
  output = TRUE
)
# Draw the Venn diagram
library(grid)
grid.draw(venn.plot)
#can also make a list of your factors and the color to give to each
#list of all your factors with the color to asign to each (replace nonCRC, CRC, H and add whatever other factors needed)

#plot(venn(sets_list),
    # fills = mycols)

#make a pcoa to evaluate each factor
phy_clr = microbiome::transform(phy, 'clr')
phy_ord= ordinate(phy_clr, "RDA", "euclidean")
plot_ordination(phy,phy_ord, color="env_medium")

sample_info_tab$geo_loc_name
#Test whether the seasons differ significantly from each other using the permutational ANOVA 
dist.uf <- phyloseq::distance(phy_clr, method = "euclidean")
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$env_medium, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
###
install.packages("gp3")

plot_heatmap(gp3, taxa.label = "Genus")

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("gp3")

citation("phyloseq")
