
library(tidyverse)
devtools::install_github("jinyizju/U.PhyloMaker")
library("U.PhyloMaker")
library(Rphylopars)
library(stringr)

# get species list
no.trees.crt1 = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/no.trees.csv", row.names = 1) # 618 data points
no.trees.drt = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/no.trees.csv", row.names = 1) # 485 data points
no.trees.crt2 = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/no.trees.csv", row.names = 1) # 507 data points


species.list.crt1 = as.data.frame(unique(no.trees.crt1$Taxon))
colnames(species.list.crt1)[1] = "species"
species.list.drt = as.data.frame(unique(no.trees.drt$Taxon))
colnames(species.list.drt)[1] = "species"
species.list.crt2 = as.data.frame(unique(no.trees.crt2$Taxon))
colnames(species.list.crt2)[1] = "species"

full.list = merge(species.list.crt1,species.list.drt, by="species", all = TRUE)
full.list.2 = merge(full.list,species.list.crt2, by="species", all = TRUE)

species.list = full.list.2
species.list$genus = word(species.list$species, 1)
colnames(species.list)[1] = "Taxon"

raw.data = read.csv("./Raw.Data/IDE_cover_2023-01-02.csv")
raw.data.2 = raw.data[,c(12,13)]
species.list.2 = left_join(species.list, raw.data.2, by = "Taxon", multiple = "first")

colnames(species.list.2)[1] = "species"
species.list.2$species = str_to_sentence(species.list.2$species)
species.list.2$genus = str_to_sentence(species.list.2$genus)

# https://github.com/jinyizju/U.PhyloMaker

# Download the plant megatree

megatree <- read.tree("plant_megatree.tre")

# Download genus-family relationship file

gen.list <- read.csv("plant_genus_list.csv")

# run the phylomaker program

tree.result <- phylo.maker(species.list.2, megatree, gen.list, nodes.type = 1, scenario = 3)

# Note: 5 species fail to be binded to the tree
# "Heliopsis_helianthoides""Leucanthemopsis_alpina""Millotia_tenuifolia""Ratibida_columnifera""Rudbeckia_hirta" 

# format new species list that adds _ to match tip labels
species.list.3 = species.list.2
trait.data.crt1 = no.trees.crt1[,c(1,11:18)]
trait.data.drt = no.trees.drt[,c(1,6:13)]
trait.data.crt2 = no.trees.crt2[,c(1,12:19)]


trait.data = rbind(trait.data.crt1, trait.data.crt2, trait.data.drt)
trait.data.2 = distinct(trait.data)
trait.data.2$Taxon = str_to_sentence(trait.data.2$Taxon)

species.list.4 = full_join(species.list.3, trait.data.2, by = join_by(species == Taxon))
species.list.4$species <- sub(" ", "_", species.list.4$species)
colnames(species.list.4) = c("species","genus","family","leafN","height","rootN",
                             "SLA","root_depth","RTD","SRL","root_diam")

species.list.5 = species.list.4[,c(1,4:11)]
species.list.6 = species.list.5[-c(220,252,276,345,348),]

trait.sim = phylopars(trait_data = species.list.6, tree = tree.result$phylo, pheno_error = FALSE,
                      phylo_correlated = TRUE, pheno_correlated = FALSE)

write.csv(trait.sim$anc_recon, file = "./Formatted.Data/trait.imputation.csv")
