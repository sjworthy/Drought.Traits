# Script for analysis of Drought Net data
# Comparing mean change in % cover between control and drought in Year 2

library(dismo)
library(gbm)
library(ggBRT)
library(multcomp)

#### read in all the data frames needs for the analyses ####

# all data
all.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/all.data.year2.csv", row.names = 1) # 579 data points
# all data without functional group WOODY
no.trees = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/no.trees.csv", row.names = 1) # 507 data points
# all data without local lifeform = TREE
trees = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/trees.csv", row.names = 1) # 576 data points
# all annual data
annual.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/annual.data.csv", row.names = 1) # 96 data points
# all perennial data
perennial.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.data.csv", row.names = 1) # 465 data points
# perennial data without functional group WOODY
perennial.tree = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.tree.csv", row.names = 1) # 395 data points
# perennial data without local lifeform = TREE
perennial.no.tree = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.no.tree.csv", row.names = 1) # 462 data points
# grasses
grass = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/grass.csv", row.names = 1) # 193 data points
# forbs
forb = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/forb.csv", row.names = 1) # 254 data points

##### imputed trait dataframes ####

impute.traits = read.csv("./Formatted.Data/trait.imputation.2.csv")
colnames(impute.traits)[1] = "Taxon"
impute.traits$Taxon = sub("_", " ", impute.traits$Taxon)

no.trees$Taxon = str_to_sentence(no.trees$Taxon)
no.trees.impute = left_join(no.trees,impute.traits)
no.trees.impute = no.trees.impute[,c(1,4,11,20:30)]

annual.data$Taxon = str_to_sentence(annual.data$Taxon)
annual.data.impute = left_join(annual.data,impute.traits)
annual.data.impute = annual.data.impute[,c(1,4,11,20:30)]

perennial.tree$Taxon = str_to_sentence(perennial.tree$Taxon)
perennial.tree.impute = left_join(perennial.tree,impute.traits)
perennial.tree.impute = perennial.tree.impute[,c(1,4,11,20:30)]

grass$Taxon = str_to_sentence(grass$Taxon)
grass.impute = left_join(grass,impute.traits)
grass.impute = grass.impute[,c(1,4,11,20:30)]

forb$Taxon = str_to_sentence(forb$Taxon)
forb.impute = left_join(forb,impute.traits)
forb.impute = forb.impute[,c(1,4,11,20:30)]

write.csv(no.trees.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/no.trees.impute.csv")
write.csv(annual.data.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/annual.data.impute.csv")
write.csv(perennial.tree.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.tree.impute.csv")
write.csv(grass.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/grass.impute.csv")
write.csv(forb.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/forb.impute.csv")


#### change site code to numeric, continuous vector ####
all.data$site.id = as.numeric(as.factor(all.data$site_code))
no.trees$site.id = as.numeric(as.factor(no.trees$site_code))
trees$site.id = as.numeric(as.factor(trees$site_code))
annual.data$site.id = as.numeric(as.factor(annual.data$site_code))
perennial.data$site.id = as.numeric(as.factor(perennial.data$site_code))
perennial.tree$site.id = as.numeric(as.factor(perennial.tree$site_code))
perennial.no.tree$site.id = as.numeric(as.factor(perennial.no.tree$site_code))
grass$site.id = as.numeric(as.factor(grass$site_code))
forb$site.id = as.numeric(as.factor(forb$site_code))

#### merge with mean annual precipitation data (MAP) ####

Site.info = read.csv("./Raw.Data/Site_Elev-Disturb.csv")
site.info.map = Site.info[,c(2,13)]

all.data = merge(all.data, site.info.map, by="site_code")
no.trees = merge(no.trees, site.info.map, by="site_code")
trees = merge(trees, site.info.map, by="site_code")
annual.data = merge(annual.data, site.info.map, by="site_code")
perennial.data = merge(perennial.data, site.info.map, by="site_code")
perennial.tree = merge(perennial.tree, site.info.map, by="site_code")
perennial.no.tree = merge(perennial.no.tree, site.info.map, by="site_code")
grass = merge(grass, site.info.map, by="site_code")
forb = merge(forb, site.info.map, by="site_code")


#### determining best parameter combination to generate 1000 trees ####
# bag fraction of 0.50 and 0.75, step.size of 25 and 50, tc = 10

set.seed(2023)
all.brt.1=gbm.step(data=all.data, gbm.x = c(12:19), gbm.y=11,
                   family = "gaussian", tree.complexity = 10, learning.rate = 0.00000001,
                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                   site.weights = all.data$site.id)

ggPerformance(all.brt.1)

set.seed(2023)
all.brt.1.no.site=gbm.step(data=all.data, gbm.x = c(12:19,24), gbm.y=11,
                           family = "gaussian", tree.complexity = 1, learning.rate = 0.0000001,
                           bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(all.brt.1.no.site)

set.seed(2023)
tree.brt.1=gbm.step(data=no.trees, gbm.x = c(12:19,24), gbm.y=11,
                    family = "gaussian", tree.complexity = 10, learning.rate = 0.0000000001,
                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                    site.weights = no.trees$site.id)

ggPerformance(tree.brt.1)


set.seed(2023)

tree.brt.1.no.site=gbm.step(data=no.trees, gbm.x = c(12:19), gbm.y=11,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0000001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(tree.brt.1.no.site)


set.seed(2023)
annual.brt.1=gbm.step(data=annual.data, gbm.x = c(12:19), gbm.y=11,
                      family = "laplace", tree.complexity = 10, learning.rate = 0.0001,
                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50,
                      site.weights = annual.data$site.id)
ggPerformance(annual.brt.1)

set.seed(2023)

annual.brt.1.no.site=gbm.step(data=annual.data, gbm.x = c(12,14,25), gbm.y=11,
                      family = "gaussian", tree.complexity = 3, learning.rate = 0.00000000001,
                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
                      
ggPerformance(annual.brt.1.no.site)

set.seed(2023)
perennial.brt.1=gbm.step(data=perennial.data, gbm.x = c(12:19), gbm.y=11,
                         family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                         bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25,
                         site.weights = perennial.data$site.id)
ggPerformance(perennial.brt.1)

set.seed(2023)
perennial.brt.1.no.site=gbm.step(data=perennial.data, gbm.x = c(12:19), gbm.y=11,
                                 family = "gaussian", tree.complexity = 10, learning.rate = 0.000001,
                                 bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(perennial.brt.1.no.site)




set.seed(2023)
perennial.tree.brt.1=gbm.step(data=perennial.tree, gbm.x = c(12:19,24), gbm.y=11,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                              bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50,
                              site.weights = perennial.tree$site.id)
ggPerformance(perennial.tree.brt.1)

set.seed(2023)

perennial.tree.brt.1.no.site=gbm.step(data=perennial.tree, gbm.x = c(12:19,24), gbm.y=11,
                                      family = "gaussian", tree.complexity = 10, learning.rate = 0.00000000001,
                                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.tree.brt.1.no.site)

set.seed(2023)
grass.brt.1=gbm.step(data=grass, gbm.x = c(12:19,24), gbm.y=11,
                     family = "gaussian", tree.complexity = 10, learning.rate = 0.00005,
                     bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50, 
                     site.weights = grass$site.id)
ggPerformance(grass.brt.1)

set.seed(2023)

grass.brt.1.no.site=gbm.step(data=grass, gbm.x = c(12:19,24), gbm.y=11,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.00000001,
                             bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(grass.brt.1.no.site)

set.seed(2023)
forb.brt.1=gbm.step(data=forb, gbm.x = c(12:19), gbm.y=11,
                    family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                    site.weights = forb$site.id)
ggPerformance(forb.brt.1)

set.seed(2023)

forb.brt.1.no.site=gbm.step(data=forb, gbm.x = c(12:19,24), gbm.y=11,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0000001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.brt.1.no.site)

set.seed(2023)
annual.grass.brt.1=gbm.step(data=annual.grass, gbm.x = c(11:18,23), gbm.y=10,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                            site.weights = annual.grass$site.id)

ggPerformance(annual.grass.brt.1)

set.seed(2023)
annual.grass.brt.1.no.site=gbm.step(data=annual.grass, gbm.x = c(11:18,23), gbm.y=10,
                                    family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.grass.brt.1.no.site)

set.seed(2023)
annual.forb.brt.1=gbm.step(data=annual.forb, gbm.x = c(11:18,23), gbm.y=10,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                           bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25,
                           site.weights = annual.forb$site.id)
ggPerformance(annual.forb.brt.1)

set.seed(2023)
annual.forb.brt.1.no.site=gbm.step(data=annual.forb, gbm.x = c(11:18,23), gbm.y=10,
                                   family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(annual.forb.brt.1.no.site)

set.seed(2023)
perennial.grass.brt.1=gbm.step(data=perennial.grass, gbm.x = c(11:18,23), gbm.y=10,
                               family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                               bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                               site.weights = perennial.grass$site.id)
ggPerformance(perennial.grass.brt.1)

set.seed(2023)
perennial.grass.brt.1.no.site=gbm.step(data=perennial.grass, gbm.x = c(11:18,23), gbm.y=10,
                                       family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                       bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(perennial.grass.brt.1.no.site)

set.seed(2023)
perennial.forb.brt.1=gbm.step(data=perennial.forb, gbm.x = c(11:18,23), gbm.y=10,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                              bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25, 
                              site.weights = perennial.forb$site.id)

ggPerformance(perennial.forb.brt.1)

set.seed(2023)
perennial.forb.brt.1.no.site=gbm.step(data=perennial.forb, gbm.x = c(11:18,23), gbm.y=10,
                                      family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(perennial.forb.brt.1.no.site)

set.seed(2023)
trees.brt.1=gbm.step(data=trees, gbm.x = c(11:18,23), gbm.y=10,
                     family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                     bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25, 
                     site.weights = trees$site.id)

ggPerformance(trees.brt.1)

set.seed(2023)
trees.brt.1.no.site=gbm.step(data=trees, gbm.x = c(11:18,23), gbm.y=10,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25) 

ggPerformance(trees.brt.1.no.site)

set.seed(2023)
perennial.no.tree.brt.1=gbm.step(data=perennial.no.tree, gbm.x = c(11:18,23), gbm.y=10,
                                 family = "gaussian", tree.complexity = 10, learning.rate = 0.0000001,
                                 bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                                 site.weights = perennial.no.tree$site.id)
ggPerformance(perennial.no.tree.brt.1)

set.seed(2023)
perennial.no.tree.brt.1.no.site=gbm.step(data=perennial.no.tree, gbm.x = c(11:18,23), gbm.y=10,
                                         family = "gaussian", tree.complexity = 10, learning.rate = 0.000001,
                                         bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.no.tree.brt.1.no.site)




#### outliers ####

# all outliers removed manually from trait.species.trt.yr1.final.new and made new file trait.species.trt.yr1.outlier
mean = mean(annual.data$rootDiam.mm, na.rm = TRUE)
std = sd(annual.data$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(2*std)
Tmax = mean+(2*std)
annual.data$rootDiam.mm[which(annual.data$rootDiam.mm <Tmin | annual.data$rootDiam.mm > Tmax)]



