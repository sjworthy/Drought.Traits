# Script for analysis of Drought Net data
# https://github.com/JBjouffray/Hawaii_RegimesPredictors for visuals and plotting

library(dismo)
library(gbm)
library(ggBRT)

trait.data.new = read.csv("./Formatted.Data/trait.species.trt.yr1.outlier.csv")

# subset traits so they must have SLA
trait.data.2 = trait.data.new[,c(1,7,8,10,12,14,15,18,20,27:29)]
trait.data.3 = subset(trait.data.2, trait.data.2$SLA_m2.kg > 0 ) # 645 data points, 

# read in cover data
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")

# remove individuals where 0 in control or drought

cover.2 = subset(cover.data, cover.data$mean.ctrl.cover > 0)
cover.3 = subset(cover.data, cover.data$mean.drt.cover > 0)

# merge new cover data with traits

all.data = merge(cover.3, trait.data.3, by="Taxon") # 856 data points with site

# merge with mean annual precipitation data (MAP)

Site.info = read.csv("./Raw.Data/Site_Elev-Disturb.csv")
site.info.map = Site.info[,c(2,13)]

all.data = merge(all.data, site.info.map, by="site_code")

#### test for correlation ####
cor.test(all.data$leafN.mg.g, all.data$height.m)
cor.test(all.data$leafN.mg.g, all.data$rootN.mg.g) # correlated r = 0.52
cor.test(all.data$leafN.mg.g, all.data$SLA_m2.kg) # correlated r = 0.35
cor.test(all.data$leafN.mg.g, all.data$root.depth_m)
cor.test(all.data$leafN.mg.g, all.data$RTD.g.cm3) # correlated r = -0.20
cor.test(all.data$leafN.mg.g, all.data$SRL_m.g) # correlated r = 0.27
cor.test(all.data$leafN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$height.m, all.data$rootN.mg.g)
cor.test(all.data$height.m, all.data$SLA_m2.kg)
cor.test(all.data$height.m, all.data$root.depth_m)
cor.test(all.data$height.m, all.data$RTD.g.cm3)
cor.test(all.data$height.m, all.data$SRL_m.g)
cor.test(all.data$height.m, all.data$rootDiam.mm) # correlated r = 0.17
cor.test(all.data$rootN.mg.g, all.data$SLA_m2.kg) # correlated r = 0.24
cor.test(all.data$rootN.mg.g, all.data$root.depth_m)
cor.test(all.data$rootN.mg.g, all.data$RTD.g.cm3) # correlated r = -0.15
cor.test(all.data$rootN.mg.g, all.data$SRL_m.g)
cor.test(all.data$rootN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$SLA_m2.kg, all.data$root.depth_m) # correlated r = -0.13
cor.test(all.data$SLA_m2.kg, all.data$RTD.g.cm3) # correlated r = -0.26
cor.test(all.data$SLA_m2.kg, all.data$SRL_m.g) # correlated r = 0.42
cor.test(all.data$SLA_m2.kg, all.data$rootDiam.mm)
cor.test(all.data$root.depth_m, all.data$RTD.g.cm3) # correlated r = 0.15
cor.test(all.data$root.depth_m, all.data$SRL_m.g) # correlated r =-0.13
cor.test(all.data$root.depth_m, all.data$rootDiam.mm) #  correlated r = 0.19
cor.test(all.data$RTD.g.cm3, all.data$SRL_m.g) # correlated r = -0.29
cor.test(all.data$RTD.g.cm3, all.data$rootDiam.mm)
cor.test(all.data$SRL_m.g, all.data$rootDiam.mm)

#### correlation between response and predictors ####

cor.test(all.data$mean.cover.response, all.data$leafN.mg.g) # r = 0.08
cor.test(all.data$mean.cover.response, all.data$height.m) # r = 0.003
cor.test(all.data$mean.cover.response, all.data$rootN.mg.g) # r = -0.02
cor.test(all.data$mean.cover.response, all.data$SLA_m2.kg) # r = 0.01
cor.test(all.data$mean.cover.response, all.data$root.depth_m) # r = 0.04
cor.test(all.data$mean.cover.response, all.data$RTD.g.cm3) # r = 0.08
cor.test(all.data$mean.cover.response, all.data$SRL_m.g) # r = -0.08
cor.test(all.data$mean.cover.response, all.data$rootDiam.mm) # r = 0.07

shapiro.test(all.data$mean.cover.response)
hist(all.data$mean.cover.response)

#### data set with trees and shrubs removed ####

no.trees = subset(all.data, !all.data$functional_group == "WOODY")

#### change site code to numeric, continuous vector ####
all.data$site.id = as.numeric(as.factor(all.data$site_code))
no.trees$site.id = as.numeric(as.factor(no.trees$site_code))

#### data set split by lifespan ####

# need to read out data and fix the lifespan to particular sites since some species have different lifespan at different sites

# write.csv(all.data, file="./Formatted.Data/all.data.response.0.csv")

all.data.ls = read.csv("./Formatted.Data/all.data.lifespan.response.0.csv", row.names = 1)
table(all.data.ls$local_lifespan)

#merge with MAP
all.data.ls = merge(all.data.ls, site.info.map, by = "site_code")

annual.data = subset(all.data.ls, all.data.ls$local_lifespan == "ANNUAL")
perennial.data = subset(all.data.ls, all.data.ls$local_lifespan == "PERENNIAL")

annual.data$site.id = as.numeric(as.factor(annual.data$site_code))
perennial.data$site.id = as.numeric(as.factor(perennial.data$site_code))

perennial.tree = subset(perennial.data, !perennial.data$functional_group == "WOODY")
perennial.tree$site.id = as.numeric(as.factor(perennial.tree$site_code))

#### split data by functional group ####

grass = subset(all.data, all.data$functional_group == "GRASS")
forb = subset(all.data, all.data$functional_group == "FORB")

grass$site.id = as.numeric(as.factor(grass$site_code))
forb$site.id = as.numeric(as.factor(forb$site_code))

annual.grass = subset(annual.data, annual.data$functional_group == "GRASS")
annual.forb = subset(annual.data, annual.data$functional_group == "FORB")

perennial.grass = subset(perennial.data, perennial.data$functional_group == "GRASS")
perennial.forb = subset(perennial.data, perennial.data$functional_group == "FORB")

annual.grass$site.id = as.numeric(as.factor(annual.grass$site_code))
annual.forb$site.id = as.numeric(as.factor(annual.forb$site_code))

perennial.grass$site.id = as.numeric(as.factor(perennial.grass$site_code))
perennial.forb$site.id = as.numeric(as.factor(perennial.forb$site_code))

#### determining best learning rate to generate 1000 trees ####

set.seed(2023)
all.brt.1=gbm.step(data=all.data, gbm.x = c(11:18,22), gbm.y=10,
                   family = "gaussian", tree.complexity = 10, learning.rate = 0.00005,
                   bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50,
                   site.weights = all.data$site.id)

ggPerformance(all.brt.1)

set.seed(2023)
all.brt.1.no.site=gbm.step(data=all.data, gbm.x = c(11:18,22), gbm.y=10,
                   family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                   bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(all.brt.1.no.site)

set.seed(2023)
tree.brt.1=gbm.step(data=no.trees, gbm.x = c(11:18,22), gbm.y=10,
                    family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                    site.weights = no.trees$site.id)

ggPerformance(tree.brt.1)


set.seed(2023)
tree.brt.1.no.site=gbm.step(data=no.trees, gbm.x = c(11:18,22), gbm.y=10,
                    family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(tree.brt.1.no.site)


set.seed(2023)
annual.brt.1=gbm.step(data=annual.data, gbm.x = c(11:18,23), gbm.y=10,
                      family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                      bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25,
                      site.weights = annual.data$site.id)
ggPerformance(annual.brt.1)


set.seed(2023)
annual.brt.1.no.site=gbm.step(data=annual.data, gbm.x = c(11:18,23), gbm.y=10,
                      family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                      bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(annual.brt.1.no.site)

set.seed(2023)
perennial.brt.1=gbm.step(data=perennial.data, gbm.x = c(11:18,23), gbm.y=10,
                         family = "gaussian", tree.complexity = 10, learning.rate = 0.000001,
                         bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25,
                         site.weights = perennial.data$site.id)
ggPerformance(perennial.brt.1)

set.seed(2023)
perennial.brt.1.no.site=gbm.step(data=perennial.data, gbm.x = c(11:18,23), gbm.y=10,
                         family = "gaussian", tree.complexity = 10, learning.rate = 0.000005,
                         bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.brt.1.no.site)

set.seed(2023)
perennial.tree.brt.1=gbm.step(data=perennial.tree, gbm.x = c(11:18,23), gbm.y=10,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                              bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25,
                              site.weights = perennial.tree$site.id)
ggPerformance(perennial.tree.brt.1)

set.seed(2023)
perennial.tree.brt.1.no.site=gbm.step(data=perennial.tree, gbm.x = c(11:18,23), gbm.y=10,
                              family = "gaussian", tree.complexity = 11, learning.rate = 0.0001,
                              bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.tree.brt.1.no.site)

set.seed(2023)
grass.brt.1=gbm.step(data=grass, gbm.x = c(11:18,22), gbm.y=10,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                              bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25, 
                     site.weights = grass$site.id)
ggPerformance(grass.brt.1)

set.seed(2023)
grass.brt.1.no.site=gbm.step(data=grass, gbm.x = c(11:18,22), gbm.y=10,
                     family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                     bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(grass.brt.1.no.site)

set.seed(2023)
forb.brt.1=gbm.step(data=forb, gbm.x = c(11:18,22), gbm.y=10,
                     family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                     bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                    site.weights = forb$site.id)
ggPerformance(forb.brt.1)

set.seed(2023)
forb.brt.1.no.site=gbm.step(data=forb, gbm.x = c(11:18,22), gbm.y=10,
                    family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)

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

#### Best Models ####

set.seed(2023)
all.brt.no.site.map=gbm.step(data=all.data, gbm.x = c(11:18,22), gbm.y=10,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                           bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(all.brt.no.site.map)
# 1150 trees Per.Expl = 9.68788975
1-(all.brt.no.site.map$self.statistics$mean.resid/all.brt.no.site.map$self.statistics$mean.null) # R2
ggInfluence(all.brt.no.site.map)

ggPD(all.brt.no.site.map, rug = T) # partial dependency plots
gbm.plot(all.brt.no.site.map, common.scale = FALSE)
gbm.plot.fits(all.brt.no.site.map)

# investigation of interactions
gbm.interactions(all.brt.no.site.map)$rank.list
ggInteract_list(all.brt.no.site.map)
# RTD x height 4.46
# height x leafN 3.29
# SRL x leafN 2.56
# precip x SRL 1.86

ggInteract_3D(all.brt.no.site.map, x = 6, y = 2, z.range = c(-1.5, 1))
ggInteract_3D(all.brt.no.site.map, x = 2, y = 1,z.range = c(-1, 1.1))
ggInteract_3D(all.brt.no.site.map, x = 7, y = 1, z.range = c(-0.5, 1.5))
ggInteract_3D(all.brt.no.site.map, x = 9, y = 7, z.range = c(0, 1.2))

set.seed(2023)
tree.no.site.map=gbm.step(data=no.trees, gbm.x = c(11:18,22), gbm.y=10,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(tree.no.site.map)
# 1650 trees Per.Expl = 14.71%
1-(tree.no.site.map$self.statistics$mean.resid/tree.no.site.map$self.statistics$mean.null) # R2
ggInfluence(tree.no.site.map)

ggPD(tree.no.site.map, rug = T) # partial dependency plots
gbm.plot(tree.no.site.map, common.scale = FALSE)
gbm.plot.fits(tree.no.site.map)

# investigation of interactions
gbm.interactions(tree.no.site.map)$rank.list
ggInteract_list(tree.no.site.map)
# RTD x height 15.35
# height x leafN 8.64
# precip x SLA 4.05
# SRL x leafN 3.67

ggInteract_3D(tree.no.site.map, x = 6, y = 2, z.range = c(-2.5, 1.2))
ggInteract_3D(tree.no.site.map, x = 2, y = 1,z.range = c(-2.5, 1.75))
ggInteract_3D(tree.no.site.map, x = 9, y = 4, z.range = c(0, 1.2))
ggInteract_3D(tree.no.site.map, x = 7, y = 1, z.range = c(-1, 1.2))

set.seed(2023)
annual.no.site.map=gbm.step(data=annual.data, gbm.x = c(11:18,23), gbm.y=10,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                              bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.no.site.map)
# 2050 trees Per.Expl = 37.13%
1-(annual.no.site.map$self.statistics$mean.resid/annual.no.site.map$self.statistics$mean.null) # R2
ggInfluence(annual.no.site.map)

ggPD(annual.no.site.map, rug = T) # partial dependency plots
gbm.plot(annual.no.site.map, common.scale = FALSE)
gbm.plot.fits(annual.no.site.map)

# investigation of interactions
gbm.interactions(annual.no.site.map)$rank.list
ggInteract_list(annual.no.site.map)
# SRL x leafN 40.56
# precip x leafN 28.86
# diam x SRL 14.01
# precip x root.depth 10.55

ggInteract_3D(annual.no.site.map, x = 7, y = 1, z.range = c(-2.5, 1.5))
ggInteract_3D(annual.no.site.map, x = 9, y = 1,z.range = c(-1.5, 1.5))
ggInteract_3D(annual.no.site.map, x = 8, y = 7, z.range = c(-1.5, 1))
ggInteract_3D(annual.no.site.map, x = 9, y = 5, z.range = c(-1, 0.8))

set.seed(2023)
perennial.brt.1.no.site=gbm.step(data=perennial.data, gbm.x = c(11:18,23), gbm.y=10,
                                 family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                                 bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)
# NA

set.seed(2023)
perennial.tree.no.site.map=gbm.step(data=perennial.tree, gbm.x = c(11:18,23), gbm.y=10,
                                      family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.tree.no.site.map)
# 1000 trees Per.Expl = 2.40%
1-(perennial.tree.no.site.map$self.statistics$mean.resid/perennial.tree.no.site.map$self.statistics$mean.null) # R2
ggInfluence(perennial.tree.no.site.map)

ggPD(perennial.tree.no.site.map, rug = T) # partial dependency plots
gbm.plot(perennial.tree.no.site.map, common.scale = FALSE)
gbm.plot.fits(perennial.tree.no.site.map)

# investigation of interactions
gbm.interactions(perennial.tree.no.site.map)$rank.list
ggInteract_list(perennial.tree.no.site.map)
# height x leafN 0.37
# SRL x leafN 0.35
# SRL x height 0.27
# precip x SLA 0.19

ggInteract_3D(perennial.tree.no.site.map, x = 2, y = 1, z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 7, y = 1,z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 7, y = 2, z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 9, y = 4, z.range = c(0.2, 0.75))

set.seed(2023)
grass.no.site.map=gbm.step(data=grass, gbm.x = c(11:18,22), gbm.y=10,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                             bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(grass.no.site.map)
# 1000 trees Per.Expl = 1.84%
1-(grass.no.site.map$self.statistics$mean.resid/grass.no.site.map$self.statistics$mean.null) # R2
ggInfluence(grass.no.site.map)

ggPD(grass.no.site.map, rug = T) # partial dependency plots
gbm.plot(grass.no.site.map, common.scale = FALSE)
gbm.plot.fits(grass.no.site.map)

# investigation of interactions
gbm.interactions(grass.no.site.map)$rank.list
ggInteract_list(grass.no.site.map)
# precip x diam 0.09
# precip x RTD 0.08
# diam x RTD 0.08
# diam x rootN 0.05

ggInteract_3D(grass.no.site.map, x = 9, y = 8, z.range = c(-0.05, 0.20))
ggInteract_3D(grass.no.site.map, x = 9, y = 6,z.range = c(-0.2, 0.20))
ggInteract_3D(grass.no.site.map, x = 8, y = 6, z.range = c(-0.2, 0.30))
ggInteract_3D(grass.no.site.map, x = 8, y = 3, z.range = c(-0.2, 0.20))


set.seed(2023)
forb.no.site.map=gbm.step(data=forb, gbm.x = c(11:18,22), gbm.y=10,
                            family = "gaussian", tree.complexity = 9, learning.rate = 0.001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.no.site.map)
# 2850 trees Per.Expl = 42.57%
1-(forb.no.site.map$self.statistics$mean.resid/forb.no.site.map$self.statistics$mean.null) # R2
ggInfluence(forb.no.site.map)

ggPD(forb.no.site.map, rug = T) # partial dependency plots
gbm.plot(forb.no.site.map, common.scale = FALSE)
gbm.plot.fits(forb.no.site.map)

# investigation of interactions
gbm.interactions(forb.no.site.map)$rank.list
ggInteract_list(forb.no.site.map)
# RTD x leafN 19.22
# depth x SLA 14.26
# precip x leafN 12.26
# height x leafN 9.15

ggInteract_3D(forb.no.site.map, x = 6, y = 1, z.range = c(-2.0, 7))
ggInteract_3D(forb.no.site.map, x = 5, y = 4,z.range = c(3, 5))
ggInteract_3D(forb.no.site.map, x = 9, y = 1, z.range = c(-0.5, 7))
ggInteract_3D(forb.no.site.map, x = 2, y = 1, z.range = c(0, 7))

set.seed(2023)
annual.grass.no.site.map=gbm.step(data=annual.grass, gbm.x = c(11:18,23), gbm.y=10,
                                    family = "gaussian", tree.complexity = 2, learning.rate = 0.001,
                                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.grass.no.site.map)
# 1950 trees Per.Expl = 34.72%
1-(annual.grass.no.site.map$self.statistics$mean.resid/annual.grass.no.site.map$self.statistics$mean.null) # R2
ggInfluence(annual.grass.no.site.map)

ggPD(annual.grass.no.site.map, rug = T) # partial dependency plots
gbm.plot(annual.grass.no.site.map, common.scale = FALSE)
gbm.plot.fits(annual.grass.no.site.map)

# investigation of interactions
gbm.interactions(annual.grass.no.site.map)$rank.list
ggInteract_list(annual.grass.no.site.map)
# height x leafN 5.82
# precip x height 5.33
# precip x SLA 0.06
# precip x depth 0.04

ggInteract_3D(annual.grass.no.site.map, x = 2, y = 1, z.range = c(-0.1, 4))
ggInteract_3D(annual.grass.no.site.map, x = 9, y = 2,z.range = c(-1, 4))
ggInteract_3D(annual.grass.no.site.map, x = 9, y = 4, z.range = c(-2, 3))
ggInteract_3D(annual.grass.no.site.map, x = 9, y = 5, z.range = c(-2, 3))

set.seed(2023)
annual.forb.no.site.map=gbm.step(data=annual.forb, gbm.x = c(11:18,23), gbm.y=10,
                                   family = "gaussian", tree.complexity = 6, learning.rate = 0.001,
                                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.forb.no.site.map)
# 1950 trees Per.Expl = 34.72%
1-(annual.forb.no.site.map$self.statistics$mean.resid/annual.forb.no.site.map$self.statistics$mean.null) # R2
ggInfluence(annual.forb.no.site.map)

ggPD(annual.forb.no.site.map, rug = T) # partial dependency plots
gbm.plot(annual.forb.no.site.map, common.scale = FALSE)
gbm.plot.fits(annual.forb.no.site.map)

# investigation of interactions
gbm.interactions(annual.forb.no.site.map)$rank.list
ggInteract_list(annual.forb.no.site.map)
# SLA x leafN 30.11
# precip x SLA 15.70
# precip x height 14.08
# SLA x height 13.23

ggInteract_3D(annual.forb.no.site.map, x = 4, y = 1, z.range = c(-2.5, 3.5))
ggInteract_3D(annual.forb.no.site.map, x = 9, y = 4,z.range = c(-4, 2))
ggInteract_3D(annual.forb.no.site.map, x = 9, y = 2, z.range = c(-2, 5.5))
ggInteract_3D(annual.forb.no.site.map, x = 4, y = 2, z.range = c(-3.5, 5))

set.seed(2023)
perennial.grass.brt.1.no.site=gbm.step(data=perennial.grass, gbm.x = c(11:18,23), gbm.y=10,
                                       family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                       bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(perennial.grass.brt.1.no.site)

set.seed(2023)
perennial.forb.brt.1.no.site=gbm.step(data=perennial.forb, gbm.x = c(11:18,23), gbm.y=10,
                                      family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(perennial.forb.brt.1.no.site)


#### initial model evaluation ####
# plots of predicted versus observed, estimate R2

# which predictors are most important
# all data without site
summary(all.brt.1)

all.predict.brt=predict(all.brt.1, n.trees = all.brt.1$n.trees)
observed = all.data$mean.cover.response
plot(all.predict.brt~observed)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.brt.1$self.statistics$mean.resid/all.brt.1$self.statistics$mean.null)
# 0.0004

# which predictors are most important
# all data without site
summary(all.brt.1.no.site)

all.predict.brt=predict(all.brt.1.no.site, n.trees = all.brt.1.no.site$n.trees)
observed = all.data$mean.cover.response
plot(all.predict.brt~observed)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.brt.1.no.site$self.statistics$mean.resid/all.brt.1.no.site$self.statistics$mean.null)
# 0.003

# which predictors are most important
# all data without tree
summary(tree.brt.1)

tree.predict.brt=predict(tree.brt.1, n.trees = tree.brt.1$n.trees)
observed = no.trees$mean.cover.response
plot(tree.predict.brt~observed)
abline(0, 1, col = 2)

R2.tree.brt = 1-(tree.brt.1$self.statistics$mean.resid/tree.brt.1$self.statistics$mean.null)
# 0.0005

# which predictors are most important
# all data without tree without site
summary(tree.brt.1.no.site)

tree.predict.brt=predict(tree.brt.1.no.site, n.trees = tree.brt.1.no.site$n.trees)
observed = no.trees$mean.cover.response
plot(tree.predict.brt~observed)
abline(0, 1, col = 2)

R2.tree.brt = 1-(tree.brt.1.no.site$self.statistics$mean.resid/tree.brt.1.no.site$self.statistics$mean.null)
# 0.04

# which predictors are most important
# only annuals
summary(annual.brt.1)

annual.predict.brt=predict(annual.brt.1, n.trees = annual.brt.1$n.trees)
observed = annual.data$mean.cover.response
plot(annual.predict.brt~observed)
abline(0, 1, col = 2)

R2.annual.brt = 1-(annual.brt.1$self.statistics$mean.resid/annual.brt.1$self.statistics$mean.null)
# 0.19

# which predictors are most important
# annuals without site
summary(annual.brt.1.no.site)

annual.predict.brt=predict(annual.brt.1.no.site, n.trees = annual.brt.1.no.site$n.trees)
observed = annual.data$mean.cover.response
plot(annual.predict.brt~observed)
abline(0, 1, col = 2)

R2.annual.brt = 1-(annual.brt.1.no.site$self.statistics$mean.resid/annual.brt.1.no.site$self.statistics$mean.null)
# 0.19

# which predictors are most important
# only perennials
summary(perennial.brt.1)

perennial.predict.brt=predict(perennial.brt.1, n.trees = perennial.brt.1$n.trees)
observed = perennial.data$mean.cover.response
plot(perennial.predict.brt~observed)
abline(0, 1, col = 2)

R2.perennial.brt = 1-(perennial.brt.1$self.statistics$mean.resid/perennial.brt.1$self.statistics$mean.null)
# 0.0006

# which predictors are most important
# only grasses
summary(grass.brt.1)

grass.predict.brt=predict(grass.brt.1, n.trees = grass.brt.1$n.trees)
observed = grass$mean.cover.response
plot(grass.predict.brt~observed)
abline(0, 1, col = 2)

R2.grass.brt = 1-(grass.brt.1$self.statistics$mean.resid/grass.brt.1$self.statistics$mean.null)
# 0.01

# which predictors are most important
# only grasses
summary(grass.brt.1.no.site)

grass.predict.brt=predict(grass.brt.1.no.site, n.trees = grass.brt.1.no.site$n.trees)
observed = grass$mean.cover.response
plot(grass.predict.brt~observed)
abline(0, 1, col = 2)

R2.grass.brt = 1-(grass.brt.1.no.site$self.statistics$mean.resid/grass.brt.1.no.site$self.statistics$mean.null)
# 0.008

# which predictors are most important
# only forbs
summary(forb.brt.1)
# RTD 89%

forb.predict.brt=predict(forb.brt.1, n.trees = forb.brt.1$n.trees)
observed = forb$mean.cover.response
plot(forb.predict.brt~observed)
abline(0, 1, col = 2)

R2.forb.brt = 1-(forb.brt.1$self.statistics$mean.resid/forb.brt.1$self.statistics$mean.null)
# 0.07

# which predictors are most important
# only forbs
summary(forb.brt.1.no.site)
# RTD 89%

forb.predict.brt=predict(forb.brt.1.no.site, n.trees = forb.brt.1.no.site$n.trees)
observed = forb$mean.cover.response
plot(forb.predict.brt~observed)
abline(0, 1, col = 2)

R2.forb.brt = 1-(forb.brt.1.no.site$self.statistics$mean.resid/forb.brt.1.no.site$self.statistics$mean.null)
# 0.14

# which predictors are most important
# only annual.grasses
summary(annual.grass.brt.1)

annual.grass.predict.brt=predict(annual.grass.brt.1, n.trees = annual.grass.brt.1$n.trees)
observed = annual.grass$mean.cover.response
plot(annual.grass.predict.brt~observed)
abline(0, 1, col = 2)

R2.annual.grass.brt = 1-(annual.grass.brt.1$self.statistics$mean.resid/annual.grass.brt.1$self.statistics$mean.null)
# 0.16

# which predictors are most important
# only annual.grasses
summary(annual.grass.brt.1.no.site)

annual.grass.predict.brt=predict(annual.grass.brt.1.no.site, n.trees = annual.grass.brt.1.no.site$n.trees)
observed = annual.grass$mean.cover.response
plot(annual.grass.predict.brt~observed)
abline(0, 1, col = 2)

R2.annual.grass.brt = 1-(annual.grass.brt.1.no.site$self.statistics$mean.resid/annual.grass.brt.1.no.site$self.statistics$mean.null)
# 0.19

# which predictors are most important
# only annual.forbes
summary(annual.forb.brt.1)
# RTD 89%

annual.forb.predict.brt=predict(annual.forb.brt.1, n.trees = annual.forb.brt.1$n.trees)
observed = annual.forb$mean.cover.response
plot(annual.forb.predict.brt~observed)
abline(0, 1, col = 2)

R2.annual.forb.brt = 1-(annual.forb.brt.1$self.statistics$mean.resid/annual.forb.brt.1$self.statistics$mean.null)
# 0.16

# which predictors are most important
# only annual.forbes
summary(annual.forb.brt.1.no.site)
# RTD 89%

annual.forb.predict.brt=predict(annual.forb.brt.1.no.site, n.trees = annual.forb.brt.1.no.site$n.trees)
observed = annual.forb$mean.cover.response
plot(annual.forb.predict.brt~observed)
abline(0, 1, col = 2)

R2.annual.forb.brt = 1-(annual.forb.brt.1.no.site$self.statistics$mean.resid/annual.forb.brt.1.no.site$self.statistics$mean.null)
# 0.13

# which predictors are most important
# only perennial.grasses
summary(perennial.grass.brt.1)

perennial.grass.predict.brt=predict(perennial.grass.brt.1, n.trees = perennial.grass.brt.1$n.trees)
observed = perennial.grass$mean.cover.response
plot(perennial.grass.predict.brt~observed)
abline(0, 1, col = 2)

R2.perennial.grass.brt = 1-(perennial.grass.brt.1$self.statistics$mean.resid/perennial.grass.brt.1$self.statistics$mean.null)
# 0.001

# which predictors are most important
# only perennial.grasses
summary(perennial.grass.brt.1.no.site)

perennial.grass.predict.brt=predict(perennial.grass.brt.1.no.site, n.trees = perennial.grass.brt.1.no.site$n.trees)
observed = perennial.grass$mean.cover.response
plot(perennial.grass.predict.brt~observed)
abline(0, 1, col = 2)

R2.perennial.grass.brt = 1-(perennial.grass.brt.1.no.site$self.statistics$mean.resid/perennial.grass.brt.1.no.site$self.statistics$mean.null)
# 0.008


# which predictors are most important
# only perennial.forbes
summary(perennial.forb.brt.1.no.site)
# RTD 89%

perennial.forb.predict.brt=predict(perennial.forb.brt.1.no.site, n.trees = perennial.forb.brt.1.no.site$n.trees)
observed = perennial.forb$mean.cover.response
plot(perennial.forb.predict.brt~observed)
abline(0, 1, col = 2)

R2.perennial.forb.brt = 1-(perennial.forb.brt.1.no.site$self.statistics$mean.resid/perennial.forb.brt.1.no.site$self.statistics$mean.null)
# 0.05

#### determining best tree complexity ####
# learning rate 0.001 only works with tree complexity = 1

# all.data without site

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:9),
                                                             ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    set.seed(2023)
    BRT.all.variables <- gbm.step(data=perennial.forb,
                                  gbm.x = c(11:18,23),
                                  gbm.y = 10,
                                  family = "gaussian",
                                  tree.complexity = tcomp,
                                  learning.rate = 0.0005,
                                  bag.fraction = 0.50,
                                  n.trees = 50,
                                  step.size = 50,
                                  plot.main=F, plot.folds=F)
                                  
                                  
    #R2 adj:
    R2Obs.all.variables[[tcomp]][i] <- 1 - (BRT.all.variables$self.statistics$mean.resid /
                                              BRT.all.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.all.variables[[tcomp]]) <- sort(rownames(summary(BRT.all.variables)))
    }
    importancePred.all.variables[[tcomp]][, i] <-
      summary(BRT.all.variables)[rownames(importancePred.all.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.all.variables, mean)
sds <- sapply(R2Obs.all.variables, sd)
plot(1:length(R2Obs.all.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.all.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

tcFactor <- as.factor(rep(1:10, each=nreps))
R2Vector <- unlist(R2Obs.all.variables)
model <- lm(R2Vector~tcFactor)
library(multcomp)
TukeyModel<-glht(model, linfct = mcp(tcFactor="Tukey"))
TukeyLetters <- cld(TukeyModel)$mcletters$Letters

plot(1:length(R2Obs.all.variables), means, ylim=c(0.63, 0.78))
for (i in 1:length(R2Obs.all.variables)){
  arrows(x0=i, x1=i, y0=means[i]-sds[i], y1=means[i]+sds[i], angle=90, code=3,
         length=0.1)
}
text(x= 1:length(R2Obs.all.variables), y= 0.77, labels=TukeyLetters)

# can't get past tc = 1

# tree variables with site

R2Obs.tree.variables <- list()
importancePred.tree.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.tree.variables[[tcomp]] <- numeric(nreps)
  importancePred.tree.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                              ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.tree.variables <- gbm.step(data=no.trees,
                                   gbm.x = c(11:18),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree..complexity = tcomp,
                                   learning.rate = 0.0001,
                                   bag.fraction = 0.75,
                                   n.tree.s = 50,
                                   plot.main=F, plot.folds=F,
                                   site.weights = no.trees$site.id)
    #R2 adj:
    R2Obs.tree.variables[[tcomp]][i] <- 1 - (BRT.tree.variables$self.statistics$mean.resid /
                                               BRT.tree.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.tree.variables[[tcomp]]) <- sort(rownames(summary(BRT.tree.variables)))
    }
    importancePred.tree.variables[[tcomp]][, i] <-
      summary(BRT.tree.variables)[rownames(importancePred.tree.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree. complexity, pag 105

means <- sapply(R2Obs.tree.variables, mean)
sds <- sapply(R2Obs.tree.variables, sd)
plot(1:length(R2Obs.tree.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.tree.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# tree variables without site

R2Obs.tree.variables <- list()
importancePred.tree.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.tree.variables[[tcomp]] <- numeric(nreps)
  importancePred.tree.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                              ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.tree.variables <- gbm.step(data=no.trees,
                                   gbm.x = c(11:18),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree..complexity = tcomp,
                                   learning.rate = 0.001,
                                   bag.fraction = 0.75,
                                   n.tree.s = 50,
                                   plot.main=F, plot.folds=F)
                                   
    #R2 adj:
    R2Obs.tree.variables[[tcomp]][i] <- 1 - (BRT.tree.variables$self.statistics$mean.resid /
                                               BRT.tree.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.tree.variables[[tcomp]]) <- sort(rownames(summary(BRT.tree.variables)))
    }
    importancePred.tree.variables[[tcomp]][, i] <-
      summary(BRT.tree.variables)[rownames(importancePred.tree.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree. complexity, pag 105

means <- sapply(R2Obs.tree.variables, mean)
sds <- sapply(R2Obs.tree.variables, sd)
plot(1:length(R2Obs.tree.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.tree.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# annual.data with site

R2Obs.annual.variables <- list()
importancePred.annual.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.annual.variables[[tcomp]] <- numeric(nreps)
  importancePred.annual.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.annual.variables <- gbm.step(data=annual.data,
                                     gbm.x = c(11:18),
                                     gbm.y = 10,
                                     family = "gaussian",
                                     tree.complexity = tcomp,
                                     learning.rate = 0.001,
                                     bag.fraction = 0.75,
                                     n.trees = 50,
                                     plot.main=F, plot.folds=F,
                                     site.weights = annual.data$site.id)
    #R2 adj:
    R2Obs.annual.variables[[tcomp]][i] <- 1 - (BRT.annual.variables$self.statistics$mean.resid /
                                                 BRT.annual.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.annual.variables[[tcomp]]) <- sort(rownames(summary(BRT.annual.variables)))
    }
    importancePred.annual.variables[[tcomp]][, i] <-
      summary(BRT.annual.variables)[rownames(importancePred.annual.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.annual.variables, mean)
sds <- sapply(R2Obs.annual.variables, sd)
plot(1:length(R2Obs.annual.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.annual.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# 9 complexity is best for all variables model with lr = 0.001

# annual.data without site

R2Obs.annual.variables <- list()
importancePred.annual.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.annual.variables[[tcomp]] <- numeric(nreps)
  importancePred.annual.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.annual.variables <- gbm.step(data=annual.data,
                                     gbm.x = c(11:18),
                                     gbm.y = 10,
                                     family = "gaussian",
                                     tree.complexity = tcomp,
                                     learning.rate = 0.0001,
                                     bag.fraction = 0.75,
                                     n.trees = 50,
                                     plot.main=F, plot.folds=F)
                                     
    #R2 adj:
    R2Obs.annual.variables[[tcomp]][i] <- 1 - (BRT.annual.variables$self.statistics$mean.resid /
                                                 BRT.annual.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.annual.variables[[tcomp]]) <- sort(rownames(summary(BRT.annual.variables)))
    }
    importancePred.annual.variables[[tcomp]][, i] <-
      summary(BRT.annual.variables)[rownames(importancePred.annual.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.annual.variables, mean)
sds <- sapply(R2Obs.annual.variables, sd)
plot(1:length(R2Obs.annual.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.annual.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# 6, lr = 0.0001

# perennial variables with site

R2Obs.perennial.variables <- list()
importancePred.perennial.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.perennial.variables[[tcomp]] <- numeric(nreps)
  importancePred.perennial.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                   ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.perennial.variables <- gbm.step(data=perennial.data,
                                        gbm.x = c(11:18),
                                        gbm.y = 10,
                                        family = "gaussian",
                                        tree.complexity = tcomp,
                                        learning.rate = 0.000001,
                                        bag.fraction = 0.75,
                                        n.trees = 50,
                                        plot.main=F, plot.folds=F,
                                        site.weights = perennial.data$site.id)
    #R2 adj:
    R2Obs.perennial.variables[[tcomp]][i] <- 1 - (BRT.perennial.variables$self.statistics$mean.resid /
                                                    BRT.perennial.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.perennial.variables[[tcomp]]) <- sort(rownames(summary(BRT.perennial.variables)))
    }
    importancePred.perennial.variables[[tcomp]][, i] <-
      summary(BRT.perennial.variables)[rownames(importancePred.perennial.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.perennial.variables, mean)
sds <- sapply(R2Obs.perennial.variables, sd)
plot(1:length(R2Obs.perennial.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.perennial.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# perennial variables without site

R2Obs.perennial.variables <- list()
importancePred.perennial.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.perennial.variables[[tcomp]] <- numeric(nreps)
  importancePred.perennial.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                   ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.perennial.variables <- gbm.step(data=perennial.data,
                                        gbm.x = c(11:18),
                                        gbm.y = 10,
                                        family = "gaussian",
                                        tree.complexity = tcomp,
                                        learning.rate = 0.0000001,
                                        bag.fraction = 0.75,
                                        n.trees = 50,
                                        plot.main=F, plot.folds=F)
                                        
    #R2 adj:
    R2Obs.perennial.variables[[tcomp]][i] <- 1 - (BRT.perennial.variables$self.statistics$mean.resid /
                                                    BRT.perennial.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.perennial.variables[[tcomp]]) <- sort(rownames(summary(BRT.perennial.variables)))
    }
    importancePred.perennial.variables[[tcomp]][, i] <-
      summary(BRT.perennial.variables)[rownames(importancePred.perennial.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.perennial.variables, mean)
sds <- sapply(R2Obs.perennial.variables, sd)
plot(1:length(R2Obs.perennial.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.perennial.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# grass variables with site

R2Obs.grass.variables <- list()
importancePred.grass.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.grass.variables[[tcomp]] <- numeric(nreps)
  importancePred.grass.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                               ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.grass.variables <- gbm.step(data=grass,
                                    gbm.x = c(11:18),
                                    gbm.y = 10,
                                    family = "gaussian",
                                    tree.complexity = tcomp,
                                    learning.rate = 0.0000001,
                                    bag.fraction = 0.75,
                                    n.trees = 50,
                                    plot.main=F, plot.folds=F,
                                    site.weights = grass$site.id)
    #R2 adj:
    R2Obs.grass.variables[[tcomp]][i] <- 1 - (BRT.grass.variables$self.statistics$mean.resid /
                                                BRT.grass.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.grass.variables[[tcomp]]) <- sort(rownames(summary(BRT.grass.variables)))
    }
    importancePred.grass.variables[[tcomp]][, i] <-
      summary(BRT.grass.variables)[rownames(importancePred.grass.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.grass.variables, mean)
sds <- sapply(R2Obs.grass.variables, sd)
plot(1:length(R2Obs.grass.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.grass.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# grass variables without site

R2Obs.grass.variables <- list()
importancePred.grass.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.grass.variables[[tcomp]] <- numeric(nreps)
  importancePred.grass.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                               ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.grass.variables <- gbm.step(data=grass,
                                    gbm.x = c(11:18),
                                    gbm.y = 10,
                                    family = "gaussian",
                                    tree.complexity = tcomp,
                                    learning.rate = 0.000001,
                                    bag.fraction = 0.75,
                                    n.trees = 50,
                                    plot.main=F, plot.folds=F)
                                    
    #R2 adj:
    R2Obs.grass.variables[[tcomp]][i] <- 1 - (BRT.grass.variables$self.statistics$mean.resid /
                                                BRT.grass.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.grass.variables[[tcomp]]) <- sort(rownames(summary(BRT.grass.variables)))
    }
    importancePred.grass.variables[[tcomp]][, i] <-
      summary(BRT.grass.variables)[rownames(importancePred.grass.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.grass.variables, mean)
sds <- sapply(R2Obs.grass.variables, sd)
plot(1:length(R2Obs.grass.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.grass.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# forb variables with site

R2Obs.forb.variables <- list()
importancePred.forb.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.forb.variables[[tcomp]] <- numeric(nreps)
  importancePred.forb.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                              ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.forb.variables <- gbm.step(data=forb,
                                   gbm.x = c(11:18),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = tcomp,
                                   learning.rate = 0.001,
                                   bag.fraction = 0.75,
                                   n.trees = 50,
                                   plot.main=F, plot.folds=F,
                                   site.weights = forb$site.id)
    #R2 adj:
    R2Obs.forb.variables[[tcomp]][i] <- 1 - (BRT.forb.variables$self.statistics$mean.resid /
                                               BRT.forb.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.forb.variables[[tcomp]]) <- sort(rownames(summary(BRT.forb.variables)))
    }
    importancePred.forb.variables[[tcomp]][, i] <-
      summary(BRT.forb.variables)[rownames(importancePred.forb.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.forb.variables, mean)
sds <- sapply(R2Obs.forb.variables, sd)
plot(1:length(R2Obs.forb.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.forb.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# tc of 7 or 9

# forb variables without site

R2Obs.forb.variables <- list()
importancePred.forb.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.forb.variables[[tcomp]] <- numeric(nreps)
  importancePred.forb.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                              ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.forb.variables <- gbm.step(data=forb,
                                   gbm.x = c(11:18),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = tcomp,
                                   learning.rate = 0.001,
                                   bag.fraction = 0.75,
                                   n.trees = 50,
                                   plot.main=F, plot.folds=F)
                                   
    #R2 adj:
    R2Obs.forb.variables[[tcomp]][i] <- 1 - (BRT.forb.variables$self.statistics$mean.resid /
                                               BRT.forb.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.forb.variables[[tcomp]]) <- sort(rownames(summary(BRT.forb.variables)))
    }
    importancePred.forb.variables[[tcomp]][, i] <-
      summary(BRT.forb.variables)[rownames(importancePred.forb.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.forb.variables, mean)
sds <- sapply(R2Obs.forb.variables, sd)
plot(1:length(R2Obs.forb.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.forb.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# tc of 10 or 6 or 8

# annual.grass variables with site

R2Obs.annual.grass.variables <- list()
importancePred.annual.grass.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.annual.grass.variables[[tcomp]] <- numeric(nreps)
  importancePred.annual.grass.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                      ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.annual.grass.variables <- gbm.step(data=annual.grass,
                                           gbm.x = c(11:18),
                                           gbm.y = 10,
                                           family = "gaussian",
                                           tree.complexity = tcomp,
                                           learning.rate = 0.001,
                                           bag.fraction = 0.75,
                                           n.trees = 50,
                                           plot.main=F, plot.folds=F,
                                           site.weights = annual.grass$site.id)
    #R2 adj:
    R2Obs.annual.grass.variables[[tcomp]][i] <- 1 - (BRT.annual.grass.variables$self.statistics$mean.resid /
                                                       BRT.annual.grass.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.annual.grass.variables[[tcomp]]) <- sort(rownames(summary(BRT.annual.grass.variables)))
    }
    importancePred.annual.grass.variables[[tcomp]][, i] <-
      summary(BRT.annual.grass.variables)[rownames(importancePred.annual.grass.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.annual.grass.variables, mean)
sds <- sapply(R2Obs.annual.grass.variables, sd)
plot(1:length(R2Obs.annual.grass.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.annual.grass.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 3

# annual.grass variables without site

R2Obs.annual.grass.variables <- list()
importancePred.annual.grass.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.annual.grass.variables[[tcomp]] <- numeric(nreps)
  importancePred.annual.grass.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                      ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.annual.grass.variables <- gbm.step(data=annual.grass,
                                           gbm.x = c(11:18),
                                           gbm.y = 10,
                                           family = "gaussian",
                                           tree.complexity = tcomp,
                                           learning.rate = 0.0001,
                                           bag.fraction = 0.75,
                                           n.trees = 50,
                                           plot.main=F, plot.folds=F)
                                           
    #R2 adj:
    R2Obs.annual.grass.variables[[tcomp]][i] <- 1 - (BRT.annual.grass.variables$self.statistics$mean.resid /
                                                       BRT.annual.grass.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.annual.grass.variables[[tcomp]]) <- sort(rownames(summary(BRT.annual.grass.variables)))
    }
    importancePred.annual.grass.variables[[tcomp]][, i] <-
      summary(BRT.annual.grass.variables)[rownames(importancePred.annual.grass.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.annual.grass.variables, mean)
sds <- sapply(R2Obs.annual.grass.variables, sd)
plot(1:length(R2Obs.annual.grass.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.annual.grass.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 4

# annual.forb variables with site

R2Obs.annual.forb.variables <- list()
importancePred.annual.forb.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.annual.forb.variables[[tcomp]] <- numeric(nreps)
  importancePred.annual.forb.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                     ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.annual.forb.variables <- gbm.step(data=annual.forb,
                                          gbm.x = c(11:18),
                                          gbm.y = 10,
                                          family = "gaussian",
                                          tree.complexity = tcomp,
                                          learning.rate = 0.0001,
                                          bag.fraction = 0.75,
                                          n.trees = 50,
                                          plot.main=F, plot.folds=F,
                                          site.weights = annual.forb$site.id)
    #R2 adj:
    R2Obs.annual.forb.variables[[tcomp]][i] <- 1 - (BRT.annual.forb.variables$self.statistics$mean.resid /
                                                      BRT.annual.forb.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.annual.forb.variables[[tcomp]]) <- sort(rownames(summary(BRT.annual.forb.variables)))
    }
    importancePred.annual.forb.variables[[tcomp]][, i] <-
      summary(BRT.annual.forb.variables)[rownames(importancePred.annual.forb.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.annual.forb.variables, mean)
sds <- sapply(R2Obs.annual.forb.variables, sd)
plot(1:length(R2Obs.annual.forb.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.annual.forb.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# tc = 9 or 10

# annual.forb variables without site

R2Obs.annual.forb.variables <- list()
importancePred.annual.forb.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.annual.forb.variables[[tcomp]] <- numeric(nreps)
  importancePred.annual.forb.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                     ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.annual.forb.variables <- gbm.step(data=annual.forb,
                                          gbm.x = c(11:18),
                                          gbm.y = 10,
                                          family = "gaussian",
                                          tree.complexity = tcomp,
                                          learning.rate = 0.0001,
                                          bag.fraction = 0.75,
                                          n.trees = 50,
                                          plot.main=F, plot.folds=F)
                                          
    #R2 adj:
    R2Obs.annual.forb.variables[[tcomp]][i] <- 1 - (BRT.annual.forb.variables$self.statistics$mean.resid /
                                                      BRT.annual.forb.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.annual.forb.variables[[tcomp]]) <- sort(rownames(summary(BRT.annual.forb.variables)))
    }
    importancePred.annual.forb.variables[[tcomp]][, i] <-
      summary(BRT.annual.forb.variables)[rownames(importancePred.annual.forb.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.annual.forb.variables, mean)
sds <- sapply(R2Obs.annual.forb.variables, sd)
plot(1:length(R2Obs.annual.forb.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.annual.forb.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# tc = 6

# perennial.grass variables with site

R2Obs.perennial.grass.variables <- list()
importancePred.perennial.grass.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.perennial.grass.variables[[tcomp]] <- numeric(nreps)
  importancePred.perennial.grass.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                         ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.perennial.grass.variables <- gbm.step(data=perennial.grass,
                                              gbm.x = c(11:18),
                                              gbm.y = 10,
                                              family = "gaussian",
                                              tree.complexity = tcomp,
                                              learning.rate = 0.00000001,
                                              bag.fraction = 0.75,
                                              n.trees = 50,
                                              plot.main=F, plot.folds=F,
                                              site.weights = perennial.grass$site.id)
    #R2 adj:
    R2Obs.perennial.grass.variables[[tcomp]][i] <- 1 - (BRT.perennial.grass.variables$self.statistics$mean.resid /
                                                          BRT.perennial.grass.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.perennial.grass.variables[[tcomp]]) <- sort(rownames(summary(BRT.perennial.grass.variables)))
    }
    importancePred.perennial.grass.variables[[tcomp]][, i] <-
      summary(BRT.perennial.grass.variables)[rownames(importancePred.perennial.grass.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.perennial.grass.variables, mean)
sds <- sapply(R2Obs.perennial.grass.variables, sd)
plot(1:length(R2Obs.perennial.grass.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.perennial.grass.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# perennial.grass variables with site

R2Obs.perennial.grass.variables <- list()
importancePred.perennial.grass.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.perennial.grass.variables[[tcomp]] <- numeric(nreps)
  importancePred.perennial.grass.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                         ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.perennial.grass.variables <- gbm.step(data=perennial.grass,
                                              gbm.x = c(11:18),
                                              gbm.y = 10,
                                              family = "gaussian",
                                              tree.complexity = tcomp,
                                              learning.rate = 0.00000001,
                                              bag.fraction = 0.75,
                                              n.trees = 50,
                                              plot.main=F, plot.folds=F)
                                              
    #R2 adj:
    R2Obs.perennial.grass.variables[[tcomp]][i] <- 1 - (BRT.perennial.grass.variables$self.statistics$mean.resid /
                                                          BRT.perennial.grass.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.perennial.grass.variables[[tcomp]]) <- sort(rownames(summary(BRT.perennial.grass.variables)))
    }
    importancePred.perennial.grass.variables[[tcomp]][, i] <-
      summary(BRT.perennial.grass.variables)[rownames(importancePred.perennial.grass.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.perennial.grass.variables, mean)
sds <- sapply(R2Obs.perennial.grass.variables, sd)
plot(1:length(R2Obs.perennial.grass.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.perennial.grass.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# perennial.forb variables with site

R2Obs.perennial.forb.variables <- list()
importancePred.perennial.forb.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.perennial.forb.variables[[tcomp]] <- numeric(nreps)
  importancePred.perennial.forb.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                        ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.perennial.forb.variables <- gbm.step(data=perennial.forb,
                                             gbm.x = c(11:18),
                                             gbm.y = 10,
                                             family = "gaussian",
                                             tree.complexity = tcomp,
                                             learning.rate = 0.000001,
                                             bag.fraction = 0.75,
                                             n.trees = 50,
                                             plot.main=F, plot.folds=F,
                                             site.weights = perennial.forb$site.id)
                                             
    #R2 adj:
    R2Obs.perennial.forb.variables[[tcomp]][i] <- 1 - (BRT.perennial.forb.variables$self.statistics$mean.resid /
                                                         BRT.perennial.forb.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.perennial.forb.variables[[tcomp]]) <- sort(rownames(summary(BRT.perennial.forb.variables)))
    }
    importancePred.perennial.forb.variables[[tcomp]][, i] <-
      summary(BRT.perennial.forb.variables)[rownames(importancePred.perennial.forb.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.perennial.forb.variables, mean)
sds <- sapply(R2Obs.perennial.forb.variables, sd)
plot(1:length(R2Obs.perennial.forb.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.perennial.forb.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

# perennial.forb variables without site

R2Obs.perennial.forb.variables <- list()
importancePred.perennial.forb.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.perennial.forb.variables[[tcomp]] <- numeric(nreps)
  importancePred.perennial.forb.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                        ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.perennial.forb.variables <- gbm.step(data=perennial.forb,
                                             gbm.x = c(11:18),
                                             gbm.y = 10,
                                             family = "gaussian",
                                             tree.complexity = tcomp,
                                             learning.rate = 0.000001,
                                             bag.fraction = 0.75,
                                             n.trees = 50,
                                             plot.main=F, plot.folds=F)
                                            
    
    #R2 adj:
    R2Obs.perennial.forb.variables[[tcomp]][i] <- 1 - (BRT.perennial.forb.variables$self.statistics$mean.resid /
                                                         BRT.perennial.forb.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.perennial.forb.variables[[tcomp]]) <- sort(rownames(summary(BRT.perennial.forb.variables)))
    }
    importancePred.perennial.forb.variables[[tcomp]][, i] <-
      summary(BRT.perennial.forb.variables)[rownames(importancePred.perennial.forb.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.perennial.forb.variables, mean)
sds <- sapply(R2Obs.perennial.forb.variables, sd)
plot(1:length(R2Obs.perennial.forb.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.perennial.forb.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# can't get past tc = 1

#### final model with all data ####

set.seed(2023)
all.final.brt.1 <- gbm.step(data=all.data,
                          gbm.x = c(11:18),
                          gbm.y = 10,
                          family = "gaussian",
                          tree.complexity = 1,
                          learning.rate = 0.00001,
                          bag.fraction = 0.50,
                          n.trees = 50,
                          plot.main=T, plot.folds=T, site.weights = all.data$site.id)
set.seed(2023)
all.final.brt.10 <- gbm.step(data=all.data,
                            gbm.x = c(11:18),
                            gbm.y = 10,
                            family = "gaussian",
                            tree.complexity = 10,
                            learning.rate = 0.00001,
                            bag.fraction = 0.75,
                            n.trees = 50,
                            plot.main=T, plot.folds=T, site.weights = all.data$site.id)

ggPerformance(all.1 = all.final.brt.1,all.10 = all.final.brt.10)

1-(all.final.brt.1$self.statistics$mean.resid/all.final.brt.1$self.statistics$mean.null)
# 0.0004
1-(all.final.brt.10$self.statistics$mean.resid/all.final.brt.10$self.statistics$mean.null)
# 0.002

summary(all.final.brt.10)
ggInfluence(all.final.brt.10)

gbm.plot(all.final.brt.10, common.scale = FALSE)
gbm.plot.fits(all.final.brt.10)

plot.gbm(all.final.brt.10, i.var = c("RTD.g.cm3"))
plot.gbm(all.final.brt.10, i.var = c("root.depth_m"))
plot.gbm(all.final.brt.10, i.var = c("leafN.mg.g"))
plot.gbm(all.final.brt.10, i.var = c("SLA_m2.kg"))
plot.gbm(all.final.brt.10, i.var = c("height.m"))
plot.gbm(all.final.brt.10, i.var = c("rootDiam.mm"))
plot.gbm(all.final.brt.10, i.var = c("SRL_m.g"))
plot.gbm(all.final.brt.10, i.var = c("rootN.mg.g"))

all.predict.brt=predict(all.final.brt.10, n.trees = all.final.brt.10$n.trees)
observed.all = all.data$mean.cover.response
plot(all.predict.brt~observed.all)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.final.brt.10$self.statistics$mean.resid/all.final.brt.10$self.statistics$mean.null)
# 0.002

# investigation of interactions
gbm.interactions(all.final.brt.10)$rank.list
# no interactions


#### final model with all data without site ####
# not done
set.seed(2023)
all.final.brt.no.site <- gbm.step(data=all.data,
                          gbm.x = c(11:18),
                          gbm.y = 10,
                          family = "gaussian",
                          tree.complexity = 6,
                          learning.rate = 0.0001,
                          bag.fraction = 0.75,
                          n.trees = 50,
                          plot.main=T, plot.folds=T)
summary(all.final.brt.no.site)

gbm.plot(all.final.brt.no.site, common.scale = FALSE)
gbm.plot.fits(all.final.brt.no.site)

plot.gbm(all.final.brt.no.site, i.var = c("RTD.g.cm3"))
plot.gbm(all.final.brt.no.site, i.var = c("root.depth_m"))
plot.gbm(all.final.brt.no.site, i.var = c("leafN.mg.g"))
plot.gbm(all.final.brt.no.site, i.var = c("SLA_m2.kg"))
plot.gbm(all.final.brt.no.site, i.var = c("height.m"))
plot.gbm(all.final.brt.no.site, i.var = c("rootDiam.mm"))
plot.gbm(all.final.brt.no.site, i.var = c("SRL_m.g"))
plot.gbm(all.final.brt.no.site, i.var = c("rootN.mg.g"))

all.predict.brt=predict(all.final.brt.no.site, n.trees = all.final.brt.no.site$n.trees)
observed.all = all.data$mean.cover.response
plot(all.predict.brt~observed.all)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.final.brt.no.site$self.statistics$mean.resid/all.final.brt.no.site$self.statistics$mean.null)
# 0.001

# investigation of interactions
gbm.interactions(all.final.brt.no.site)$rank.list

all.brt.simple = gbm.simplify(all.final.brt.no.site)
# keep only leafN and height

#### model without trees ####
set.seed(2023)
tree.final.brt <- gbm.step(data=no.trees,
                           gbm.x = c(11:18),
                           gbm.y = 10,
                           family = "gaussian",
                           tree.complexity = 10,
                           learning.rate = 0.0001,
                           bag.fraction = 0.75,
                           n.trees = 50,
                           plot.main=T, plot.folds=T, site.weights = no.trees$site.id)

summary(tree.final.brt)
ggInfluence(tree.final.brt)

gbm.plot(tree.final.brt, common.scale = FALSE)
gbm.plot.fits(tree.final.brt)

plot.gbm(tree.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(tree.final.brt, i.var = c("root.depth_m"))
plot.gbm(tree.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(tree.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(tree.final.brt, i.var = c("height.m"))
plot.gbm(tree.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(tree.final.brt, i.var = c("SRL_m.g"))
plot.gbm(tree.final.brt, i.var = c("rootN.mg.g"))

tree.predict.brt=predict(tree.final.brt, n.trees = tree.final.brt$n.trees)
observed.tree = no.trees$mean.cover.response
plot(tree.predict.brt~observed.tree)
abline(0, 1, col = 2)

R2.tree.brt = 1-(tree.final.brt$self.statistics$mean.resid/tree.final.brt$self.statistics$mean.null)
# 0.06

# investigation of interactions
gbm.interactions(tree.final.brt)$rank.list

# significant interactions
# RTD x height size = 16.55
# RTD x SLA = 2.29
# diam x RTD size = 2.14

ggInteract_2D(gbm.object = tree.final.brt,x="RTD.g.cm3",y="height.m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = tree.final.brt,x="RTD.g.cm3",y="SLA_m2.kg",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = tree.final.brt,x="rootDiam.mm",y="RTD.g.cm3",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(tree.final.brt,x="RTD.g.cm3",y="height.m")
ggInteract_3D(tree.final.brt,x="RTD.g.cm3",y="SLA_m2.kg")
ggInteract_3D(tree.final.brt,x="rootDiam.mm",y="RTD.g.cm3")



#### model without trees without site ####
# not done
set.seed(2023)
tree.final.brt.no.site <- gbm.step(data=no.trees,
                           gbm.x = c(11:18),
                           gbm.y = 10,
                           family = "gaussian",
                           tree.complexity = 10,
                           learning.rate = 0.0001,
                           bag.fraction = 0.75,
                           n.trees = 50,
                           plot.main=T, plot.folds=T)

summary(tree.final.brt.no.site)

gbm.plot(tree.final.brt.no.site, common.scale = FALSE)
gbm.plot.fits(tree.final.brt.no.site)

plot.gbm(tree.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(tree.final.brt, i.var = c("root.depth_m"))
plot.gbm(tree.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(tree.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(tree.final.brt, i.var = c("height.m"))
plot.gbm(tree.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(tree.final.brt, i.var = c("SRL_m.g"))
plot.gbm(tree.final.brt, i.var = c("rootN.mg.g"))

tree.predict.brt=predict(tree.final.brt.no.site, n.trees = tree.final.brt.no.site$n.trees)
observed.tree = no.trees$mean.cover.response
plot(tree.predict.brt~observed.tree)
abline(0, 1, col = 2)

R2.tree.brt = 1-(tree.final.brt.no.site$self.statistics$mean.resid/tree.final.brt.no.site$self.statistics$mean.null)
# 0.08

# investigation of interactions
gbm.interactions(tree.final.brt.no.site)$rank.list

# significant interactions
# none

#### annual final model ####
set.seed(2023)
annual.final.brt <- gbm.step(data=annual.data,
                             gbm.x = c(11:18),
                             gbm.y = 10,
                             family = "gaussian",
                             tree.complexity = 9,
                             learning.rate = 0.001,
                             bag.fraction = 0.75,
                             n.trees = 50,
                             plot.main=T, plot.folds=T, site.weights = annual.data$site.id)

summary(annual.final.brt)
ggInfluence(annual.final.brt)

gbm.plot(annual.final.brt, common.scale = FALSE)
gbm.plot.fits(annual.final.brt)

plot.gbm(annual.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(annual.final.brt, i.var = c("root.depth_m"))
plot.gbm(annual.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(annual.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(annual.final.brt, i.var = c("height.m"))
plot.gbm(annual.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(annual.final.brt, i.var = c("SRL_m.g"))
plot.gbm(annual.final.brt, i.var = c("rootN.mg.g"))

annual.predict.brt=predict(annual.final.brt, n.trees = annual.final.brt$n.trees)
observed.annaul = annual.data$mean.cover.response
plot(annual.predict.brt~observed.annaul)
abline(0, 1, col = 2)

R2.annual.brt = 1-(annual.final.brt$self.statistics$mean.resid/annual.final.brt$self.statistics$mean.null)
# 0.24

# investigation of interactions
gbm.interactions(annual.final.brt)$rank.list

# significant interactions
# SRL x leafN, size = 44.94
# height x leafN size = 8.43
# SRL x depth, size = 3.44

ggInteract_2D(gbm.object = annual.final.brt,x="SRL_m.g",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = annual.final.brt,x="height.m",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = annual.final.brt,x="SRL_m.g",y="root.depth_m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(annual.final.brt,x="SRL_m.g",y="leafN.mg.g")
ggInteract_3D(annual.final.brt,x="height.m",y="leafN.mg.g")
ggInteract_3D(annual.final.brt,x="SRL_m.g",y="root.depth_m")


#### annual final model without site ####
# not done
set.seed(2023)
annual.final.brt.no.site <- gbm.step(data=annual.data,
                             gbm.x = c(11:18),
                             gbm.y = 10,
                             family = "gaussian",
                             tree.complexity = 6,
                             learning.rate = 0.0001,
                             bag.fraction = 0.75,
                             n.trees = 50,
                             plot.main=T, plot.folds=T)
summary(annual.final.brt.no.site)


gbm.plot(annual.final.brt.no.site, common.scale = FALSE)
gbm.plot.fits(annual.final.brt.no.site)

plot.gbm(annual.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(annual.final.brt, i.var = c("root.depth_m"))
plot.gbm(annual.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(annual.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(annual.final.brt, i.var = c("height.m"))
plot.gbm(annual.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(annual.final.brt, i.var = c("SRL_m.g"))
plot.gbm(annual.final.brt, i.var = c("rootN.mg.g"))

annual.predict.brt=predict(annual.final.brt.no.site, n.trees = annual.final.brt.no.site$n.trees)
observed.annaul = annual.data$mean.cover.response
plot(annual.predict.brt~observed.annaul)
abline(0, 1, col = 2)

R2.annual.brt = 1-(annual.final.brt.no.site$self.statistics$mean.resid/annual.final.brt.no.site$self.statistics$mean.null)
# 0.18

# investigation of interactions
gbm.interactions(annual.final.brt.no.site)$rank.list

# significant interactions
# SRL x leafN, size = 42.91
# rootDiam x SRL size = 3.97
# SRL x depth, size = 2.58

#### perennial final model without site####
# perennial with site wouldn't fit
# perennial without site wouldn't fit
set.seed(2023)
perennial.final.brt <- gbm.step(data=perennial.data,
                                gbm.x = c(11:18),
                                gbm.y = 10,
                                family = "gaussian",
                                tree.complexity = 10,
                                learning.rate = 0.0001,
                                bag.fraction = 0.75,
                                n.trees = 50,
                                plot.main=T, plot.folds=T)
                                
summary(perennial.final.brt)
ggInfluence(perennial.final.brt)

gbm.plot(perennial.final.brt, common.scale = FALSE)
gbm.plot.fits(perennial.final.brt)

plot.gbm(perennial.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(perennial.final.brt, i.var = c("root.depth_m"))
plot.gbm(perennial.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(perennial.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(perennial.final.brt, i.var = c("height.m"))
plot.gbm(perennial.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(perennial.final.brt, i.var = c("SRL_m.g"))
plot.gbm(perennial.final.brt, i.var = c("rootN.mg.g"))

perennial.predict.brt=predict(perennial.final.brt, n.trees = perennial.final.brt$n.trees)
observed.perennial = perennial.data$mean.cover.response
plot(annual.predict.brt~observed.perennial)
abline(0, 1, col = 2)

R2.perennail.brt = 1-(perennial.final.brt$self.statistics$mean.resid/perennial.final.brt$self.statistics$mean.null)
# 0.01

# investigation of interactions
gbm.interactions(perennial.final.brt)$rank.list
# none

ggInteract_2D(gbm.object = perennial.final.brt,x="SRL_m.g",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = perennial.final.brt,x="SLA_m2.kg",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = perennial.final.brt,x="height.m",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(perennial.final.brt,x="SRL_m.g",y="leafN.mg.g")
ggInteract_3D(perennial.final.brt,x="SLA_m2.kg",y="leafN.mg.g")
ggInteract_3D(perennial.final.brt,x="height.m",y="leafN.mg.g")


#### perennial without tree final model ####
# won't fit with site.id
set.seed(2023)
perennial.tree.final.brt <- gbm.step(data=perennial.tree,
                                     gbm.x = c(11:18),
                                     gbm.y = 10,
                                     family = "gaussian",
                                     tree.complexity = 10,
                                     learning.rate = 0.0001,
                                     bag.fraction = 0.75,
                                     n.trees = 50,
                                     plot.main=T, plot.folds=T)
                                     
summary(perennial.tree.final.brt)
ggInfluence(perennial.tree.final.brt)

gbm.plot(perennial.tree.final.brt, common.scale = FALSE)
gbm.plot.fits(perennial.tree.final.brt)

plot.gbm(perennial.tree.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(perennial.tree.final.brt, i.var = c("root.depth_m"))
plot.gbm(perennial.tree.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(perennial.tree.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(perennial.tree.final.brt, i.var = c("height.m"))
plot.gbm(perennial.tree.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(perennial.tree.final.brt, i.var = c("SRL_m.g"))
plot.gbm(perennial.tree.final.brt, i.var = c("rootN.mg.g"))

perennial.tree.predict.brt=predict(perennial.tree.final.brt, n.trees = perennial.tree.final.brt$n.trees)
observed.perennial.tree = perennial.tree$mean.cover.response
plot(perennial.tree.predict.brt~observed.perennial.tree)
abline(0, 1, col = 2)

R2.perennail.brt.tree = 1-(perennial.tree.final.brt$self.statistics$mean.resid/perennial.tree.final.brt$self.statistics$mean.null)
# 0.02

# investigation of interactions
gbm.interactions(perennial.tree.final.brt)$rank.list
# three weak interactions

ggInteract_2D(gbm.object = perennial.tree.final.brt,x="SRL_m.g",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = perennial.tree.final.brt,x="height.m",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = perennial.tree.final.brt,x="SRL_m.g",y="height.m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(perennial.tree.final.brt,x="SRL_m.g",y="leafN.mg.g")
ggInteract_3D(perennial.tree.final.brt,x="height.m",y="leafN.mg.g")
ggInteract_3D(perennial.tree.final.brt,x="SRL_m.g",y="height.m")


#### grass final model with site ####

set.seed(2023)
grass.final.brt <- gbm.step(data=grass,
                            gbm.x = c(11:18),
                            gbm.y = 10,
                            family = "gaussian",
                            tree.complexity = 10,
                            learning.rate = 0.0001,
                            bag.fraction = 0.75,
                            n.trees = 50,
                            plot.main=T, plot.folds=T, 
                            site.weights = grass$site.id)
summary(grass.final.brt)
ggInfluence(grass.final.brt)

gbm.plot(grass.final.brt, common.scale = FALSE)
gbm.plot.fits(grass.final.brt)

plot.gbm(grass.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(grass.final.brt, i.var = c("root.depth_m"))
plot.gbm(grass.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(grass.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(grass.final.brt, i.var = c("height.m"))
plot.gbm(grass.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(grass.final.brt, i.var = c("SRL_m.g"))
plot.gbm(grass.final.brt, i.var = c("rootN.mg.g"))

grass.predict.brt=predict(grass.final.brt, n.trees = grass.final.brt$n.trees)
observed.grass = grass$mean.cover.response
plot(grass.predict.brt~observed.grass)
abline(0, 1, col = 2)

R2.grass.brt = 1-(grass.final.brt$self.statistics$mean.resid/grass.final.brt$self.statistics$mean.null)
# 0.04

# investigation of interactions
gbm.interactions(grass.final.brt)$rank.list
# three weak interactions

ggInteract_2D(gbm.object = grass.final.brt,x="RTD.g.cm3",y="SLA_m2.kg",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = grass.final.brt,x="rootDiam.mm",y="RTD.g.cm3",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = grass.final.brt,x="RTD.g.cm3",y="rootN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(grass.final.brt,x="RTD.g.cm3",y="SLA_m2.kg")
ggInteract_3D(grass.final.brt,x="rootDiam.mm",y="RTD.g.cm3")
ggInteract_3D(grass.final.brt,x="RTD.g.cm3",y="rootN.mg.g")


#### forb final model with site ####

set.seed(2023)
forb.final.brt <- gbm.step(data=forb,
                           gbm.x = c(11:18),
                           gbm.y = 10,
                           family = "gaussian",
                           tree.complexity = 10,
                           learning.rate = 0.001,
                           bag.fraction = 0.75,
                           n.trees = 50,
                           plot.main=T, plot.folds=T, 
                           site.weights = forb$site.id)
summary(forb.final.brt)
ggInfluence(forb.final.brt)

gbm.plot(forb.final.brt, common.scale = FALSE)
gbm.plot.fits(forb.final.brt)

plot.gbm(forb.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(forb.final.brt, i.var = c("root.depth_m"))
plot.gbm(forb.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(forb.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(forb.final.brt, i.var = c("height.m"))
plot.gbm(forb.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(forb.final.brt, i.var = c("SRL_m.g"))
plot.gbm(forb.final.brt, i.var = c("rootN.mg.g"))

forb.predict.brt=predict(forb.final.brt, n.trees = forb.final.brt$n.trees)
observed.forb = forb$mean.cover.response
plot(forb.predict.brt~observed.forb)
abline(0, 1, col = 2)

R2.forb.brt = 1-(forb.final.brt$self.statistics$mean.resid/forb.final.brt$self.statistics$mean.null)
# 0.19

# investigation of interactions
gbm.interactions(forb.final.brt)$rank.list
# depth x SLA 80.92
# height x leafN 38.51
# RTD x height 8.30

ggInteract_2D(gbm.object = forb.final.brt,x="root.depth_m",y="SLA_m2.kg",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = forb.final.brt,x="height.m",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = forb.final.brt,x="RTD.g.cm3",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(forb.final.brt,x="root.depth_m",y="SLA_m2.kg")
ggInteract_3D(forb.final.brt,x="height.m",y="leafN.mg.g")
ggInteract_3D(forb.final.brt,x="RTD.g.cm3",y="leafN.mg.g")


#### annual.grass final model with site ####

set.seed(2023)
annual.grass.final.brt <- gbm.step(data=annual.grass,
                                   gbm.x = c(11:18),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = 3,
                                   learning.rate = 0.001,
                                   bag.fraction = 0.75,
                                   n.trees = 50,
                                   plot.main=T, plot.folds=T, 
                                   site.weights = annual.grass$site.id)
summary(annual.grass.final.brt)
ggInfluence(annual.grass.final.brt)

gbm.plot(annual.grass.final.brt, common.scale = FALSE)
gbm.plot.fits(annual.grass.final.brt)

plot.gbm(annual.grass.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(annual.grass.final.brt, i.var = c("root.depth_m"))
plot.gbm(annual.grass.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(annual.grass.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(annual.grass.final.brt, i.var = c("height.m"))
plot.gbm(annual.grass.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(annual.grass.final.brt, i.var = c("SRL_m.g"))
plot.gbm(annual.grass.final.brt, i.var = c("rootN.mg.g"))

annual.grass.predict.brt=predict(annual.grass.final.brt, n.trees = annual.grass.final.brt$n.trees)
observed.annual.grass = annual.grass$mean.cover.response
plot(annual.grass.predict.brt~observed.annual.grass)
abline(0, 1, col = 2)

R2.annual.grass.brt = 1-(annual.grass.final.brt$self.statistics$mean.resid/annual.grass.final.brt$self.statistics$mean.null)
# 0.17

# investigation of interactions
gbm.interactions(annual.grass.final.brt)$rank.list
# three weak interactions

ggInteract_2D(gbm.object = annual.grass.final.brt,x="SLA_m2.kg",y="height.m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = annual.grass.final.brt,x="height.m",y="leafN.mg.g",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = annual.grass.final.brt,x="root.depth_m",y="height.m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(annual.grass.final.brt,x="SLA_m2.kg",y="height.m")
ggInteract_3D(annual.grass.final.brt,x="height.m",y="leafN.mg.g")
ggInteract_3D(annual.grass.final.brt,x="root.depth_m",y="height.m")


#### annual.forb final model with site ####

set.seed(2023)
annual.forb.final.brt <- gbm.step(data=annual.forb,
                                  gbm.x = c(11:18),
                                  gbm.y = 10,
                                  family = "gaussian",
                                  tree.complexity = 4,
                                  learning.rate = 0.0001,
                                  bag.fraction = 0.75,
                                  n.trees = 50,
                                  plot.main=T, plot.folds=T, 
                                  site.weights = annual.forb$site.id)
summary(annual.forb.final.brt)
ggInfluence(annual.forb.final.brt)

gbm.plot(annual.forb.final.brt, common.scale = FALSE)
gbm.plot.fits(annual.forb.final.brt)

plot.gbm(annual.forb.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(annual.forb.final.brt, i.var = c("root.depth_m"))
plot.gbm(annual.forb.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(annual.forb.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(annual.forb.final.brt, i.var = c("height.m"))
plot.gbm(annual.forb.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(annual.forb.final.brt, i.var = c("SRL_m.g"))
plot.gbm(annual.forb.final.brt, i.var = c("rootN.mg.g"))

annual.forb.predict.brt=predict(annual.forb.final.brt, n.trees = annual.forb.final.brt$n.trees)
observed.annual.forb = annual.forb$mean.cover.response
plot(annual.forb.predict.brt~observed.annual.forb)
abline(0, 1, col = 2)

R2.annual.forb.brt = 1-(annual.forb.final.brt$self.statistics$mean.resid/annual.forb.final.brt$self.statistics$mean.null)
# 0.20

# investigation of interactions
gbm.interactions(annual.forb.final.brt)$rank.list
# SLA x height 3.69
# RTD x height 3.20
# SLA x leafN 1.54 

ggInteract_2D(gbm.object = annual.forb.final.brt,x="SLA_m2.kg",y="height.m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = annual.forb.final.brt,x="RTD.g.cm3",y="height.m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)
ggInteract_2D(gbm.object = annual.forb.final.brt,x="SLA_m2.kg",y="height.m",
              col.gradient = c("white","#5576AF"),show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,col.contour = "#254376",show.axis = T,legend = T)

ggInteract_3D(annual.forb.final.brt,x="SLA_m2.kg",y="height.m")
ggInteract_3D(annual.forb.final.brt,x="RTD.g.cm3",y="leafN.mg.g")
ggInteract_3D(annual.forb.final.brt,x="SLA_m2.kg",y="leafN.mg.g")


#### perennial.grass final model with site ####

set.seed(2023)
perennial.grass.final.brt <- gbm.step(data=perennial.grass,
                                      gbm.x = c(11:18),
                                      gbm.y = 10,
                                      family = "gaussian",
                                      tree.complexity = 10,
                                      learning.rate = 0.00001,
                                      bag.fraction = 0.75,
                                      n.trees = 50,
                                      plot.main=T, plot.folds=T, 
                                      site.weights = perennial.grass$site.id)
summary(perennial.grass.final.brt)
ggInfluence(perennial.grass.final.brt)

gbm.plot(perennial.grass.final.brt, common.scale = FALSE)
gbm.plot.fits(perennial.grass.final.brt)

plot.gbm(perennial.grass.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(perennial.grass.final.brt, i.var = c("root.depth_m"))
plot.gbm(perennial.grass.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(perennial.grass.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(perennial.grass.final.brt, i.var = c("height.m"))
plot.gbm(perennial.grass.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(perennial.grass.final.brt, i.var = c("SRL_m.g"))
plot.gbm(perennial.grass.final.brt, i.var = c("rootN.mg.g"))

perennial.grass.predict.brt=predict(perennial.grass.final.brt, n.trees = perennial.grass.final.brt$n.trees)
observed.perennial.grass = perennial.grass$mean.cover.response
plot(perennial.grass.predict.brt~observed.perennial.grass)
abline(0, 1, col = 2)

R2.perennial.grass.brt = 1-(perennial.grass.final.brt$self.statistics$mean.resid/perennial.grass.final.brt$self.statistics$mean.null)
# 0.004

# investigation of interactions
gbm.interactions(perennial.grass.final.brt)$rank.list
# two weak interactions

#### perennial.forb final model with site ####

set.seed(2023)
perennial.forb.final.brt <- gbm.step(data=perennial.forb,
                                     gbm.x = c(11:18),
                                     gbm.y = 10,
                                     family = "gaussian",
                                     tree.complexity = 10,
                                     learning.rate = 0.0001,
                                     bag.fraction = 0.50,
                                     n.trees = 50,
                                     plot.main=T, plot.folds=T, 
                                     site.weights = perennial.forb$site.id)
summary(perennial.forb.final.brt)
ggInfluence(perennial.forb.final.brt)

gbm.plot(perennial.forb.final.brt, common.scale = FALSE)
gbm.plot.fits(perennial.forb.final.brt)

plot.gbm(perennial.forb.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(perennial.forb.final.brt, i.var = c("root.depth_m"))
plot.gbm(perennial.forb.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(perennial.forb.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(perennial.forb.final.brt, i.var = c("height.m"))
plot.gbm(perennial.forb.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(perennial.forb.final.brt, i.var = c("SRL_m.g"))
plot.gbm(perennial.forb.final.brt, i.var = c("rootN.mg.g"))

perennial.forb.predict.brt=predict(perennial.forb.final.brt, n.trees = perennial.forb.final.brt$n.trees)
observed.perennial.forb = perennial.forb$mean.cover.response
plot(perennial.forb.predict.brt~observed.perennial.forb)
abline(0, 1, col = 2)

R2.perennial.forb.brt = 1-(perennial.forb.final.brt$self.statistics$mean.resid/perennial.forb.final.brt$self.statistics$mean.null)
# 0.03

# investigation of interactions
gbm.interactions(perennial.forb.final.brt)$rank.list
# three weak interactions


