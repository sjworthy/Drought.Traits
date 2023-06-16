# Script for analysis of Drought Net data
# Comparing mean change in % cover between drought year 1 and drought year 2
# models of year 2 - year 1 cover

library(dismo)
library(gbm)
library(ggBRT)
library(multcomp)

#### read in all the data frames needs for the analyses ####

# all data
all.data = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/all.data.year2.csv", row.names = 1) # 563 data points
# all data without functional group WOODY
no.trees = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/no.trees.csv", row.names = 1) # 485 data points
# all data without local lifeform = TREE
trees = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/trees.csv", row.names = 1) # 559 data points
# all annual data
annual.data = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/annual.data.csv", row.names = 1) # 83 data points
annual.data = subset(annual.data, !annual.data$functional_group == "WOODY") # 82 data points
# all perennial data
perennial.data = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.data.csv", row.names = 1) # 463 data points
# perennial data without functional group WOODY
perennial.tree = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.tree.csv", row.names = 1) # 388 data points
# perennial data without local lifeform = TREE
perennial.no.tree = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.no.tree.csv", row.names = 1) # 459 data points
# grasses
grass = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/grass.csv", row.names = 1) # 175 data points
# forbs
forb = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/forb.csv", row.names = 1) # 252 data points

##### imputed trait dataframes ####

impute.traits = read.csv("./Formatted.Data/trait.imputation.2.csv")
colnames(impute.traits)[1] = "Taxon"
impute.traits$Taxon = sub("_", " ", impute.traits$Taxon)

# make negative imputed traits equal to NA

impute.traits = impute.traits %>% mutate(height = replace(height, which(height<0), NA)) %>%
  mutate(RTD = replace(RTD, which(RTD<0), NA)) %>%
  mutate(SRL = replace(SRL, which(SRL<0), NA))

no.trees$Taxon = str_to_sentence(no.trees$Taxon)
no.trees.impute = left_join(no.trees,impute.traits)
no.trees.impute = no.trees.impute[,c(1,2,5,14:24)]

annual.data$Taxon = str_to_sentence(annual.data$Taxon)
annual.data.impute = left_join(annual.data,impute.traits)
annual.data.impute = annual.data.impute[,c(1,2,5,14:24)]

perennial.tree$Taxon = str_to_sentence(perennial.tree$Taxon)
perennial.tree.impute = left_join(perennial.tree,impute.traits)
perennial.tree.impute = perennial.tree.impute[,c(1,2,5,14:24)]

grass$Taxon = str_to_sentence(grass$Taxon)
grass.impute = left_join(grass,impute.traits)
grass.impute = grass.impute[,c(1,2,5,14:24)]

forb$Taxon = str_to_sentence(forb$Taxon)
forb.impute = left_join(forb,impute.traits)
forb.impute = forb.impute[,c(1,2,5,14:24)]

#write.csv(no.trees.impute, file = "./Formatted.Data/Drt.yr1.v.drt.yr2.data/no.trees.impute.csv")
#write.csv(annual.data.impute, file = "./Formatted.Data/Drt.yr1.v.drt.yr2.data/annual.data.impute.csv")
#write.csv(perennial.tree.impute, file = "./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.tree.impute.csv")
#write.csv(grass.impute, file = "./Formatted.Data/Drt.yr1.v.drt.yr2.data/grass.impute.csv")
#write.csv(forb.impute, file = "./Formatted.Data/Drt.yr1.v.drt.yr2.data/forb.impute.csv")

no.trees.impute = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/no.trees.impute.csv") # 485 data points
annual.data.impute = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/annual.data.impute.csv")
annual.data.impute = subset(annual.data.impute, !annual.data.impute$functional_group == "WOODY") # 82 data points
perennial.tree.impute = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.tree.impute.csv") # 388 data points
grass.impute = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/grass.impute.csv") # 175 data points
forb.impute = read.csv("./Formatted.Data/Drt.yr1.v.drt.yr2.data/forb.impute.csv") # 252 data points


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

no.trees.impute$site.id = as.numeric(as.factor(no.trees.impute$site_code))
annual.data.impute$site.id = as.numeric(as.factor(annual.data.impute$site_code))
perennial.tree.impute$site.id = as.numeric(as.factor(perennial.tree.impute$site_code))
grass.impute$site.id = as.numeric(as.factor(grass.impute$site_code))
forb.impute$site.id = as.numeric(as.factor(forb.impute$site_code))


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

no.trees.impute = merge(no.trees.impute, site.info.map, by="site_code")
annual.data.impute = merge(annual.data.impute, site.info.map, by="site_code")
perennial.tree.impute = merge(perennial.tree.impute, site.info.map, by="site_code")
grass.impute = merge(grass.impute, site.info.map, by="site_code")
forb.impute = merge(forb.impute, site.info.map, by="site_code")


#### determining best parameter combination to generate 1000 trees ####
# bag fraction of 0.50 and 0.75, step.size of 25 and 50, tc = 10

set.seed(2023)
all.brt.1=gbm.step(data=all.data, gbm.x = c(6:13,18), gbm.y=5,
                   family = "gaussian", tree.complexity = 10, learning.rate = 0.00005,
                   bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50,
                   site.weights = all.data$site.id)

ggPerformance(all.brt.1)

set.seed(2023)
all.brt.1.no.site=gbm.step(data=all.data, gbm.x = c(6:13,18), gbm.y=5,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                           bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(all.brt.1.no.site)

set.seed(2023)

tree.brt.1=gbm.step(data=no.trees, gbm.x = c(6:13,18), gbm.y=5,
                    family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25, 
                    site.weights = no.trees$site.id)

ggPerformance(tree.brt.1)

set.seed(2023)

tree.brt.1.no.site=gbm.step(data=no.trees, gbm.x = c(9), gbm.y=5,
                            family = "gaussian", tree.complexity = 1, learning.rate = 0.00000000005,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(tree.brt.1.no.site)


tree.impute.no.site=gbm.step(data=no.trees.impute, gbm.x = c(8:15,17), gbm.y=4,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.00000001,
                             bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
                             
ggPerformance(tree.impute.no.site)


set.seed(2023)

annual.brt.1=gbm.step(data=annual.data, gbm.x = c(6:13,18), gbm.y=5,
                      family = "gaussian", tree.complexity = 10, learning.rate = 0.005,
                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50,
                      site.weights = annual.data$site.id)

ggPerformance(annual.brt.1) 

set.seed(2023)

annual.brt.1.no.site=gbm.step(data=annual.data, gbm.x = c(6:13,18), gbm.y=5,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                              bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(annual.brt.1.no.site) 

annual.impute.no.site=gbm.step(data=annual.data.impute, gbm.x = c(8:15,17), gbm.y=4,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                              bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(annual.impute.no.site) 



set.seed(2023)
perennial.brt.1=gbm.step(data=perennial.data, gbm.x = c(6:13,18), gbm.y=5,
                         family = "gaussian", tree.complexity = 10, learning.rate = 0.00005,
                         bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50, 
                         site.weights = perennial.data$site.id)

ggPerformance(perennial.brt.1) 

set.seed(2023)
perennial.brt.1.no.site=gbm.step(data=perennial.data, gbm.x = c(6:13,18), gbm.y=5,
                                 family = "gaussian", tree.complexity = 10, learning.rate = 0.00005,
                                 bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(perennial.brt.1.no.site) 

set.seed(2023)

perennial.tree.brt.1=gbm.step(data=perennial.tree, gbm.x = c(6:13), gbm.y=5,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                              bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25,
                              site.weights = perennial.tree$site.id)

ggPerformance(perennial.tree.brt.1) 

set.seed(2023)

perennial.tree.brt.1.no.site=gbm.step(data=perennial.tree, gbm.x = c(6:13), gbm.y=5,
                                      family = "gaussian", tree.complexity = 10, learning.rate = 0.000001,
                                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(perennial.tree.brt.1.no.site) 

perennial.tree.impute.no.site=gbm.step(data=perennial.tree.impute, gbm.x = c(8:15,17), gbm.y=4,
                                      family = "gaussian", tree.complexity = 10, learning.rate = 0.0000005,
                                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(perennial.tree.impute.no.site) 


set.seed(2023)
grass.brt.1=gbm.step(data=grass, gbm.x = c(6:13,18), gbm.y=5,
                     family = "gaussian", tree.complexity = 10, learning.rate = 0.00000000000001,
                     bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25,
                     site.weights = grass$site.id)

ggPerformance(grass.brt.1) 

set.seed(2023)

grass.brt.1.no.site=gbm.step(data=grass, gbm.x = c(6:13), gbm.y=5,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.brt.1.no.site) 

grass.impute.no.site=gbm.step(data=grass.impute, gbm.x = c(8:15,17), gbm.y=4,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.impute.no.site) 

set.seed(2023)
forb.brt.1=gbm.step(data=forb, gbm.x = c(6:13,18), gbm.y=5,
                    family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE,step.size = 50,
                    site.weights = forb$site.id)

ggPerformance(forb.brt.1) 

set.seed(2023)

forb.brt.1.no.site=gbm.step(data=forb, gbm.x = c(6:13,18), gbm.y=5,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.brt.1.no.site) 

forb.impute.no.site=gbm.step(data=forb.impute, gbm.x = c(8:15,17), gbm.y=4,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.impute.no.site) 


set.seed(2023)
annual.grass.brt.1=gbm.step(data=annual.grass, gbm.x = c(6:13), gbm.y=5,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                            bag.fraction = 0.25, n.trees = 50, verbose = TRUE, step.size = 25,
                            site.weights = annual.grass$site.id)

ggPerformance(annual.grass.brt.1) 

set.seed(2023)
annual.grass.brt.1.no.site=gbm.step(data=annual.grass, gbm.x = c(6:13), gbm.y=5,
                                    family = "gaussian", tree.complexity = 10, learning.rate = 0.005,
                                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.grass.brt.1.no.site) 


set.seed(2023)
annual.forb.brt.1=gbm.step(data=annual.forb, gbm.x = c(6:13,17), gbm.y=5,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                           bag.fraction = 0.75, n.trees = 50, verbose = TRUE,step.size = 25,
                           site.weights = annual.forb$site.id)
ggPerformance(annual.forb.brt.1) 

set.seed(2023)
annual.forb.brt.1.no.site=gbm.step(data=annual.forb, gbm.x = c(6:13,17), gbm.y=5,
                                   family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.forb.brt.1.no.site) 

set.seed(2023)
perennial.grass.brt.1=gbm.step(data=perennial.grass, gbm.x = c(6:13,17), gbm.y=5,
                               family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                               bag.fraction = 0.75, n.trees = 50, verbose = TRUE,step.size = 50,
                               site.weights = perennial.grass$site.id)
ggPerformance(perennial.grass.brt.1)

set.seed(2023)
perennial.grass.brt.1.no.site=gbm.step(data=perennial.grass, gbm.x = c(6:13,17), gbm.y=5,
                                       family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                       bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.grass.brt.1.no.site)

set.seed(2023)
perennial.forb.brt.1=gbm.step(data=perennial.forb, gbm.x = c(6:13,17), gbm.y=5,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                              bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25,
                              site.weights = perennial.forb$site.id)
ggPerformance(perennial.forb.brt.1)

set.seed(2023)
perennial.forb.brt.1.no.site=gbm.step(data=perennial.forb, gbm.x = c(6:13,17), gbm.y=5,
                                      family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.forb.brt.1.no.site)


set.seed(2023)
trees.brt.1=gbm.step(data=trees, gbm.x = c(6:13,18), gbm.y=5,
                     family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                     bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25, 
                     site.weights = trees$site.id)

ggPerformance(trees.brt.1)

set.seed(2023)
trees.brt.1.no.site=gbm.step(data=trees, gbm.x = c(6:13,18), gbm.y=5,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50) 

ggPerformance(trees.brt.1.no.site)

set.seed(2023)
perennial.no.tree.brt.1=gbm.step(data=perennial.no.tree, gbm.x = c(6:13,18), gbm.y=5,
                                 family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                                 bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25,
                                 site.weights = perennial.no.tree$site.id)
ggPerformance(perennial.no.tree.brt.1)

set.seed(2023)
perennial.no.tree.brt.1.no.site=gbm.step(data=perennial.no.tree, gbm.x = c(6:13,18), gbm.y=5,
                                         family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                                         bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(perennial.no.tree.brt.1.no.site)


#### determining best tree complexity ####

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:9),
                                                             ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.all.variables <- gbm.step(data=forb.impute,
                                  gbm.x = c(8:15,17),
                                  gbm.y = 4,
                                  family = "gaussian",
                                  tree.complexity = tcomp,
                                  learning.rate = 0.0001,
                                  bag.fraction = 0.75,
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

# examine how R2 improves with tree complexity

means <- sapply(R2Obs.all.variables, mean)
sds <- sapply(R2Obs.all.variables, sd)
plot(1:length(R2Obs.all.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.all.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

tcFactor <- as.factor(rep(1:8, each=nreps))
R2Vector <- unlist(R2Obs.all.variables)
model <- lm(R2Vector~tcFactor)
TukeyModel<-glht(model, linfct = mcp(tcFactor="Tukey"))
TukeyLetters <- cld(TukeyModel)$mcletters$Letters

plot(1:length(R2Obs.all.variables), means)
for (i in 1:length(R2Obs.all.variables)){
  arrows(x0=i, x1=i, y0=means[i]-sds[i], y1=means[i]+sds[i], angle=90, code=3,
         length=0.1)
}
text(x= 1:length(R2Obs.all.variables),labels=TukeyLetters)

# elected the lowest tc value that did not show significant differences 
# compared to the largest tc value
# no trees without map tc = 1, won't fit past tc = 1
# annual tc = 1,  won't fit past tc = 1
# perennial tree tc = NA
# grass tc = NA
# forb tc = 3

# annual.impute tc = 1
# forb impute tc = 3

#### Best Models ####

set.seed(2023)
all.brt.no.site.map=gbm.step(data=all.data, gbm.x = c(6:13,18), gbm.y=5,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)


ggPerformance(all.brt.no.site.map)
# 1100 trees Per.Expl = 11.19586552
1-(all.brt.no.site.map$self.statistics$mean.resid/all.brt.no.site.map$self.statistics$mean.null) # R2
ggInfluence(all.brt.no.site.map)

ggPD(all.brt.no.site.map, rug = T) # partial dependency plots
gbm.plot(all.brt.no.site.map, common.scale = FALSE)
gbm.plot.fits(all.brt.no.site.map)

# investigation of interactions
gbm.interactions(all.brt.no.site.map)$rank.list
ggInteract_list(all.brt.no.site.map)
# diam x leafN 18.98
# rootN x leafN 14.12
# MAP x leafN 7.58
# SLA x rootN 2.43

ggInteract_3D(all.brt.no.site.map, x = 6, y = 2, z.range = c(-1.5, 0.45))
ggInteract_3D(all.brt.no.site.map, x = 6, y = 1,z.range = c(-1, 1.1))
ggInteract_3D(all.brt.no.site.map, x = 7, y = 1, z.range = c(-0.5, 1.5))
ggInteract_3D(all.brt.no.site.map, x = 6, y = 4, z.range = c(0, 1.2))


#### all data without woody ####
tree.no.site.map=gbm.step(data=no.trees, gbm.x = c(6:13,18), gbm.y=5,
                          family = "gaussian", tree.complexity = 1, learning.rate = 0.00001,
                          bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(tree.no.site.map)
# 1500 trees Per.Expl = 14.41%
1-(tree.no.site.map$self.statistics$mean.resid/tree.no.site.map$self.statistics$mean.null) # R2
ggInfluence(tree.no.site.map)

ggPD(tree.no.site.map, rug = T) # partial dependency plots
gbm.plot(tree.no.site.map, common.scale = FALSE)
gbm.plot.fits(tree.no.site.map)

# investigation of interactions
gbm.interactions(tree.no.site.map)$rank.list
ggInteract_list(tree.no.site.map)
# height x leafN 35.06
# RTD x height 27.94
# rootDiam x height 6.23
# SRL x height 5.45

ggInteract_3D(tree.no.site.map, x = 2, y = 1, z.range = c(-2.5, 1.2))
ggInteract_3D(tree.no.site.map, x = 6, y = 2,z.range = c(-2.5, 1.75))
ggInteract_3D(tree.no.site.map, x = 8, y = 2, z.range = c(0, 1.2))
ggInteract_3D(tree.no.site.map, x = 7, y = 2, z.range = c(-1, 1.2))

#### annual ####
annual.no.site.map=gbm.step(data=annual.data, gbm.x = c(6:13,18), gbm.y=5,
                            family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.no.site.map)
# 1550 trees Per.Expl = 23.06%
1-(annual.no.site.map$self.statistics$mean.resid/annual.no.site.map$self.statistics$mean.null) # R2
ggInfluence(annual.no.site.map)

ggPD(annual.no.site.map, rug = T) # partial dependency plots
ggPDfit(annual.no.site.map)
gbm.plot(annual.no.site.map, common.scale = FALSE)
gbm.plot.fits(annual.no.site.map)

annual.prerun<- plot.gbm.4list(annual.no.site.map)
annual.boot <- gbm.bootstrap.functions(annual.no.site.map, list.predictors=annual.prerun, n.reps=100)
ggPD_boot(annual.no.site.map, predictor="leafN.mg.g", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map, predictor="root.depth_m", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map, predictor="rootN", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)

save(annual.no.site.map, annual.prerun, annual.boot, file = "./Results/drtyr1.v.drtyr2/annual.output.RData")



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
# 2000 trees Per.Expl = 4.92%
1-(perennial.tree.no.site.map$self.statistics$mean.resid/perennial.tree.no.site.map$self.statistics$mean.null) # R2
ggInfluence(perennial.tree.no.site.map)

ggPD(perennial.tree.no.site.map, rug = T) # partial dependency plots
gbm.plot(perennial.tree.no.site.map, common.scale = FALSE)
gbm.plot.fits(perennial.tree.no.site.map)

# investigation of interactions
gbm.interactions(perennial.tree.no.site.map)$rank.list
ggInteract_list(perennial.tree.no.site.map)
# height x leafN 4.33
# SRL x height 1.63
# SRL x leafN 3.30
# SRL x SLA 0.33

ggInteract_3D(perennial.tree.no.site.map, x = 2, y = 1, z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 7, y = 1,z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 7, y = 2, z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 9, y = 4, z.range = c(0.2, 0.75))

set.seed(2023)
grass.no.site.map=gbm.step(data=grass, gbm.x = c(11:18,23), gbm.y=10,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                           bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

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

#### annual impute ####

annual.no.site.map.impute=gbm.step(data=annual.data.impute, gbm.x = c(8:15,17), gbm.y=4,
                            family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.no.site.map.impute)
# 2950 trees Per.Expl = 30.21 %
1-(annual.no.site.map.impute$self.statistics$mean.resid/annual.no.site.map.impute$self.statistics$mean.null) # R2
ggInfluence(annual.no.site.map.impute)

ggPD(annual.no.site.map.impute, rug = T) # partial dependency plots
ggPDfit(annual.no.site.map.impute)
gbm.plot(annual.no.site.map.impute, common.scale = FALSE)
gbm.plot.fits(annual.no.site.map.impute)

annual.prerun<- plot.gbm.4list(annual.no.site.map.impute)
annual.boot <- gbm.bootstrap.functions(annual.no.site.map.impute, list.predictors=annual.prerun, n.reps=100)
ggPD_boot(annual.no.site.map.impute, predictor="rootN", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map.impute, predictor="leafN", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map.impute, predictor="height", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map.impute, predictor="SRL", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map.impute, predictor="root_diam", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)


save(annual.no.site.map.impute, annual.prerun, annual.boot, file = "./Results/drtyr1.v.drtyr2/annual.impute.RData")


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
# 2000 trees Per.Expl = 4.92%
1-(perennial.tree.no.site.map$self.statistics$mean.resid/perennial.tree.no.site.map$self.statistics$mean.null) # R2
ggInfluence(perennial.tree.no.site.map)

ggPD(perennial.tree.no.site.map, rug = T) # partial dependency plots
gbm.plot(perennial.tree.no.site.map, common.scale = FALSE)
gbm.plot.fits(perennial.tree.no.site.map)

# investigation of interactions
gbm.interactions(perennial.tree.no.site.map)$rank.list
ggInteract_list(perennial.tree.no.site.map)
# height x leafN 4.33
# SRL x height 1.63
# SRL x leafN 3.30
# SRL x SLA 0.33

ggInteract_3D(perennial.tree.no.site.map, x = 2, y = 1, z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 7, y = 1,z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 7, y = 2, z.range = c(0, 0.75))
ggInteract_3D(perennial.tree.no.site.map, x = 9, y = 4, z.range = c(0.2, 0.75))

set.seed(2023)
grass.no.site.map=gbm.step(data=grass, gbm.x = c(11:18,23), gbm.y=10,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                           bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

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

#### Forb ####
forb.no.site.map=gbm.step(data=forb, gbm.x = c(6:13,18), gbm.y=5,
                          family = "gaussian", tree.complexity = 3, learning.rate = 0.001,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.no.site.map)
# 1050 trees Per.Expl = 13.67%
1-(forb.no.site.map$self.statistics$mean.resid/forb.no.site.map$self.statistics$mean.null) # R2
ggInfluence(forb.no.site.map)

ggPD(forb.no.site.map, rug = T) # partial dependency plots
ggPDfit(forb.no.site.map)
gbm.plot(forb.no.site.map, common.scale = FALSE)
gbm.plot.fits(forb.no.site.map)

forb.prerun<- plot.gbm.4list(forb.no.site.map)
forb.boot <- gbm.bootstrap.functions(forb.no.site.map, list.predictors=forb.prerun, n.reps=100)
ggPD_boot(forb.no.site.map, predictor="rootDiam.mm", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.no.site.map, predictor="rootN.mg.g", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.no.site.map, predictor="height.m", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.no.site.map, predictor="leafN.mg.g", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)

# investigation of interactions
gbm.interactions(forb.no.site.map)$rank.list
ggInteract_list(forb.no.site.map)
# MAP x rootN 1.70
# depth x rootN 1.00
# MAP x height 0.93
# Diam x rootN 0.84

ggInteract_3D(forb.no.site.map, x = 6, y = 1, z.range = c(-2.0, 7))
ggInteract_3D(forb.no.site.map, x = 5, y = 4,z.range = c(3, 5))
ggInteract_3D(forb.no.site.map, x = 9, y = 1, z.range = c(-0.5, 7))
ggInteract_3D(forb.no.site.map, x = 2, y = 1, z.range = c(0, 7))

save(forb.no.site.map, forb.prerun, forb.boot, file = "./Results/drtyr1.v.drtyr2/forb.output.RData")


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



#### Forb Impute ####
forb.no.site.map.impute=gbm.step(data=forb.impute, gbm.x = c(8:15,17), gbm.y=4,
                          family = "gaussian", tree.complexity = 4, learning.rate = 0.0001,
                          bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.no.site.map.impute)
# 3750 trees Per.Expl = 8.69%
1-(forb.no.site.map.impute$self.statistics$mean.resid/forb.no.site.map.impute$self.statistics$mean.null) # R2
ggInfluence(forb.no.site.map.impute)

ggPD(forb.no.site.map.impute, rug = T) # partial dependency plots
ggPDfit(forb.no.site.map.impute)
gbm.plot(forb.no.site.map.impute, common.scale = FALSE)
gbm.plot.fits(forb.no.site.map.impute)

forb.prerun<- plot.gbm.4list(forb.no.site.map.impute)
forb.boot <- gbm.bootstrap.functions(forb.no.site.map.impute, list.predictors=forb.prerun, n.reps=100)
ggPD_boot(forb.no.site.map.impute, predictor="precip", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.no.site.map.impute, predictor="root_diam", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.no.site.map.impute, predictor="height", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)

# investigation of interactions
gbm.interactions(forb.no.site.map.impute)$rank.list
ggInteract_list(forb.no.site.map.impute)
# MAP x height 15.07
# RTD x leafN 0.71
# RTD x height 0.52
# Diam x SLA 0.21

ggInteract_3D(forb.no.site.map.impute, x = 6, y = 1, z.range = c(-2.0, 7))
ggInteract_3D(forb.no.site.map.impute, x = 5, y = 4,z.range = c(3, 5))
ggInteract_3D(forb.no.site.map.impute, x = 9, y = 1, z.range = c(-0.5, 7))
ggInteract_3D(forb.no.site.map.impute, x = 2, y = 1, z.range = c(0, 7))

save(forb.no.site.map.impute, forb.prerun, forb.boot, file = "./Results/drtyr1.v.drtyr2/forb.impute.RData")


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


