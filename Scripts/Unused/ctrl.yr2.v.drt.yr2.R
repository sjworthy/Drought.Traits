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
annual.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/annual.data.csv", row.names = 1)
annual.data = subset(annual.data, !annual.data$functional_group == "WOODY") # 95 data points
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

# make negative imputed traits equal to NA

impute.traits = impute.traits %>% mutate(height = replace(height, which(height<0), NA)) %>%
  mutate(RTD = replace(RTD, which(RTD<0), NA)) %>%
  mutate(SRL = replace(SRL, which(SRL<0), NA))

no.trees$Taxon = str_to_sentence(no.trees$Taxon)
no.trees.impute = left_join(no.trees,impute.traits)
no.trees.impute = no.trees.impute[,c(1,4,11,20:30)]

annual.data$Taxon = str_to_sentence(annual.data$Taxon)
annual.data.impute = left_join(annual.data,impute.traits)
annual.data.impute = annual.data.impute[,c(1,4,11,20:22,24:31)]

perennial.tree$Taxon = str_to_sentence(perennial.tree$Taxon)
perennial.tree.impute = left_join(perennial.tree,impute.traits)
perennial.tree.impute = perennial.tree.impute[,c(1,4,11,20:30)]

grass$Taxon = str_to_sentence(grass$Taxon)
grass.impute = left_join(grass,impute.traits)
grass.impute = grass.impute[,c(1,4,11,20:30)]

forb$Taxon = str_to_sentence(forb$Taxon)
forb.impute = left_join(forb,impute.traits)
forb.impute = forb.impute[,c(1,4,11,20:30)]

#write.csv(no.trees.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/no.trees.impute.csv")
#write.csv(annual.data.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/annual.data.impute.csv")
#write.csv(perennial.tree.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.tree.impute.csv")
#write.csv(grass.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/grass.impute.csv")
#write.csv(forb.impute, file = "./Formatted.Data/Ctrl.v.drt.yr2.data/forb.impute.csv")

no.trees.impute = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/no.trees.impute.csv") # 507 data points
annual.data.impute = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/annual.data.impute.csv") # 95 data points
perennial.tree.impute = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.tree.impute.csv") # 395 data points
grass.impute = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/grass.impute.csv") # 193 data points
forb.impute = read.csv("./Formatted.Data/Ctrl.v.drt.yr2.data/forb.impute.csv") # 254 data points


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

tree.impute.no.site=gbm.step(data=no.trees.impute, gbm.x = c(8:15,17), gbm.y=4,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.00000005,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(tree.impute.no.site)


set.seed(2023)
annual.brt.1=gbm.step(data=annual.data, gbm.x = c(12:19), gbm.y=11,
                      family = "laplace", tree.complexity = 10, learning.rate = 0.0001,
                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50,
                      site.weights = annual.data$site.id)
ggPerformance(annual.brt.1)

set.seed(2023)

annual.brt.1.no.site=gbm.step(data=annual.data, gbm.x = c(12:19,25), gbm.y=11,
                      family = "gaussian", tree.complexity = 10, learning.rate = 0.00000001,
                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
                      
ggPerformance(annual.brt.1.no.site)

annual.data.impute$scale = scale(annual.data.impute$mean.cover.response)

annual.impute.no.site=gbm.step(data=annual.data.impute, gbm.x = c(8:13,15,17), gbm.y=18,
                             family = "gaussian", tree.complexity = 1, learning.rate = 0.0000005,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(annual.impute.no.site)


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

perennial.impute.no.site=gbm.step(data=perennial.tree.impute, gbm.x = c(8:15,17), gbm.y=4,
                               family = "gaussian", tree.complexity = 10, learning.rate = 0.0000005,
                               bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.impute.no.site)


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

grass.impute.no.site=gbm.step(data=grass.impute, gbm.x = c(8:15,17), gbm.y=4,
                                  family = "gaussian", tree.complexity = 10, learning.rate = 0.000005,
                                  bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.impute.no.site)

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

forb.impute.no.site=gbm.step(data=forb.impute, gbm.x = c(8:15,17), gbm.y=4,
                              family = "gaussian", tree.complexity = 10, learning.rate = 0.000005,
                              bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(forb.impute.no.site)

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




#### determining best tree complexity to use ####

# code edited to run for each of the data sets

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
    
    BRT.all.variables <- gbm.step(data=annual.data.impute,
                                  gbm.x = c(8:15,17),
                                  gbm.y = 4,
                                  family = "gaussian",
                                  tree.complexity = tcomp,
                                  learning.rate = 0.001,
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

tcFactor <- as.factor(rep(1:10, each=nreps))
R2Vector <- unlist(R2Obs.all.variables)
model <- lm(R2Vector~tcFactor)
TukeyModel<-glht(model, linfct = mcp(tcFactor="Tukey"))
TukeyLetters <- cld(TukeyModel)$mcletters$Letters

plot(1:length(R2Obs.all.variables), means)
for (i in 1:length(R2Obs.all.variables)){
  arrows(x0=i, x1=i, y0=means[i]-sds[i], y1=means[i]+sds[i], angle=90, code=3,
         length=0.1)
}
text(x= 1:length(R2Obs.all.variables), y= 0.77, labels=TukeyLetters)

# elected the lowest tc value that did not show significant differences 
# compared to the largest tc value

# annual impute tc = 3


#### annual data impute ####
annual.no.site.map.impute=gbm.step(data=annual.data.impute, gbm.x = c(9:15,17), gbm.y=4,
                                   family = "gaussian", tree.complexity = 3, learning.rate = 0.001,
                                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.no.site.map.impute)
# 1200 trees Per.Expl = 37.38%
1-(annual.no.site.map.impute$self.statistics$mean.resid/annual.no.site.map.impute$self.statistics$mean.null) # R2
ggInfluence(annual.no.site.map.impute)

ggPD(annual.no.site.map.impute, rug = T) # partial dependency plots
ggPDfit(annual.no.site.map.impute)
gbm.plot(annual.no.site.map.impute, common.scale = FALSE)
gbm.plot.fits(annual.no.site.map.impute)

annual.prerun<- plot.gbm.4list(annual.no.site.map.impute)
annual.boot <- gbm.bootstrap.functions(annual.no.site.map.impute, list.predictors=annual.prerun, n.reps=100)
ggPD_boot(annual.no.site.map.impute, predictor="RTD", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map.impute, predictor="root_diam", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map.impute, predictor="SRL", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map.impute, predictor="height", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)

# investigation of interactions
gbm.interactions(annual.no.site.map.impute)$rank.list
ggInteract_list(annual.no.site.map.impute)
# MAP x RTD 39.94
# SRL x RTD 35.58
# diam x height 34.95
# MAP x SRL 24.52

ggInteract_3D(annual.no.site.map.impute, x = 7, y = 1, z.range = c(-2.5, 1.5))
ggInteract_3D(annual.no.site.map.impute, x = 9, y = 1,z.range = c(-1.5, 1.5))
ggInteract_3D(annual.no.site.map.impute, x = 8, y = 7, z.range = c(-1.5, 1))
ggInteract_3D(annual.no.site.map.impute, x = 9, y = 5, z.range = c(-1, 0.8))

save(annual.no.site.map.impute, annual.prerun, annual.boot, file = "./Results/ctrl.v.drt.yr1/annual.impute.RData")


set.seed(2023)
perennial.brt.1.no.site=gbm.step(data=perennial.data, gbm.x = c(11:18,23), gbm.y=10,
                                 family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                                 bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)
# NA



