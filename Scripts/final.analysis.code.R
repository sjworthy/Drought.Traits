# Script for analysis of Drought Net data using BRTs
# Comparing mean change in % cover between control and drought in Year 1 

# https://github.com/JBjouffray/Hawaii_RegimesPredictors for visuals and plotting

library(dismo)
library(gbm)
library(ggBRT)
library(multcomp)
library(cowplot)

#### read in all the data frames needs for the analyses ####

# all data
all.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/all.data.csv", row.names = 1) # 632 data points
# annual data
annual.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/annual.data.csv", row.names = 1) # 127 data points
# all perennial data
perennial.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.data.csv", row.names = 1) # 490 data points
# grass
grass = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/grass.csv", row.names = 1) # 230 data points
# forbs
forb = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/forb.csv", row.names = 1) # 324 data points

#### change site code to numeric, continuous vector ####
all.data$site.id = as.numeric(as.factor(all.data$site_code))
annual.data$site.id = as.numeric(as.factor(annual.data$site_code))
perennial.data$site.id = as.numeric(as.factor(perennial.data$site_code))
grass$site.id = as.numeric(as.factor(grass$site_code))
forb$site.id = as.numeric(as.factor(forb$site_code))

#### merge with environmental data ####
# mean annual precipitation data (MAP)
# aridity (arid)

Site.info = read.csv("./Raw.Data/Site_Elev-Disturb.csv")
site.info.map = Site.info[,c(2,13,42)]

all.data = merge(all.data, site.info.map, by="site_code")
annual.data = merge(annual.data, site.info.map, by="site_code")
perennial.data = merge(perennial.data, site.info.map, by="site_code")
grass = merge(grass, site.info.map, by="site_code")
forb = merge(forb, site.info.map, by="site_code")

#### determining best parameter combination to generate 1000 trees ####
# bag fraction of 0.50 and 0.75, step.size of 25 and 50, tc = 10, 

all.data.map=gbm.step(data=all.data, gbm.x = c(12:19,24), gbm.y=10,
                      family = "gaussian", tree.complexity = 9, learning.rate = 0.0005,
                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(all.data.map)

all.data.arid=gbm.step(data=all.data, gbm.x = c(12:19,25), gbm.y=10,
                       family = "gaussian", tree.complexity = 9, learning.rate = 0.0005,
                       bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(all.data.arid)

all.data.map.arid=gbm.step(data=all.data, gbm.x = c(12:19,24,25), gbm.y=10,
                           family = "gaussian", tree.complexity = 9, learning.rate = 0.001,
                           bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(all.data.map.arid)

annual.data.map=gbm.step(data=annual.data, gbm.x = c(12:19,24), gbm.y=10,
                         family = "gaussian", tree.complexity = 9, learning.rate = 0.001,
                         bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(annual.data.map)

annual.data.arid=gbm.step(data=annual.data, gbm.x = c(12:19,25), gbm.y=10,
                          family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(annual.data.arid)

annual.data.map.arid=gbm.step(data=annual.data, gbm.x = c(12:19,24,25), gbm.y=10,
                              family = "gaussian", tree.complexity = 9, learning.rate = 0.001,
                              bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.data.map.arid)


perennial.data.map=gbm.step(data=perennial.data, gbm.x = c(12:19,24), gbm.y=10,
                            family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.data.map)

perennial.data.arid=gbm.step(data=perennial.data, gbm.x = c(12:19,25), gbm.y=10,
                             family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.data.arid)

perennial.data.map.arid=gbm.step(data=perennial.data, gbm.x = c(12:19,24,25), gbm.y=10,
                                 family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                                 bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.data.map.arid)

grass.data.map=gbm.step(data=grass, gbm.x = c(12:19,24), gbm.y=10,
                     family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                     bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.data.map)

grass.data.arid=gbm.step(data=grass, gbm.x = c(12:19,25), gbm.y=10,
                        family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                        bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.data.arid)

grass.data.map.arid=gbm.step(data=grass, gbm.x = c(12:19,24,25), gbm.y=10,
                         family = "gaussian", tree.complexity = 9, learning.rate = 0.00005,
                         bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.data.map.arid)

forb.data.map=gbm.step(data=forb, gbm.x = c(12:19,24), gbm.y=10,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.data.map)

forb.data.arid=gbm.step(data=forb, gbm.x = c(12:19,25), gbm.y=10,
                       family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                       bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.data.arid)

forb.data.map.arid=gbm.step(data=forb, gbm.x = c(12:19,24,25), gbm.y=10,
                        family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                        bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.data.map.arid)

#### determining best tree complexity to use ####

# code edited to run for each of the data sets

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 6:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:10),
                                                             ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    
    BRT.all.variables <- gbm.step(data=grass,
                                  gbm.x = c(12:19,24,25),
                                  gbm.y = 10,
                                  family = "gaussian",
                                  tree.complexity = tcomp,
                                  learning.rate = 0.00001,
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

tcFactor <- as.factor(rep(6:9, each=nreps))
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
# all.data tc = 1
# annual tc = 1
# perennial tree tc = 7
# grass tc = 1, won't fit past tc 1
# forb tc = 4

#### all data without woody ####
load("./Results/ctrl.v.drt.yr1/all.data.no.woody.output.RData") 

all.data.map=gbm.step(data=all.data, gbm.x = c(12:19,24), gbm.y=10,
                          family = "gaussian", tree.complexity = 6, learning.rate = 0.0005,
                          bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
all.data.arid=gbm.step(data=all.data, gbm.x = c(12:19,25), gbm.y=10,
                      family = "gaussian", tree.complexity = 4, learning.rate = 0.0005,
                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
all.data.map.arid=gbm.step(data=all.data, gbm.x = c(12:19,24,25), gbm.y=10,
                       family = "gaussian", tree.complexity = 6, learning.rate = 0.001,
                       bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(all.data.map)
# 1950 trees Per.Expl = 12.02%
ggPerformance(all.data.arid)
# 1325 trees Per.Expl = 6.48%
ggPerformance(all.data.map.arid)
# 900 trees Per.Expl = 9.9%

1-(all.data.map$self.statistics$mean.resid/all.data.map$self.statistics$mean.null) # R2
all.data.map$contributions$var = c("Height","SLA","LeafN","RTD","SRL","Rooting Depth","Precipitation","Root Diameter", "RootN")
all.data.influence.plot = ggInfluence(all.data.map, main = "All (n = 632)", 
                                      col.bar = "gray70", col.signif = "#B50200")

# saved as pdf from device (3.71 x 4.43)
ggInfluence(all.data.arid)
ggInfluence(all.data.map.arid)

ggPD(all.data.map, rug = T, common.scale = FALSE) # partial dependency plots
ggPD(all.data.map, predictor = "height.m",rug = T, smooth = TRUE)
ggPDfit(all.data.map)
gbm.plot(all.data.map,common.scale = FALSE, rug = TRUE)
gbm.plot.fits(all.data.map)
gg_partial_plot(all.data.map,
                var = c("height.m","SLA_m2.kg",
                        "leafN.mg.g"))
plotmo(tree.no.site.map, pmethod="partdep", all1=TRUE, all2=TRUE)

# get data to plot partial dependecy plots
plot(tree.no.site.map, "height.m",return.grid = TRUE)
plot(tree.no.site.map, "SLA_m2.kg", return.grid = TRUE)
plot(tree.no.site.map, "leafN.mg.g", return.grid = TRUE)

all.data.map.prerun<- plot.gbm.4list(all.data.map)
all.data.map$contributions$var = c("height.m","SLA_m2.kg","leafN.mg.g","RTD.g.cm3","SRL_m.g","root.depth_m","precip","rootDiam.mm", "rootN.mg.g")
all.data.map.boot <- gbm.bootstrap.functions(all.data.map, list.predictors=all.data.map.prerun, n.reps=1000)

all.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(all.data.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.map$gbm.call$dataframe)[24] = "Precipitation"

all.pd.plots = ggPD_boot(all.data.map,n.plots = 2,nrow=1,list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

all.pd.plot.all = ggPD_boot(all.data.map,list.4.preds=all.data.map.prerun, col.line="#769370",
                            booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
                            alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")


# output in landscape at 6.89 x 4.43 pdf

ggPD_boot(all.data.map, predictor="Height", list.4.preds=all.data.map.prerun, 
          booted.preds=all.data.map.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(all.data.map, predictor="SLA_m2.kg", list.4.preds=all.data.map.prerun, 
          booted.preds=all.data.map.boot$function.preds, type.ci = "ribbon",rug = T)


# investigation of interactions
gbm.interactions(all.data.map)$interactions
ggInteract_list(all.data.map, index = T)
# RTD x height 9.13
# RTD x SLA 7.65
# SLA x height 7.27
# height x leafN 4.71

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_brt2.all<-ggInteract_boot(c('RTD.g.cm3','height.m'),c('RTD.g.cm3','SLA_m2.kg'),
                                        c('SLA_m2.kg','height.m'),c('height.m','leafN.mg.g'),
                                         nboots = 100,data=all.data, predictors =  c(12:19,24), 
                                         response="mean.cover.response",
                                         family = "gaussian", tc = 6, lr = 0.0005, bf= 0.75,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 2,obs = 9.13) # significant
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 3,obs = 7.65) # non-significant
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 4,obs = 7.27) # non-significant
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 5,obs = 4.71) # non-significant

ggInteract_3D(all.data.map, x = 6, y = 2, z.range = c(-2.0, 0.35))
ggInteract_3D(all.data.map, x = 6, y = 4, z.range = c(-2.0, 0.30))
ggInteract_3D(all.data.map, x = 4, y = 2, z.range = c(-1.5, 0.75))
ggInteract_3D(all.data.map, x = 2, y = 1, z.range = c(-2, 0.85))

all.data.map$contributions$var = c("height.m","SLA_m2.kg","leafN.mg.g","RTD.g.cm3","SRL_m.g","root.depth_m","precip","rootDiam.mm", "rootN.mg.g")
all.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(all.data.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(all.data.map$gbm.call$dataframe)[24] = "precip"


ggInteract_2D(gbm.object = all.data.map,x="RTD.g.cm3",y="height.m",col.gradient = c("white","#769370"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "Height",
              z.range = c(-3.0, 0.35))

# output in landscape at 6.89 x 4.43 pdf

save(all.data.map, all.data.map.prerun, all.data.map.boot, file = "./Results/ctrl.v.drt.yr1/all.data.output.RData")

#### annual data ####
load("./Results/ctrl.v.drt.yr1/annual.output.RData") 

annual.no.site.map=gbm.step(data=annual.data, gbm.x = c(12:19,24), gbm.y=10,
                            family = "gaussian", tree.complexity = 3, learning.rate = 0.001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
annual.no.site.arid=gbm.step(data=annual.data, gbm.x = c(12:19,25), gbm.y=10,
                            family = "gaussian", tree.complexity = 4, learning.rate = 0.001,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
annual.no.site.map.arid=gbm.step(data=annual.data, gbm.x = c(12:19,24,25), gbm.y=10,
                             family = "gaussian", tree.complexity = 4, learning.rate = 0.001,
                             bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.no.site.map)
# 1175 trees Per.Expl = 25.48%
ggPerformance(annual.no.site.arid)
# 2075 trees Per.Expl = 36.50%
ggPerformance(annual.no.site.map.arid)
# 2900 trees Per.Expl = 45.57%

1-(annual.no.site.map$self.statistics$mean.resid/annual.no.site.map$self.statistics$mean.null) # R2

annual.no.site.map$contributions$var = c("Rooting Depth","SRL","LeafN","RTD","Precipitation","Height","SLA","Root Diameter", "RootN")
annual.influence.plot = ggInfluence(annual.no.site.map, main = "Annuals (n = 127)", 
                                    col.bar = "gray70", col.signif = "#B50200")

ggInfluence(annual.no.site.map)
ggInfluence(annual.no.site.arid)
ggInfluence(annual.no.site.map.arid)

ggPD(annual.no.site.map, rug = T) # partial dependency plots
ggPDfit(annual.no.site.map)
gbm.plot(annual.no.site.map, common.scale = FALSE)
gbm.plot.fits(annual.no.site.map)

annual.prerun<- plot.gbm.4list(annual.no.site.map)
annual.no.site.map$contributions$var = c("root.depth_m","SRL_m.g","leafN.mg.g","RTD.g.cm3","precip","height.m","SLA_m2.kg","rootDiam.mm","rootN.mg.g")

annual.no.site.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(annual.no.site.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(annual.no.site.map$gbm.call$dataframe)[24] = "precip"

annual.boot <- gbm.bootstrap.functions(annual.no.site.map, list.predictors=annual.prerun, n.reps=1000)

annual.no.site.map$contributions$var = c("Rooting Depth","SRL","LeafN","RTD","Precipitation","Height","SLA","Root Diameter", "RootN")
annual.no.site.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(annual.no.site.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.no.site.map$gbm.call$dataframe)[24] = "Precipitation"

annual.pd.plot = ggPD_boot(annual.no.site.map,n.plots = 4,nrow=1,list.4.preds=annual.prerun, col.line="#BDB2A7",
                         booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#BDB2A7",
                         alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

# output in landscape at 6.89 x 4.43 pdf

annual.pd.plot.all = ggPD_boot(annual.no.site.map,list.4.preds=annual.prerun, col.line="#BDB2A7",
                           booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#BDB2A7",
                           alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(annual.no.site.map, predictor="root.depth_m", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map, predictor="leafN.mg.g", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map, predictor="SRL_m.g", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(annual.no.site.map, predictor="RTD.g.cm3", list.4.preds=annual.prerun, 
          booted.preds=annual.boot$function.preds, type.ci = "ribbon",rug = T)

# investigation of interactions
gbm.interactions(annual.no.site.map)$rank.list
ggInteract_list(annual.no.site.map)
# SRL x leafN 16.43
# MAP x SRL 10.28
# SRL x height 4.29
# MAP x leafN 3.77

ggInteract_3D(annual.no.site.map, x = 7, y = 1, z.range = c(-2, 1.3))
ggInteract_3D(annual.no.site.map, x = 9, y = 7, z.range = c(-2.1, 0.75))
ggInteract_3D(annual.no.site.map, x = 7, y = 2, z.range = c(-1.2, 1.2))
ggInteract_3D(annual.no.site.map, x = 9, y = 1, z.range = c(-1, 0.85))


# bootstrap interactions

annual.no.site.map$contributions$var = c("root.depth_m","SRL_m.g","leafN.mg.g","RTD.g.cm3","precip","height.m","SLA_m2.kg","rootDiam.mm","rootN.mg.g")
annual.no.site.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(annual.no.site.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(annual.no.site.map$gbm.call$dataframe)[24] = "precip"


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_annual<-ggInteract_boot(c('SRL_m.g','leafN.mg.g'),c('precip','SRL_m.g'),
                                    c('SRL_m.g','height.m'),c('precip','leafN.mg.g'),
                                    nboots = 100,data=annual.data, predictors =  c(12:19,24), 
                                    response="mean.cover.response",
                                    family = "gaussian", tc = 3, lr = 0.001, bf= 0.75,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_annual, column = 2,obs = 16.43) # significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 3,obs = 10.28) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 4,obs = 4.29) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 5,obs = 3.77) # non-significant

ggInteract_2D(gbm.object = all.data.map,x="SRL_m.g",y="leafN.mg.g",col.gradient = c("white","#BDB2A7"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SRL", y.label = "LeafN",
              z.range = c(-1.5, 1))

# output in landscape at 6.89 x 4.43 pdf

save(annual.no.site.map, annual.prerun, annual.boot, file = "./Results/ctrl.v.drt.yr1/annual.output.RData")

#### Perennial without woody ####
load("./Results/ctrl.v.drt.yr1/perennial.no.woody.output.RData") 

perennial.data.map=gbm.step(data=perennial.data, gbm.x = c(12:19,24), gbm.y=10,
                                    family = "gaussian", tree.complexity = 8, learning.rate = 0.0001,
                                    bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

# perennial.data.arid=gbm.step(data=perennial.tree, gbm.x = c(12:19,25), gbm.y=10,
#                            family = "gaussian", tree.complexity = 7, learning.rate = 0.0001,
#                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
# perennial.data.map.arid=gbm.step(data=perennial.tree, gbm.x = c(12:19,24,25), gbm.y=10,
#                            family = "gaussian", tree.complexity = 7, learning.rate = 0.0001,
#                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(perennial.data.map)
# 1000 trees Per.Expl = 1.97%
1-(perennial.data.map$self.statistics$mean.resid/perennial.data.map$self.statistics$mean.null) # R2

perennial.data.map$contributions$var = c("Height","SLA","SRL","LeafN","Rooting Depth","RTD","Precipitation","Root Diameter", "RootN")
perennial.influence.plot = ggInfluence(perennial.data.map, main = "Perennials (n = 296)", 
                                    col.bar = "gray70", col.signif = "#B50200")

ggInfluence(perennial.data.map)

ggPD(perennial.data.map, rug = T) # partial dependency plots
ggPDfit(perennial.data.map)
gbm.plot(perennial.data.map, common.scale = FALSE)
gbm.plot.fits(perennial.data.map)

perennial.prerun<- plot.gbm.4list(perennial.data.map)
perennial.boot <- gbm.bootstrap.functions(perennial.data.map, list.predictors=perennial.prerun, n.reps=1000)

perennial.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(perennial.data.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.map$gbm.call$dataframe)[24] = "Precipitation"

perennial.pd.plot = ggPD_boot(perennial.data.map,n.plots = 3,nrow=1,list.4.preds=perennial.prerun, col.line="#F1C646",
                           booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
                           alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

# output in landscape at 6.89 x 4.43 pdf

perennial.pd.plot.all = ggPD_boot(perennial.data.map,list.4.preds=perennial.prerun, col.line="#F1C646",
                              booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
                              alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(perennial.data.map, predictor="height.m", list.4.preds=perennial.prerun, 
          booted.preds=perennial.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(perennial.data.map, predictor="SLA_m2.kg", list.4.preds=perennial.prerun, 
          booted.preds=perennial.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(perennial.data.map, predictor="precip", list.4.preds=perennial.prerun, 
          booted.preds=perennial.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(perennial.data.map, predictor="SRL_m.g", list.4.preds=perennial.prerun, 
          booted.preds=perennial.boot$function.preds, type.ci = "ribbon",rug = T)

# investigation of interactions
gbm.interactions(perennial.data.map)$rank.list
ggInteract_list(perennial.data.map)
# SRL x height 0.65
# height x leafN 0.48
# SLA x height 0.41
# SRL x leafN 0.14

ggInteract_3D(perennial.data.map, x = 7, y = 2, z.range = c(-0.5, 0.5))
ggInteract_3D(perennial.data.map, x = 2, y = 1, z.range = c(-0.5, 0))
ggInteract_3D(perennial.data.map, x = 4, y = 2, z.range = c(-0.5, 0.1))
ggInteract_3D(perennial.data.map, x = 7, y = 1, z.range = c(-0.5, 0))

# bootstrap interactions

perennial.data.map$contributions$var = c("height.m","SLA_m2.kg","SRL_m.g","leafN.mg.g","root.depth_m","RTD.g.cm3","precip","rootDiam.mm", "rootN.mg.g")
perennial.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(perennial.data.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map$gbm.call$dataframe)[24] = "precip"


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_perennial<-ggInteract_boot(c('SRL_m.g','height.m'),c('height.m','leafN.mg.g'),
                                      c('SLA_m2.kg','height.m'),c('SRL_m.g','leafN.mg.g'),
                                      nboots = 100,data=perennial.data, predictors =  c(12:19,24), 
                                      response="mean.cover.response",
                                      family = "gaussian", tc = 8, lr = 0.0001, bf= 0.75,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_perennial, column = 2,obs = 0.65) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 3,obs = 0.48) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 4,obs = 0.41) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 5,obs = 0.14) # non-significant

save(perennial.data.map, perennial.prerun, perennial.boot, file = "./Results/ctrl.v.drt.yr1/perennial.data.output.RData")

#### Grass ####
load("./Results/ctrl.v.drt.yr1/grass.output.RData") 

grass.map=gbm.step(data=grass, gbm.x = c(12:19,24), gbm.y=10,
                           family = "gaussian", tree.complexity = 2, learning.rate = 0.0001,
                           bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
grass.arid=gbm.step(data=grass, gbm.x = c(12:19,25), gbm.y=10,
                   family = "gaussian", tree.complexity = 2, learning.rate = 0.0001,
                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
grass.map.arid=gbm.step(data=grass, gbm.x = c(12:19,24,25), gbm.y=10,
                   family = "gaussian", tree.complexity = 2, learning.rate = 0.00005,
                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(grass.map)
# 1000 trees Per.Expl = 1.02%
ggPerformance(grass.arid)
# 1000 trees Per.Expl = 1.03%
ggPerformance(grass.map.arid)
# 1000 trees Per.Expl = 0.52%

1-(grass.map$self.statistics$mean.resid/grass.map$self.statistics$mean.null) # R2

grass.map$contributions$var = c("Height","RTD","SLA","LeafN","Root Diameter","Rooting Depth","RootN","SRL", "Precipitation")
grass.influence.plot = ggInfluence(grass.map, main = "Grasses (n = 230)", 
                                       col.bar = "gray70", col.signif = "#B50200")
ggInfluence(grass.map)

ggPD(grass.map, rug = T) # partial dependency plots
ggPDfit(grass.map)
gbm.plot(grass.map, common.scale = FALSE)
gbm.plot.fits(grass.map)

grass.prerun<- plot.gbm.4list(grass.map)
grass.map$contributions$var = c("height.m","RTD.g.cm3","SLA_m2.kg","leafN.mg.g","rootDiam.mm","root.depth_m","rootN.mg.g","SRL_m.g", "precip")
grass.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","Precipitation")
colnames(grass.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.map$gbm.call$dataframe)[24] = "precip"

grass.boot <- gbm.bootstrap.functions(grass.map, list.predictors=grass.prerun, n.reps=1000)

grass.map$contributions$var = c("Height","RTD","SLA","LeafN","Root Diameter","Rooting Depth","RootN","SRL", "Precipitation")
grass.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(grass.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.map$gbm.call$dataframe)[24] = "Precipitation"

grass.pd.plot = ggPD_boot(grass.map,n.plots = 2,nrow=1,list.4.preds=grass.prerun, col.line="#6E687E",
                              booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
                              alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

# output in landscape at 6.89 x 4.43 pdf

grass.pd.plot.all = ggPD_boot(grass.map,list.4.preds=grass.prerun, col.line="#6E687E",
                          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
                          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")



ggPD_boot(grass.map, predictor="height.m", list.4.preds=grass.prerun, 
          booted.preds=grass.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(grass.map, predictor="RTD.g.cm3", list.4.preds=grass.prerun, 
          booted.preds=grass.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(grass.map, predictor="leafN.mg.g", list.4.preds=grass.prerun, 
          booted.preds=grass.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(grass.map, predictor="rootDiam.mm", list.4.preds=grass.prerun, 
          booted.preds=grass.boot$function.preds, type.ci = "ribbon",rug = T)

# investigation of interactions
gbm.interactions(grass.map)$rank.list
ggInteract_list(grass.map)
# RTD x height 1.34
# SLA x height 0.06
# height x leafN 0.05

ggInteract_3D(grass.map, x = 6, y = 2, z.range = c(-0.5, 0))
ggInteract_3D(grass.map, x = 4, y = 2, z.range = c(-0.5, -0.1))
ggInteract_3D(grass.map, x = 2, y = 1, z.range = c(-0.5, -0.1))

# bootstrap interactions

grass.map$contributions$var = c("height.m","RTD.g.cm3","SLA_m2.kg","leafN.mg.g","rootDiam.mm","root.depth_m","rootN.mg.g","SRL_m.g", "precip")
grass.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(grass.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.map$gbm.call$dataframe)[24] = "precip"


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_grass<-ggInteract_boot(c('RTD.g.cm3','height.m'),c('SLA_m2.kg','height.m'),
                                         c('height.m','leafN.mg.g'),
                                         nboots = 100,data=grass, predictors =  c(12:19,24), 
                                         response="mean.cover.response",
                                         family = "gaussian", tc = 2, lr = 0.0001, bf= 0.75,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_grass, column = 2,obs = 1.34) # significant
ggInteract_boot_hist(data = Interact_boot_grass, column = 3,obs = 0.06) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass, column = 4,obs = 0.05) # non-significant

ggInteract_2D(gbm.object = grass.map, x="RTD.g.cm3",y="height.m",col.gradient = c("white","#F1C646"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "Height",
              z.range = c(-0.5, 0.5))

# output in landscape at 6.89 x 4.43 pdf

save(grass.map, grass.prerun, grass.boot, file = "./Results/ctrl.v.drt.yr1/grass.output.RData")

#### Forb ####
load("./Results/ctrl.v.drt.yr1/forb.output.RData") 

forb.map=gbm.step(data=forb, gbm.x = c(12:19,24), gbm.y=10,
                          family = "gaussian", tree.complexity = 5, learning.rate = 0.001,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
forb.arid=gbm.step(data=forb, gbm.x = c(12:19,25), gbm.y=10,
                  family = "gaussian", tree.complexity = 2, learning.rate = 0.001,
                  bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
forb.map.arid=gbm.step(data=forb, gbm.x = c(12:19,24,25), gbm.y=10,
                  family = "gaussian", tree.complexity = 3, learning.rate = 0.001,
                  bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.map)
# 2550 trees Per.Expl = 31.04%
ggPerformance(forb.arid)
# 3200 trees Per.Expl = 25.68%
ggPerformance(forb.map.arid)
# 2450 trees Per.Expl = 29.14%

1-(forb.map$self.statistics$mean.resid/forb.map$self.statistics$mean.null) # R2
forb.map$contributions$var = c("SLA","Height","LeafN","Precipitation","SRL","Rooting Depth","RTD","Root Diameter", "RootN")
forb.influence.plot = ggInfluence(forb.map, main = "Forbs (n = 324)", 
                                   col.bar = "gray70", col.signif = "#B50200")
ggInfluence(forb.map)

ggPD(forb.map, rug = T) # partial dependency plots
ggPDfit(forb.map)
gbm.plot(forb.map, common.scale = FALSE)
gbm.plot.fits(forb.map)

forb.prerun<- plot.gbm.4list(forb.map)
forb.boot <- gbm.bootstrap.functions(forb.map, list.predictors=forb.prerun, n.reps=1000)

forb.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(forb.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.map$gbm.call$dataframe)[24] = "Precipitation"

forb.pd.plot = ggPD_boot(forb.map,n.plots = 4,nrow=1,list.4.preds=forb.prerun, col.line="#F17236",
                          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
                          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

# output in landscape at 6.89 x 4.43 pdf

forb.pd.plot.all = ggPD_boot(forb.map,list.4.preds=forb.prerun, col.line="#F17236",
                         booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
                         alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(forb.map, predictor="height.m", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.map, predictor="RTD.g.cm3", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.map, predictor="leafN.mg.g", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)
ggPD_boot(forb.map, predictor="rootDiam.mm", list.4.preds=forb.prerun, 
          booted.preds=forb.boot$function.preds, type.ci = "ribbon",rug = T)

# investigation of interactions
gbm.interactions(forb.map)$rank.list
ggInteract_list(forb.map)
# depth x SLA 27.82
# SLA x height 25.01
# MAP x SLA 9.29
# SLA x leafN 7.37

ggInteract_3D(forb.map, x = 5, y = 4, z.range = c(-0.5, 3))
ggInteract_3D(forb.map, x = 4, y = 2, z.range = c(-2.2, 2.5))
ggInteract_3D(forb.map, x = 9, y = 4, z.range = c(-0.75, 2))
ggInteract_3D(forb.map, x = 4, y = 1, z.range = c(-1.5, 3.75))

# bootstrap interactions
forb.map$contributions$var = c("SLA_m2.kg","height.m","leafN.mg.g","precip","SRL_m.g","root.depth_m","RTD.g.cm3","rootDiam.mm", "rootN.mg.g")
forb.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(forb.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(forb.map$gbm.call$dataframe)[24] = "precip"


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_forb<-ggInteract_boot(c('root.depth_m','SLA_m2.kg'),c('SLA_m2.kg','height.m'),
                                     c('precip','SLA_m2.kg'),c('SLA_m2.kg','leafN.mg.g'),
                                     nboots = 100,data=forb, predictors =  c(12:19,24), 
                                     response="mean.cover.response",
                                     family = "gaussian", tc = 5, lr = 0.001, bf= 0.50,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_forb, column = 2,obs = 27.82) # significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 3,obs = 25.01) # significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 4,obs = 9.29) # significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 5,obs = 7.37) # non-significant


ggInteract_2D(gbm.object = forb.map, x="root.depth_m",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "Rooting Depth", y.label = "SLA",
              z.range = c(-1.5, 3))
ggInteract_2D(gbm.object = forb.map, x="SLA_m2.kg",y="height.m",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA", y.label = "Height",
              z.range = c(-2.8, 2.5))
ggInteract_2D(gbm.object = forb.map, x="precip",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "SLA",
              z.range = c(-2, 2))

# output in landscape at 6.89 x 4.43 pdf

save(forb.map, forb.prerun, forb.boot, file = "./Results/ctrl.v.drt.yr1/forb.output.RData")



