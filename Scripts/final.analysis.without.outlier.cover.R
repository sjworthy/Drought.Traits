# Redoing analyses removing outlier cover change 
# Use the same data files just removing outlier cover change tested for in All Data

library(dismo)
library(gbm)
library(ggBRT)
library(multcomp)
library(cowplot)
library(corrplot)
library(treezy)

load("final.analysis.wo.outlier.environment.RData")

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

#### Testing for outliers in all data ####

hist(all.data$mean.cover.response)
boxplot(all.data$mean.cover.response)
mean = mean(all.data$mean.cover.response, na.rm = TRUE)
std = sd(all.data$mean.cover.response, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(all.data$mean.cover.response[which(all.data$mean.cover.response <Tmin | all.data$mean.cover.response > Tmax)])
# removed mean.cover.response -27.500 to -19.39497 and 21.33333 to 50.50000

all.data.otl.rm = subset(all.data, all.data$mean.cover.response >-19.39496 & all.data$mean.cover.response < 21.33332)
# lose 16 data points, lose 8 species
annual.data.otl.rm = subset(annual.data, annual.data$mean.cover.response >-19.39496 & annual.data$mean.cover.response < 21.33332)
# lose 3 data points, lose 1 species
perennial.data.otl.rm = subset(perennial.data, perennial.data$mean.cover.response >-19.39496 & perennial.data$mean.cover.response < 21.33332)
# lose 13 data points, lose 7 species
grass.data.otl.rm = subset(grass, grass$mean.cover.response >-19.39496 & grass$mean.cover.response < 21.33332)
# lose 9 data points, lose 5 species
forb.data.otl.rm = subset(forb, forb$mean.cover.response >-19.39496 & forb$mean.cover.response < 21.33332)
# lose 6 data points, lose 2 species

#### change site code to numeric, continuous vector ####
all.data.otl.rm$site.id = as.numeric(as.factor(all.data.otl.rm$site_code))
annual.data.otl.rm$site.id = as.numeric(as.factor(annual.data.otl.rm$site_code))
perennial.data.otl.rm$site.id = as.numeric(as.factor(perennial.data.otl.rm$site_code))
grass.data.otl.rm$site.id = as.numeric(as.factor(grass.data.otl.rm$site_code))
forb.data.otl.rm$site.id = as.numeric(as.factor(forb.data.otl.rm$site_code))

#### merge with environmental data ####
# mean annual precipitation data (MAP)

Site.info = read.csv("./Raw.Data/Site_Elev-Disturb.csv")
site.info.map = Site.info[,c(2,13)]

# override original data names so can use same code
all.data = merge(all.data.otl.rm, site.info.map, by="site_code")
annual.data = merge(annual.data.otl.rm, site.info.map, by="site_code")
perennial.data = merge(perennial.data.otl.rm, site.info.map, by="site_code")
grass = merge(grass.data.otl.rm, site.info.map, by="site_code")
forb = merge(forb.data.otl.rm, site.info.map, by="site_code")



#### testing for correlation among traits ####

# all.data
all.data.2 = all.data
colnames(all.data.2)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.2)[24] = "Precipitation"

cor.mat.all = cor(all.data.2[,c(12:19,24)],use = "pairwise") 
corrplot(cor.mat.all, method="number",tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black",
         title = "All Species")

annual.data.2 = annual.data
colnames(annual.data.2)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.data.2)[24] = "Precipitation"

cor.mat.annual = cor(annual.data.2[,c(12:19,24)],use = "pairwise") 
corrplot(cor.mat.annual, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black",
         title = "Annual Species")

perennial.data.2 = perennial.data
colnames(perennial.data.2)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.2)[24] = "Precipitation"

cor.mat.perennial = cor(perennial.data.2[,c(12:19,24)],use = "pairwise") 
corrplot(cor.mat.perennial, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black",
         title = "Perennial Species")

grass.2 = grass
colnames(grass.2)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.2)[24] = "Precipitation"

cor.mat.grass = cor(grass.2[,c(12:19,24)],use = "pairwise") 
corrplot(cor.mat.grass, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black",
         title = "Grass Species")

forb.2 = forb
colnames(forb.2)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.2)[24] = "Precipitation"

cor.mat.forb = cor(forb.2[,c(12:19,24)],use = "pairwise") 
corrplot(cor.mat.forb, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black",
         title = "Forb Species")

#### determining best parameter combination to generate 1000 trees ####
# bag fraction of 0.50 and 0.75, step.size of 25 and 50, tc = 10, 

all.data.map=gbm.step(data=all.data, gbm.x = c(12:19,24), gbm.y=10,
                      family = "gaussian", tree.complexity = 9, learning.rate = 0.0005,
                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(all.data.map)

annual.data.map=gbm.step(data=annual.data, gbm.x = c(12:19,24), gbm.y=10,
                         family = "gaussian", tree.complexity = 9, learning.rate = 0.001,
                         bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(annual.data.map)

perennial.data.map=gbm.step(data=perennial.data, gbm.x = c(12:19,24), gbm.y=10,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.00005,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.data.map)

grass.data.map=gbm.step(data=grass, gbm.x = c(12:19,24), gbm.y=10,
                        family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                        bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.data.map)

forb.data.map=gbm.step(data=forb, gbm.x = c(12:19,24), gbm.y=10,
                       family = "gaussian", tree.complexity = 9, learning.rate = 0.0005,
                       bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(forb.data.map)

#### determining best tree complexity to use ####

# code edited to run for each of the data sets

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 2:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:9),
                                                             ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    
    BRT.all.variables <- gbm.step(data=forb,
                                  gbm.x = c(12:19,24),
                                  gbm.y = 10,
                                  family = "gaussian",
                                  tree.complexity = tcomp,
                                  learning.rate = 0.0005,
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

tcFactor <- as.factor(rep(2:10, each=nreps))
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

# selected the lowest tc value that did not show significant differences 
# compared to the largest tc value
# all.data tc = 6
# annual tc = 2
# perennial tree tc = 3
# grass tc = 2
# forb tc = 3

#### all data without woody ####
# load("./Results/ctrl.v.drt.yr1/all.data.no.woody.output.RData") 

all.data.map=gbm.step(data=all.data, gbm.x = c(12:19,24), gbm.y=10,
                      family = "gaussian", tree.complexity = 6, learning.rate = 0.0005,
                      bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(all.data.map)
# 1350 trees Per.Expl = 8.23%

1-(all.data.map$self.statistics$mean.resid/all.data.map$self.statistics$mean.null) # R2

all.data.map$contributions$var = c("Height","LeafN","SLA","Root Diameter","RTD","Rooting Depth",
                                   "SRL","Precipitation", "RootN")
all.data.influence.plot = ggInfluence(all.data.map, main = "All Species (n = 616)", 
                                      col.bar = c("#769370","#769370","#769370","#769370",
                                                  "gray70","gray70","gray70","gray70","gray70"), 
                                      col.signif = "#B50200")

# saved as pdf from device (4 x 4)
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

# get data to plot partial dependency plots
all.data.map.prerun<- plot.gbm.4list(all.data.map)
all.data.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","rootDiam.mm",
                                   "RTD.g.cm3","root.depth_m","SRL_m.g","precip", "rootN.mg.g")
all.data.map.boot <- gbm.bootstrap.functions(all.data.map, list.predictors=all.data.map.prerun, n.reps=1000)

all.data.map$contributions$var = c("Height","LeafN","SLA","Root Diameter","RTD","Rooting Depth",
                                   "SRL","Precipitation", "RootN")
all.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(all.data.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(all.data.map,predictor = "Height",list.4.preds=all.data.map.prerun, col.line="#769370",
                         booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
                         alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,
          y.label = "Percent Cover Change")
ggPD_boot(all.data.map,predictor = "LeafN",list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(all.data.map,predictor = "SLA",list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(all.data.map,predictor = "Root Diameter",list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(all.data.map,list.4.preds=all.data.map.prerun, col.line="#769370",
                            booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
                            alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(all.data.map, col.line="#769370", common.scale = FALSE, y.label = "Percent Cover Change")

# output individual plots as 3x3
# output all plots as one 6x6

# investigation of interactions
all.data.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","rootDiam.mm",
                                   "RTD.g.cm3","root.depth_m","SRL_m.g","precip", "rootN.mg.g")
all.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(all.data.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(all.data.map$gbm.call$dataframe)[24] = "precip"

gbm.interactions(all.data.map)$interactions
ggInteract_list(all.data.map, index = T)
# RTD x height 2.63
# rootDiam x leafN 1.49
# SRL x height 1.41
# precip x height 1.29

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_brt2.all<-ggInteract_boot(c('RTD.g.cm3','height.m'),c('rootDiam.mm','leafN.mg.g'),
                                        c('SRL_m.g','height.m'),c('precip','height.m'),
                                        nboots = 100,data=all.data, predictors =  c(12:19,24), 
                                        response="mean.cover.response",
                                        family = "gaussian", tc = 6, lr = 0.0005, bf= 0.50,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 2,obs = 2.63) # significant
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 3,obs = 1.49) # non-significant
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 4,obs = 1.41) # non-significant
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 5,obs = 1.29) # non-significant

ggInteract_3D(all.data.map, x = 6, y = 2, z.range = c(-0.75, 0.30))
ggInteract_3D(all.data.map, x = 6, y = 4, z.range = c(-2.0, 0.30))
ggInteract_3D(all.data.map, x = 4, y = 2, z.range = c(-1.5, 0.75))
ggInteract_3D(all.data.map, x = 2, y = 1, z.range = c(-2, 0.85))

ggInteract_2D(gbm.object = all.data.map, x="RTD.g.cm3",y="height.m",col.gradient = c("white","#769370"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "Height",
              z.range = c(-1, 0.28),z.label = "% Cover Change")

# output in landscape at 4 x 6

# save(all.data.map, all.data.map.prerun, all.data.map.boot, file = "./Results/ctrl.v.drt.yr1/all.data.output.RData")

#### annual data ####
# load("./Results/ctrl.v.drt.yr1/annual.output.RData") 

annual.no.site.map=gbm.step(data=annual.data, gbm.x = c(12:19,24), gbm.y=10,
                            family = "gaussian", tree.complexity = 2, learning.rate = 0.001,
                            bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(annual.no.site.map)
# 1100 trees Per.Expl = 14.43%

1-(annual.no.site.map$self.statistics$mean.resid/annual.no.site.map$self.statistics$mean.null) # R2

annual.no.site.map$contributions$var = c("SLA","SRL","Height","LeafN","Rooting Depth",
                                         "Precipitation","Root Diameter","RTD", "RootN")
annual.influence.plot = ggInfluence(annual.no.site.map, main = "Annuals (n = 124)", 
                                    col.bar = c("#979461","#979461","#979461","#979461",
                                                "gray70","gray70","gray70","gray70","gray70"),
                                    col.signif = "#B50200")

ggInfluence(annual.no.site.map)
ggInfluence(annual.no.site.arid)
ggInfluence(annual.no.site.map.arid)

ggPD(annual.no.site.map, rug = T) # partial dependency plots
ggPDfit(annual.no.site.map)
gbm.plot(annual.no.site.map, common.scale = FALSE)
gbm.plot.fits(annual.no.site.map)

annual.prerun<- plot.gbm.4list(annual.no.site.map)
annual.no.site.map$contributions$var = c("SLA_m2.kg","SRL_m.g","height.m","leafN.mg.g","root.depth_m",
                                         "precip","rootDiam.mm","RTD.g.cm3","rootN.mg.g")
annual.boot <- gbm.bootstrap.functions(annual.no.site.map, list.predictors=annual.prerun, n.reps=1000)

annual.no.site.map$contributions$var = c("SLA","SRL","Height","LeafN","Rooting Depth","Precipitation","Root Diameter","RTD", "RootN")
annual.no.site.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(annual.no.site.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.no.site.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(annual.no.site.map,predictor = "SLA",list.4.preds=annual.prerun, col.line="#979461",
                           booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
                           alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
ggPD_boot(annual.no.site.map,predictor = "SRL",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(annual.no.site.map,predictor = "Height",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(annual.no.site.map,predictor = "LeafN",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(annual.no.site.map,list.4.preds=annual.prerun, col.line="#979461",
                               booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
                               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(annual.no.site.map, col.line="#979461", common.scale = FALSE, y.label = "Percent Cover Change")

annual.no.site.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(annual.no.site.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(annual.no.site.map$gbm.call$dataframe)[24] = "precip"

# investigation of interactions
gbm.interactions(annual.no.site.map)$rank.list
ggInteract_list(annual.no.site.map)
# SLA x height 1.60
# MAP x leafN 0.70
# SLA x leafN 0.66
# height x leafN 0.38

ggInteract_3D(annual.no.site.map, x = 7, y = 1, z.range = c(-2, 1.3))
ggInteract_3D(annual.no.site.map, x = 9, y = 7, z.range = c(-2.1, 0.75))
ggInteract_3D(annual.no.site.map, x = 7, y = 2, z.range = c(-1.2, 1.2))
ggInteract_3D(annual.no.site.map, x = 9, y = 1, z.range = c(-1, 0.85))

# bootstrap interactions

annual.no.site.map$contributions$var = c("SLA_m2.kg","SRL_m.g","height.m","leafN.mg.g",
                                         "root.depth_m","precip","rootDiam.mm","RTD.g.cm3","rootN.mg.g")

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_annual<-ggInteract_boot(c('SLA_m2.kg','height.m'),c('precip','leafN.mg.g'),
                                      c('SLA_m2.kg','leafN.mg.g'),c('height.m','leafN.mg.g'),
                                      nboots = 100,data=annual.data, predictors =  c(12:19,24), 
                                      response="mean.cover.response",
                                      family = "gaussian", tc = 2, lr = 0.001, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_annual, column = 2,obs = 1.60) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 3,obs = 0.70) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 4,obs = 0.66) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 5,obs = 0.38) # non-significant


# save(annual.no.site.map, annual.prerun, annual.boot, file = "./Results/ctrl.v.drt.yr1/annual.output.RData")

#### Perennial without woody ####
#load("./Results/ctrl.v.drt.yr1/perennial.no.woody.output.RData") 

perennial.data.map=gbm.step(data=perennial.data, gbm.x = c(12:19,24), gbm.y=10,
                            family = "gaussian", tree.complexity = 4, learning.rate = 0.00005,
                            bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)


ggPerformance(perennial.data.map)
# 1000 trees Per.Expl = 0.69%
1-(perennial.data.map$self.statistics$mean.resid/perennial.data.map$self.statistics$mean.null) # R2

perennial.data.map$contributions$var = c("Height","SRL","SLA","Root Diameter","LeafN",
                                         "Rooting Depth","RTD","Precipitation", "RootN")
perennial.influence.plot = ggInfluence(perennial.data.map, main = "Perennials (n = 477)", 
                                       col.bar = c("#F1C646","#F1C646","#F1C646","gray70",
                                                   "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")

ggPD(perennial.data.map, rug = T) # partial dependency plots
ggPDfit(perennial.data.map)
gbm.plot(perennial.data.map, common.scale = FALSE)
gbm.plot.fits(perennial.data.map)

perennial.prerun<- plot.gbm.4list(perennial.data.map)
perennial.data.map$contributions$var = c("height.m","SRL_m.g","SLA_m2.kg","rootDiam.mm","leafN.mg.g",
                                         "root.depth_m","RTD.g.cm3","precip","rootN.mg.g")
perennial.boot <- gbm.bootstrap.functions(perennial.data.map, list.predictors=perennial.prerun, n.reps=1000)

perennial.data.map$contributions$var = c("Height","SRL","SLA","Root Diameter","LeafN",
                                         "Rooting Depth","RTD","Precipitation","RootN")
perennial.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(perennial.data.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(perennial.data.map,predictor = "Height",list.4.preds=perennial.prerun, col.line="#F1C646",
                              booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
                              alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
ggPD_boot(perennial.data.map,predictor = "SRL",list.4.preds=perennial.prerun, col.line="#F1C646",
          booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(perennial.data.map,predictor = "SLA",list.4.preds=perennial.prerun, col.line="#F1C646",
          booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")


ggPD_boot(perennial.data.map,list.4.preds=perennial.prerun, col.line="#F1C646",
                                  booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
                                  alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(perennial.data.map, col.line="#F1C646", common.scale = FALSE, y.label = "Percent Cover Change")


# investigation of interactions
gbm.interactions(perennial.data.map)$rank.list
ggInteract_list(perennial.data.map)
# SRL x height 0.21
# precip x height 0.02
# SRL x SLA 0.02
# RTD x height 0.01

ggInteract_3D(perennial.data.map, x = 7, y = 2, z.range = c(-0.5, 0.5))
ggInteract_3D(perennial.data.map, x = 2, y = 1, z.range = c(-0.5, 0))
ggInteract_3D(perennial.data.map, x = 4, y = 2, z.range = c(-0.5, 0.1))
ggInteract_3D(perennial.data.map, x = 7, y = 1, z.range = c(-0.5, 0))

# bootstrap interactions

perennial.data.map$contributions$var = c("height.m","SRL_m.g","SLA_m2.kg","rootDiam.mm",
                                         "leafN.mg.g","root.depth_m","RTD.g.cm3","precip","rootN.mg.g")
perennial.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(perennial.data.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map$gbm.call$dataframe)[24] = "precip"


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_perennial<-ggInteract_boot(c('SRL_m.g','height.m'),c('precip','height.m'),
                                         c('SRL_m.g','SLA_m2.kg'),c('RTD.g.cm3','height.m'),
                                         nboots = 100,data=perennial.data, predictors =  c(12:19,24), 
                                         response="mean.cover.response",
                                         family = "gaussian", tc = 4, lr = 0.00005, bf= 0.75,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_perennial, column = 2,obs = 0.21, bindwidth = 0.1) # significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 3,obs = 0.02, bindwidth = 0.1) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 4,obs = 0.02, bindwidth = 0.1) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 5,obs = 0.01, bindwidth = 0.1) # non-significant

ggInteract_2D(gbm.object = perennial.data.map,x="SRL_m.g",y="height.m",col.gradient = c("white","#F1C646"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SRL", y.label = "Height",
              z.range = c(-0.3, 0), z.label = "% Cover Change")

#save(perennial.data.map, perennial.prerun, perennial.boot, file = "./Results/ctrl.v.drt.yr1/perennial.data.output.RData")

#### Grass ####
#load("./Results/ctrl.v.drt.yr1/grass.output.RData") 

grass.map=gbm.step(data=grass, gbm.x = c(12:19,24), gbm.y=10,
                   family = "gaussian", tree.complexity = 2, learning.rate = 0.0005,
                   bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.map)
# 1000 trees Per.Expl = 6.26%


1-(grass.map$self.statistics$mean.resid/grass.map$self.statistics$mean.null) # R2

grass.map$contributions$var = c("Root Diameter","SLA","Rooting Depth","SRL",
                                "LeafN","Height","RTD","RootN", "Precipitation")
grass.influence.plot = ggInfluence(grass.map, main = "Grasses (n = 221)", 
                                   col.bar = c("#6E687E","#6E687E","#6E687E","gray70",
                                               "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")
ggInfluence(grass.map)

ggPD(grass.map, rug = T) # partial dependency plots
ggPDfit(grass.map)
gbm.plot(grass.map, common.scale = FALSE)
gbm.plot.fits(grass.map)

grass.prerun<- plot.gbm.4list(grass.map)
grass.map$contributions$var = c("rootDiam.mm","SLA_m2.kg","root.depth_m","SRL_m.g","leafN.mg.g","height.m","RTD.g.cm3","rootN.mg.g", "precip")
grass.boot <- gbm.bootstrap.functions(grass.map, list.predictors=grass.prerun, n.reps=1000)

grass.map$contributions$var = c("Root Diameter","SLA","Rooting Depth","SRL",
                                "LeafN","Height","RTD","RootN", "Precipitation")
grass.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(grass.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.map$gbm.call$dataframe)[24] = "Precipitation"


ggPD_boot(grass.map,predictor = "Root Diameter",list.4.preds=grass.prerun, col.line="#6E687E",
                          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
                          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
ggPD_boot(grass.map,predictor = "SLA",list.4.preds=grass.prerun, col.line="#6E687E",
          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(grass.map,predictor = "Rooting Depth",list.4.preds=grass.prerun, col.line="#6E687E",
          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(grass.map,list.4.preds=grass.prerun, col.line="#6E687E",
                              booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
                              alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(grass.map, col.line="#6E687E", common.scale = FALSE, y.label = "Percent Cover Change")

# investigation of interactions
gbm.interactions(grass.map)$rank.list
ggInteract_list(grass.map)
# rootDiam x SLA 2.11
# rootDiam x SRL 0.39
# root.depth x SLA 0.22
# rootDiam x rootN 0.21

ggInteract_3D(grass.map, x = 6, y = 2, z.range = c(-0.5, 0))
ggInteract_3D(grass.map, x = 4, y = 2, z.range = c(-0.5, -0.1))
ggInteract_3D(grass.map, x = 2, y = 1, z.range = c(-0.5, -0.1))

# bootstrap interactions

grass.map$contributions$var = c("rootDiam.mm","SLA_m2.kg","root.depth_m","SRL_m.g","leafN.mg.g","height.m","RTD.g.cm3","rootN.mg.g", "precip")
grass.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(grass.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.map$gbm.call$dataframe)[24] = "precip"


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_grass<-ggInteract_boot(c('rootDiam.mm','SLA_m2.kg'),c('rootDiam.mm','SRL_m.g'),
                                     c('height.m','SLA_m2.kg'),c('rootDiam.mm','rootN.mg.g'),
                                     nboots = 100,data=grass, predictors =  c(12:19,24), 
                                     response="mean.cover.response",
                                     family = "gaussian", tc = 2, lr = 0.0005, bf= 0.50,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_grass, column = 2,obs = 2.11) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass, column = 3,obs = 0.39) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass, column = 4,obs = 0.22) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass, column = 5,obs = 0.21) # non-significant


# save(grass.map, grass.prerun, grass.boot, file = "./Results/ctrl.v.drt.yr1/grass.output.RData")

#### Forb ####
#load("./Results/ctrl.v.drt.yr1/forb.output.RData") 

forb.map=gbm.step(data=forb, gbm.x = c(12:19,24), gbm.y=10,
                  family = "gaussian", tree.complexity = 3, learning.rate = 0.0005,
                  bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)


ggPerformance(forb.map)
# 1900 trees Per.Expl = 12.94%

1-(forb.map$self.statistics$mean.resid/forb.map$self.statistics$mean.null) # R2
forb.map$contributions$var = c("Height","LeafN","SLA","RTD","SRL","Rooting Depth","Root Diameter","Precipitation", "RootN")
forb.influence.plot = ggInfluence(forb.map, main = "Forbs (n = 318)", 
                                  col.bar = c("#F17236","#F17236","#F17236","#F17236",
                                              "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")
ggInfluence(forb.map)

ggPD(forb.map, rug = T) # partial dependency plots
ggPDfit(forb.map)
gbm.plot(forb.map, common.scale = FALSE)
gbm.plot.fits(forb.map)

forb.prerun<- plot.gbm.4list(forb.map)
forb.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","RTD.g.cm3","SRL_m.g",
                               "root.depth_m","rootDiam.mm","precip", "rootN.mg.g")
forb.boot <- gbm.bootstrap.functions(forb.map, list.predictors=forb.prerun, n.reps=1000)

forb.map$contributions$var = c("Height","LeafN","SLA","RTD","SRL","Rooting Depth","Root Diameter","Precipitation", "RootN")
forb.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(forb.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(forb.map,predictor = "Height",list.4.preds=forb.prerun, col.line="#F17236",
                         booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
                         alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
ggPD_boot(forb.map,predictor = "LeafN",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(forb.map,predictor = "SLA",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
ggPD_boot(forb.map,predictor = "RTD",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(forb.map,list.4.preds=forb.prerun, col.line="#F17236",
                             booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
                             alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(forb.map, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")


# investigation of interactions
forb.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","RTD.g.cm3","SRL_m.g",
                               "root.depth_m","rootDiam.mm","precip", "rootN.mg.g")
forb.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(forb.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(forb.map$gbm.call$dataframe)[24] = "precip"

gbm.interactions(forb.map)$rank.list
ggInteract_list(forb.map)
# height x leafN 4.20
# RTD x height 2.28
# SLA x height 1.16
# SRL x leafN 0.90

ggInteract_3D(forb.map, x = 5, y = 4, z.range = c(-0.5, 3))
ggInteract_3D(forb.map, x = 4, y = 2, z.range = c(-2.2, 2.5))
ggInteract_3D(forb.map, x = 9, y = 4, z.range = c(-0.75, 2))
ggInteract_3D(forb.map, x = 4, y = 1, z.range = c(-1.5, 3.75))

# bootstrap interactions
# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_forb<-ggInteract_boot(c('height.m','leafN.mg.g'),c('RTD.g.cm3','height.m'),
                                    c('SLA_m2.kg','height.m'),c('SRL_m.g','leafN.mg.g'),
                                    nboots = 100,data=forb, predictors =  c(12:19,24), 
                                    response="mean.cover.response",
                                    family = "gaussian", tc = 3, lr = 0.0005, bf= 0.75,global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_forb, column = 2,obs = 4.20) # significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 3,obs = 2.28) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 4,obs = 1.16) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 5,obs = 0.90) # non-significant


ggInteract_2D(gbm.object = forb.map, x="height.m",y="leafN.mg.g",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "Height", y.label = "LeafN",
              z.range = c(-1.5, 1), z.label = "% Cover Change")

#save(forb.map, forb.prerun, forb.boot, file = "./Results/ctrl.v.drt.yr1/forb.output.RData")






#### Table S1 ####
mean(all.data$leafN.mg.g, na.rm = TRUE)
mean(all.data$height.m, na.rm = TRUE)
mean(all.data$rootN.mg.g, na.rm = TRUE)
mean(all.data$SLA_m2.kg, na.rm = TRUE)
mean(all.data$root.depth_m, na.rm = TRUE)
mean(all.data$RTD.g.cm3, na.rm = TRUE)
mean(all.data$SRL_m.g, na.rm = TRUE)
mean(all.data$rootDiam.mm, na.rm = TRUE)
mean(all.data$precip, na.rm = TRUE)

sd(all.data$leafN.mg.g, na.rm = TRUE)
sd(all.data$height.m, na.rm = TRUE)
sd(all.data$rootN.mg.g, na.rm = TRUE)
sd(all.data$SLA_m2.kg, na.rm = TRUE)
sd(all.data$root.depth_m, na.rm = TRUE)
sd(all.data$RTD.g.cm3, na.rm = TRUE)
sd(all.data$SRL_m.g, na.rm = TRUE)
sd(all.data$rootDiam.mm, na.rm = TRUE)
sd(all.data$precip, na.rm = TRUE)

range(all.data$leafN.mg.g, na.rm = TRUE)
range(all.data$height.m, na.rm = TRUE)
range(all.data$rootN.mg.g, na.rm = TRUE)
range(all.data$SLA_m2.kg, na.rm = TRUE)
range(all.data$root.depth_m, na.rm = TRUE)
range(all.data$RTD.g.cm3, na.rm = TRUE)
range(all.data$SRL_m.g, na.rm = TRUE)
range(all.data$rootDiam.mm, na.rm = TRUE)
range(all.data$precip, na.rm = TRUE)

sum(is.na(all.data$leafN.mg.g))/616*100
sum(is.na(all.data$height.m))/616*100
sum(is.na(all.data$rootN.mg.g))/616*100
sum(is.na(all.data$SLA_m2.kg))/616*100
sum(is.na(all.data$root.depth_m))/616*100
sum(is.na(all.data$RTD.g.cm3))/616*100
sum(is.na(all.data$SRL_m.g))/616*100
sum(is.na(all.data$rootDiam.mm))/616*100
sum(is.na(all.data$precip))/616*100

mean(annual.data$leafN.mg.g, na.rm = TRUE)
mean(annual.data$height.m, na.rm = TRUE)
mean(annual.data$rootN.mg.g, na.rm = TRUE)
mean(annual.data$SLA_m2.kg, na.rm = TRUE)
mean(annual.data$root.depth_m, na.rm = TRUE)
mean(annual.data$RTD.g.cm3, na.rm = TRUE)
mean(annual.data$SRL_m.g, na.rm = TRUE)
mean(annual.data$rootDiam.mm, na.rm = TRUE)

sd(annual.data$leafN.mg.g, na.rm = TRUE)
sd(annual.data$height.m, na.rm = TRUE)
sd(annual.data$rootN.mg.g, na.rm = TRUE)
sd(annual.data$SLA_m2.kg, na.rm = TRUE)
sd(annual.data$root.depth_m, na.rm = TRUE)
sd(annual.data$RTD.g.cm3, na.rm = TRUE)
sd(annual.data$SRL_m.g, na.rm = TRUE)
sd(annual.data$rootDiam.mm, na.rm = TRUE)
sd(annual.data$precip, na.rm = TRUE)

range(annual.data$leafN.mg.g, na.rm = TRUE)
range(annual.data$height.m, na.rm = TRUE)
range(annual.data$rootN.mg.g, na.rm = TRUE)
range(annual.data$SLA_m2.kg, na.rm = TRUE)
range(annual.data$root.depth_m, na.rm = TRUE)
range(annual.data$RTD.g.cm3, na.rm = TRUE)
range(annual.data$SRL_m.g, na.rm = TRUE)
range(annual.data$rootDiam.mm, na.rm = TRUE)
range(annual.data$precip, na.rm = TRUE)

mean(perennial.data$leafN.mg.g, na.rm = TRUE)
mean(perennial.data$height.m, na.rm = TRUE)
mean(perennial.data$rootN.mg.g, na.rm = TRUE)
mean(perennial.data$SLA_m2.kg, na.rm = TRUE)
mean(perennial.data$root.depth_m, na.rm = TRUE)
mean(perennial.data$RTD.g.cm3, na.rm = TRUE)
mean(perennial.data$SRL_m.g, na.rm = TRUE)
mean(perennial.data$rootDiam.mm, na.rm = TRUE)

sd(perennial.data$leafN.mg.g, na.rm = TRUE)
sd(perennial.data$height.m, na.rm = TRUE)
sd(perennial.data$rootN.mg.g, na.rm = TRUE)
sd(perennial.data$SLA_m2.kg, na.rm = TRUE)
sd(perennial.data$root.depth_m, na.rm = TRUE)
sd(perennial.data$RTD.g.cm3, na.rm = TRUE)
sd(perennial.data$SRL_m.g, na.rm = TRUE)
sd(perennial.data$rootDiam.mm, na.rm = TRUE)

range(perennial.data$leafN.mg.g, na.rm = TRUE)
range(perennial.data$height.m, na.rm = TRUE)
range(perennial.data$rootN.mg.g, na.rm = TRUE)
range(perennial.data$SLA_m2.kg, na.rm = TRUE)
range(perennial.data$root.depth_m, na.rm = TRUE)
range(perennial.data$RTD.g.cm3, na.rm = TRUE)
range(perennial.data$SRL_m.g, na.rm = TRUE)
range(perennial.data$rootDiam.mm, na.rm = TRUE)

mean(grass$leafN.mg.g, na.rm = TRUE)
mean(grass$height.m, na.rm = TRUE)
mean(grass$rootN.mg.g, na.rm = TRUE)
mean(grass$SLA_m2.kg, na.rm = TRUE)
mean(grass$root.depth_m, na.rm = TRUE)
mean(grass$RTD.g.cm3, na.rm = TRUE)
mean(grass$SRL_m.g, na.rm = TRUE)
mean(grass$rootDiam.mm, na.rm = TRUE)

sd(grass$leafN.mg.g, na.rm = TRUE)
sd(grass$height.m, na.rm = TRUE)
sd(grass$rootN.mg.g, na.rm = TRUE)
sd(grass$SLA_m2.kg, na.rm = TRUE)
sd(grass$root.depth_m, na.rm = TRUE)
sd(grass$RTD.g.cm3, na.rm = TRUE)
sd(grass$SRL_m.g, na.rm = TRUE)
sd(grass$rootDiam.mm, na.rm = TRUE)

range(grass$leafN.mg.g, na.rm = TRUE)
range(grass$height.m, na.rm = TRUE)
range(grass$rootN.mg.g, na.rm = TRUE)
range(grass$SLA_m2.kg, na.rm = TRUE)
range(grass$root.depth_m, na.rm = TRUE)
range(grass$RTD.g.cm3, na.rm = TRUE)
range(grass$SRL_m.g, na.rm = TRUE)
range(grass$rootDiam.mm, na.rm = TRUE)

mean(forb$leafN.mg.g, na.rm = TRUE)
mean(forb$height.m, na.rm = TRUE)
mean(forb$rootN.mg.g, na.rm = TRUE)
mean(forb$SLA_m2.kg, na.rm = TRUE)
mean(forb$root.depth_m, na.rm = TRUE)
mean(forb$RTD.g.cm3, na.rm = TRUE)
mean(forb$SRL_m.g, na.rm = TRUE)
mean(forb$rootDiam.mm, na.rm = TRUE)

sd(forb$leafN.mg.g, na.rm = TRUE)
sd(forb$height.m, na.rm = TRUE)
sd(forb$rootN.mg.g, na.rm = TRUE)
sd(forb$SLA_m2.kg, na.rm = TRUE)
sd(forb$root.depth_m, na.rm = TRUE)
sd(forb$RTD.g.cm3, na.rm = TRUE)
sd(forb$SRL_m.g, na.rm = TRUE)
sd(forb$rootDiam.mm, na.rm = TRUE)

range(forb$leafN.mg.g, na.rm = TRUE)
range(forb$height.m, na.rm = TRUE)
range(forb$rootN.mg.g, na.rm = TRUE)
range(forb$SLA_m2.kg, na.rm = TRUE)
range(forb$root.depth_m, na.rm = TRUE)
range(forb$RTD.g.cm3, na.rm = TRUE)
range(forb$SRL_m.g, na.rm = TRUE)
range(forb$rootDiam.mm, na.rm = TRUE)
