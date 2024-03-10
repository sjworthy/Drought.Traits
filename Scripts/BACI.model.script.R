# Script to peform analyses where % cover calculated using BACI method
# drought after-drought before) - (control after - control before)

# Use the same data files just removing outlier cover change tested for in All Data

library(dismo)
library(gbm)
library(ggBRT)
library(multcomp)
library(cowplot)
library(corrplot)
library(treezy)

# Modified ggPD_boot function to make y-axis the same for all plots.
# edited function is script: ggPD_boot_edit.R

#### read in all the data frames needs for the analyses ####

# all data
all.data = read.csv("./New.dfs/BACI.all.data.csv", row.names = 1) # 1234 data points
# annual data
annual.data = read.csv("./New.dfs/BACI.annual.data.csv", row.names = 1) # 292 data points
# all perennial data
perennial.data = read.csv("./New.dfs/BACI.perennial.data.csv", row.names = 1) # 890 data points
# grass
grass = read.csv("./New.dfs/BACI.grass.csv", row.names = 1) # 392 data points
# forbs
forb = read.csv("./New.dfs/BACI.forb.csv", row.names = 1) # 675 data points
# grass perennial
grass.perennial = read.csv("./New.dfs/BACI.grass.perennial.csv", row.names = 1) # 310 data points
# grass annual
grass.annual = read.csv("./New.dfs/BACI.grass.annual.csv", row.names = 1) # 73 data points

#### Testing for outliers in all data ####

hist(all.data$cover.change)
boxplot(all.data$cover.change)
mean = mean(all.data$cover.change, na.rm = TRUE)
std = sd(all.data$cover.change, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(all.data$cover.change[which(all.data$cover.change <Tmin | all.data$cover.change > Tmax)])
# removed cover.change -79.47619 to -20.58333
# removed cover.change 20.00000 to 43.00000

all.data.otl.rm = subset(all.data, all.data$cover.change >-20.58333 & all.data$cover.change < 20.00000)
# lose 30 data points, lose 13 species
annual.data.otl.rm = subset(annual.data, annual.data$cover.change >-20.58333 & annual.data$cover.change < 20.00000)
# lose 8 data points, lose 1 species
perennial.data.otl.rm = subset(perennial.data, perennial.data$cover.change >-20.58333 & perennial.data$cover.change < 20.00000)
# lose 21 data points, lose 11 species
grass.data.otl.rm = subset(grass, grass$cover.change >-20.58333 & grass$cover.change < 20.00000)
# lose 19 data points, lose 7 species
forb.data.otl.rm = subset(forb, forb$cover.change >-20.58333 & forb$cover.change < 20.00000)
# lose 9 data points, lose 4 species
grass.perennial.data.otl.rm = subset(grass.perennial, grass.perennial$cover.change >-20.58333 & grass.perennial$cover.change < 20.00000)
# lose 13 data points, lose 6 species
grass.annual.data.otl.rm = subset(grass.annual, grass.annual$cover.change >-20.58333 & grass.annual$cover.change < 20.00000)
# lose 6 data points, lose 1 species

#### change site code to numeric, continuous vector ####
all.data.otl.rm$site.id = as.numeric(as.factor(all.data.otl.rm$site_code))
annual.data.otl.rm$site.id = as.numeric(as.factor(annual.data.otl.rm$site_code))
perennial.data.otl.rm$site.id = as.numeric(as.factor(perennial.data.otl.rm$site_code))
grass.data.otl.rm$site.id = as.numeric(as.factor(grass.data.otl.rm$site_code))
forb.data.otl.rm$site.id = as.numeric(as.factor(forb.data.otl.rm$site_code))
grass.perennial.data.otl.rm$site.id = as.numeric(as.factor(grass.perennial.data.otl.rm$site_code))
grass.annual.data.otl.rm$site.id = as.numeric(as.factor(grass.annual.data.otl.rm$site_code))

#### merge with environmental data ####
# mean annual precipitation data (MAP)

Site.info = read.csv("./Raw.Data/Site_Elev-Disturb.csv")
site.info.map = Site.info[,c(2,13)]

# drought severity index

dsi = read.csv("./Formatted.Data/site.drt.dev.index.csv", row.names = 1)

# merge MAP and DSI

enviro.data = left_join(site.info.map,dsi, by = "site_code")

# override original data names so can use same code
all.data = left_join(all.data.otl.rm, enviro.data, by="site_code") # 83 unique sites
annual.data = left_join(annual.data.otl.rm, enviro.data, by="site_code")
perennial.data = left_join(perennial.data.otl.rm, enviro.data, by="site_code")
grass = left_join(grass.data.otl.rm, enviro.data, by="site_code")
forb = left_join(forb.data.otl.rm, enviro.data, by="site_code")
grass.perennial = left_join(grass.perennial.data.otl.rm, enviro.data, by="site_code")
grass.annual = left_join(grass.annual.data.otl.rm, enviro.data, by="site_code")

#### testing for correlation among traits ####

# all.data
all.data.2 = all.data
colnames(all.data.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.2)[22] = "Precipitation"
colnames(all.data.2)[23] = "DSI"

cor.mat.all = cor(all.data.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.all, method="number",tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

annual.data.2 = annual.data
colnames(annual.data.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.data.2)[22] = "Precipitation"
colnames(annual.data.2)[23] = "DSI"

cor.mat.annual = cor(annual.data.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.annual, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

perennial.data.2 = perennial.data
colnames(perennial.data.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.2)[22] = "Precipitation"
colnames(perennial.data.2)[23] = "DSI"

cor.mat.perennial = cor(perennial.data.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.perennial, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

grass.2 = grass
colnames(grass.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.2)[22] = "Precipitation"
colnames(grass.2)[23] = "DSI"

cor.mat.grass = cor(grass.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.grass, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

forb.2 = forb
colnames(forb.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.2)[22] = "Precipitation"
colnames(forb.2)[23] = "DSI"

cor.mat.forb = cor(forb.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.forb, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

grass.perennial.2 = grass.perennial
colnames(grass.perennial.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.perennial.2)[22] = "Precipitation"
colnames(grass.perennial.2)[23] = "DSI"

cor.mat.grass.perennial = cor(grass.perennial.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.grass.perennial, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

grass.annual.2 = grass.annual
colnames(grass.annual.2)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.annual.2)[22] = "Precipitation"
colnames(grass.annual.2)[23] = "DSI"

cor.mat.grass.annual = cor(grass.annual.2[,c(10:17,22,23)],use = "pairwise") 
corrplot(cor.mat.grass.annual, method="number", tl.col = "black", is.corr = TRUE,
         col.lim = c(-1,1), col = COL2('BrBG', 200), addgrid.col = "black")

#### determining best parameter combination to generate 1000 trees ####
# bag fraction of 0.50 and 0.75, step.size of 25 and 50, tc = 10, 

all.data.map.dsi=gbm.step(data=all.data, gbm.x = c(10:17,22,23), gbm.y=9,
                          family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(all.data.map.dsi)

all.data.map=gbm.step(data=all.data, gbm.x = c(10:17,22), gbm.y=9,
                      family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(all.data.map)

# not used
all.data.dsi=gbm.step(data=all.data, gbm.x = c(10:17,23), gbm.y=9,
                      family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                      bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(all.data.dsi)

annual.data.map.dsi=gbm.step(data=annual.data, gbm.x = c(10:17,22,23), gbm.y=9,
                             family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                             bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.data.map.dsi)

annual.data.map=gbm.step(data=annual.data, gbm.x = c(10:17,22), gbm.y=9,
                         family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                         bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.data.map)

# not used
annual.data.dsi=gbm.step(data=annual.data, gbm.x = c(10:17,23), gbm.y=9,
                         family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                         bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(annual.data.dsi)

perennial.data.map.dsi=gbm.step(data=perennial.data, gbm.x = c(10:17,22,23), gbm.y=9,
                                family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(perennial.data.map.dsi)

perennial.data.map=gbm.step(data=perennial.data, gbm.x = c(10:17,22), gbm.y=9,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(perennial.data.map)

# not used
perennial.data.dsi=gbm.step(data=perennial.data, gbm.x = c(10:17,23), gbm.y=9,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(perennial.data.dsi)

grass.data.map.dsi=gbm.step(data=grass, gbm.x = c(10:17,22,23), gbm.y=9,
                            family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                            bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(grass.data.map.dsi)

grass.data.map=gbm.step(data=grass, gbm.x = c(10:17,22), gbm.y=9,
                        family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                        bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(grass.data.map)

# not used
grass.data.dsi=gbm.step(data=grass, gbm.x = c(10:17,23), gbm.y=9,
                        family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                        bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.data.dsi)

forb.data.map.dsi=gbm.step(data=forb, gbm.x = c(10:17,22,23), gbm.y=9,
                           family = "gaussian", tree.complexity = 10, learning.rate = 0.00001,
                           bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.data.map.dsi)

forb.data.map=gbm.step(data=forb, gbm.x = c(10:17,22), gbm.y=9,
                       family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                       bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(forb.data.map)

# not used
forb.data.dsi=gbm.step(data=forb, gbm.x = c(10:17,23), gbm.y=9,
                       family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                       bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(forb.data.dsi)

grass.perennial.data.map.dsi=gbm.step(data=grass.perennial, gbm.x = c(10:17,22,23), gbm.y=9,
                                      family = "gaussian", tree.complexity = 10, learning.rate = 0.005,
                                      bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.perennial.data.map.dsi)

grass.perennial.data.map=gbm.step(data=grass.perennial, gbm.x = c(10:17,22), gbm.y=9,
                                  family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                                  bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.perennial.data.map)

# not used
grass.perennial.data.dsi=gbm.step(data=grass.perennial, gbm.x = c(10:17,23), gbm.y=9,
                                  family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                                  bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(grass.perennial.data.dsi)

grass.annual.data.map.dsi=gbm.step(data=grass.annual, gbm.x = c(10:17,22,23), gbm.y=9,
                                   family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                                   bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(grass.annual.data.map.dsi)

grass.annual.data.map=gbm.step(data=grass.annual, gbm.x = c(10:17,22), gbm.y=9,
                               family = "gaussian", tree.complexity = 10, learning.rate = 0.0005,
                               bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.annual.data.map)

# not used
grass.annual.data.dsi=gbm.step(data=grass.annual, gbm.x = c(10:17,23), gbm.y=9,
                               family = "gaussian", tree.complexity = 10, learning.rate = 0.001,
                               bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 25)
ggPerformance(grass.annual.data.dsi)

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
    
    BRT.all.variables <- gbm.step(data=grass.perennial,
                                  gbm.x = c(10:17,22),
                                  gbm.y = 9,
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
text(x= 1:length(R2Obs.all.variables), y=77, labels=TukeyLetters)

# selected the lowest tc value that did not show significant differences 
# compared to the largest tc value
# with just map
# all.data tc = 8
# annual tc = 2
# perennial tree tc = 8
# grass tc = 5
# forb tc = 3
# grass.perennial = 2

# with map and dsi
# all.data tc = 5
# annual tc = 1
# perennial tree tc = 9
# grass tc = 6
# forb tc = 7
# grass.perennial = 4

#### all data without woody ####

all.data.map.dsi=gbm.step(data=all.data, gbm.x = c(10:17,22,23), gbm.y=9,
                          family = "gaussian", tree.complexity = 5, learning.rate = 0.0005,
                          bag.fraction = 0.5, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(all.data.map.dsi)
# 1650 trees Per.Expl = 6.44%

all.data.map=gbm.step(data=all.data, gbm.x = c(10:17,22), gbm.y=9,
                      family = "gaussian", tree.complexity = 8, learning.rate = 0.0005,
                      bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(all.data.map)
# 1200 trees Per.Expl = 6.76%

all.data.map$contributions$var = c("LeafN","Height","Root Diameter","Precipitation", 
                                   "RTD","SLA","Rooting Depth","SRL","RootN")
all.data.influence.plot = ggInfluence(all.data.map, main = "All Species (n = 612)", 
                                      col.bar = c("#769370","#769370","gray70","gray70",
                                                  "gray70","gray70","gray70","gray70","gray70"), 
                                      col.signif = "#B50200")

all.data.influence.plot = ggInfluence(all.data.map.dsi, main = "All Species (n = 612)", 
                                      col.bar = c("#769370","#769370","gray70","gray70",
                                                  "gray70","gray70","gray70","gray70","gray70","gray70"), 
                                      col.signif = "#B50200")

# saved as pdf from device (4 x 4)

ggPD(all.data.map, rug = T, common.scale = FALSE) # partial dependency plots
ggPD(all.data.map, predictor = "height.m",rug = T, smooth = TRUE)
ggPDfit(all.data.map)
gbm.plot(all.data.map,common.scale = FALSE, rug = TRUE)
gbm.plot(all.data.map.dsi,common.scale = FALSE, rug = TRUE)

gbm.plot.fits(all.data.map)

gg_partial_plot(all.data.map,
                var = c("height.m","SLA_m2.kg",
                        "leafN.mg.g"))

# get data to plot partial dependency plots
all.data.map.prerun<- plot.gbm.4list(all.data.map)
all.data.map$contributions$var = c("leafN.mg.g","height.m","rootDiam.mm","precip",
                                   "RTD.g.cm3","SLA_m2.kg","root.depth_m","SRL_m.g","rootN.mg.g")
all.data.map.boot <- gbm.bootstrap.functions(all.data.map, list.predictors=all.data.map.prerun, n.reps=1000)

all.data.map$contributions$var = c("LeafN","Height","Root Diameter","Precipitation","RTD",
                                   "SLA","Rooting Depth","SRL", "RootN")
all.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(all.data.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(all.data.map,predictor = "Height",list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,
          y.label = "")
# same y-limits
ggPD_boot_test(all.data.map,predictor = "Height",list.4.preds=all.data.map.prerun, col.line="#769370",
               booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,
               y.label = "")

ggPD_boot(all.data.map,predictor = "LeafN",list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-limits
ggPD_boot_test(all.data.map,predictor = "LeafN",list.4.preds=all.data.map.prerun, col.line="#769370",
               booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(all.data.map,predictor = "Root Diameter",list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-limits
ggPD_boot_test(all.data.map,predictor = "Root Diameter",list.4.preds=all.data.map.prerun, col.line="#769370",
               booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(all.data.map,list.4.preds=all.data.map.prerun, col.line="#769370",
          booted.preds=all.data.map.boot$function.preds, cex.line=1, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(all.data.map, col.line="#769370", common.scale = FALSE, y.label = "Percent Cover Change")

# output individual plots as 3x3
# output all plots as one 6x6

# investigation of interactions
all.data.map$contributions$var = c("leafN.mg.g","height.m","rootDiam.mm","precip",
                                   "RTD.g.cm3","SLA_m2.kg","root.depth_m","SRL_m.g","rootN.mg.g")
all.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(all.data.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(all.data.map$gbm.call$dataframe)[24] = "precip"

gbm.interactions(all.data.map)$interactions
ggInteract_list(all.data.map, index = T)
# rootDiam x height.m 125.19
# SLA x height.m 99.20
# root.depth x SLA 47.97
# RTD x SLA 32.43

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_brt2.all<-ggInteract_boot(c('rootDiam.mm','height.m'),c('SLA_m2.kg','height.m'),
                                        c('root.depth_m','SLA_m2.kg'),c('RTD.g.cm3','SLA_m2.kg'),
                                        nboots = 500, data=all.data, predictors =  c(10:17,22), 
                                        response="cover.change",
                                        family = "gaussian", tc = 3, lr = 0.001, bf= 0.75, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 2,obs = 125.19) # significant
# low root diameter x taller height
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 3,obs = 99.20) # significant
# larger SLA x taller height
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 4,obs = 47.97) # significant
# higher SLA x lower rooting depth
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 5,obs = 32.43) # significant
# higher SLA x higher RTD

ggInteract_2D(gbm.object = all.data.map, x="rootDiam.mm",y="height.m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "rootDiam.mm", y.label = "height.m",
              z.range = c(-3.5, 2), z.label = "% Cover Change")

ggInteract_2D(gbm.object = all.data.map, x="SLA_m2.kg",y="height.m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "height.m",
              z.range = c(-1.5, 7.5), z.label = "% Cover Change")

ggInteract_2D(gbm.object = all.data.map, x="SLA_m2.kg",y="root.depth_m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "root.depth_m",
              z.range = c(-2, 3.5), z.label = "% Cover Change")

ggInteract_2D(gbm.object = all.data.map, x="SLA_m2.kg",y="RTD.g.cm3",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "RTD.g.cm3",
              z.range = c(-2.5, 4), z.label = "% Cover Change")

ggInteract_list(all.data.map.dsi, index = T)
# SLA x height.m 46.59
# root.depth x SLA 24.01
# rootDiam x height.m 16.24
# height x leafN 8.24
# SLA x leafN 6.95

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_brt2.all<-ggInteract_boot(c('SLA_m2.kg','height.m'),c('root.depth_m','SLA_m2.kg'),
                                        c('rootDiam.mm','height.m'),c('height.m','leafN.mg.g'),
                                        c('SLA_m2.kg','leafN.mg.g'),
                                        nboots = 500, data=all.data, predictors =  c(10:17,22,23), 
                                        response="cover.change",
                                        family = "gaussian", tc = 4, lr = 0.001, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 2,obs = 46.59) # significant
# larger SLA x taller height
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 3,obs = 24.01) # significant
# higher SLA x lower rooting depth
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 4,obs = 16.24) # significant
# higher height x lower rooting depth
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 5,obs = 8.24) # significant
# higher leafN x taller height
ggInteract_boot_hist(data = Interact_boot_brt2.all, column = 6,obs = 6.95) # significant
# higher SLA x higher RTD

ggInteract_2D(gbm.object = all.data.map.dsi, x="SLA_m2.kg",y="height.m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "height.m",
              z.range = c(-2, 4), z.label = "% Cover Change")
ggInteract_2D(gbm.object = all.data.map.dsi, x="SLA_m2.kg",y="root.depth_m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "root.depth_m",
              z.range = c(-2, 2.5), z.label = "% Cover Change")
ggInteract_2D(gbm.object = all.data.map.dsi, x="rootDiam.mm",y="height.m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "rootDiam.mm", y.label = "height.m",
              z.range = c(-3, 0), z.label = "% Cover Change")
ggInteract_2D(gbm.object = all.data.map.dsi, x="leafN.mg.g",y="height.m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "leafN.mg.g", y.label = "height.m",
              z.range = c(-2.5, 1), z.label = "% Cover Change")
ggInteract_2D(gbm.object = all.data.map.dsi, x="leafN.mg.g",y="SLA_m2.kg",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "leafN.mg.g", y.label = "SLA_m2.kg",
              z.range = c(-2, 3), z.label = "% Cover Change")

#### annual data ####
annual.no.site.map.dsi=gbm.step(data=annual.data, gbm.x = c(10:17,22,23), gbm.y=9,
                                family = "gaussian", tree.complexity = 1, learning.rate = 0.0005,
                                bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.no.site.map.dsi)
# 1350 trees Per.Expl = 4.08%

annual.no.site.map=gbm.step(data=annual.data, gbm.x = c(10:17,22), gbm.y=9,
                            family = "gaussian", tree.complexity = 2, learning.rate = 0.0005,
                            bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(annual.no.site.map)
# 1400 trees Per.Expl = 6.64%

annual.no.site.map$contributions$var = c("SLA","Height","SRL","Precipitation","LeafN","Rooting Depth",
                                         "RTD","Root Diameter", "RootN")
annual.influence.plot = ggInfluence(annual.no.site.map, main = "Annuals (n = 124)", 
                                    col.bar = c("#979461","#979461","#979461","#979461",
                                                "gray70","gray70","gray70","gray70","gray70"),
                                    col.signif = "#B50200")

annual.influence.plot = ggInfluence(annual.no.site.map.dsi, main = "Annuals (n = 124)", 
                                    col.bar = c("#979461","#979461","#979461","#979461",
                                                "gray70","gray70","gray70","gray70","gray70","gray70"),
                                    col.signif = "#B50200")

ggPD(annual.no.site.map, rug = T) # partial dependency plots
ggPDfit(annual.no.site.map)
gbm.plot(annual.no.site.map, common.scale = FALSE)
gbm.plot(annual.no.site.map.dsi, common.scale = FALSE)

gbm.plot.fits(annual.no.site.map)

annual.prerun<- plot.gbm.4list(annual.no.site.map)
annual.no.site.map$contributions$var = c("SLA_m2.kg","height.m","SRL_m.g","precip","leafN.mg.g","root.depth_m",
                                         "RTD.g.cm3","rootDiam.mm","rootN.mg.g")
annual.boot <- gbm.bootstrap.functions(annual.no.site.map, list.predictors=annual.prerun, n.reps=1000)

annual.no.site.map$contributions$var = c("SLA","Height","SRL","Precipitation","LeafN","Rooting Depth","RTD", "Root Diameter","RootN")
annual.no.site.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(annual.no.site.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.no.site.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(annual.no.site.map,predictor = "SLA",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
# same y-axis
ggPD_boot_test(annual.no.site.map,predictor = "SLA",list.4.preds=annual.prerun, col.line="#979461",
               booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(annual.no.site.map,predictor = "SRL",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(annual.no.site.map,predictor = "SRL",list.4.preds=annual.prerun, col.line="#979461",
               booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(annual.no.site.map,predictor = "Height",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(annual.no.site.map,predictor = "Height",list.4.preds=annual.prerun, col.line="#979461",
               booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(annual.no.site.map,predictor = "LeafN",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(annual.no.site.map,predictor = "LeafN",list.4.preds=annual.prerun, col.line="#979461",
               booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(annual.no.site.map,predictor = "Precipitation",list.4.preds=annual.prerun, col.line="#979461",
          booted.preds=annual.boot$function.preds, cex.line=1, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(annual.no.site.map,predictor = "Precipitation",list.4.preds=annual.prerun, col.line="#979461",
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
# SLA x height 4.34
# MAP x leafN 3.56
# MAP x height 1.62
# MAP x SRL 0.93

ggInteract_3D(annual.no.site.map, x = 7, y = 1, z.range = c(-2, 1.3))
ggInteract_3D(annual.no.site.map, x = 9, y = 7, z.range = c(-2.1, 0.75))
ggInteract_3D(annual.no.site.map, x = 7, y = 2, z.range = c(-1.2, 1.2))
ggInteract_3D(annual.no.site.map, x = 9, y = 1, z.range = c(-1, 0.85))

# bootstrap interactions

annual.no.site.map$contributions$var = c("SLA_m2.kg","height.m","SRL_m.g","precip","leafN.mg.g",
                                         "root.depth_m","RTD.g.cm3","rootDiam.mm","rootN.mg.g")

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_annual<-ggInteract_boot(c('SLA_m2.kg','height.m'),c('precip','leafN.mg.g'),
                                      c('precip','height.m'),c('precip','SRL_m.g'),
                                      nboots = 500,data=annual.data, predictors =  c(12:19,24), 
                                      response="mean.cover.response",
                                      family = "gaussian", tc = 3, lr = 0.001, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_annual, column = 2,obs = 4.34) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 3,obs = 3.56) # significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 4,obs = 1.62) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 5,obs = 0.93) # non-significant

ggInteract_2D(gbm.object = annual.no.site.map, x="precip",y="leafN.mg.g",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "LeafN",
              z.range = c(-1.8, 0.25), z.label = "% Cover Change")

# save(annual.no.site.map, annual.prerun, annual.boot, file = "./Results/ctrl.v.drt.yr1/annual.output.RData")

#### Perennial without woody ####
perennial.data.map.dsi=gbm.step(data=perennial.data, gbm.x = c(10:17,22,23), gbm.y=9,
                                family = "gaussian", tree.complexity = 9, learning.rate = 0.0001,
                                bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 50)


ggPerformance(perennial.data.map.dsi)
# 1000 trees Per.Expl = 1.70%

perennial.data.map=gbm.step(data=perennial.data, gbm.x = c(10:17,22), gbm.y=9,
                            family = "gaussian", tree.complexity = 8, learning.rate = 0.0005,
                            bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(perennial.data.map)
# 1000 trees Per.Expl = 6.55%

perennial.data.map$contributions$var = c("Height","LeafN","SLA","Precipitation","Root Diameter",
                                         "SRL","Rooting Depth","RTD","RootN")
perennial.influence.plot = ggInfluence(perennial.data.map, main = "Perennials (n = 474)", 
                                       col.bar = c("#F1C646","#F1C646","#F1C646","gray70",
                                                   "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")


perennial.influence.plot = ggInfluence(perennial.data.map.dsi, main = "Perennials (n = 474)", 
                                       col.bar = c("#F1C646","#F1C646","gray70","gray70",
                                                   "gray70","gray70","gray70","gray70","gray70","gray70"
                                       ), col.signif = "#B50200")


ggPD(perennial.data.map, rug = T) # partial dependency plots
ggPDfit(perennial.data.map)
gbm.plot(perennial.data.map, common.scale = FALSE)
gbm.plot(perennial.data.map.dsi, common.scale = FALSE)

gbm.plot.fits(perennial.data.map)

perennial.prerun<- plot.gbm.4list(perennial.data.map)
perennial.data.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","precip","rootDiam.mm","SRL_m.g",
                                         "root.depth_m","RTD.g.cm3","rootN.mg.g")
perennial.boot <- gbm.bootstrap.functions(perennial.data.map, list.predictors=perennial.prerun, n.reps=1000)

perennial.data.map$contributions$var = c("Height","LeafN","SLA","Precipitation","Root Diameter",
                                         "SRL","Rooting Depth","RTD","RootN")
perennial.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(perennial.data.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(perennial.data.map,predictor = "Height",list.4.preds=perennial.prerun, col.line="#F1C646",
          booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
# same y-axis
ggPD_boot_test(perennial.data.map,predictor = "Height",list.4.preds=perennial.prerun, col.line="#F1C646",
               booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(perennial.data.map,predictor = "SLA",list.4.preds=perennial.prerun, col.line="#F1C646",
          booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(perennial.data.map,predictor = "SLA",list.4.preds=perennial.prerun, col.line="#F1C646",
               booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(perennial.data.map,predictor = "LeafN",list.4.preds=perennial.prerun, col.line="#F1C646",
          booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(perennial.data.map,predictor = "LeafN",list.4.preds=perennial.prerun, col.line="#F1C646",
               booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(perennial.data.map,list.4.preds=perennial.prerun, col.line="#F1C646",
          booted.preds=perennial.boot$function.preds, cex.line=1, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(perennial.data.map, col.line="#F1C646", common.scale = FALSE, y.label = "Percent Cover Change")


perennial.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(perennial.data.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map$gbm.call$dataframe)[24] = "precip"


# investigation of interactions
gbm.interactions(perennial.data.map)$rank.list
ggInteract_list(perennial.data.map)
# SRL x SLA 0.05
# SRL x height 0.05
# precip x height 0.03
# SLA x height 0.02

ggInteract_3D(perennial.data.map, x = 7, y = 2, z.range = c(-0.5, 0.5))
ggInteract_3D(perennial.data.map, x = 2, y = 1, z.range = c(-0.5, 0))
ggInteract_3D(perennial.data.map, x = 4, y = 2, z.range = c(-0.5, 0.1))
ggInteract_3D(perennial.data.map, x = 7, y = 1, z.range = c(-0.5, 0))

# bootstrap interactions

perennial.data.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","precip","rootDiam.mm",
                                         "SRL_m.g","root.depth_m","RTD.g.cm3","rootN.mg.g")
perennial.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(perennial.data.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map$gbm.call$dataframe)[24] = "precip"


# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_perennial<-ggInteract_boot(c('SRL_m.g','SLA_m2.kg'),c('SRL_m.g','height.m'),
                                         c('precip','height.m'),c('SLA_m2.kg','height.m'),
                                         nboots = 500,data=perennial.data, predictors =  c(12:19,24), 
                                         response="mean.cover.response",
                                         family = "gaussian", tc = 8, lr = 0.00005, bf= 0.75, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_perennial, column = 2,obs = 0.05) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 3,obs = 0.05) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 4,obs = 0.03) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial, column = 5,obs = 0.02) # non-significant

#### Grass ####
grass.map.dsi=gbm.step(data=grass, gbm.x = c(10:17,22,23), gbm.y=9,
                       family = "gaussian", tree.complexity = 6, learning.rate = 0.0001,
                       bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.map.dsi)
# 1850 trees Per.Expl = 3.51%

grass.map=gbm.step(data=grass, gbm.x = c(10:17,22), gbm.y=9,
                   family = "gaussian", tree.complexity = 5, learning.rate = 0.0001,
                   bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)
ggPerformance(grass.map)
# 1000 trees Per.Expl = 1.77%

grass.map$contributions$var = c("Root Diameter","Rooting Depth","SLA",
                                "LeafN","Height","SRL","RTD","RootN", "Precipitation")
grass.influence.plot = ggInfluence(grass.map, main = "Grasses (n = 221)", 
                                   col.bar = c("#6E687E","#6E687E","#6E687E","gray70",
                                               "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")

grass.influence.plot = ggInfluence(grass.map.dsi, main = "Grasses (n = 221)", 
                                   col.bar = c("#6E687E","#6E687E","#6E687E","gray70",
                                               "gray70","gray70","gray70","gray70","gray70","gray70"
                                   ), col.signif = "#B50200")


ggInfluence(grass.map)

ggPD(grass.map, rug = T) # partial dependency plots
ggPDfit(grass.map)
gbm.plot(grass.map, common.scale = FALSE)
gbm.plot(grass.map.dsi, common.scale = FALSE)

gbm.plot.fits(grass.map)

grass.prerun<- plot.gbm.4list(grass.map)
grass.map$contributions$var = c("rootDiam.mm","root.depth_m","SLA_m2.kg","leafN.mg.g","height.m","SRL_m.g","RTD.g.cm3","rootN.mg.g", "precip")
grass.boot <- gbm.bootstrap.functions(grass.map, list.predictors=grass.prerun, n.reps=1000)

grass.map$contributions$var = c("Root Diameter","Rooting Depth","SLA",
                                "LeafN","Height","SRL","RTD","RootN", "Precipitation")
grass.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(grass.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.map$gbm.call$dataframe)[24] = "Precipitation"


ggPD_boot(grass.map,predictor = "Root Diameter",list.4.preds=grass.prerun, col.line="#6E687E",
          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
# same y-axis
ggPD_boot_test(grass.map,predictor = "Root Diameter",list.4.preds=grass.prerun, col.line="#6E687E",
               booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(grass.map,predictor = "SLA",list.4.preds=grass.prerun, col.line="#6E687E",
          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(grass.map,predictor = "SLA",list.4.preds=grass.prerun, col.line="#6E687E",
               booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(grass.map,predictor = "Rooting Depth",list.4.preds=grass.prerun, col.line="#6E687E",
          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(grass.map,predictor = "Rooting Depth",list.4.preds=grass.prerun, col.line="#6E687E",
               booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(grass.map,list.4.preds=grass.prerun, col.line="#6E687E",
          booted.preds=grass.boot$function.preds, cex.line=1, col.ci="#6E687E",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(grass.map, col.line="#6E687E", common.scale = FALSE, y.label = "Percent Cover Change")

ggInteract_list(grass.map)
# SLA x height 5.64
# precip x SLA 3.83
# rootN x height 3.02
# RTD x SLA 1.32

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_brt2.grass<-ggInteract_boot(c('SLA_m2.kg','height.m'),c('precip','SLA_m2.kg'),
                                          c('rootN.mg.g','height.m '),c('RTD.g.cm3','SLA_m2.kg'),
                                          nboots = 500, data=grass, predictors =  c(10:17,22), 
                                          response="cover.change",
                                          family = "gaussian", tc = 2, lr = 0.001, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 2,obs = 5.64) # significant
# larger SLA x taller height also for ALL data
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 3,obs = 3.83) # significant
# larger SLA x more precip
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 4,obs = 3.02) # significant
# higher rootN x taller plants
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 5,obs = 1.32) # non-significant


ggInteract_2D(gbm.object = grass.map, x="SLA_m2.kg",y="height.m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "height.m",
              z.range = c(-1, 3), z.label = "% Cover Change")

ggInteract_2D(gbm.object = grass.map, x="SLA_m2.kg",y="precip",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "precip",
              z.range = c(-1, 2.5), z.label = "% Cover Change")

ggInteract_2D(gbm.object = grass.map, x="rootN.mg.g",y="height.m",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "rootN.mg.g", y.label = "height.m",
              z.range = c(-1.5, 1.5), z.label = "% Cover Change")


ggInteract_list(grass.map.dsi)
# precip x SLA 6.13
# rootN x height 4.54
# SLA x height 3.39
# RTD x SLA 1.80
# DSI x precip 1.72

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_brt2.grass<-ggInteract_boot(c('precip','SLA_m2.kg'),c('rootN.mg.g ','SLA_m2.kg'),
                                          c('SLA_m2.kg','height.m '),c('RTD.g.cm3','SLA_m2.kg'),
                                          c('mean.drt.sev.index','precip'),
                                          nboots = 500, data=grass, predictors =  c(10:17,22,23), 
                                          response="cover.change",
                                          family = "gaussian", tc = 3, lr = 0.001, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 2,obs = 6.13) # significant
# larger SLA x more precip
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 3,obs = 4.54) # non-significant
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 4,obs = 3.39) # non-significant
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 5,obs = 1.80) # non-significant
ggInteract_boot_hist(data = Interact_boot_brt2.grass, column = 6,obs = 1.80) # non-significant



ggInteract_2D(gbm.object = grass.map.dsi, x="SLA_m2.kg",y="precip",col.gradient = c("white","#979461"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "SLA_m2.kg", y.label = "precip",
              z.range = c(-1, 2.5), z.label = "% Cover Change")

#### Forb ####
forb.map.dsi=gbm.step(data=forb, gbm.x = c(10:17,22,23), gbm.y=9,
                      family = "gaussian", tree.complexity = 7, learning.rate = 0.0001,
                      bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.map.dsi)
# 1000 trees Per.Expl = 1.48%

forb.map=gbm.step(data=forb, gbm.x = c(10:17,22), gbm.y=9,
                  family = "gaussian", tree.complexity = 3, learning.rate = 0.0005,
                  bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(forb.map)
# 1200 trees Per.Expl = 4.15%

forb.map$contributions$var = c("Height","LeafN","SLA","Root Diameter","Precipitation","SRL","RTD","Rooting Depth", "RootN")
forb.influence.plot = ggInfluence(forb.map, main = "Forbs (n = 314)", 
                                  col.bar = c("#F17236","#F17236","#F17236","#F17236",
                                              "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")
forb.influence.plot = ggInfluence(forb.map.dsi, main = "Forbs (n = 314)", 
                                  col.bar = c("#F17236","#F17236","#F17236","gray70",
                                              "gray70","gray70","gray70","gray70","gray70","gray70"
                                  ), col.signif = "#B50200")



ggInfluence(forb.map)

ggPD(forb.map, rug = T) # partial dependency plots
ggPDfit(forb.map)
gbm.plot(forb.map, common.scale = FALSE)
gbm.plot(forb.map.dsi, common.scale = FALSE)
gbm.plot.fits(forb.map)

forb.prerun<- plot.gbm.4list(forb.map)
forb.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","rootDiam.mm","precip", "SRL_m.g",
                               "RTD.g.cm3", "root.depth_m","rootN.mg.g")
forb.boot <- gbm.bootstrap.functions(forb.map, list.predictors=forb.prerun, n.reps=1000)

forb.map$contributions$var = c("Height","LeafN","SLA","Root Diameter","Precipitation","SRL","RTD","Rooting Depth","RootN")
forb.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(forb.map$gbm.call$dataframe)[12:19] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.map$gbm.call$dataframe)[24] = "Precipitation"

ggPD_boot(forb.map,predictor = "Height",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
# same y-axis
ggPD_boot_test(forb.map,predictor = "Height",list.4.preds=forb.prerun, col.line="#F17236",
               booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(forb.map,predictor = "LeafN",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(forb.map,predictor = "LeafN",list.4.preds=forb.prerun, col.line="#F17236",
               booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(forb.map,predictor = "SLA",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(forb.map,predictor = "SLA",list.4.preds=forb.prerun, col.line="#F17236",
               booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(forb.map,list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(forb.map, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")


# investigation of interactions
forb.map$contributions$var = c("height.m","leafN.mg.g","SLA_m2.kg","rootDiam.mm","precip", "SRL_m.g",
                               "RTD.g.cm3","root.depth_m","rootN.mg.g")
forb.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(forb.map$gbm.call$dataframe)[12:19] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(forb.map$gbm.call$dataframe)[24] = "precip"

gbm.interactions(forb.map)$rank.list
ggInteract_list(forb.map)
# rootDiam x SLA 27.19
# root.depth x SLA 11.20
# rootDiam x leafN 6.18
# precip x rootDiam 6.07

ggInteract_3D(forb.map, x = 5, y = 4, z.range = c(-0.5, 3))
ggInteract_3D(forb.map, x = 4, y = 2, z.range = c(-2.2, 2.5))
ggInteract_3D(forb.map, x = 9, y = 4, z.range = c(-0.75, 2))
ggInteract_3D(forb.map, x = 4, y = 1, z.range = c(-1.5, 3.75))

# bootstrap interactions
# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_forb<-ggInteract_boot(c('rootDiam.mm','SLA_m2.kg'),c('root.depth_m ','SLA_m2.kg'),
                                    c('rootDiam.mm','leafN.mg.g'),c('precip','rootDiam.mm'),
                                    nboots = 500, data=forb, predictors =  c(10:17,22), 
                                    response="cover.change",
                                    family = "gaussian", tc = 4, lr = 0.001, bf= 0.75, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_forb, column = 2,obs = 27.19) # significant
# higher SLA x lower root Diam
ggInteract_boot_hist(data = Interact_boot_forb, column = 3,obs = 11.20) # significant
# higher SLA x higher rootong depth
ggInteract_boot_hist(data = Interact_boot_forb, column = 4,obs = 6.18) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 5,obs = 6.07) # non-significant

ggInteract_2D(gbm.object = forb.map, x="rootDiam.mm",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "rootDiam.mm", y.label = "SLA_m2.kg",
              z.range = c(-3.5, 1), z.label = "% Cover Change")
ggInteract_2D(gbm.object = forb.map, x="root.depth_m",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "root.depth_m", y.label = "SLA_m2.kg",
              z.range = c(-3, 0.5), z.label = "% Cover Change")

ggInteract_list(forb.map.dsi)
# rootDiam x SLA 38.12
# root.depth x SLA 9.96
# SLA x leafN 9.42
# DSI x rootDiam 5.34
# precip x rootDiam 5.01

# bootstrap interactions
# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_forb<-ggInteract_boot(c('rootDiam.mm','SLA_m2.kg'),c('root.depth_m ','SLA_m2.kg'),
                                    c('SLA_m2.kg','leafN.mg.g'),c('mean.drt.sev.index','rootDiam.mm'),
                                    c('precip','rootDiam.mm'),
                                    nboots = 500, data=forb, predictors =  c(10:17,22,23), 
                                    response="cover.change",
                                    family = "gaussian", tc = 4, lr = 0.001, bf= 0.75, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_forb, column = 2,obs = 38.12) # significant
# higher SLA x lower root Diam
ggInteract_boot_hist(data = Interact_boot_forb, column = 3,obs = 9.96) # significant
# higher SLA x higher rootong depth
ggInteract_boot_hist(data = Interact_boot_forb, column = 4,obs = 9.42) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 5,obs = 5.34) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 6,obs = 5.01) # non-significant


ggInteract_2D(gbm.object = forb.map, x="rootDiam.mm",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "rootDiam.mm", y.label = "SLA_m2.kg",
              z.range = c(-3.5, 1), z.label = "% Cover Change")
ggInteract_2D(gbm.object = forb.map, x="root.depth_m",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "root.depth_m", y.label = "SLA_m2.kg",
              z.range = c(-3, 0.5), z.label = "% Cover Change")



ggInteract_list(forb.map.dsi)
# rootDiam x SLA 27.19
# root.depth x SLA 11.20
# rootDiam x leafN 6.18
# precip x rootDiam 6.07

# bootstrap interactions
# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_forb<-ggInteract_boot(c('rootDiam.mm','SLA_m2.kg'),c('root.depth_m ','SLA_m2.kg'),
                                    c('rootDiam.mm','leafN.mg.g'),c('precip','rootDiam.mm'),
                                    nboots = 500, data=forb, predictors =  c(10:17,22,23), 
                                    response="cover.change",
                                    family = "gaussian", tc = 4, lr = 0.001, bf= 0.75, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_forb, column = 2,obs = 27.19) # significant
# higher SLA x lower root Diam
ggInteract_boot_hist(data = Interact_boot_forb, column = 3,obs = 11.20) # significant
# higher SLA x higher rootong depth
ggInteract_boot_hist(data = Interact_boot_forb, column = 4,obs = 6.18) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb, column = 5,obs = 6.07) # non-significant

ggInteract_2D(gbm.object = forb.map.dsi, x="rootDiam.mm",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "rootDiam.mm", y.label = "SLA_m2.kg",
              z.range = c(-3.5, 1), z.label = "% Cover Change")

ggInteract_2D(gbm.object = forb.map.dsi, x="root.depth_m",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = T,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "root.depth_m", y.label = "SLA_m2.kg",
              z.range = c(-3.5, 0), z.label = "% Cover Change")


#### Grass.Perennial ####
grass.perennial.map.dsi=gbm.step(data=grass.perennial, gbm.x = c(10:17,22,23), gbm.y=9,
                                 family = "gaussian", tree.complexity = 4, learning.rate = 0.0005,
                                 bag.fraction = 0.75, n.trees = 50, verbose = TRUE, step.size = 25)

ggPerformance(grass.perennial.map.dsi)
# 1200 trees Per.Expl = 10.07%

grass.perennial.map=gbm.step(data=grass.perennial, gbm.x = c(10:17,22), gbm.y=9,
                             family = "gaussian", tree.complexity = 2, learning.rate = 0.0005,
                             bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(grass.perennial.map)
# 1550 trees Per.Expl = 6.97%

grass.perennial.map$contributions$var = c("Root Diameter","Rooting Depth","SLA","LeafN","Height","RTD","SRL","RootN","Precipitation")
grass.influence.plot = ggInfluence(grass.perennial.map, main = "Graminoids (n = 244)", 
                                   col.bar = c("#F17236","#F17236","gray70","gray70",
                                               "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")

grass.influence.plot = ggInfluence(grass.perennial.map.dsi, main = "Graminoids (n = 244)", 
                                   col.bar = c("#F17236","#F17236","gray70","gray70","gray70",
                                               "gray70","gray70","gray70","gray70","gray70"),
                                   col.signif = "#B50200")

ggInfluence(grass.perennial.map)

ggPD(graminoid.map, rug = T) # partial dependency plots
ggPDfit(grass.perennial.map)
gbm.plot(grass.perennial.map, common.scale = FALSE)

gbm.plot(grass.perennial.map.dsi, common.scale = FALSE)

gbm.plot.fits(grass.perennial.map)

grass.perennial.prerun<- plot.gbm.4list(grass.perennial.map)
grass.perennial.map$contributions$var = c("rootDiam.mm","root.depth_m","SLA_m2.kg","leafN.mg.g","height.m","RTD.g.cm3","SRL_m.g",
                                          "rootN.mg.g","precip")
grass.perennial.boot <- gbm.bootstrap.functions(grass.perennial.map, list.predictors=grass.perennial.prerun, n.reps=1000)

graminoid.map$contributions$var = c("Root Diameter","Rooting Depth","SLA","LeafN","Height","RTD","SRL","RootN","Precipitation")
graminoid.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(graminoid.map$gbm.call$dataframe)[8:15] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(graminoid.map$gbm.call$dataframe)[20] = "Precipitation"

ggPD_boot(graminoid.map,predictor = "Height",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")
# same y-axis
ggPD_boot_test(forb.map,predictor = "Height",list.4.preds=forb.prerun, col.line="#F17236",
               booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD_boot(forb.map,predictor = "LeafN",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(forb.map,predictor = "LeafN",list.4.preds=forb.prerun, col.line="#F17236",
               booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(forb.map,predictor = "SLA",list.4.preds=forb.prerun, col.line="#F17236",
          booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")
# same y-axis
ggPD_boot_test(forb.map,predictor = "SLA",list.4.preds=forb.prerun, col.line="#F17236",
               booted.preds=forb.boot$function.preds, cex.line=1, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "")

ggPD_boot(grass.perennial.map,list.4.preds=grass.perennial.prerun, col.line="#F17236",
          booted.preds=grass.perennial.boot$function.preds, cex.line=1, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,y.label = "Percent Cover Change")

ggPD(grass.perennial.map, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")

ggInteract_list(grass.perennial.map)
# rootDiam x SLA 3.47
# precip x SLA 1.46
# SLA x height 0.74
# height x leafN 0.72

# bootstrap interactions
# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_grass_perennial<-ggInteract_boot(c('rootDiam.mm','SLA_m2.kg'),c('precip ','SLA_m2.kg'),
                                               c('SLA_m2.kg','height.m'),c('height.m','leafN.mg.g'),
                                               nboots = 500, data=grass.perennial, predictors =  c(10:17,22), 
                                               response="cover.change",
                                               family = "gaussian", tc = 2, lr = 0.001, bf= 0.5, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_grass_perennial, column = 2,obs = 3.47) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial, column = 3,obs = 1.46) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial, column = 4,obs = 0.74) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial, column = 5,obs = 0.72) # non-significant


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
mean(all.data$mean.cover.response)

sd(all.data$leafN.mg.g, na.rm = TRUE)
sd(all.data$height.m, na.rm = TRUE)
sd(all.data$rootN.mg.g, na.rm = TRUE)
sd(all.data$SLA_m2.kg, na.rm = TRUE)
sd(all.data$root.depth_m, na.rm = TRUE)
sd(all.data$RTD.g.cm3, na.rm = TRUE)
sd(all.data$SRL_m.g, na.rm = TRUE)
sd(all.data$rootDiam.mm, na.rm = TRUE)
sd(all.data$precip, na.rm = TRUE)
sd(all.data$mean.cover.response)

range(all.data$leafN.mg.g, na.rm = TRUE)
range(all.data$height.m, na.rm = TRUE)
range(all.data$rootN.mg.g, na.rm = TRUE)
range(all.data$SLA_m2.kg, na.rm = TRUE)
range(all.data$root.depth_m, na.rm = TRUE)
range(all.data$RTD.g.cm3, na.rm = TRUE)
range(all.data$SRL_m.g, na.rm = TRUE)
range(all.data$rootDiam.mm, na.rm = TRUE)
range(all.data$precip, na.rm = TRUE)
range(all.data$mean.cover.response)

sum(is.na(all.data$leafN.mg.g))/612*100
sum(is.na(all.data$height.m))/612*100
sum(is.na(all.data$rootN.mg.g))/612*100
sum(is.na(all.data$SLA_m2.kg))/612*100
sum(is.na(all.data$root.depth_m))/612*100
sum(is.na(all.data$RTD.g.cm3))/612*100
sum(is.na(all.data$SRL_m.g))/612*100
sum(is.na(all.data$rootDiam.mm))/612*100
sum(is.na(all.data$precip))/612*100
sum(is.na(all.data$mean.cover.response))/612*100

mean(annual.data$leafN.mg.g, na.rm = TRUE)
mean(annual.data$height.m, na.rm = TRUE)
mean(annual.data$rootN.mg.g, na.rm = TRUE)
mean(annual.data$SLA_m2.kg, na.rm = TRUE)
mean(annual.data$root.depth_m, na.rm = TRUE)
mean(annual.data$RTD.g.cm3, na.rm = TRUE)
mean(annual.data$SRL_m.g, na.rm = TRUE)
mean(annual.data$rootDiam.mm, na.rm = TRUE)
mean(annual.data$precip, na.rm = TRUE)
mean(annual.data$mean.cover.response, na.rm = TRUE)

sd(annual.data$leafN.mg.g, na.rm = TRUE)
sd(annual.data$height.m, na.rm = TRUE)
sd(annual.data$rootN.mg.g, na.rm = TRUE)
sd(annual.data$SLA_m2.kg, na.rm = TRUE)
sd(annual.data$root.depth_m, na.rm = TRUE)
sd(annual.data$RTD.g.cm3, na.rm = TRUE)
sd(annual.data$SRL_m.g, na.rm = TRUE)
sd(annual.data$rootDiam.mm, na.rm = TRUE)
sd(annual.data$precip, na.rm = TRUE)
sd(annual.data$mean.cover.response, na.rm = TRUE)

range(annual.data$leafN.mg.g, na.rm = TRUE)
range(annual.data$height.m, na.rm = TRUE)
range(annual.data$rootN.mg.g, na.rm = TRUE)
range(annual.data$SLA_m2.kg, na.rm = TRUE)
range(annual.data$root.depth_m, na.rm = TRUE)
range(annual.data$RTD.g.cm3, na.rm = TRUE)
range(annual.data$SRL_m.g, na.rm = TRUE)
range(annual.data$rootDiam.mm, na.rm = TRUE)
range(annual.data$precip, na.rm = TRUE)
range(annual.data$mean.cover.response, na.rm = TRUE)

mean(perennial.data$leafN.mg.g, na.rm = TRUE)
mean(perennial.data$height.m, na.rm = TRUE)
mean(perennial.data$rootN.mg.g, na.rm = TRUE)
mean(perennial.data$SLA_m2.kg, na.rm = TRUE)
mean(perennial.data$root.depth_m, na.rm = TRUE)
mean(perennial.data$RTD.g.cm3, na.rm = TRUE)
mean(perennial.data$SRL_m.g, na.rm = TRUE)
mean(perennial.data$rootDiam.mm, na.rm = TRUE)
mean(perennial.data$precip, na.rm = TRUE)
mean(perennial.data$mean.cover.response, na.rm = TRUE)

sd(perennial.data$leafN.mg.g, na.rm = TRUE)
sd(perennial.data$height.m, na.rm = TRUE)
sd(perennial.data$rootN.mg.g, na.rm = TRUE)
sd(perennial.data$SLA_m2.kg, na.rm = TRUE)
sd(perennial.data$root.depth_m, na.rm = TRUE)
sd(perennial.data$RTD.g.cm3, na.rm = TRUE)
sd(perennial.data$SRL_m.g, na.rm = TRUE)
sd(perennial.data$rootDiam.mm, na.rm = TRUE)
sd(perennial.data$precip, na.rm = TRUE)
sd(perennial.data$mean.cover.response, na.rm = TRUE)

range(perennial.data$leafN.mg.g, na.rm = TRUE)
range(perennial.data$height.m, na.rm = TRUE)
range(perennial.data$rootN.mg.g, na.rm = TRUE)
range(perennial.data$SLA_m2.kg, na.rm = TRUE)
range(perennial.data$root.depth_m, na.rm = TRUE)
range(perennial.data$RTD.g.cm3, na.rm = TRUE)
range(perennial.data$SRL_m.g, na.rm = TRUE)
range(perennial.data$rootDiam.mm, na.rm = TRUE)
range(perennial.data$precip, na.rm = TRUE)
range(perennial.data$mean.cover.response, na.rm = TRUE)

mean(grass$leafN.mg.g, na.rm = TRUE)
mean(grass$height.m, na.rm = TRUE)
mean(grass$rootN.mg.g, na.rm = TRUE)
mean(grass$SLA_m2.kg, na.rm = TRUE)
mean(grass$root.depth_m, na.rm = TRUE)
mean(grass$RTD.g.cm3, na.rm = TRUE)
mean(grass$SRL_m.g, na.rm = TRUE)
mean(grass$rootDiam.mm, na.rm = TRUE)
mean(grass$precip, na.rm = TRUE)
mean(grass$mean.cover.response, na.rm = TRUE)

sd(grass$leafN.mg.g, na.rm = TRUE)
sd(grass$height.m, na.rm = TRUE)
sd(grass$rootN.mg.g, na.rm = TRUE)
sd(grass$SLA_m2.kg, na.rm = TRUE)
sd(grass$root.depth_m, na.rm = TRUE)
sd(grass$RTD.g.cm3, na.rm = TRUE)
sd(grass$SRL_m.g, na.rm = TRUE)
sd(grass$rootDiam.mm, na.rm = TRUE)
sd(grass$precip, na.rm = TRUE)
sd(grass$mean.cover.response, na.rm = TRUE)

range(grass$leafN.mg.g, na.rm = TRUE)
range(grass$height.m, na.rm = TRUE)
range(grass$rootN.mg.g, na.rm = TRUE)
range(grass$SLA_m2.kg, na.rm = TRUE)
range(grass$root.depth_m, na.rm = TRUE)
range(grass$RTD.g.cm3, na.rm = TRUE)
range(grass$SRL_m.g, na.rm = TRUE)
range(grass$rootDiam.mm, na.rm = TRUE)
range(grass$precip, na.rm = TRUE)
range(grass$mean.cover.response, na.rm = TRUE)

mean(forb$leafN.mg.g, na.rm = TRUE)
mean(forb$height.m, na.rm = TRUE)
mean(forb$rootN.mg.g, na.rm = TRUE)
mean(forb$SLA_m2.kg, na.rm = TRUE)
mean(forb$root.depth_m, na.rm = TRUE)
mean(forb$RTD.g.cm3, na.rm = TRUE)
mean(forb$SRL_m.g, na.rm = TRUE)
mean(forb$rootDiam.mm, na.rm = TRUE)
mean(forb$precip, na.rm = TRUE)
mean(forb$mean.cover.response, na.rm = TRUE)

sd(forb$leafN.mg.g, na.rm = TRUE)
sd(forb$height.m, na.rm = TRUE)
sd(forb$rootN.mg.g, na.rm = TRUE)
sd(forb$SLA_m2.kg, na.rm = TRUE)
sd(forb$root.depth_m, na.rm = TRUE)
sd(forb$RTD.g.cm3, na.rm = TRUE)
sd(forb$SRL_m.g, na.rm = TRUE)
sd(forb$rootDiam.mm, na.rm = TRUE)
sd(forb$precip, na.rm = TRUE)
sd(forb$mean.cover.response, na.rm = TRUE)

range(forb$leafN.mg.g, na.rm = TRUE)
range(forb$height.m, na.rm = TRUE)
range(forb$rootN.mg.g, na.rm = TRUE)
range(forb$SLA_m2.kg, na.rm = TRUE)
range(forb$root.depth_m, na.rm = TRUE)
range(forb$RTD.g.cm3, na.rm = TRUE)
range(forb$SRL_m.g, na.rm = TRUE)
range(forb$rootDiam.mm, na.rm = TRUE)
range(forb$precip, na.rm = TRUE)
range(forb$mean.cover.response, na.rm = TRUE)

table(perennial.data$functional_group)
# 230 forbs, 180 grasses
table(annual.data$functional_group)
# 76 forbs, 39 grasses
