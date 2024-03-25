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

# Modified ggInfluence function to make text size larger on plots 
# edited function is script: ggInfluence_edit.R

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
                                   family = "gaussian", tree.complexity = 10, learning.rate = 0.0001,
                                   bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 25)
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
# 2000 trees Per.Expl = 7.54%

all.data.map=gbm.step(data=all.data, gbm.x = c(10:17,22), gbm.y=9,
                      family = "gaussian", tree.complexity = 8, learning.rate = 0.0005,
                      bag.fraction = 0.50, n.trees = 50, verbose = TRUE, step.size = 50)

ggPerformance(all.data.map)
# 1500 trees Per.Expl = 8.12%


ggInfluence(all.data.map.dsi)
all.data.map.dsi$contributions$var = c("LeafN","Height","SLA","Precipitation",
                                   "RTD","RootN","DSI",
                                   "Rooting Depth","SRL","Root Diameter")
all.data.dsi.influence.plot = ggInfluence_test(all.data.map.dsi, main = expression("All Species (n = 1204), R"^2*" = 7.54%"), 
            col.bar = c("#769370","#769370","#769370","gray70",
                        "gray70","gray70","gray70","gray70","gray70","gray70"), 
            col.signif = "#B50200")


all.data.map$contributions$var = c("LeafN","Height","SLA","Precipitation",
                                   "Rooting Depth","SRL","RTD",
                                   "Root Diameter","RootN")
all.data.influence.plot = ggInfluence_test(all.data.map, main = expression("All Species (n = 1204), R"^2*" = 8.12%"), 
                                      col.bar = c("#769370","#769370","#769370","gray70",
                                                  "gray70","gray70","gray70","gray70","gray70"), 
                                      col.signif = "#B50200")

# saved as pdf from device (4 x 5)

# get data to plot partial dependency plots

all.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg","precip",
                                       "RTD.g.cm3","rootN.mg.g","mean.drt.sev.index",
                                       "root.depth_m","SRL_m.g","rootDiam.mm")

all.data.map$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg","precip",
                                   "root.depth_m","SRL_m.g","RTD.g.cm3",
                                   "rootDiam.mm","rootN.mg.g")

all.data.dsi.prerun<- plot.gbm.4list(all.data.map.dsi)
all.data.map.prerun<- plot.gbm.4list(all.data.map)

all.data.dsi.boot <- gbm.bootstrap.functions(all.data.map.dsi, list.predictors=all.data.dsi.prerun, n.reps=1000)
all.data.map.boot <- gbm.bootstrap.functions(all.data.map, list.predictors=all.data.map.prerun, n.reps=1000)

all.data.map.dsi$contributions$var = c("LeafN","Height","SLA","Precipitation",
                                       "RTD","RootN","DSI",
                                       "Rooting Depth","SRL","Root Diameter")
all.data.map$contributions$var = c("LeafN","Height","SLA","Precipitation",
                                   "Rooting Depth","SRL","RTD",
                                   "Root Diameter","RootN")

all.data.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(all.data.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")


all.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(all.data.map$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(all.data.map$gbm.call$dataframe)[22] = "Precipitation"

# significant for both: leafN, height, SLA

ggPD_boot_test(all.data.map.dsi,predictor = "LeafN",list.4.preds=all.data.dsi.prerun, col.line="#769370",
          booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
          y.label = "Percent Cover Change")
ggPD_boot_test(all.data.map.dsi,predictor = "Height",list.4.preds=all.data.dsi.prerun, col.line="#769370",
               booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(all.data.map.dsi,predictor = "SLA",list.4.preds=all.data.dsi.prerun, col.line="#769370",
               booted.preds=all.data.dsi.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test(all.data.map,predictor = "LeafN",list.4.preds=all.data.map.prerun, col.line="#769370",
               booted.preds=all.data.map.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(all.data.map,predictor = "Height",list.4.preds=all.data.map.prerun, col.line="#769370",
               booted.preds=all.data.map.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(all.data.map,predictor = "SLA",list.4.preds=all.data.map.prerun, col.line="#769370",
               booted.preds=all.data.map.boot$function.preds, cex.line=2, col.ci="#769370",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")


ggPD(all.data.map.dsi, col.line="#769370", common.scale = FALSE, y.label = "Percent Cover Change")

ggPD(all.data.map, col.line="#769370", common.scale = FALSE, y.label = "Percent Cover Change")

# output individual plots as 3x3
# output all plots as one 6x6

# investigation of interactions

all.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg","precip",
                                       "RTD.g.cm3","rootN.mg.g","mean.drt.sev.index",
                                       "root.depth_m","SRL_m.g","rootDiam.mm")

all.data.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip","mean.drt.sev.index")
colnames(all.data.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(all.data.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

all.data.map$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg","precip",
                                   "root.depth_m","SRL_m.g","RTD.g.cm3",
                                   "rootDiam.mm","rootN.mg.g")

all.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(all.data.map$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(all.data.map$gbm.call$dataframe)[22] = "precip"

gbm.interactions(all.data.map.dsi)$interactions
ggInteract_list(all.data.map.dsi, index = T)
# 6 RTD x 2 height 6.46
# 6 RTD x 1 leafN 2.88
# 9 precip x 2 height 2.56
# 10 DSI x 2 height 2.10
# 9 precip x 1 leafN 1.92

gbm.interactions(all.data.map)$interactions
ggInteract_list(all.data.map, index = T)
# 6 RTD x 2 height 6.09
# 9 precip x 2 height 3.11
# 6 RTD x 1 leafN 2.52
# 7 SRL x 2 height 2.19

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi<-ggInteract_boot(c('RTD.g.cm3','height.m'),c('RTD.g.cm3','leafN.mg.g'),
                                        c('precip','height.m'),c('mean.drt.sev.inde','height.m'),
                                   c('precip','leafN.mg.g'),
                                        nboots = 500, data=all.data, predictors =  c(10:17,22,23), 
                                        response="cover.change",
                                        family = "gaussian", tc = 5, lr = 0.0005, bf= 0.50, global.env=F)

Interact_boot_map<-ggInteract_boot(c('RTD.g.cm3','height.m'),c('precip','height.m'),
                                        c('RTD.g.cm3','leafN.mg.g'),c('SRL_m.g','height.m'),
                                        nboots = 500, data=all.data, predictors =  c(10:17,22), 
                                        response="cover.change",
                                        family = "gaussian", tc = 8, lr = 0.0005, bf= 0.50, global.env=F)

# Significance histogram p-value< 0.05
ggInteract_boot_hist(data = Interact_boot_dsi, column = 2,obs = 6.46) # significant
# lower RTD x higher height
ggInteract_boot_hist(data = Interact_boot_dsi, column = 3,obs = 2.88) # significant
# lower RTD x lower leafN 
ggInteract_boot_hist(data = Interact_boot_dsi, column = 4,obs = 2.56) # significant
# higher precip x higher height
ggInteract_boot_hist(data = Interact_boot_dsi, column = 5,obs = 2.10) # significant
# higher DSI x higher height
ggInteract_boot_hist(data = Interact_boot_dsi, column = 6,obs = 1.92) # significant
# higher precip x higher leafN

ggInteract_2D_test(gbm.object = all.data.map.dsi, x="RTD.g.cm3",y="height.m",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "Height",
                   z.range = c(-0.63, 0.95), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="RTD.g.cm3",y="leafN.mg.g",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "LeafN",
                   z.range = c(-1, 0.14), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="precip",y="height.m",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "Height",
                   z.range = c(-0.51, 0.90), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="mean.drt.sev.index",y="height.m",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "DSI", y.label = "Height",
                   z.range = c(-0.23, 1.25), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map.dsi, x="precip",y="leafN.mg.g",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "LeafN",
                   z.range = c(-1.01, 0.12), z.label = "% Cover Change", smooth = "average")

ggInteract_boot_hist(data = Interact_boot_map, column = 2,obs = 6.09) # significant
# lower RTD x higher height
ggInteract_boot_hist(data = Interact_boot_map, column = 3,obs = 3.11) # significant
# higher precip x higher height
ggInteract_boot_hist(data = Interact_boot_map, column = 4,obs = 2.52) # significant
# lower RTD x lower leafN
ggInteract_boot_hist(data = Interact_boot_map, column = 5,obs = 2.19) # significant
# higher SRL x higher height

ggInteract_2D_test(gbm.object = all.data.map, x="RTD.g.cm3",y="height.m",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "Height",
                   z.range = c(-0.57, 0.97), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map, x="precip",y="height.m",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Precipitation", y.label = "Height",
                   z.range = c(-0.48, 0.90), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map, x="RTD.g.cm3",y="leafN.mg.g",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "RTD", y.label = "LeafN",
                   z.range = c(-0.81, 0.18), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = all.data.map, x="SRL_m.g",y="height.m",col.gradient = c("white","#979461"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "SRL", y.label = "Height",
                   z.range = c(-0.34, 0.91), z.label = "% Cover Change", smooth = "average")

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

ggInfluence(annual.no.site.map)

annual.no.site.map.dsi$contributions$var = c( "RootN","RTD","DSI","Precipitation",
                                              "LeafN","SRL","Rooting Depth",
                                              "Height", "SLA","Root Diameter")
annual.influence.plot = ggInfluence_test(annual.no.site.map.dsi, main = expression("Annuals (n = 284), R"^2*" = 4.08%"), 
                                    col.bar = c("#979461","#979461","#979461","gray70",
                                                "gray70","gray70","gray70","gray70","gray70","gray70"),
                                    col.signif = "#B50200")

annual.no.site.map$contributions$var = c("RootN","RTD","Precipitation",
                                          "Height","SRL","LeafN","Rooting Depth",
                                           "SLA","Root Diameter")
annual.influence.plot = ggInfluence_test(annual.no.site.map, main = expression("Annuals (n = 284), R"^2*" = 6.64%"), 
                                    col.bar = c("#979461","#979461","#979461","#979461",
                                                "gray70","gray70","gray70","gray70","gray70"),
                                    col.signif = "#B50200")

# saved as pdf from device (4 x 5)

# get data to plot partial dependency plots
annual.no.site.map.dsi$contributions$var = c("rootN.mg.g","RTD.g.cm3","mean.drt.sev.index",
                                             "precip","leafN.mg.g","SRL_m.g","root.depth_m",
                                             "height.m","SLA_m2.kg","rootDiam.mm")

annual.no.site.map$contributions$var = c("rootN.mg.g","RTD.g.cm3","precip",
                                         "height.m","SRL_m.g","leafN.mg.g",
                                         "root.depth_m","SLA_m2.kg","rootDiam.mm")

annual.dsi.prerun<- plot.gbm.4list(annual.no.site.map.dsi)
annual.map.prerun<- plot.gbm.4list(annual.no.site.map)

annual.dsi.boot <- gbm.bootstrap.functions(annual.no.site.map.dsi, list.predictors=annual.dsi.prerun, n.reps=1000)
annual.map.boot <- gbm.bootstrap.functions(annual.no.site.map, list.predictors=annual.map.prerun, n.reps=1000)

annual.no.site.map.dsi$contributions$var = c( "RootN","RTD","DSI","Precipitation",
                                              "LeafN","SRL","Rooting Depth",
                                              "Height", "SLA","Root Diameter")

annual.no.site.map$contributions$var = c("RootN","RTD","Precipitation",
                                         "Height","SRL","LeafN","Rooting Depth",
                                         "SLA","Root Diameter")

annual.no.site.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(annual.no.site.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.no.site.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

annual.no.site.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(annual.no.site.map$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(annual.no.site.map$gbm.call$dataframe)[22] = "Precipitation"

# significant for DSI: RootN, RTD, DSI
# significant for MAP: RootN, RTD, Precip, Height

ggPD_boot_test(annual.no.site.map.dsi,predictor = "RootN",list.4.preds=annual.dsi.prerun, col.line="#979461",
          booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
          y.label = "Percent Cover Change")
ggPD_boot_test(annual.no.site.map.dsi,predictor = "RTD",list.4.preds=annual.dsi.prerun, col.line="#979461",
               booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(annual.no.site.map.dsi,predictor = "DSI",list.4.preds=annual.dsi.prerun, col.line="#979461",
               booted.preds=annual.dsi.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")

ggPD_boot_test(annual.no.site.map,predictor = "RootN",list.4.preds=annual.map.prerun, col.line="#979461",
               booted.preds=annual.map.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(annual.no.site.map,predictor = "RTD",list.4.preds=annual.map.prerun, col.line="#979461",
               booted.preds=annual.map.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(annual.no.site.map,predictor = "Precipitation",list.4.preds=annual.map.prerun, col.line="#979461",
               booted.preds=annual.map.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(annual.no.site.map,predictor = "Height",list.4.preds=annual.map.prerun, col.line="#979461",
               booted.preds=annual.map.boot$function.preds, cex.line=2, col.ci="#979461",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")

ggPD(annual.no.site.map.dsi, col.line="#979461", common.scale = FALSE, y.label = "Percent Cover Change")
ggPD(annual.no.site.map, col.line="#979461", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions: only for MAP b/c tc = 1 for DSI

annual.no.site.map$contributions$var = c("rootN.mg.g","RTD.g.cm3","precip",
                                         "height.m","SRL_m.g","leafN.mg.g",
                                         "root.depth_m","SLA_m2.kg","rootDiam.mm")

annual.no.site.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(annual.no.site.map$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(annual.no.site.map$gbm.call$dataframe)[22] = "precip"


gbm.interactions(annual.no.site.map)$rank.list
ggInteract_list(annual.no.site.map)
# 9 precip x 2 height 0.21
# 9 precip x 1 leafN 0.04
# 7 SRL x 1 leafN 0.04
# 5 root.depth x 1 leafN 0.04

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_annual<-ggInteract_boot(c('precip','height.m'),c('precip','leafN.mg.g'),
                                      c('SRL_m.g','leafN.mg.g'),c('"root.depth_m"','leafN.mg.g'),
                                      nboots = 500, data=annual.data, predictors = c(10:17,22), 
                                      response="cover.change",
                                      family = "gaussian", tc = 2, lr = 0.0005, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_annual, column = 2,obs = 0.21) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 3,obs = 0.04) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 4,obs = 0.04) # non-significant
ggInteract_boot_hist(data = Interact_boot_annual, column = 5,obs = 0.04) # non-significant

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

ggInfluence(perennial.data.map.dsi)

perennial.data.map.dsi$contributions$var = c("LeafN","Height","SLA",
                                         "SRL","RootN","Precipitation",
                                         "Rooting Depth","Root Diameter","RTD","DSI")
perennial.influence.plot = ggInfluence_test(perennial.data.map.dsi, main = expression("Perennials (n = 869), R"^2*" = 1.70%"), 
                                       col.bar = c("#F1C646","#F1C646","#F1C646","gray70",
                                                   "gray70","gray70","gray70","gray70","gray70","gray70"
                                       ), col.signif = "#B50200")

perennial.data.map$contributions$var = c("LeafN","SLA","Height",
                                         "SRL","Rooting Depth","Precipitation",
                                         "RootN","RTD","Root Diameter")
perennial.influence.plot = ggInfluence_test(perennial.data.map, main = expression("Perennials (n = 869), R"^2*" = 6.55%"), 
                                       col.bar = c("#F1C646","#F1C646","#F1C646","gray70",
                                                   "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")

# get data to plot partial dependency plots

perennial.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg",
                                             "SRL_m.g","rootN.mg.g","precip",
                                             "root.depth_m","rootDiam.mm","RTD.g.cm3","mean.drt.sev.index")

perennial.data.map$contributions$var = c("leafN.mg.g","SLA_m2.kg","height.m",
                                         "SRL_m.g","root.depth_m","precip",
                                         "rootN.mg.g","RTD.g.cm3","rootDiam.mm")

perennial.dsi.prerun<- plot.gbm.4list(perennial.data.map.dsi)
perennial.map.prerun<- plot.gbm.4list(perennial.data.map)

perennial.dsi.boot <- gbm.bootstrap.functions(perennial.data.map.dsi, list.predictors=perennial.dsi.prerun, n.reps=1000)
perennial.map.boot <- gbm.bootstrap.functions(perennial.data.map, list.predictors=perennial.map.prerun, n.reps=1000)

perennial.data.map.dsi$contributions$var = c("LeafN","Height","SLA",
                                             "SRL","RootN","Precipitation",
                                             "Rooting Depth","Root Diameter","RTD","DSI")

perennial.data.map$contributions$var = c("LeafN","SLA","Height",
                                         "SRL","Rooting Depth","Precipitation",
                                         "RootN","RTD","Root Diameter")

perennial.data.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

perennial.data.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(perennial.data.map$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(perennial.data.map$gbm.call$dataframe)[22] = "Precipitation"

# significant for both: leafN, SLA, height

ggPD_boot_test(perennial.data.map.dsi,predictor = "LeafN",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
          booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
          y.label = "Percent Cover Change")
ggPD_boot_test(perennial.data.map.dsi,predictor = "SLA",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
               booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(perennial.data.map.dsi,predictor = "Height",list.4.preds=perennial.dsi.prerun, col.line="#F1C646",
               booted.preds=perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test(perennial.data.map,predictor = "LeafN",list.4.preds=perennial.map.prerun, col.line="#F1C646",
               booted.preds=perennial.map.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(perennial.data.map,predictor = "SLA",list.4.preds=perennial.map.prerun, col.line="#F1C646",
               booted.preds=perennial.map.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(perennial.data.map,predictor = "Height",list.4.preds=perennial.map.prerun, col.line="#F1C646",
               booted.preds=perennial.map.boot$function.preds, cex.line=2, col.ci="#F1C646",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")

ggPD(perennial.data.map.dsi, col.line="#F1C646", common.scale = FALSE, y.label = "Percent Cover Change")
ggPD(perennial.data.map, col.line="#F1C646", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions

perennial.data.map.dsi$contributions$var = c("leafN.mg.g","height.m","SLA_m2.kg",
                                             "SRL_m.g","rootN.mg.g","precip",
                                             "root.depth_m","rootDiam.mm","RTD.g.cm3","mean.drt.sev.index")

perennial.data.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip","mean.drt.sev.index")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

perennial.data.map$contributions$var = c("leafN.mg.g","SLA_m2.kg","height.m",
                                         "SRL_m.g","root.depth_m","precip",
                                         "rootN.mg.g","RTD.g.cm3","rootDiam.mm")

perennial.data.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm", "precip")
colnames(perennial.data.map$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(perennial.data.map$gbm.call$dataframe)[22] = "precip"

gbm.interactions(perennial.data.map.dsi)$rank.list
ggInteract_list(perennial.data.map.dsi)
# 7 SRL x 2 height 0.17
# 6 RTD x 1 leafN 0.09
# 7 SRL x 3 rootN 0.06
# 9 precip x 1 leafN 0.05
# 6 RTD x 2 height 0.04

gbm.interactions(perennial.data.map)$rank.list
ggInteract_list(perennial.data.map)
# 6 RTD x 2 height 1.33
# 2 height x 1 leafN 1.18
# 6 RTD x 1 leafN 1.15
# 9 precip x 1 leafN 0.92

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_perennial_dsi<-ggInteract_boot(c('SRL_m.g','height.m'),c('RTD.g.cm3','leafN.mg.g'),
                                         c('SRL_m.g','rootN.mg.g'),c('precip','leafN.mg.g'),
                                         c('RTD.g.cm3','height.m'),
                                         nboots = 500,data=perennial.data, predictors =  c(10:17,22,23), 
                                         response="cover.change",
                                         family = "gaussian", tc = 9, lr = 0.0001, bf= 0.75, global.env=F)

Interact_boot_perennial_map<-ggInteract_boot(c('RTD.g.cm3','height.m'),c('height.m','leafN.mg.g'),
                                             c('RTD.g.cm3','leafN.mg.g'),c('precip','leafN.mg.g'),
                                             nboots = 500,data=perennial.data, predictors =  c(10:17,22), 
                                             response="cover.change",
                                             family = "gaussian", tc = 8, lr = 0.0005, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 2,obs = 0.17) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 3,obs = 0.09) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 4,obs = 0.06) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 5,obs = 0.05) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_dsi, column = 6,obs = 0.04) # non-significant

ggInteract_boot_hist(data = Interact_boot_perennial_map, column = 2,obs = 1.33) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_map, column = 3,obs = 1.18) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_map, column = 4,obs = 1.15) # non-significant
ggInteract_boot_hist(data = Interact_boot_perennial_map, column = 5,obs = 0.92) # non-significant


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

ggInfluence(grass.map)

grass.map.dsi$contributions$var = c("LeafN", "Precipitation","Height",
                                    "SLA","DSI","Root Diameter","RootN",
                                    "RTD","SRL","Rooting Depth")
                                
grass.influence.plot = ggInfluence_test(grass.map.dsi, main = expression("Grasses (n = 297), R"^2*" = 3.51%"), 
                                   col.bar = c("#6E687E","#6E687E","#6E687E","#6E687E",
                                               "gray70","gray70","gray70","gray70","gray70","gray70"
                                   ), col.signif = "#B50200")

grass.map$contributions$var =  c("LeafN", "Precipitation","Height",
                                 "SLA","Root Diameter","RootN",
                                 "RTD","SRL","Rooting Depth")
grass.influence.plot = ggInfluence_test(grass.map, main = expression("Grasses (n = 297), R"^2*" = 1.77%"), 
                                   col.bar = c("#6E687E","#6E687E","#6E687E","#6E687E",
                                               "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")

# get data to plot partial dependency plots

grass.map.dsi$contributions$var = c("leafN.mg.g", "precip","height.m",
                                    "SLA_m2.kg","mean.drt.sev.index","rootDiam.mm","rootN.mg.g",
                                    "RTD.g.cm3","SRL_m.g","root.depth_m")

grass.map$contributions$var =  c("leafN.mg.g", "precip","height.m",
                                 "SLA_m2.kg","rootDiam.mm","rootN.mg.g",
                                 "RTD.g.cm3","SRL_m.g","root.depth_m")


grass.dsi.prerun<- plot.gbm.4list(grass.map.dsi)
grass.map.prerun<- plot.gbm.4list(grass.map)

grass.dsi.boot <- gbm.bootstrap.functions(grass.map.dsi, list.predictors=grass.dsi.prerun, n.reps=1000)
grass.map.boot <- gbm.bootstrap.functions(grass.map, list.predictors=grass.map.prerun, n.reps=1000)

grass.map.dsi$contributions$var = c("LeafN", "Precipitation","Height",
                                    "SLA","DSI","Root Diameter","RootN",
                                    "RTD","SRL","Rooting Depth")

grass.map$contributions$var =  c("LeafN", "Precipitation","Height",
                                 "SLA","Root Diameter","RootN",
                                 "RTD","SRL","Rooting Depth")

grass.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(grass.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

grass.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(grass.map$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.map$gbm.call$dataframe)[22] = "Precipitation"

# significant for both: leafN, precip, height, SLA

ggPD_boot_test(grass.map.dsi,predictor = "LeafN",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map.dsi,predictor = "Precipitation",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map.dsi,predictor = "Height",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map.dsi,predictor = "SLA",list.4.preds=grass.dsi.prerun, col.line="#6E687E",
               booted.preds=grass.dsi.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test(grass.map,predictor = "LeafN",list.4.preds=grass.map.prerun, col.line="#6E687E",
               booted.preds=grass.map.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map,predictor = "Precipitation",list.4.preds=grass.map.prerun, col.line="#6E687E",
               booted.preds=grass.map.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map,predictor = "Height",list.4.preds=grass.map.prerun, col.line="#6E687E",
               booted.preds=grass.map.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.map,predictor = "SLA",list.4.preds=grass.map.prerun, col.line="#6E687E",
               booted.preds=grass.map.boot$function.preds, cex.line=2, col.ci="#6E687E",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")

ggPD(grass.map.dsi, col.line="#6E687E", common.scale = FALSE, y.label = "Percent Cover Change")
ggPD(grass.map, col.line="#6E687E", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions

grass.map.dsi$contributions$var = c("leafN.mg.g", "precip","height.m",
                                    "SLA_m2.kg","mean.drt.sev.index","rootDiam.mm","rootN.mg.g",
                                    "RTD.g.cm3","SRL_m.g","root.depth_m")

grass.map$contributions$var =  c("leafN.mg.g", "precip","height.m",
                                 "SLA_m2.kg","rootDiam.mm","rootN.mg.g",
                                 "RTD.g.cm3","SRL_m.g","root.depth_m")

grass.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip","mean.drt.sev.index")
colnames(grass.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

grass.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(grass.map$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.map$gbm.call$dataframe)[22] = "precip"

gbm.interactions(grass.map.dsi)$interactions
ggInteract_list(grass.map.dsi, index = T)
# 9 precip x 2 height 0.29
# 10 DSI x 2 height 0.23
# 9 precip x 8 rootDiam 0.12
# 8 rootDiam x 1 leafN 0.11
# 8 rootDiam x 3 rootN 0.08

gbm.interactions(grass.map)$interactions
ggInteract_list(grass.map, index = T)
# 9 precip x 2 height 0.08
# 8 rootDiam x 1 leafN 0.06
# 9 precip x 8 rootDiam 0.04
# 9 precip x 4 SLA 0.02

# bootstrap interactions

# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model
Interact_boot_dsi.grass<-ggInteract_boot(c('precip','height.m'),c('mean.drt.sev.index','height.m'),
                                          c('precip','rootDiam.mm  '),c('rootDiam.mm ','leafN.mg.g'),
                                         c('rootDiam.mm ','rootN.mg.g'),
                                          nboots = 500, data=grass, predictors =  c(10:17,22,23), 
                                          response="cover.change",
                                          family = "gaussian", tc = 6, lr = 0.0001, bf= 0.50, global.env=F)

Interact_boot_map.grass<-ggInteract_boot(c('precip','height.m'),c('rootDiam.mm','leafN.mg.g'),
                                         c('precip','rootDiam.mm  '),c('precip ','SLA_m2.kg'),
                                         nboots = 500, data=grass, predictors =  c(10:17,22), 
                                         response="cover.change",
                                         family = "gaussian", tc = 5, lr = 0.0001, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 2,obs = 0.29) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 3,obs = 0.23) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 4,obs = 0.12) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 5,obs = 0.11) # non-significant
ggInteract_boot_hist(data = Interact_boot_dsi.grass, column = 6,obs = 0.08) # non-significant

ggInteract_boot_hist(data = Interact_boot_map.grass, column = 2,obs = 0.08) # non-significant
ggInteract_boot_hist(data = Interact_boot_map.grass, column = 3,obs = 0.06) # non-significant
ggInteract_boot_hist(data = Interact_boot_map.grass, column = 4,obs = 0.04) # non-significant
ggInteract_boot_hist(data = Interact_boot_map.grass, column = 5,obs = 0.02) # non-significant

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


ggInfluence(forb.map)

forb.map.dsi$contributions$var = c("LeafN","SLA","Rooting Depth",
                                   "DSI","Height","SRL",
                                   "Precipitation","Root Diameter",
                                   "RTD", "RootN")

forb.influence.plot = ggInfluence_test(forb.map.dsi, main = expression("Forbs (n = 666), R"^2*" = 1.48%"),
                                  col.bar = c("#F17236","#F17236","#F17236","#F17236",
                                              "gray70","gray70","gray70","gray70","gray70","gray70"
                                  ), col.signif = "#B50200")

forb.map$contributions$var = c("LeafN","SLA","Rooting Depth",
                               "Root Diameter","SRL",
                               "RTD","Height","Precipitation",
                                "RootN")
forb.influence.plot = ggInfluence_test(forb.map, main = expression("Forbs (n = 666), R"^2*" = 4.15%"), 
                                  col.bar = c("#F17236","#F17236","#F17236","gray70",
                                              "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")

# get data to plot partial dependency plots

forb.map.dsi$contributions$var = c("leafN.mg.g","SLA_m2.kg","root.depth_m",
                                   "mean.drt.sev.index","height.m","SRL_m.g",
                                   "precip","rootDiam.mm",
                                   "RTD.g.cm3", "rootN.mg.g")

forb.map$contributions$var = c("leafN.mg.g","SLA_m2.kg","root.depth_m",
                               "rootDiam.mm","SRL_m.g",
                               "RTD.g.cm3","height.m","precip",
                               "rootN.mg.g")

forb.dsi.prerun<- plot.gbm.4list(forb.map.dsi)
forb.map.prerun<- plot.gbm.4list(forb.map)

forb.dsi.boot <- gbm.bootstrap.functions(forb.map.dsi, list.predictors=forb.dsi.prerun, n.reps=1000)
forb.map.boot <- gbm.bootstrap.functions(forb.map, list.predictors=forb.map.prerun, n.reps=1000)

forb.map.dsi$contributions$var = c("LeafN","SLA","Rooting Depth",
                                   "DSI","Height","SRL",
                                   "Precipitation","Root Diameter",
                                   "RTD", "RootN")

forb.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(forb.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

forb.map$contributions$var = c("LeafN","SLA","Rooting Depth",
                               "Root Diameter","SRL",
                               "RTD","Height","Precipitation",
                               "RootN")

forb.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(forb.map$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(forb.map$gbm.call$dataframe)[22] = "Precipitation"

# significant for DSI: leafN, SLA, depth, DSI
# significant for MAP: leafN, SLA, depth

ggPD_boot_test(forb.map.dsi,predictor = "LeafN",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map.dsi,predictor = "SLA",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map.dsi,predictor = "Rooting Depth",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map.dsi,predictor = "DSI",list.4.preds=forb.dsi.prerun, col.line="#F17236",
               booted.preds=forb.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test(forb.map,predictor = "LeafN",list.4.preds=forb.map.prerun, col.line="#F17236",
               booted.preds=forb.map.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map,predictor = "SLA",list.4.preds=forb.map.prerun, col.line="#F17236",
               booted.preds=forb.map.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(forb.map,predictor = "Rooting Depth",list.4.preds=forb.map.prerun, col.line="#F17236",
               booted.preds=forb.map.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T, common.scale = FALSE,
               y.label = "Percent Cover Change")

ggPD(forb.map.dsi, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")
ggPD(forb.map, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions

forb.map.dsi$contributions$var = c("leafN.mg.g","SLA_m2.kg","root.depth_m",
                                   "mean.drt.sev.index","height.m","SRL_m.g",
                                   "precip","rootDiam.mm",
                                   "RTD.g.cm3", "rootN.mg.g")

forb.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip","mean.drt.sev.index")
colnames(forb.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(forb.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")

forb.map$contributions$var = c("leafN.mg.g","SLA_m2.kg","root.depth_m",
                               "rootDiam.mm","SRL_m.g",
                               "RTD.g.cm3","height.m","precip",
                               "rootN.mg.g")

forb.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(forb.map$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(forb.map$gbm.call$dataframe)[22] = "precip"

gbm.interactions(forb.map.dsi)$interactions
ggInteract_list(forb.map.dsi, index = T)
# 4 SLA x 1 leafN 0.06
# 5 depth x 4 SLA 0.04
# 10 DSI x 5 depth 0.01
# 10 DSI x 4 SLA 0.01
# 10 DSI x 2 height 0.01

gbm.interactions(forb.map)$interactions
ggInteract_list(forb.map, index = T)
# 4 SLA x 1 leafN 1.06
# 5 depth x 4 SLA 0.93
# 6 RTD x 1 leafN 0.09
# 8 Diam x 4 SLA 0.08

# bootstrap interactions
# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_forb_dsi<-ggInteract_boot(c('SLA_m2.kg','leafN.mg.g'),c('root.depth_m ','SLA_m2.kg'),
                                    c('mean.drt.sev.index','root.depth_m'),c('mean.drt.sev.index','SLA_m2.kg'),
                                    c('mean.drt.sev.index','height.m'),
                                    nboots = 500, data=forb, predictors =  c(10:17,22,23), 
                                    response="cover.change",
                                    family = "gaussian", tc = 7, lr = 0.0001, bf= 0.50, global.env=F)

Interact_boot_forb_map<-ggInteract_boot(c('SLA_m2.kg','leafN.mg.g'),c('root.depth_m ','SLA_m2.kg'),
                                        c('RTD.g.cm3','leafN.mg.g'),c('rootDiam.mm','SLA_m2.kg'),
                                        nboots = 500, data=forb, predictors =  c(10:17,22), 
                                        response="cover.change",
                                        family = "gaussian", tc = 3, lr = 0.0005, bf= 0.50, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 2,obs = 0.06) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 3,obs = 0.04) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 4,obs = 0.01) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 5,obs = 0.01) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_dsi, column = 6,obs = 0.01) # non-significant

ggInteract_boot_hist(data = Interact_boot_forb_map, column = 2,obs = 1.06) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_map, column = 3,obs = 0.93) # significant
# 5 depth x 4 SLA 0.93, higher rooting depth x higher SLA
ggInteract_boot_hist(data = Interact_boot_forb_map, column = 4,obs = 0.09) # non-significant
ggInteract_boot_hist(data = Interact_boot_forb_map, column = 5,obs = 0.08) # non-significant

ggInteract_2D_test(gbm.object = forb.map, x="root.depth_m",y="SLA_m2.kg",col.gradient = c("white","#F17236"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
              col.contour = "#254376",show.axis = T,legend = T, x.label = "Rooting Depth", y.label = "SLA",
              z.range = c(-0.19, 0.44), z.label = "% Cover Change", smooth = "average")

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

ggInfluence(grass.perennial.map)

grass.perennial.map.dsi$contributions$var = c("LeafN","Height","Root Diameter",
                                              "SLA","Precipitation","DSI",
                                              "SRL","RootN","RTD","Rooting Depth")
                                              
grass.influence.plot = ggInfluence_test(grass.perennial.map.dsi, main = expression("Perennial Grasses (n = 297), R"^2*" = 10.07%"), 
                                   col.bar = c("#F17236","#F17236","#F17236","#F17236","gray70",
                                               "gray70","gray70","gray70","gray70","gray70"),
                                   col.signif = "#B50200")

grass.perennial.map$contributions$var = c("SLA","LeafN","Height","Root Diameter",
                                          "Precipitation",
                                          "SRL","RootN","RTD","Rooting Depth")
grass.influence.plot = ggInfluence_test(grass.perennial.map, main = expression("Perennial Grasses (n = 297), R"^2*" = 6.97%"), 
                                   col.bar = c("#F17236","#F17236","#F17236","#F17236",
                                               "gray70","gray70","gray70","gray70","gray70"), col.signif = "#B50200")


# exported as 4 x 6 instead of 4 x 5 like the others

# get data to plot partial dependency plots

grass.perennial.map.dsi$contributions$var = c("leafN.mg.g","height.m","rootDiam.mm",
                                              "SLA_m2.kg","precip","mean.drt.sev.index",
                                              "SRL_m.g","rootN.mg.g","RTD.g.cm3","root.depth_m")

grass.perennial.map$contributions$var = c("SLA_m2.kg","leafN.mg.g","height.m","rootDiam.mm",
                                          "precip",
                                          "SRL_m.g","rootN.mg.g","RTD.g.cm3","root.depth_m")

grass.perennial.dsi.prerun<- plot.gbm.4list(grass.perennial.map.dsi)
grass.perennial.map.prerun<- plot.gbm.4list(grass.perennial.map)

grass.perennial.dsi.boot <- gbm.bootstrap.functions(grass.perennial.map.dsi, list.predictors=grass.perennial.dsi.prerun, n.reps=1000)
grass.perennial.map.boot <- gbm.bootstrap.functions(grass.perennial.map, list.predictors=grass.perennial.map.prerun, n.reps=1000)

grass.perennial.map.dsi$contributions$var = c("LeafN","Height","Root Diameter",
                                              "SLA","Precipitation","DSI",
                                              "SRL","RootN","RTD","Rooting Depth")

grass.perennial.map$contributions$var = c("SLA","LeafN","Height","Root Diameter",
                                          "Precipitation",
                                          "SRL","RootN","RTD","Rooting Depth")

grass.perennial.map.dsi$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation","DSI")
colnames(grass.perennial.map.dsi$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.perennial.map.dsi$gbm.call$dataframe)[22:23] = c("Precipitation","DSI")

grass.perennial.map$gbm.call$predictor.names = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter","Precipitation")
colnames(grass.perennial.map$gbm.call$dataframe)[10:17] = c("LeafN","Height","RootN","SLA","Rooting Depth","RTD","SRL","Root Diameter")
colnames(grass.perennial.map$gbm.call$dataframe)[22] = "Precipitation"

# significant for both: leafN, height, Root diameter SLA

ggPD_boot_test(grass.perennial.map.dsi,predictor = "LeafN",list.4.preds=grass.perennial.dsi.prerun, col.line="#F17236",
          booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
          alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
          y.label = "Percent Cover Change")
ggPD_boot_test(grass.perennial.map.dsi,predictor = "Height",list.4.preds=grass.perennial.dsi.prerun, col.line="#F17236",
               booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.perennial.map.dsi,predictor = "Root Diameter",list.4.preds=grass.perennial.dsi.prerun, col.line="#F17236",
               booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.perennial.map.dsi,predictor = "SLA",list.4.preds=grass.perennial.dsi.prerun, col.line="#F17236",
               booted.preds=grass.perennial.dsi.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = TRUE,
               y.label = "Percent Cover Change")

ggPD_boot_test(grass.perennial.map,predictor = "LeafN",list.4.preds=grass.perennial.map.prerun, col.line="#F17236",
               booted.preds=grass.perennial.map.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.perennial.map,predictor = "Height",list.4.preds=grass.perennial.map.prerun, col.line="#F17236",
               booted.preds=grass.perennial.map.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.perennial.map,predictor = "Root Diameter",list.4.preds=grass.perennial.map.prerun, col.line="#F17236",
               booted.preds=grass.perennial.map.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")
ggPD_boot_test(grass.perennial.map,predictor = "SLA",list.4.preds=grass.perennial.map.prerun, col.line="#F17236",
               booted.preds=grass.perennial.map.boot$function.preds, cex.line=2, col.ci="#F17236",
               alpha.dot=0.2,type.ci = "ribbon",alpha.ci= 0.3,rug = T,common.scale = FALSE,
               y.label = "Percent Cover Change")

ggPD(grass.perennial.map.dsi, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")
ggPD(grass.perennial.map, col.line="#F17236", common.scale = FALSE, y.label = "Percent Cover Change")

# output all plots as one 6x6

# investigation of interactions

grass.perennial.map.dsi$contributions$var = c("leafN.mg.g","height.m","rootDiam.mm",
                                              "SLA_m2.kg","precip","mean.drt.sev.index",
                                              "SRL_m.g","rootN.mg.g","RTD.g.cm3","root.depth_m")

grass.perennial.map$contributions$var = c("SLA_m2.kg","leafN.mg.g","height.m","rootDiam.mm",
                                          "precip",
                                          "SRL_m.g","rootN.mg.g","RTD.g.cm3","root.depth_m")

grass.perennial.map.dsi$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip","mean.drt.sev.index")
colnames(grass.perennial.map.dsi$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.perennial.map.dsi$gbm.call$dataframe)[22:23] = c("precip","mean.drt.sev.index")


grass.perennial.map$gbm.call$predictor.names = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm","precip")
colnames(grass.perennial.map$gbm.call$dataframe)[10:17] = c("leafN.mg.g","height.m","rootN.mg.g","SLA_m2.kg","root.depth_m","RTD.g.cm3","SRL_m.g","rootDiam.mm")
colnames(grass.perennial.map$gbm.call$dataframe)[22] = "precip"

gbm.interactions(grass.perennial.map.dsi)$interactions
ggInteract_list(grass.perennial.map.dsi, index = T)
# 8 rootDiam x 1 leafN 10.56
# 10 DSI x 2 height 4.11
# 2 height x 1 leafN 2.60
# 7 SRL x 1 leafN 0.94
# 9 precip x 6 RTD 0.82

gbm.interactions(grass.perennial.map)$interactions
ggInteract_list(grass.perennial.map, index = T)
# 8 diam x 1 leafN 0.61
# 8 diam x 4 SLA 0.20
# 8 diam x 2 height 0.12
# 7 SRL x 1 leafN 0.11


# bootstrap interactions
# Randomization of the response to test significance of the three strongest interactions. Note, the parameters need to be the same as your BRT model

Interact_boot_grass_perennial_dsi<-ggInteract_boot(c('rootDiam.mm','leafN.mg.g'),c('mean.drt.sev.index ','height.m'),
                                               c('height.m','leafN.mg.g'),c('SRL_m.g','leafN.mg.g'),
                                               c('precip','RTD.g.cm3'),
                                               nboots = 500, data=grass.perennial, predictors =  c(10:17,22,23), 
                                               response="cover.change",
                                               family = "gaussian", tc = 4, lr = 0.0005, bf= 0.75, global.env=F)

Interact_boot_grass_perennial_map<-ggInteract_boot(c('rootDiam.mm','leafN.mg.g'),c('rootDiam.mm ','SLA_m2.kg'),
                                                   c('rootDiam.mm','height.m'),c('SRL_m.g','leafN.mg.g'),
                                                   nboots = 500, data=grass.perennial, predictors =  c(10:17,22,23), 
                                                   response="cover.change",
                                                   family = "gaussian", tc = 2, lr = 0.0005, bf= 0.5, global.env=F)

# Significance histogram p-value<0.05)
ggInteract_boot_hist(data = Interact_boot_grass_perennial_dsi, column = 2,obs = 10.56) # significant
# highter Root diameter x higher leafN
ggInteract_boot_hist(data = Interact_boot_grass_perennial_dsi, column = 3,obs = 4.11) # significant
# higher drought, taller height
ggInteract_boot_hist(data = Interact_boot_grass_perennial_dsi, column = 4,obs = 2.60) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial_dsi, column = 5,obs = 0.94) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial_dsi, column = 6,obs = 0.82) # non-significant

ggInteract_boot_hist(data = Interact_boot_grass_perennial_map, column = 2,obs = 0.61) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial_map, column = 3,obs = 0.20) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial_map, column = 4,obs = 0.12) # non-significant
ggInteract_boot_hist(data = Interact_boot_grass_perennial_map, column = 5,obs = 0.11) # non-significant


ggInteract_2D_test(gbm.object = grass.perennial.map.dsi, x="rootDiam.mm",y="leafN.mg.g",col.gradient = c("white","#F17236"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "Root Diameter", y.label = "LeafN",
                   z.range = c(-1.11, 0.18), z.label = "% Cover Change", smooth = "average")
ggInteract_2D_test(gbm.object = grass.perennial.map.dsi, x="mean.drt.sev.index",y="height.m",col.gradient = c("white","#F17236"),
                   show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,label.contour = F,
                   col.contour = "#254376",show.axis = T,legend = T, x.label = "DSI", y.label = "Height",
                   z.range = c(-0.49, 0.74), z.label = "% Cover Change", smooth = "average")


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
