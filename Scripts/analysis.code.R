# Script for analysis of Drought Net data

library(rpart)
library(rpart.plot)
library(dismo)
library(gbm)

# read in test data framex.var.pred.frame


test.data = read.csv("./Formatted.Data/analysis.test.data.csv")
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")

# subset traits to only those of interest

trait.data = test.data[,c(1,7,8,10,12,14,15,18,20)]

# subset traits so they must have SLA

trait.data.2 = subset(trait.data, trait.data$SLA_m2.kg > 0 ) # 640 species

# merge trait data and cover data

all.data = merge(cover.data, trait.data.2, by="Taxon") # 1076 data points with site

# prior to running below code, we need to:
# test for correlation between predictors
# test for normality of data, if non-normal transform
# scale data for comparisons

# test for correlation
cor.test(all.data$leafN.mg.g, all.data$height.m)
cor.test(all.data$leafN.mg.g, all.data$rootN.mg.g) # correlated
cor.test(all.data$leafN.mg.g, all.data$SLA_m2.kg) # correlated
cor.test(all.data$leafN.mg.g, all.data$root.depth_m)
cor.test(all.data$leafN.mg.g, all.data$RTD.g.cm3) # correlated
cor.test(all.data$leafN.mg.g, all.data$SRL_m.g) # correlated
cor.test(all.data$leafN.mg.g, all.data$rootDiam.mm)

cor.test(all.data$height.m, all.data$rootN.mg.g)
cor.test(all.data$height.m, all.data$SLA_m2.kg)
cor.test(all.data$height.m, all.data$root.depth_m) # correlated
cor.test(all.data$height.m, all.data$RTD.g.cm3) # correlated
cor.test(all.data$height.m, all.data$SRL_m.g) # correlated
cor.test(all.data$height.m, all.data$rootDiam.mm) # correlated

cor.test(all.data$rootN.mg.g, all.data$SLA_m2.kg)
cor.test(all.data$rootN.mg.g, all.data$root.depth_m)
cor.test(all.data$rootN.mg.g, all.data$RTD.g.cm3) # correlated
cor.test(all.data$rootN.mg.g, all.data$SRL_m.g)
cor.test(all.data$rootN.mg.g, all.data$rootDiam.mm) # correlated

cor.test(all.data$SLA_m2.kg, all.data$root.depth_m)
cor.test(all.data$SLA_m2.kg, all.data$RTD.g.cm3) # correlated
cor.test(all.data$SLA_m2.kg, all.data$SRL_m.g) # correlated
cor.test(all.data$SLA_m2.kg, all.data$rootDiam.mm) # correlated

cor.test(all.data$root.depth_m, all.data$RTD.g.cm3) # correlated
cor.test(all.data$root.depth_m, all.data$SRL_m.g) # correlated
cor.test(all.data$root.depth_m, all.data$rootDiam.mm)

cor.test(all.data$RTD.g.cm3, all.data$SRL_m.g) # correlated
cor.test(all.data$RTD.g.cm3, all.data$rootDiam.mm)

cor.test(all.data$SRL_m.g, all.data$rootDiam.mm) # correlated

# uncorrelated trait combinations
# leafN, root.depth, rootDiam
# height, rootN, SLA
# rootN, SLA, root.depth, 

# test for normality
shapiro.test(all.data$leafN.mg.g)
hist(all.data$leafN.mg.g)
shapiro.test(all.data$height.m)
hist(all.data$height.m)
shapiro.test(all.data$rootN.mg.g)
hist(all.data$rootN.mg.g)
shapiro.test(all.data$SLA_m2.kg)
hist(all.data$SLA_m2.kg)
shapiro.test(all.data$root.depth_m)
hist(all.data$root.depth_m)
shapiro.test(all.data$RTD.g.cm3)
hist(all.data$RTD.g.cm3)
shapiro.test(all.data$SRL_m.g)
hist(all.data$SRL_m.g)
shapiro.test(all.data$rootDiam.mm)
hist(all.data$rootDiam.mm)

shapiro.test(all.data$mean.cover.response)


# log and scale values to transform if non-normal, , mean of 0, sd = 1
all.data$log.leafN = as.vector(scale(log(all.data$leafN.mg.g)))
all.data$log.height = as.vector(scale(log(all.data$height.m)))
all.data$log.rootN = as.vector(scale(log(all.data$rootN.mg.g)))
all.data$log.SLA = as.vector(scale(log(all.data$SLA_m2.kg)))
all.data$log.root.depth = as.vector(scale(log(all.data$root.depth_m)))
all.data$log.RTD = as.vector(scale(log(all.data$RTD.g.cm3)))
all.data$log.SRL = as.vector(scale(log(all.data$SRL_m.g)))
all.data$log.root.diam = as.vector(scale(log(all.data$rootDiam.mm)))

# get data set of only complete trait cases

all.complete = all.data[complete.cases(all.data),]

all.data$site.id = as.numeric(as.factor(all.data$site_code))
all.complete$site.id = as.numeric(as.factor(all.complete$site_code))

#### basic regression tree analysis ####
# page 84

# all variables
reg.tree.1=rpart(mean.cover.response ~ log.leafN+log.height+log.rootN+
                   log.SLA+log.root.depth+log.RTD+log.SRL+log.root.diam,data=all.data)
rpart.plot(reg.tree.1)
summary(reg.tree.1)
# RTD only important variable

reg.tree.1.site=rpart(mean.cover.response ~ log.leafN+log.height+log.rootN+
                   log.SLA+log.root.depth+log.RTD+log.SRL+log.root.diam,
                   weights = all.data$site.id, data=all.data)
rpart.plot(reg.tree.1.site)
summary(reg.tree.1.site)

# variable subset
reg.tree.2=rpart(mean.cover.response ~ log.leafN+log.root.depth+log.root.diam, data=all.data)
rpart.plot(reg.tree.2)
summary(reg.tree.2)
# no information

reg.tree.2.site=rpart(mean.cover.response ~ log.leafN+log.root.depth+log.root.diam, 
                      weights = all.data$site.id, data=all.data)
rpart.plot(reg.tree.2.site)
summary(reg.tree.2.site)
# no information

# variable subset
reg.tree.3=rpart(mean.cover.response ~ log.height+log.rootN+log.SLA, data=all.data)
rpart.plot(reg.tree.3)
summary(reg.tree.3)
# no information

reg.tree.3.site=rpart(mean.cover.response ~ log.height+log.rootN+log.SLA, 
                 weights = all.data$site.id, data=all.data)
rpart.plot(reg.tree.3.site)
summary(reg.tree.3.site)
# SLA

# variable subset
reg.tree.4=rpart(mean.cover.response ~ log.rootN+log.SLA+log.root.depth, data=all.data)
rpart.plot(reg.tree.4)
summary(reg.tree.4)
# no information

reg.tree.4.site=rpart(mean.cover.response ~ log.rootN+log.SLA+log.root.depth,
                      weights = all.data$site.id, data=all.data)
rpart.plot(reg.tree.4.site)
summary(reg.tree.4.site)
# root depth

# complete cases
reg.tree.5=rpart(mean.cover.response ~ log.leafN+log.height+log.rootN+
                   log.SLA+log.root.depth+log.RTD+log.SRL+log.root.diam, data=all.complete)
rpart.plot(reg.tree.5)
summary(reg.tree.5)

reg.tree.5.site=rpart(mean.cover.response ~ log.leafN+log.height+log.rootN+
                   log.SLA+log.root.depth+log.RTD+log.SRL+log.root.diam, 
                   weights = all.complete$site.id, data=all.complete)
rpart.plot(reg.tree.5.site)
summary(reg.tree.5)
# root diameter


#### boosted regression trees ####
# to avoid overfitting, regression trees are pruned by means of cross validation

brt.1=rpart(cover.response ~ Leaf.area, data=test.data,
            control = rpart.control(cp=0.05))

# add in the parameter cp which tells the tree not to incorporate a new partition unless the 
# R2 of the model increases at least cp units, 0.05 in this case.

rpart.plot(brt.1)
summary(brt.1)

# using maxdepth control, sets maximum depth of any node of the tree with 1 meaning there is only 1 node (1 decision)

brt.2=rpart(cover.change ~ total.leaf.area + total.leaf.mass + 
              stem.mass + fine.root.mass + total.Length.cm., data=test.data,
            control = rpart.control(maxdepth = 1))
rpart.plot(brt.2)

# changing maxdepth = 5

brt.3=rpart(cover.change ~ total.leaf.area + total.leaf.mass + 
              stem.mass + fine.root.mass + total.Length.cm., data=test.data,
            control = rpart.control(maxdepth = 5))
rpart.plot(brt.3)

# need to make decisions for performing BRT, page 99
# tree complexity = max number of nodes per tree
# learning rate = scalar that multiplies the contribution of each tree
# bag fraction = proportion of data used to train the model
# number of trees to adjust, minimum of 1000

# page 100
# bag fraction = 75%, use bigger proportion of data to fit model and less observations to test it
# learning rate = keep a value that results in fitting at least 1000 trees
# n.trees = number of initial trees to fit

# start with learning rate 0.1
#### Model to use ####
# all variables
brt.all.1=gbm.step(data=all.data, gbm.x = c(19:26), gbm.y=10,
               family = "gaussian", tree.complexity = 1, learning.rate = 0.1,
               bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = all.data$site.id)

brt.sub1.1=gbm.step(data=all.data, gbm.x = c(19,23,26), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.1,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub2.1=gbm.step(data=all.data, gbm.x = c(20,21,22), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.1,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub3.1=gbm.step(data=all.data, gbm.x = c(21,22,23), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.1,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.complete.1=gbm.step(data=all.complete, gbm.x = c(19:26), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.1,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = all.complete$site.id)

# change learning rate to 0.01

brt.all.2=gbm.step(data=all.data, gbm.x = c(19:26), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.01,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub1.2=gbm.step(data=all.data, gbm.x = c(19,23,26), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.01,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub2.2=gbm.step(data=all.data, gbm.x = c(20,21,22), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.01,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub3.2=gbm.step(data=all.data, gbm.x = c(21,22,23), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.01,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.complete.2=gbm.step(data=all.complete, gbm.x = c(19:26), gbm.y=10,
                        family = "gaussian", tree.complexity = 1, learning.rate = 0.01,
                        bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.complete$site.id)

# change learning rate to 0.001

brt.all.3=gbm.step(data=all.data, gbm.x = c(19:26), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)
# worked 650 trees without sites

brt.sub1.3=gbm.step(data=all.data, gbm.x = c(19,23,26), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)
# worked 300 trees

brt.sub2.3=gbm.step(data=all.data, gbm.x = c(20,21,22), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub3.3=gbm.step(data=all.data, gbm.x = c(21,22,23), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)
# worked 950 trees

brt.complete.3=gbm.step(data=all.complete, gbm.x = c(19:26), gbm.y=10,
                        family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                        bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.complete$site.id)

# change learning rate to 0.0001

brt.all.4=gbm.step(data=all.data, gbm.x = c(19:26), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)
# worked 1000 trees

brt.sub1.4=gbm.step(data=all.data, gbm.x = c(19,23,26), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub2.4=gbm.step(data=all.data, gbm.x = c(20,21,22), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)

brt.sub3.4=gbm.step(data=all.data, gbm.x = c(21,22,23), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)
# worked 1000 trees

brt.complete.4=gbm.step(data=all.complete, gbm.x = c(19:26), gbm.y=10,
                        family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                        bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.complete$site.id)

# change learning rate to 0.00001

brt.sub1.5=gbm.step(data=all.data, gbm.x = c(19,23,26), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.00001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)
# worked 1000 trees

brt.sub2.5=gbm.step(data=all.data, gbm.x = c(19,23,26), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.00001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.data$site.id)
# worked, 1000 trees

brt.complete.5=gbm.step(data=all.complete, gbm.x = c(19:26), gbm.y=10,
                        family = "gaussian", tree.complexity = 1, learning.rate = 0.00001,
                        bag.fraction = 0.75, n.trees = 50, verbose = FALSE,site.weights = all.complete$site.id)
# worked 1000 trees


# model evaluation
# plots of predicted versus observed, estimate R2

# All model with learning rate = 0.0001
predicted.brt.all.4=predict(brt.all.4, n.trees = brt.all.4$n.trees)
observed = all.data$mean.cover.response
plot(predicted.brt.all.4~observed)
abline(0, 1, col = 2)

R2.brt.all.4 = 1-(brt.all.4$self.statistics$mean.resid/brt.all.4$self.statistics$mean.null)
# 0.0045

# which predictors are most important
summary(brt.all.4)
# RTD

# Subset model 1 with learning rate = 0.00001
predicted.brt.sub1.5=predict(brt.sub1.5, n.trees = brt.sub1.5$n.trees)
observed = all.data$mean.cover.response
plot(predicted.brt.sub1.5~observed)
abline(0, 1, col = 2)

R2.brt.sub1.5 = 1-(brt.sub1.5$self.statistics$mean.resid/brt.sub1.5$self.statistics$mean.null)
# 0.00023

# which predictors are most important
summary(brt.sub1.5)
# root depth

# Subset model 2 with learning rate = 0.00001
predicted.brt.sub2.5=predict(brt.sub2.5, n.trees = brt.sub2.5$n.trees)
observed = all.data$mean.cover.response
plot(predicted.brt.sub2.5~observed)
abline(0, 1, col = 2)

R2.brt.sub2.5 = 1-(brt.sub2.5$self.statistics$mean.resid/brt.sub2.5$self.statistics$mean.null)
# 0.00023

# which predictors are most important
summary(brt.sub2.5)
# root depth

# Subset model 3 with learning rate = 0.0001
predicted.brt.sub3.4=predict(brt.sub3.4, n.trees = brt.sub3.4$n.trees)
observed = all.data$mean.cover.response
plot(predicted.brt.sub3.4~observed)
abline(0, 1, col = 2)

R2.brt.sub3.4 = 1-(brt.sub3.4$self.statistics$mean.resid/brt.sub3.4$self.statistics$mean.null)
# 0.002

# which predictors are most important
summary(brt.sub3.4)
# root depth

# Complete cases model 5 with learning rate = 0.00001
predicted.brt.complete.5=predict(brt.complete.5, n.trees = brt.complete.5$n.trees)
observed = all.complete$mean.cover.response
plot(predicted.brt.complete.5~observed)
abline(0, 1, col = 2)

R2.brt.complete.5 = 1-(brt.complete.5$self.statistics$mean.resid/brt.complete.5$self.statistics$mean.null)
# 0.0008349757

# which predictors are most important
summary(brt.complete.5)
# root diameter

# determining best tree complexity, pages 104-105
# add silent = T so it doesn't print everything
# removed silent = T because wasn't working and error told me to change learning rate
# learning rate 0.001 only works with tree complexity = 1

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                       ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.all.variables <- gbm.step(data=all.data,
                            gbm.x = c(19,20,21,22,23,24,25,26),
                            gbm.y = 10,
                            family = "gaussian",
                            tree.complexity = tcomp,
                            learning.rate = 0.0001,
                            bag.fraction = 0.75,
                            n.trees = 50,
                            plot.main=F, plot.folds=F,
                            site.weights = all.data$site.id)
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

# model complexity = 8

# subset of variables 1
R2Obs.sub.var.1 <- list()
importancePred.sub.var.1 <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.sub.var.1[[tcomp]] <- numeric(nreps)
  importancePred.sub.var.1[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:3),
                                                         ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.sub.var.1 <- gbm.step(data=all.data,
                              gbm.x = c(19,23,26),
                              gbm.y = 10,
                              family = "gaussian",
                              tree.complexity = tcomp,
                              learning.rate = 0.0000000001,
                              bag.fraction = 0.75,
                              n.trees = 50,
                              plot.main=F, plot.folds=F, site.weights = all.data$site.id)
    #R2 adj:
    R2Obs.sub.var.1[[tcomp]][i] <- 1 - (BRT.sub.var.1$self.statistics$mean.resid /
                                          BRT.sub.var.1$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.sub.var.1[[tcomp]]) <- sort(rownames(summary(BRT.sub.var.1)))
    }
    importancePred.sub.var.1[[tcomp]][, i] <-
      summary(BRT.sub.var.1)[rownames(importancePred.sub.var.1[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.sub.var.1, mean)
sds <- sapply(R2Obs.sub.var.1, sd)
plot(1:length(R2Obs.sub.var.1), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.sub.var.1)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# subset of variables 2
R2Obs.sub.var.2 <- list()
importancePred.sub.var.2 <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.sub.var.2[[tcomp]] <- numeric(nreps)
  importancePred.sub.var.2[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:3),
                                                         ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.sub.var.2 <- gbm.step(data=all.data,
                              gbm.x = c(20,21,22),
                              gbm.y = 10,
                              family = "gaussian",
                              tree.complexity = tcomp,
                              learning.rate = 0.000000001,
                              bag.fraction = 0.75,
                              n.trees = 50,
                              plot.main=F, plot.folds=F, site.weights = all.data$site.id)
    #R2 adj:
    R2Obs.sub.var.2[[tcomp]][i] <- 1 - (BRT.sub.var.2$self.statistics$mean.resid /
                                          BRT.sub.var.2$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.sub.var.2[[tcomp]]) <- sort(rownames(summary(BRT.sub.var.2)))
    }
    importancePred.sub.var.2[[tcomp]][, i] <-
      summary(BRT.sub.var.2)[rownames(importancePred.sub.var.2[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.sub.var.2, mean)
sds <- sapply(R2Obs.sub.var.2, sd)
plot(1:length(R2Obs.sub.var.2), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.sub.var.2)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# subset of variables 3
R2Obs.sub.var.3 <- list()
importancePred.sub.var.3 <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.sub.var.3[[tcomp]] <- numeric(nreps)
  importancePred.sub.var.3[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:3),
                                                         ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.sub.var.3 <- gbm.step(data=all.data,
                              gbm.x = c(21,22,23),
                              gbm.y = 10,
                              family = "gaussian",
                              tree.complexity = tcomp,
                              learning.rate = 0.0001,
                              bag.fraction = 0.75,
                              n.trees = 50,
                              plot.main=F, plot.folds=F, site.weights = all.data$site.id)
    #R2 adj:
    R2Obs.sub.var.3[[tcomp]][i] <- 1 - (BRT.sub.var.3$self.statistics$mean.resid /
                                          BRT.sub.var.3$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.sub.var.3[[tcomp]]) <- sort(rownames(summary(BRT.sub.var.3)))
    }
    importancePred.sub.var.3[[tcomp]][, i] <-
      summary(BRT.sub.var.3)[rownames(importancePred.sub.var.3[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.sub.var.3, mean)
sds <- sapply(R2Obs.sub.var.3, sd)
plot(1:length(R2Obs.sub.var.3), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.sub.var.3)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# complexity = 4

# complete cases variables
R2Obs.complete <- list()
importancePred.complete <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.complete[[tcomp]] <- numeric(nreps)
  importancePred.complete[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                        ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.complete <- gbm.step(data=all.complete,
                             gbm.x = c(19,20,21,22,23,24,25,26),
                             gbm.y = 10,
                             family = "gaussian",
                             tree.complexity = tcomp,
                             learning.rate = 0.00001,
                             bag.fraction = 0.75,
                             n.trees = 50,
                             plot.main=F, plot.folds=F, site.weights = all.complete$site.id)
    #R2 adj:
    R2Obs.complete[[tcomp]][i] <- 1 - (BRT.complete$self.statistics$mean.resid /
                                         BRT.complete$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.complete[[tcomp]]) <- sort(rownames(summary(BRT.complete)))
    }
    importancePred.complete[[tcomp]][, i] <-
      summary(BRT.complete)[rownames(importancePred.complete[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.complete, mean)
sds <- sapply(R2Obs.complete, sd)
plot(1:length(R2Obs.complete), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.complete)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}



# to select tree complexity using ANOVA on R2 with Tukey test

nreps <- 10
tcFactor <- as.factor(rep(1:3, each = nreps))
R2Vector <- unlist(R2Obs.all.variables)
model <- lm(R2Vector ~ tcFactor)
library(multcomp)
TukeyModel <- glht(model, linfct = mcp(tcFactor = "Tukey"))
TukeyLetters <- cld(TukeyModel)$mcletters$Letters

plot(1:length(R2Obs.all.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.all.variables)){
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}
text(x = 1:length(R2Obs.all.variables), y = 0.77, labels = TukeyLetters)

# best tree complexity is 2

# how important are the different predictors across the simulations
sort(rowMeans(importancePred.all.variables[[3]]), decreasing = T)
# always root tissue density

# examine shape of cover.change to root tissue density using partial dependence plot
# all other explanatory variables are fixed

#### Final models ####
BRT.final.all.variable <- gbm.step(data=all.data,
                        gbm.x = c(19,20,21,22,23,24,25,26),
                        gbm.y = 10,
                        family = "gaussian",
                        tree.complexity = 8,
                        learning.rate = 0.0001,
                        bag.fraction = 0.75,
                        n.trees = 50,
                        plot.main=T, plot.folds=F, site.weights = all.data$site.id)

summary(BRT.final.all.variable)

gbm.plot(BRT.final.all.variable)
gbm.plot.fits(BRT.final.all.variable)

plot.gbm(BRT.final.all.variable, i.var = c("log.SLA"))
# cover change increases with root tissue density

# investigation of interactions
gbm.interactions(BRT.final.all.variable)$rank.list

# significant interactions
# RTD x root depth
# root diam x height
# RTD x height

gbm.perspec(BRT.final.all.variable,
            6, # log.RTD
            5, # log.root.depth
            y.label = "root depth",
            x.label = "RTD",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-0.70,0.66),
            theta=10)
#theta=10, perspective = TRUE

gbm.perspec(BRT.final.all.variable,
            8, # root diam
            2, # height
            y.label = "height",
            x.label = "root diam",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.75,-0.28))

gbm.perspec(BRT.final.all.variable,
            6, # RTD
            2, # height
            y.label = "height",
            x.label = "RTD",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.8,0.56))

# assesses the potential to remove predictors using k-fold cross validation

brt.simple.all.variables = gbm.simplify(BRT.final.all.variable)

BRT.final.sub.1.5 <- gbm.step(data=all.data,
                                   gbm.x = c(19,23,26),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = 8,
                                   learning.rate = 0.00001,
                                   bag.fraction = 0.75,
                                   n.trees = 50,
                                   plot.main=T, plot.folds=F, site.weights = all.data$site.id)

summary(BRT.final.sub.1.5)

gbm.plot(BRT.final.sub.1.5)
gbm.plot.fits(BRT.final.sub.1.5)

plot.gbm(BRT.final.sub.1.5, i.var = c("log.root.depth"))
# cover change increases with root tissue density

# investigation of interactions
gbm.interactions(BRT.final.sub.1.5)$rank.list

BRT.final.sub.2.5 <- gbm.step(data=all.data,
                              gbm.x = c(20,21,22),
                              gbm.y = 10,
                              family = "gaussian",
                              tree.complexity = 8,
                              learning.rate = 0.00001,
                              bag.fraction = 0.75,
                              n.trees = 50,
                              plot.main=T, plot.folds=F, site.weights = all.data$site.id)

summary(BRT.final.sub.2.5)

gbm.plot(BRT.final.sub.2.5)
gbm.plot.fits(BRT.final.sub.2.5)

plot.gbm(BRT.final.sub.2.5, i.var = c("log.SLA"))
# cover change increases with root tissue density

# investigation of interactions
gbm.interactions(BRT.final.sub.2.5)$rank.list

BRT.final.sub.3.4 <- gbm.step(data=all.data,
                              gbm.x = c(21,22,23),
                              gbm.y = 10,
                              family = "gaussian",
                              tree.complexity = 8,
                              learning.rate = 0.0001,
                              bag.fraction = 0.75,
                              n.trees = 50,
                              plot.main=T, plot.folds=F, site.weights = all.data$site.id)

summary(BRT.final.sub.3.4)

gbm.plot(BRT.final.sub.3.4)
gbm.plot.fits(BRT.final.sub.3.4)

plot.gbm(BRT.final.sub.3.4, i.var = c("log.SLA"))
# cover change increases with root tissue density

# investigation of interactions
gbm.interactions(BRT.final.sub.3.4)$rank.list

# significant interactions
# root depth x rootN
# root depth x SLA

gbm.perspec(BRT.final.sub.3.4,
            3, # log.root.depth
            1, # log.rootN
            y.label = "root depth",
            x.label = "RootN",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-1,0.08))

gbm.perspec(BRT.final.sub.3.4,
            3, # log.root.depth
            2, # log.SLA
            y.label = "root depth",
            x.label = "SLA",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.6,0.15))

BRT.final.complete.5<- gbm.step(data=all.complete,
                              gbm.x = c(19,20,21,22,23,24,25,26),
                              gbm.y = 10,
                              family = "gaussian",
                              tree.complexity = 8,
                              learning.rate = 0.00001,
                              bag.fraction = 0.75,
                              n.trees = 50,
                              plot.main=T, plot.folds=F, site.weights = all.complete$site.id,
                              interaction.depth = 3)

summary(BRT.final.complete.5)

gbm.plot(BRT.final.complete.5)
gbm.plot.fits(BRT.final.complete.5)

plot.gbm(BRT.final.complete.5, i.var = c("log.root.diam"))
# cover change increases with root tissue density

# investigation of interactions
gbm.interactions(BRT.final.complete.5)$rank.list

gbm.perspec(BRT.final.complete.5,
            8, # log.root.depth
            4, # log.SLA
            y.label = "root depth",
            x.label = "SLA",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-1.22,-1.15))


gbm.perspec(BRT.final.complete.5,
            8, # log.root.depth
            4, # log.SLA
            y.label = "root depth",
            x.label = "SLA",
            perspective = TRUE,
            z.range = c(-1.22,-1.15))





