# Script for analysis of Drought Net data

library(dismo)
library(gbm)

trait.data = read.csv("./Formatted.Data/trait.species.trt.yr1.final.csv")
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")


# subset traits to only those of interest

trait.data.2 = trait.data[,c(1,7,8,10,12,14,15,18,20,27:29)]

# subset traits so they must have SLA

trait.data.3 = subset(trait.data.2, trait.data.2$SLA_m2.kg > 0 ) # 657 species

table(is.na(trait.data.3$leafN.mg.g))
table(is.na(trait.data.3$height.m))
table(is.na(trait.data.3$rootN.mg.g))
table(is.na(trait.data.3$root.depth_m))
table(is.na(trait.data.3$RTD.g.cm3))
table(is.na(trait.data.3$SRL_m.g))
table(is.na(trait.data.3$rootDiam.mm))

# merge trait data and cover data

all.data = merge(cover.data, trait.data.3, by="Taxon") # 1096 data points with site

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
cor.test(all.data$height.m, all.data$root.depth_m)
cor.test(all.data$height.m, all.data$RTD.g.cm3)
cor.test(all.data$height.m, all.data$SRL_m.g) # correlated
cor.test(all.data$height.m, all.data$rootDiam.mm)
cor.test(all.data$rootN.mg.g, all.data$SLA_m2.kg)
cor.test(all.data$rootN.mg.g, all.data$root.depth_m)
cor.test(all.data$rootN.mg.g, all.data$RTD.g.cm3) # correlated
cor.test(all.data$rootN.mg.g, all.data$SRL_m.g) # correlated
cor.test(all.data$rootN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$SLA_m2.kg, all.data$root.depth_m)
cor.test(all.data$SLA_m2.kg, all.data$RTD.g.cm3) # correlated
cor.test(all.data$SLA_m2.kg, all.data$SRL_m.g) # correlated
cor.test(all.data$SLA_m2.kg, all.data$rootDiam.mm)
cor.test(all.data$root.depth_m, all.data$RTD.g.cm3) # correlated
cor.test(all.data$root.depth_m, all.data$SRL_m.g) # correlated
cor.test(all.data$root.depth_m, all.data$rootDiam.mm)
cor.test(all.data$RTD.g.cm3, all.data$SRL_m.g) # correlated
cor.test(all.data$RTD.g.cm3, all.data$rootDiam.mm)
cor.test(all.data$SRL_m.g, all.data$rootDiam.mm)

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
hist(all.data$mean.cover.response)

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

# log and scale values to transform if non-normal, , mean of 0, sd = 1
all.complete$log.leafN = as.vector(scale(log(all.complete$leafN.mg.g)))
all.complete$log.height = as.vector(scale(log(all.complete$height.m)))
all.complete$log.rootN = as.vector(scale(log(all.complete$rootN.mg.g)))
all.complete$log.SLA = as.vector(scale(log(all.complete$SLA_m2.kg)))
all.complete$log.root.depth = as.vector(scale(log(all.complete$root.depth_m)))
all.complete$log.RTD = as.vector(scale(log(all.complete$RTD.g.cm3)))
all.complete$log.SRL = as.vector(scale(log(all.complete$SRL_m.g)))
all.complete$log.root.diam = as.vector(scale(log(all.complete$rootDiam.mm)))

# data set with trees and shrubs removed

no.trees = subset(all.data, !all.data$functional_group == "WOODY")

no.trees$log.leafN = as.vector(scale(log(no.trees$leafN.mg.g)))
no.trees$log.height = as.vector(scale(log(no.trees$height.m)))
no.trees$log.rootN = as.vector(scale(log(no.trees$rootN.mg.g)))
no.trees$log.SLA = as.vector(scale(log(no.trees$SLA_m2.kg)))
no.trees$log.root.depth = as.vector(scale(log(no.trees$root.depth_m)))
no.trees$log.RTD = as.vector(scale(log(no.trees$RTD.g.cm3)))
no.trees$log.SRL = as.vector(scale(log(no.trees$SRL_m.g)))
no.trees$log.root.diam = as.vector(scale(log(no.trees$rootDiam.mm)))

# change site code to numeric, continuous vector
all.data$site.id = as.numeric(as.factor(all.data$site_code))
all.complete$site.id = as.numeric(as.factor(all.complete$site_code))
no.trees$site.id = as.numeric(as.factor(no.trees$site_code))

# determining best learning rate to generate 1000 trees
# start learning rate at 0.1
all.brt.1=gbm.step(data=all.data, gbm.x = c(22:29), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = all.data$site.id)
# 0.0001 worked

comp.brt.1=gbm.step(data=all.complete, gbm.x = c(22:29), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.00000000000001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = all.complete$site.id)
# 0.00000000000001 work, will probably need to alter step size

tree.brt.1=gbm.step(data=no.trees, gbm.x = c(22:29), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = no.trees$site.id)
# 0.0001 worked 

# initial model evaluation
# plots of predicted versus observed, estimate R2

# which predictors are most important
summary(all.brt.1)
# RTD 96%

all.predict.brt=predict(all.brt.1, n.trees = all.brt.1$n.trees)
observed = all.data$mean.cover.response
plot(all.predict.brt~observed)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.brt.1$self.statistics$mean.resid/all.brt.1$self.statistics$mean.null)
# 0.0033

# which predictors are most important
summary(comp.brt.1)
# height 41%, root diam: 40%

comp.predict.brt=predict(comp.brt.1, n.trees = comp.brt.1$n.trees)
observed = all.complete$mean.cover.response
plot(comp.predict.brt~observed)
abline(0, 1, col = 2)

R2.comp.brt = 1-(comp.brt.1$self.statistics$mean.resid/comp.brt.1$self.statistics$mean.null)
# essentially nothing

# which predictors are most important
summary(tree.brt.1)
# RTD 82%

tree.predict.brt=predict(tree.brt.1, n.trees = tree.brt.1$n.trees)
observed = no.trees$mean.cover.response
plot(tree.predict.brt~observed)
abline(0, 1, col = 2)

R2.tree.brt = 1-(tree.brt.1$self.statistics$mean.resid/tree.brt.1$self.statistics$mean.null)
# 0.004

# determining best tree complexity
# learning rate 0.001 only works with tree complexity = 1

# all.data

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
                                  gbm.x = c(22:29),
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

# 4 or 5 complexity is best for all variables model with lr = 0.0001

# complete.data

# comp variables
R2Obs.comp.variables <- list()
importancePred.comp.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.comp.variables[[tcomp]] <- numeric(nreps)
  importancePred.comp.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                              ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.comp.variables <- gbm.step(data=all.complete,
                                   gbm.x = c(22:29),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = tcomp,
                                   learning.rate = 0.0000000000000001,
                                   bag.fraction = 0.75,
                                   n.trees = 50,
                                   plot.main=F, plot.folds=F,
                                   site.weights = all.complete$site.id)
    #R2 adj:
    R2Obs.comp.variables[[tcomp]][i] <- 1 - (BRT.comp.variables$self.statistics$mean.resid /
                                               BRT.comp.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.comp.variables[[tcomp]]) <- sort(rownames(summary(BRT.comp.variables)))
    }
    importancePred.comp.variables[[tcomp]][, i] <-
      summary(BRT.comp.variables)[rownames(importancePred.comp.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.comp.variables, mean)
sds <- sapply(R2Obs.comp.variables, sd)
plot(1:length(R2Obs.comp.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.comp.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# won't run as is because need to tweak step size

# no trees

# comp variables
R2Obs.tree.variables <- list()
importancePred.tree.variables <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.tree.variables[[tcomp]] <- numeric(nreps)
  importancePred.tree.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                              ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.tree.variables <- gbm.step(data=no.trees,
                                   gbm.x = c(22:29),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree..complexity = tcomp,
                                   learning.rate = 0.000001,
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

# model complexity of 1, lr = 0.000001

# final models

# all.final.brt <- gbm.step(data=all.data,
                                   gbm.x = c(22:29),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = 5,
                                   learning.rate = 0.0001,
                                   bag.fraction = 0.75,
                                   n.trees = 50,
                                   plot.main=T, plot.folds=T, site.weights = all.data$site.id)
summary(all.final.brt)
# RTD 37%

gbm.plot(all.final.brt)
gbm.plot.fits(all.final.brt)

plot.gbm(all.final.brt, i.var = c("log.RTD"))

# investigation of interactions
gbm.interactions(all.final.brt)$rank.list

# significant interactions
# RTD x height
# RTD x leafN
# root diam x RTD

gbm.perspec(all.final.brt,
            6, # log.RTD
            2, # log.height
            y.label = "height",
            x.label = "RTD",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-0.7,0.57),
            theta=5)

gbm.perspec(all.final.brt,
            6, # log.RTD
            2, # log.height
            y.label = "height",
            x.label = "RTD",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.7,0.57))

gbm.perspec(all.final.brt,
            6, # log.RTD
            1, # log.leafN
            y.label = "leafN",
            x.label = "RTD",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-0.7,0.57),
            theta=5)

gbm.perspec(all.final.brt,
            6, # log.RTD
            1, # log.leafN
            y.label = "leafN",
            x.label = "RTD",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.7,0.57))

gbm.perspec(all.final.brt,
            8, # log.root.diam
            6, # log.RTD
            y.label = "RTD",
            x.label = "root.diam",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-0.7,0.57),
            theta=50)

gbm.perspec(all.final.brt,
            8, # log.root.diam
            6, # log.RTD
            y.label = "RTD",
            x.label = "root.diam",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.7,0.57))












