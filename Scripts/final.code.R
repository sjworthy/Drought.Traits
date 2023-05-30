# Script for analysis of Drought Net data

library(dismo)
library(gbm)

trait.data = read.csv("./Formatted.Data/trait.species.trt.yr1.final.new.csv")

# subset traits to only those of interest

trait.data.2 = trait.data[,c(1,7,8,10,12,14,15,18,20,27:29)]

#### checking for outliers ####
# all outliers removed manually from trait.species.trt.yr1.final.new and made new file trait.species.trt.yr1.outlier
hist(trait.data.2$leafN.mg.g)
boxplot(trait.data.2$leafN.mg.g)
mean = mean(trait.data.2$leafN.mg.g, na.rm = TRUE)
std = sd(trait.data.2$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$leafN.mg.g[which(trait.data.2$leafN.mg.g <Tmin | trait.data.2$leafN.mg.g > Tmax)]
# removed leafN 58.30000 and 67.58333

hist(trait.data.2$height.m)
boxplot(trait.data.2$height.m)
mean = mean(trait.data.2$height.m, na.rm = TRUE)
std = sd(trait.data.2$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$height.m[which(trait.data.2$height.m <Tmin | trait.data.2$height.m > Tmax)]
# removed height > 14.37587 16.68757 19.05000 20.17067 20.42708 20.46877 20.50100 20.83333 
# 21.00000 21.58000 22.64920 23.39467 23.60517 24.86933 25.75100 27.10186 28.99781 32.57366 46.30267

hist(trait.data.2$rootN.mg.g)
boxplot(trait.data.2$rootN.mg.g)
mean = mean(trait.data.2$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data.2$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$rootN.mg.g[which(trait.data.2$rootN.mg.g <Tmin | trait.data.2$rootN.mg.g > Tmax)]
# remove rootN 33.34667 34.39000 39.26274

hist(trait.data.2$SLA_m2.kg)
boxplot(trait.data.2$SLA_m2.kg)
mean = mean(trait.data.2$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data.2$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$SLA_m2.kg[which(trait.data.2$SLA_m2.kg <Tmin | trait.data.2$SLA_m2.kg > Tmax)]
# remove SLA 53.50000 57.70000 61.28000 63.13000 64.32508 66.41054 68.90000 75.90000 86.00000 96.66667 98.20000

hist(trait.data.2$root.depth_m)
boxplot(trait.data.2$root.depth_m)
mean = mean(trait.data.2$root.depth_m, na.rm = TRUE)
std = sd(trait.data.2$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$root.depth_m[which(trait.data.2$root.depth_m <Tmin | trait.data.2$root.depth_m > Tmax)]
# remove depth 3.352800  3.500000  3.709333  3.828571  6.451600 11.625000

hist(trait.data.2$RTD.g.cm3)
boxplot(trait.data.2$RTD.g.cm3)
mean = mean(trait.data.2$RTD.g.cm3, na.rm = TRUE)
std = sd(trait.data.2$RTD.g.cm3, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$RTD.g.cm3[which(trait.data.2$RTD.g.cm3 <Tmin | trait.data.2$RTD.g.cm3 > Tmax)]
# remove RTD 0.9167417 1.1945504

hist(trait.data.2$SRL_m.g)
boxplot(trait.data.2$SRL_m.g)
mean = mean(trait.data.2$SRL_m.g, na.rm = TRUE)
std = sd(trait.data.2$SRL_m.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$SRL_m.g[which(trait.data.2$SRL_m.g <Tmin | trait.data.2$SRL_m.g > Tmax)]
# remove SRL 601.5858 627.5050 630.3000 657.7263 760.8500 833.5418 929.9185

hist(trait.data.2$rootDiam.mm)
boxplot(trait.data.2$rootDiam.mm)
mean = mean(trait.data.2$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data.2$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data.2$rootDiam.mm[which(trait.data.2$rootDiam.mm <Tmin | trait.data.2$rootDiam.mm > Tmax)]
# remove diam 81.81050 27.13865

#### read in new trait data without outliers and cover data ####

trait.data.new = read.csv("./Formatted.Data/trait.species.trt.yr1.outlier.csv")
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")


# subset traits so they must have SLA
trait.data.2 = trait.data.new[,c(1,7,8,10,12,14,15,18,20,27:29)]
trait.data.3 = subset(trait.data.2, trait.data.2$SLA_m2.kg > 0 ) # 646 data points, 

table(is.na(trait.data.3$leafN.mg.g))
table(is.na(trait.data.3$height.m))
table(is.na(trait.data.3$rootN.mg.g))
table(is.na(trait.data.3$root.depth_m))
table(is.na(trait.data.3$RTD.g.cm3))
table(is.na(trait.data.3$SRL_m.g))
table(is.na(trait.data.3$rootDiam.mm))

# merge trait data and cover data

all.data = merge(cover.data, trait.data.3, by="Taxon") # 1077 data points with site

#### test for correlation ####
cor.test(all.data$leafN.mg.g, all.data$height.m)
cor.test(all.data$leafN.mg.g, all.data$rootN.mg.g) # correlated r = 0.49
cor.test(all.data$leafN.mg.g, all.data$SLA_m2.kg) # correlated r = 0.34
cor.test(all.data$leafN.mg.g, all.data$root.depth_m)
cor.test(all.data$leafN.mg.g, all.data$RTD.g.cm3) # correlated r = -0.20
cor.test(all.data$leafN.mg.g, all.data$SRL_m.g) # correlated r = 0.24
cor.test(all.data$leafN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$height.m, all.data$rootN.mg.g)
cor.test(all.data$height.m, all.data$SLA_m2.kg) # correlated r = -0.08
cor.test(all.data$height.m, all.data$root.depth_m)
cor.test(all.data$height.m, all.data$RTD.g.cm3)
cor.test(all.data$height.m, all.data$SRL_m.g) # correlated r = -0.10
cor.test(all.data$height.m, all.data$rootDiam.mm) # correlated r = 0.14
cor.test(all.data$rootN.mg.g, all.data$SLA_m2.kg) # correlated r = 0.22
cor.test(all.data$rootN.mg.g, all.data$root.depth_m)
cor.test(all.data$rootN.mg.g, all.data$RTD.g.cm3) # correlated r = -0.16
cor.test(all.data$rootN.mg.g, all.data$SRL_m.g)
cor.test(all.data$rootN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$SLA_m2.kg, all.data$root.depth_m) # correlated r = -0.15
cor.test(all.data$SLA_m2.kg, all.data$RTD.g.cm3) # correlated r = -0.24
cor.test(all.data$SLA_m2.kg, all.data$SRL_m.g) # correlated r = 0.42
cor.test(all.data$SLA_m2.kg, all.data$rootDiam.mm)
cor.test(all.data$root.depth_m, all.data$RTD.g.cm3) # correlated r = 0.13
cor.test(all.data$root.depth_m, all.data$SRL_m.g) # correlated r =-0.15
cor.test(all.data$root.depth_m, all.data$rootDiam.mm) #  correlated r = 0.15
cor.test(all.data$RTD.g.cm3, all.data$SRL_m.g) # correlated r = -0.27
cor.test(all.data$RTD.g.cm3, all.data$rootDiam.mm) # correlated r = 0.13
cor.test(all.data$SRL_m.g, all.data$rootDiam.mm) #  correlated r = -0.09

#### correlation between response and predictors ####

cor.test(all.data$mean.cover.response, all.data$leafN.mg.g) # r = 0.07
cor.test(all.data$mean.cover.response, all.data$height.m) # r = -0.01
cor.test(all.data$mean.cover.response, all.data$rootN.mg.g) # r = -0.01
cor.test(all.data$mean.cover.response, all.data$SLA_m2.kg) # r = 0.04
cor.test(all.data$mean.cover.response, all.data$root.depth_m) # r = 0.04
cor.test(all.data$mean.cover.response, all.data$RTD.g.cm3) # r = 0.06
cor.test(all.data$mean.cover.response, all.data$SRL_m.g) # r = -0.04
cor.test(all.data$mean.cover.response, all.data$rootDiam.mm) # r = 0.06

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

#### get data set of only complete trait cases ####

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

#### data set with trees and shrubs removed ####

no.trees = subset(all.data, !all.data$functional_group == "WOODY")

no.trees$log.leafN = as.vector(scale(log(no.trees$leafN.mg.g)))
no.trees$log.height = as.vector(scale(log(no.trees$height.m)))
no.trees$log.rootN = as.vector(scale(log(no.trees$rootN.mg.g)))
no.trees$log.SLA = as.vector(scale(log(no.trees$SLA_m2.kg)))
no.trees$log.root.depth = as.vector(scale(log(no.trees$root.depth_m)))
no.trees$log.RTD = as.vector(scale(log(no.trees$RTD.g.cm3)))
no.trees$log.SRL = as.vector(scale(log(no.trees$SRL_m.g)))
no.trees$log.root.diam = as.vector(scale(log(no.trees$rootDiam.mm)))

#### change site code to numeric, continuous vector ####
all.data$site.id = as.numeric(as.factor(all.data$site_code))
all.complete$site.id = as.numeric(as.factor(all.complete$site_code))
no.trees$site.id = as.numeric(as.factor(no.trees$site_code))

# data set split by lifespan

# need to read out data anf fix the lifespan to particular sites since some species have different lifespan at different sites

# write.csv(all.data, file="./Formatted.Data/all.data.response.0.csv")

all.data.ls = read.csv("./Formatted.Data/all.data.lifespan.csv", row.names = 1)
table(all.data.ls$local_lifespan)

annual.data = subset(all.data.ls, all.data.ls$local_lifespan == "ANNUAL")

annual.data$log.leafN = as.vector(scale(log(annual.data$leafN.mg.g)))
annual.data$log.height = as.vector(scale(log(annual.data$height.m)))
annual.data$log.rootN = as.vector(scale(log(annual.data$rootN.mg.g)))
annual.data$log.SLA = as.vector(scale(log(annual.data$SLA_m2.kg)))
annual.data$log.root.depth = as.vector(scale(log(annual.data$root.depth_m)))
annual.data$log.RTD = as.vector(scale(log(annual.data$RTD.g.cm3)))
annual.data$log.SRL = as.vector(scale(log(annual.data$SRL_m.g)))
annual.data$log.root.diam = as.vector(scale(log(annual.data$rootDiam.mm)))

perennial.data = subset(all.data.ls, all.data.ls$local_lifespan == "PERENNIAL")

perennial.data$log.leafN = as.vector(scale(log(perennial.data$leafN.mg.g)))
perennial.data$log.height = as.vector(scale(log(perennial.data$height.m)))
perennial.data$log.rootN = as.vector(scale(log(perennial.data$rootN.mg.g)))
perennial.data$log.SLA = as.vector(scale(log(perennial.data$SLA_m2.kg)))
perennial.data$log.root.depth = as.vector(scale(log(perennial.data$root.depth_m)))
perennial.data$log.RTD = as.vector(scale(log(perennial.data$RTD.g.cm3)))
perennial.data$log.SRL = as.vector(scale(log(perennial.data$SRL_m.g)))
perennial.data$log.root.diam = as.vector(scale(log(perennial.data$rootDiam.mm)))

annual.data$site.id = as.numeric(as.factor(annual.data$site_code))
perennial.data$site.id = as.numeric(as.factor(perennial.data$site_code))

# perennial data without woody functional group

perennial.tree = subset(perennial.data, !perennial.data$functional_group == "WOODY")

perennial.tree$log.leafN = as.vector(scale(log(perennial.tree$leafN.mg.g)))
perennial.tree$log.height = as.vector(scale(log(perennial.tree$height.m)))
perennial.tree$log.rootN = as.vector(scale(log(perennial.tree$rootN.mg.g)))
perennial.tree$log.SLA = as.vector(scale(log(perennial.tree$SLA_m2.kg)))
perennial.tree$log.root.depth = as.vector(scale(log(perennial.tree$root.depth_m)))
perennial.tree$log.RTD = as.vector(scale(log(perennial.tree$RTD.g.cm3)))
perennial.tree$log.SRL = as.vector(scale(log(perennial.tree$SRL_m.g)))
perennial.tree$log.root.diam = as.vector(scale(log(perennial.tree$rootDiam.mm)))

perennial.tree$site.id = as.numeric(as.factor(perennial.tree$site_code))

# determining best learning rate to generate 1000 trees
# start learning rate at 0.1
set.seed(2023)
all.brt.1=gbm.step(data=all.data, gbm.x = c(11:18), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = all.data$site.id)
# 0.0001 worked
# tolerance 0.0407

set.seed(2023)
comp.brt.1=gbm.step(data=all.complete, gbm.x = c(11:18), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.00000000000000001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = all.complete$site.id)
# the standard deviation is zero

set.seed(2023)
tree.brt.1=gbm.step(data=no.trees, gbm.x = c(11:18), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.0001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = no.trees$site.id)
# 0.0001 worked 
# tolerance 0.0362

set.seed(2023)
annual.brt.1=gbm.step(data=annual.data, gbm.x = c(11:18), gbm.y=10,
                    family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
                    bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = annual.data$site.id)
# 0.001 worked 
# tolerance 0.0256

set.seed(2023)
perennial.brt.1=gbm.step(data=perennial.data, gbm.x = c(11:18), gbm.y=10,
                      family = "gaussian", tree.complexity = 1, learning.rate = 0.000000000000000001,
                      bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = perennial.data$site.id)
# the standard deviation is zero

set.seed(2023)
perennial.tree.brt.1=gbm.step(data=perennial.tree, gbm.x = c(11:18), gbm.y=10,
                         family = "gaussian", tree.complexity = 1, learning.rate = 0.000000000000000001,
                         bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = perennial.tree$site.id)
# the standard deviation is zero

# initial model evaluation
# plots of predicted versus observed, estimate R2

# which predictors are most important
summary(all.brt.1)
# RTD 95%

all.predict.brt=predict(all.brt.1, n.trees = all.brt.1$n.trees)
observed = all.data$mean.cover.response
plot(all.predict.brt~observed)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.brt.1$self.statistics$mean.resid/all.brt.1$self.statistics$mean.null)
# 0.0033

# which predictors are most important
summary(comp.brt.1)
# rootDiam 46%, height: 34%

comp.predict.brt=predict(comp.brt.1, n.trees = comp.brt.1$n.trees)
observed = all.complete$mean.cover.response
plot(comp.predict.brt~observed)
abline(0, 1, col = 2)

R2.comp.brt = 1-(comp.brt.1$self.statistics$mean.resid/comp.brt.1$self.statistics$mean.null)
# essentially nothing

# which predictors are most important
summary(tree.brt.1)
# RTD 89%

tree.predict.brt=predict(tree.brt.1, n.trees = tree.brt.1$n.trees)
observed = no.trees$mean.cover.response
plot(tree.predict.brt~observed)
abline(0, 1, col = 2)

R2.tree.brt = 1-(tree.brt.1$self.statistics$mean.resid/tree.brt.1$self.statistics$mean.null)
# 0.004

# which predictors are most important
summary(annual.brt.1)
# RTD 29%, leafN 25 %

annual.predict.brt=predict(annual.brt.1, n.trees = annual.brt.1$n.trees)
observed = annual.data$mean.cover.response
plot(annual.predict.brt~observed)
abline(0, 1, col = 2)

R2.annual.brt = 1-(annual.brt.1$self.statistics$mean.resid/annual.brt.1$self.statistics$mean.null)
# 0.12

# which predictors are most important
summary(perennial.brt.1)
# RTD 87%

perennial.predict.brt=predict(perennial.brt.1, n.trees = perennial.brt.1$n.trees)
observed = perennial.data$mean.cover.response
plot(perennial.predict.brt~observed)
abline(0, 1, col = 2)

R2.perennial.brt = 1-(perennial.brt.1$self.statistics$mean.resid/perennial.brt.1$self.statistics$mean.null)
# essentially nothing

# which predictors are most important
summary(perennial.tree.brt.1)
# RTD 87%

perennial.tree.predict.brt=predict(perennial.tree.brt.1, n.trees = perennial.tree.brt.1$n.trees)
observed = perennial.tree$mean.cover.response
plot(perennial.tree.predict.brt~observed)
abline(0, 1, col = 2)

R2.perennial.brt = 1-(perennial.tree.brt.1$self.statistics$mean.resid/perennial.tree.brt.1$self.statistics$mean.null)
# essentially nothing

# determining best tree complexity
# learning rate 0.001 only works with tree complexity = 1

# all.data

# all variables
R2Obs.all.variables <- list()
importancePred.all.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.all.variables[[tcomp]] <- numeric(nreps)
  importancePred.all.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                             ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.all.variables <- gbm.step(data=all.data,
                                  gbm.x = c(11:18),
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

# 6 complexity is best for all variables model with lr = 0.0001

# complete.data

# comp variables
R2Obs.comp.variables <- list()
importancePred.comp.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.comp.variables[[tcomp]] <- numeric(nreps)
  importancePred.comp.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                              ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.comp.variables <- gbm.step(data=all.complete,
                                   gbm.x = c(11:18),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = tcomp,
                                   learning.rate = 0.00000000000000001,
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

# complexity of 2 with lr = 0.00000000000000001


# tree variables
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

# model complexity of 3, lr = 0.0001

# annual.data

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

# 6 or 10 complexity is best for all variables model with lr = 0.001

# perennial variables

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
                                        learning.rate = 0.000000000000000001,
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

# 1 or 2 complexity is best for all variables model with lr = 0.00000000000000001


# perennial without tree variables

R2Obs.perennial.tree.variables <- list()
importancePred.perennial.tree.variables <- list()
nreps <- 10 #number of simulations
set.seed(2023)
for (tcomp in 1:10) {
  R2Obs.perennial.tree.variables[[tcomp]] <- numeric(nreps)
  importancePred.perennial.tree.variables[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:8),
                                                                        ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.perennial.tree.variables <- gbm.step(data=perennial.tree,
                                             gbm.x = c(11:18),
                                             gbm.y = 10,
                                             family = "gaussian",
                                             tree.complexity = tcomp,
                                             learning.rate = 0.000000000000000001,
                                             bag.fraction = 0.75,
                                             n.trees = 50,
                                             plot.main=F, plot.folds=F,
                                             site.weights = perennial.tree$site.id)
    #R2 adj:
    R2Obs.perennial.tree.variables[[tcomp]][i] <- 1 - (BRT.perennial.tree.variables$self.statistics$mean.resid /
                                                         BRT.perennial.tree.variables$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.perennial.tree.variables[[tcomp]]) <- sort(rownames(summary(BRT.perennial.tree.variables)))
    }
    importancePred.perennial.tree.variables[[tcomp]][, i] <-
      summary(BRT.perennial.tree.variables)[rownames(importancePred.perennial.tree.variables[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.perennial.tree.variables, mean)
sds <- sapply(R2Obs.perennial.tree.variables, sd)
plot(1:length(R2Obs.perennial.tree.variables), means, ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.perennial.tree.variables)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# 1 complexity is best for all variables model with lr = 0.000000000000000001


#### final models ####

set.seed(2023)
all.final.brt <- gbm.step(data=all.data,
                                   gbm.x = c(11:18),
                                   gbm.y = 10,
                                   family = "gaussian",
                                   tree.complexity = 6,
                                   learning.rate = 0.0001,
                                   bag.fraction = 0.75,
                                   n.trees = 50,
                                   plot.main=T, plot.folds=T, site.weights = all.data$site.id)
summary(all.final.brt)
# RTD 30%

gbm.plot(all.final.brt, common.scale = FALSE)
gbm.plot.fits(all.final.brt)

plot.gbm(all.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(all.final.brt, i.var = c("root.depth_m"))
plot.gbm(all.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(all.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(all.final.brt, i.var = c("height.m"))
plot.gbm(all.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(all.final.brt, i.var = c("SRL_m.g"))
plot.gbm(all.final.brt, i.var = c("rootN.mg.g"))

all.predict.brt=predict(all.final.brt, n.trees = all.final.brt$n.trees)
observed.all = all.data$mean.cover.response
plot(all.predict.brt~observed.all)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.final.brt$self.statistics$mean.resid/all.final.brt$self.statistics$mean.null)
# 0.01

# investigation of interactions
gbm.interactions(all.final.brt)$rank.list

all.brt.simple = gbm.simplify(all.final.brt)
# keep only SLA and RTD



#### final complete cases model ####
set.seed(2023)
comp.final.brt <- gbm.step(data=all.complete,
gbm.x = c(11:18),
gbm.y = 10,
family = "gaussian",
tree.complexity = 1,
learning.rate = 0.0000000000000001,
bag.fraction = 0.75,
n.trees = 50,
plot.main=T, plot.folds=T, site.weights = all.complete$site.id)

summary(comp.final.brt)
# root diam 49%

gbm.plot(comp.final.brt)
gbm.plot.fits(comp.final.brt)

plot.gbm(comp.final.brt, i.var = c("RTD.g.cm3"))
plot.gbm(comp.final.brt, i.var = c("root.depth_m"))
plot.gbm(comp.final.brt, i.var = c("leafN.mg.g"))
plot.gbm(comp.final.brt, i.var = c("SLA_m2.kg"))
plot.gbm(comp.final.brt, i.var = c("height.m"))
plot.gbm(comp.final.brt, i.var = c("rootDiam.mm"))
plot.gbm(comp.final.brt, i.var = c("SRL_m.g"))
plot.gbm(comp.final.brt, i.var = c("rootN.mg.g"))

comp.predict.brt=predict(comp.final.brt, n.trees = comp.final.brt$n.trees)
observed.comp = all.complete$mean.cover.response
plot(comp.predict.brt~observed.comp)
abline(0, 1, col = 2)

R2.comp.brt = 1-(comp.final.brt$self.statistics$mean.resid/comp.final.brt$self.statistics$mean.null)
# explains little of the variation

# investigation of interactions
gbm.interactions(comp.final.brt)$rank.list

#### model without trees ####
set.seed(2023)
tree.final.brt <- gbm.step(data=no.trees,
gbm.x = c(11:18),
gbm.y = 10,
family = "gaussian",
tree.complexity = 3,
learning.rate = 0.0001,
bag.fraction = 0.75,
n.trees = 50,
plot.main=T, plot.folds=T, site.weights = no.trees$site.id)

summary(tree.final.brt)
# RTD 42%

gbm.plot(tree.final.brt)
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
# 0.02

# investigation of interactions
gbm.interactions(tree.final.brt)$rank.list

# significant interactions
# RTD x height


gbm.perspec(tree.final.brt,
            6, # log.RTD
            2, # log.height
            y.label = "height",
            x.label = "RTD",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-0.9,0.82),
            theta = 320)

source("gbm.perspec2.R")

gbm.perspec(tree.final.brt,
            6, # log.RTD
            2, # log.height
            y.label = "height",
            x.label = "RTD",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.9,0.82))


#### annual final model ####
set.seed(2023)
annual.final.brt <- gbm.step(data=annual.data,
                          gbm.x = c(11:18),
                          gbm.y = 10,
                          family = "gaussian",
                          tree.complexity = 10,
                          learning.rate = 0.0001,
                          bag.fraction = 0.75,
                          n.trees = 50,
                          plot.main=T, plot.folds=T, site.weights = annual.data$site.id)
summary(annual.final.brt)
# leafN 26%
# RTD 18%

gbm.plot(annual.final.brt)
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
# 0.16

# investigation of interactions
gbm.interactions(annual.final.brt)$rank.list

# significant interactions
# SRL x leafN
# SRL x depth
# diam x SRL

gbm.perspec(annual.final.brt,
            7, # SRL
            1, # leafN
            y.label = "leafN",
            x.label = "SRL",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-1.5,0.93),
            theta=50)

gbm.perspec(annual.final.brt,
            7, # SRL
            1, # leafN
            y.label = "leafN",
            x.label = "SRL",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-1.5,0.93))

gbm.perspec(annual.final.brt,
            7, # SRL
            5, # depth
            y.label = "depth",
            x.label = "SRL",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-0.9,0.90),
            theta=30)

gbm.perspec(annual.final.brt,
            7, # SRL
            5, # depth
            y.label = "depth",
            x.label = "SRL",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.9,0.90))

gbm.perspec(annual.final.brt,
            8, # diam
            7, # SRL
            y.label = "SRL",
            x.label = "diam",
            smooth = "average",
            perspective = TRUE,
            z.range = c(-0.9,0.22),
            theta=100)

gbm.perspec(annual.final.brt,
            8, # diam
            7, # SRL
            y.label = "SRL",
            x.label = "diam",
            smooth = "average",
            perspective = FALSE,
            z.range = c(-0.9,0.22))

#### perennial final model ####
set.seed(2023)
perennial.final.brt <- gbm.step(data=perennial.data,
                             gbm.x = c(11:18),
                             gbm.y = 10,
                             family = "gaussian",
                             tree.complexity = 1,
                             learning.rate = 0.00000000000000001,
                             bag.fraction = 0.75,
                             n.trees = 50,
                             plot.main=T, plot.folds=T, site.weights = perennial.data$site.id)
summary(perennial.final.brt)
# leafN 26%
# RTD 18 %

gbm.plot(perennial.final.brt)
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
# 0.0

# investigation of interactions
gbm.interactions(perennial.final.brt)$rank.list

#### perennial without tree final model ####
set.seed(2023)
perennial.tree.final.brt <- gbm.step(data=perennial.tree,
                                gbm.x = c(11:18),
                                gbm.y = 10,
                                family = "gaussian",
                                tree.complexity = 1,
                                learning.rate = 0.00000000000000001,
                                bag.fraction = 0.75,
                                n.trees = 50,
                                plot.main=T, plot.folds=T)
summary(perennial.tree.final.brt)
# SRL 27%
# leafN 20 %
# SLA 20%
# height 17%

gbm.plot(perennial.tree.final.brt)
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

R2.perennail.brt = 1-(perennial.tree.final.brt$self.statistics$mean.resid/perennial.tree.final.brt$self.statistics$mean.null)
# 0.0

# investigation of interactions
gbm.interactions(perennial.tree.final.brt)$rank.list

#### ggBRT ####

install.packages("devtools") # in case "devtools" has not already been installed
devtools::install_github("JBjouffray/ggBRT") # will take several minutes to install
library(ggBRT)
# https://jbjouffray.github.io/ggBRT/ggBRT_Tutorial.html

ggPerformance(all=all.final.brt,complete=comp.final.brt,no.trees.all=tree.final.brt,annual=annual.final.brt,
              perennial=perennial.final.brt, perennial.no.tree=perennial.tree.final.brt)

ggInfluence(all.final.brt, signif = TRUE)
ggInfluence(comp.final.brt, signif = TRUE)
ggInfluence(tree.final.brt, signif = TRUE)
ggInfluence(annual.final.brt, signif = TRUE)
ggInfluence(perennial.final.brt, signif = TRUE)
ggInfluence(perennial.tree.final.brt, signif = TRUE)


ggMultiInfluence(all=all.final.brt,complete=comp.final.brt,no.trees.all=tree.final.brt,annual=annual.final.brt,
                 perennial=perennial.final.brt, perennial.no.tree=perennial.tree.final.brt)

ggPD(all.final.brt,rug = T)
ggPD(annual.final.brt,rug = T)
ggPDfit(all.final.brt)
ggPDfit(annual.final.brt)

ggInteract_list(all.final.brt,index = F)
ggInteract_list(comp.final.brt,index = F)
ggInteract_list(tree.final.brt,index = F)
ggInteract_list(annual.final.brt,index = F)
ggInteract_list(perennial.final.brt,index = F)
ggInteract_list(perennial.tree.final.brt,index = F)


ggInteract_2D(gbm.object = annual.final.brt,x="SRL_m.g",y="leafN.mg.g",col.gradient = c("white","#5576AF"),
              show.dot = T,col.dot = "grey20",alpha.dot = 0.5,cex.dot = 0.2,
              label.contour = T,col.contour = "#254376",show.axis = T,legend = T)


ggInteract_3D(annual.final.brt,x="SRL_m.g",y="leafN.mg.g")
ggInteract_3D(annual.final.brt,x="leafN.mg.g",y="root.depth_m")
ggInteract_3D(annual.final.brt,x="rootDiam.mm",y="SRL_m.g")

