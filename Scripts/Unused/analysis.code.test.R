# Script for analysis of Drought Net data

install.packages("rpart")
library(rpart)
install.packages("rpart.plot")
library(rpart.plot)
install.packages("dismo")
library(dismo)
install.packages("gbm")
library(gbm)

# read in test data frame
# don't forget to set your working directory to where the file is stored

test.data=read.csv("test.data.csv")

# prior to running below code, we need to:
# test for correlation between predictors
# test for normality of data, if non-normal transform
# scale data for comparisons

# test for correlation
cor.test()

# test for normality
shapiro.test()

# log transform if non-normal
log()

# scale all data, mean of 0, sd = 1
scale()

#### basic regression tree analysis ####
# page 84

reg.tree.1=rpart(cover.change ~ total.leaf.area + total.leaf.mass + 
                   stem.mass + fine.root.mass + total.Length.cm., data=test.data)
rpart.plot(reg.tree.1)
summary(reg.tree.1)
# very deep tree that is overfitted, page 86

#### boosted regression trees ####
# to avoid overfitting, regression trees are pruned by means of cross validation

brt.1=rpart(cover.change ~ total.leaf.area + total.leaf.mass + 
                   stem.mass + fine.root.mass + total.Length.cm., data=test.data,
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

brt.4=gbm.step(data=test.data, gbm.x = c(4,5,8,11,12), gbm.y=3,
               family = "gaussian", tree.complexity = 1, learning.rate = 0.1,
               bag.fraction = 0.75, n.trees = 50, verbose = FALSE)

# change learning rate to 0.01

brt.5=gbm.step(data=test.data, gbm.x = c(4,5,8,11,12), gbm.y=3,
               family = "gaussian", tree.complexity = 1, learning.rate = 0.01,
               bag.fraction = 0.75, n.trees = 50, verbose = FALSE)

# change learning rate to 0.001

brt.6=gbm.step(data=test.data, gbm.x = c(4,5,8,11,12), gbm.y=3,
               family = "gaussian", tree.complexity = 1, learning.rate = 0.001,
               bag.fraction = 0.75, n.trees = 50, verbose = FALSE)


# model evaluation
# plots of predicted versus observed, estimate R2

# model with learning rate = 0.1
predicted.brt.4=predict(brt.4, n.trees = brt.4$n.trees)
observed = test.data$cover.change
plot(predicted.brt.4~observed)
abline(0, 1, col = 2)

R2.brt.4 = 1-(brt.4$self.statistics$mean.resid/brt.4$self.statistics$mean.null)
# 0.5553155 # best model

# model with learning rate = 0.01
predicted.brt.5=predict(brt.5, n.trees = brt.5$n.trees)
observed = test.data$cover.change
plot(predicted.brt.5~observed)
abline(0, 1, col = 2)

R2.brt.5 = 1-(brt.5$self.statistics$mean.resid/brt.5$self.statistics$mean.null)
# 0.4902875

# model with learning rate = 0.001
predicted.brt.6=predict(brt.6, n.trees = brt.6$n.trees)
observed = test.data$cover.change
plot(predicted.brt.6~observed)
abline(0, 1, col = 2)

R2.brt.6 = 1-(brt.6$self.statistics$mean.resid/brt.6$self.statistics$mean.null)
# 0.3987002

# which predictors are most important
summary(brt.4)
# total leaf area most important

# determining best tree complexity, pag 104-105
# add silent = T so it doesn't print everything
# removed silent = T because wasn't working and error told me to change learning rate
# learning rate 0.1 only works with tree complexity = 1
# learnign rate 0.01 works for tree complexity = 1-9

R2Obs.meadows <- list()
importancePred.meadows <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.meadows[[tcomp]] <- numeric(nreps)
  importancePred.meadows[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:5),
                                                       ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.meadows <- gbm.step(data=test.data,
                            gbm.x = c(4,5,8,11,12),
                            gbm.y = 3,
                            family = "gaussian",
                            tree.complexity = tcomp,
                            learning.rate = 0.01,
                            bag.fraction = 0.75,
                            n.trees = 50,
                            plot.main=F, plot.folds=F)
    #R2 adj:
    R2Obs.meadows[[tcomp]][i] <- 1 - (BRT.meadows$self.statistics$mean.resid /
                                        BRT.meadows$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.meadows[[tcomp]]) <- sort(rownames(summary(BRT.meadows)))
    }
    importancePred.meadows[[tcomp]][, i] <-
      summary(BRT.meadows)[rownames(importancePred.meadows[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.meadows, mean)
sds <- sapply(R2Obs.meadows, sd)
plot(1:length(R2Obs.meadows), means, ylim = c(0.40, 0.80), ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.meadows)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# lowering learning rate to try tree complexity = 10

R2Obs.meadows <- list()
importancePred.meadows <- list()
nreps <- 10 #number of simulations
for (tcomp in 1:10) {
  R2Obs.meadows[[tcomp]] <- numeric(nreps)
  importancePred.meadows[[tcomp]] <- data.frame(matrix(NA, nrow = length(1:5),
                                                       ncol = nreps))
  for(i in 1:nreps){
    if (i == 1) {
      cat(paste("Starting tc =", tcomp, "\n"))
    }
    BRT.meadows <- gbm.step(data=test.data,
                            gbm.x = c(4,5,8,11,12),
                            gbm.y = 3,
                            family = "gaussian",
                            tree.complexity = tcomp,
                            learning.rate = 0.001,
                            bag.fraction = 0.75,
                            n.trees = 50,
                            plot.main=F, plot.folds=F)
    #R2 adj:
    R2Obs.meadows[[tcomp]][i] <- 1 - (BRT.meadows$self.statistics$mean.resid /
                                        BRT.meadows$self.statistics$mean.null)
    if (i == 1) {
      rownames(importancePred.meadows[[tcomp]]) <- sort(rownames(summary(BRT.meadows)))
    }
    importancePred.meadows[[tcomp]][, i] <-
      summary(BRT.meadows)[rownames(importancePred.meadows[[tcomp]]), ]$rel.inf
  }
}

# examine how R2 improves with tree complexity, pag 105

means <- sapply(R2Obs.meadows, mean)
sds <- sapply(R2Obs.meadows, sd)
plot(1:length(R2Obs.meadows), means, ylim = c(0.40, 0.80), ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.meadows)) {
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}

# to select tree complexity using ANOVA on R2 with Tukey test

nreps <- 10
tcFactor <- as.factor(rep(1:10, each = nreps))
R2Vector <- unlist(R2Obs.meadows)
model <- lm(R2Vector ~ tcFactor)
library(multcomp)
TukeyModel <- glht(model, linfct = mcp(tcFactor = "Tukey"))
TukeyLetters <- cld(TukeyModel)$mcletters$Letters

plot(1:length(R2Obs.meadows), means, ylim=c(0.30, 0.80), ylab = "R squared",
     xlab = "Model complexity")
for (i in 1:length(R2Obs.meadows)){
  arrows(x0 = i, x1 = i, y0 = means[i] - sds[i], y1 = means[i] + sds[i],
         angle = 90, code = 3, length = 0.1)
}
text(x = 1:length(R2Obs.meadows), y = 0.77, labels = TukeyLetters)

# best tree complexity is 9? R2 = 0.75

# how important are the different predictors across the simulations
sort(rowMeans(importancePred.meadows[[9]]), decreasing = T)
# total leaf area most important

# examine shape of cover.change to total leaf area using partial dependence plot
# all other explanatory variables are fixed

BRT.final <- gbm.step(data=test.data,
                        gbm.x = c(4,5,8,11,12),
                        gbm.y = 3,
                        family = "gaussian",
                        tree.complexity = 9,
                        learning.rate = 0.001,
                        bag.fraction = 0.75,
                        n.trees = 50,
                        plot.main=F, plot.folds=F)

plot.gbm(BRT.final, i.var = c("total.leaf.area"))
# cover change increases with total.leaf.area following a saturating relationship

# investigation of interactions
gbm.interactions(BRT.final)$rank.list

# significant interactions
# total.length.cm x total.leaf.mass
# stem.mass x total.leaf.mass

# plot of interaction
# max values = 0.22

gbm.perspec(BRT.final,
            5, # total.length
            1, #total.leaf.mass
            y.label = "mass",
            x.label = "root length",
            smooth = "average")

gbm.perspec(BRT.final,
            3, # stem mass
            1, #total.leaf.mass
            y.label = "stem mass",
            x.label = "root length",
            smooth = "average")

# in place of random effects in linear models,
# include site as fixed effect
# not best option, investigate other options
# Package ‘gpboost’ 
# https://github.com/fabsig/GPBoost/tree/master/R-package


