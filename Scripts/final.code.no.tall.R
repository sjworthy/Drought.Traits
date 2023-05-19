# Script for analysis of Drought Net data
# excluding species with outlier heights
# heights > 12.2, removes 20 species

library(dismo)
library(gbm)

trait.data = read.csv("./Formatted.Data/trait.species.trt.yr1.final.no.woody.csv")

# checking for outliers
# all outliers removed manually from trait.species.trt.yr1.final
shapiro.test(trait.data$leafN.mg.g)
hist(trait.data$leafN.mg.g)
boxplot(trait.data$leafN.mg.g)
mean = mean(trait.data$leafN.mg.g, na.rm = TRUE)
std = sd(trait.data$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$leafN.mg.g[which(trait.data$leafN.mg.g <Tmin | trait.data$leafN.mg.g > Tmax)]
# removed leafN > 58.3

shapiro.test(trait.data$height.m)
hist(trait.data$height.m)
boxplot(trait.data$height.m)
mean = mean(trait.data$height.m, na.rm = TRUE)
std = sd(trait.data$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$height.m[which(trait.data$height.m <Tmin | trait.data$height.m > Tmax)]
# removed height > 4.5

shapiro.test(trait.data$rootN.mg.g)
hist(trait.data$rootN.mg.g)
boxplot(trait.data$rootN.mg.g)
mean = mean(trait.data$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$rootN.mg.g[which(trait.data$rootN.mg.g <Tmin | trait.data$rootN.mg.g > Tmax)]
# remove rootN > 33.33

shapiro.test(trait.data$SLA_m2.kg)
hist(trait.data$SLA_m2.kg)
boxplot(trait.data$SLA_m2.kg)
mean = mean(trait.data$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$SLA_m2.kg[which(trait.data$SLA_m2.kg <Tmin | trait.data$SLA_m2.kg > Tmax)]
# remove SLA > 53.4

shapiro.test(trait.data$root.depth_m)
hist(trait.data$root.depth_m)
boxplot(trait.data$root.depth_m)
mean = mean(trait.data$root.depth_m, na.rm = TRUE)
std = sd(trait.data$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$root.depth_m[which(trait.data$root.depth_m <Tmin | trait.data$root.depth_m > Tmax)]
# remove depth > 3.5

shapiro.test(trait.data$RTD.g.cm3)
hist(trait.data$RTD.g.cm3)
boxplot(trait.data$RTD.g.cm3)
mean = mean(trait.data$RTD.g.cm3, na.rm = TRUE)
std = sd(trait.data$RTD.g.cm3, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$RTD.g.cm3[which(trait.data$RTD.g.cm3 <Tmin | trait.data$RTD.g.cm3 > Tmax)]
# remove RTD > 0.92

shapiro.test(trait.data$SRL_m.g)
hist(trait.data$SRL_m.g)
boxplot(trait.data$SRL_m.g)
mean = mean(trait.data$SRL_m.g, na.rm = TRUE)
std = sd(trait.data$SRL_m.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$SRL_m.g[which(trait.data$SRL_m.g <Tmin | trait.data$SRL_m.g > Tmax)]
# remove SRL > 601

shapiro.test(trait.data$rootDiam.mm)
hist(trait.data$rootDiam.mm)
boxplot(trait.data$rootDiam.mm)
mean = mean(trait.data$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
trait.data$rootDiam.mm[which(trait.data$rootDiam.mm <Tmin | trait.data$rootDiam.mm > Tmax)]
# remove diam > 27

trait.data = read.csv("./Formatted.Data/trait.species.trt.yr1.final.no.woody.csv")
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

all.data = merge(cover.data, trait.data.3, by="Taxon") # 1055 data points with site

# test for correlation
cor.test(all.data$leafN.mg.g, all.data$height.m)
cor.test(all.data$leafN.mg.g, all.data$rootN.mg.g) # correlated
cor.test(all.data$leafN.mg.g, all.data$SLA_m2.kg) # correlated
cor.test(all.data$leafN.mg.g, all.data$root.depth_m)
cor.test(all.data$leafN.mg.g, all.data$RTD.g.cm3) # correlated
cor.test(all.data$leafN.mg.g, all.data$SRL_m.g) # correlated
cor.test(all.data$leafN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$height.m, all.data$rootN.mg.g)
cor.test(all.data$height.m, all.data$SLA_m2.kg) # correlated
cor.test(all.data$height.m, all.data$root.depth_m)
cor.test(all.data$height.m, all.data$RTD.g.cm3) # correlated
cor.test(all.data$height.m, all.data$SRL_m.g) # correlated
cor.test(all.data$height.m, all.data$rootDiam.mm)
cor.test(all.data$rootN.mg.g, all.data$SLA_m2.kg) # correlated
cor.test(all.data$rootN.mg.g, all.data$root.depth_m)
cor.test(all.data$rootN.mg.g, all.data$RTD.g.cm3) # correlated
cor.test(all.data$rootN.mg.g, all.data$SRL_m.g)
cor.test(all.data$rootN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$SLA_m2.kg, all.data$root.depth_m) # correlated
cor.test(all.data$SLA_m2.kg, all.data$RTD.g.cm3) # correlated
cor.test(all.data$SLA_m2.kg, all.data$SRL_m.g) # correlated
cor.test(all.data$SLA_m2.kg, all.data$rootDiam.mm)
cor.test(all.data$root.depth_m, all.data$RTD.g.cm3) # correlated
cor.test(all.data$root.depth_m, all.data$SRL_m.g) # correlated
cor.test(all.data$root.depth_m, all.data$rootDiam.mm) # correlated
cor.test(all.data$RTD.g.cm3, all.data$SRL_m.g) # correlated
cor.test(all.data$RTD.g.cm3, all.data$rootDiam.mm) # correlated
cor.test(all.data$SRL_m.g, all.data$rootDiam.mm)

# checking for normality

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

# consider log transforming
shapiro.test(all.data$mean.cover.response)
hist(all.data$mean.cover.response)
range(all.data$mean.cover.response)

all.data$log.mean.ctrl.cover = log(all.data$mean.ctrl.cover + 0.1)
all.data$log.mean.drt.cover = log(all.data$mean.drt.cover + 0.1)
all.data$log.mean.cover.response = all.data$log.mean.drt.cover - all.data$log.mean.ctrl.cover
all.data$scale.mean.cover.response = as.vector(scale(all.data$log.mean.cover.response))


# log and scale values to transform if non-normal, , mean of 0, sd = 1
all.data$log.leafN = log(all.data$leafN.mg.g)
all.data$log.height = log(all.data$height.m)
all.data$log.rootN = log(all.data$rootN.mg.g)
all.data$log.SLA = log(all.data$SLA_m2.kg)
all.data$log.root.depth = log(all.data$root.depth_m)
all.data$log.RTD = log(all.data$RTD.g.cm3)
all.data$log.SRL = log(all.data$SRL_m.g)
all.data$log.root.diam = log(all.data$rootDiam.mm)

all.data$log.leafN = as.vector(scale(log(all.data$leafN.mg.g)))
all.data$log.height = as.vector(scale(log(all.data$height.m)))
all.data$log.rootN = as.vector(scale(log(all.data$rootN.mg.g)))
all.data$log.SLA = as.vector(scale(log(all.data$SLA_m2.kg)))
all.data$log.root.depth = as.vector(scale(log(all.data$root.depth_m)))
all.data$log.RTD = as.vector(scale(log(all.data$RTD.g.cm3)))
all.data$log.SRL = as.vector(scale(log(all.data$SRL_m.g)))
all.data$log.root.diam = as.vector(scale(log(all.data$rootDiam.mm)))

# change site code to numeric, continuous vector
all.data$site.id = as.numeric(as.factor(all.data$site_code))

non.woody=subset(all.data, !all.data$functional_group == "WOODY")
non.woody$log.leafN = as.vector(scale(log(non.woody$leafN.mg.g)))
non.woody$log.height = as.vector(scale(log(non.woody$height.m)))
non.woody$log.rootN = as.vector(scale(log(non.woody$rootN.mg.g)))
non.woody$log.SLA = as.vector(scale(log(non.woody$SLA_m2.kg)))
non.woody$log.root.depth = as.vector(scale(log(non.woody$root.depth_m)))
non.woody$log.RTD = as.vector(scale(log(non.woody$RTD.g.cm3)))
non.woody$log.SRL = as.vector(scale(log(non.woody$SRL_m.g)))
non.woody$log.root.diam = as.vector(scale(log(non.woody$rootDiam.mm)))

non.woody$log.mean.ctrl.cover = log(non.woody$mean.ctrl.cover + 0.1)
non.woody$log.mean.drt.cover = log(non.woody$mean.drt.cover + 0.1)
non.woody$log.mean.cover.response = non.woody$log.mean.drt.cover - non.woody$log.mean.ctrl.cover
non.woody$scale.mean.cover.response = as.vector(scale(non.woody$log.mean.cover.response))


non.woody$site.id = as.numeric(as.factor(non.woody$site_code))

# determining best learning rate to generate 1000 trees
# start learning rate at 0.1

all.brt.1=gbm.step(data=all.data, gbm.x = c(22:29), gbm.y=33,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.00001, tolerance.method = "fixed",
                   tolerance = 0.0001, site.weights = all.data$site.id,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE)

all.brt.1=gbm.step(data=non.woody, gbm.x = c(22:29), gbm.y=34,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.00001,tolerance.method = "fixed",
                   tolerance = 0.0001,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE, site.weights = non.woody$site.id)
# 0.000001 worked

set.seed(2023)
all.brt.1=gbm.step(data=all.data, gbm.x = c(11:18), gbm.y=10,
                   family = "gaussian", tree.complexity = 1, learning.rate = 0.00001, site.weights = all.data$site.id,
                   bag.fraction = 0.75, n.trees = 50, verbose = FALSE)


# initial model evaluation
# plots of predicted versus observed, estimate R2

# which predictors are most important
summary(all.brt.1)
# RTD 92%

all.predict.brt=predict(all.brt.1, n.trees = all.brt.1$n.trees)
observed = all.data$mean.cover.response
plot(all.predict.brt~observed)
abline(0, 1, col = 2)

R2.all.brt = 1-(all.brt.1$self.statistics$mean.resid/all.brt.1$self.statistics$mean.null)
# 0.0033





