# Script for analysis of Drought Net data

library(dismo)
library(gbm)
library(corrplot)

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
sort(trait.data.2$leafN.mg.g[which(trait.data.2$leafN.mg.g <Tmin | trait.data.2$leafN.mg.g > Tmax)])
# removed leafN 3.60, 38.6 - 67.58333

hist(trait.data.2$height.m)
boxplot(trait.data.2$height.m)
mean = mean(trait.data.2$height.m, na.rm = TRUE)
std = sd(trait.data.2$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.2$height.m[which(trait.data.2$height.m <Tmin | trait.data.2$height.m > Tmax)])
# removed height > 10 - 46.30267

hist(trait.data.2$rootN.mg.g)
boxplot(trait.data.2$rootN.mg.g)
mean = mean(trait.data.2$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data.2$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.2$rootN.mg.g[which(trait.data.2$rootN.mg.g <Tmin | trait.data.2$rootN.mg.g > Tmax)])
# remove rootN 26.12779 - 39.26274

hist(trait.data.2$SLA_m2.kg)
boxplot(trait.data.2$SLA_m2.kg)
mean = mean(trait.data.2$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data.2$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.2$SLA_m2.kg[which(trait.data.2$SLA_m2.kg <Tmin | trait.data.2$SLA_m2.kg > Tmax)])
# remove SLA 41.82567 - 98.20000

hist(trait.data.2$root.depth_m)
boxplot(trait.data.2$root.depth_m)
mean = mean(trait.data.2$root.depth_m, na.rm = TRUE)
std = sd(trait.data.2$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.2$root.depth_m[which(trait.data.2$root.depth_m <Tmin | trait.data.2$root.depth_m > Tmax)])
# remove depth 2.5 - 11.625000

hist(trait.data.2$RTD.g.cm3)
boxplot(trait.data.2$RTD.g.cm3)
mean = mean(trait.data.2$RTD.g.cm3, na.rm = TRUE)
std = sd(trait.data.2$RTD.g.cm3, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.2$RTD.g.cm3[which(trait.data.2$RTD.g.cm3 <Tmin | trait.data.2$RTD.g.cm3 > Tmax)])
# remove RTD 0.6342500 - 1.1945504

hist(trait.data.2$SRL_m.g)
boxplot(trait.data.2$SRL_m.g)
mean = mean(trait.data.2$SRL_m.g, na.rm = TRUE)
std = sd(trait.data.2$SRL_m.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.2$SRL_m.g[which(trait.data.2$SRL_m.g <Tmin | trait.data.2$SRL_m.g > Tmax)])
# remove SRL 437.3523 - 929.9185

hist(trait.data.2$rootDiam.mm)
boxplot(trait.data.2$rootDiam.mm)
mean = mean(trait.data.2$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data.2$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.2$rootDiam.mm[which(trait.data.2$rootDiam.mm <Tmin | trait.data.2$rootDiam.mm > Tmax)])
# remove diam 27.13865 - 81.81050 

#### read in new trait data without outliers and cover data ####

trait.data.new = read.csv("./Formatted.Data/trait.species.trt.yr1.outlier.std2.csv")
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")

# subset traits so they must have SLA
trait.data.2 = trait.data.new[,c(1,7,8,10,12,14,15,18,20,27:29)]
trait.data.3 = subset(trait.data.2, trait.data.2$SLA_m2.kg > 0 ) # 646 data points, 

table(is.na(trait.data.3$leafN.mg.g)) # 204 missing
table(is.na(trait.data.3$height.m)) # 93 missing
table(is.na(trait.data.3$rootN.mg.g)) # 433 missing
table(is.na(trait.data.3$root.depth_m)) # 270 missing
table(is.na(trait.data.3$RTD.g.cm3)) # 386 missing
table(is.na(trait.data.3$SRL_m.g)) # 334 missing
table(is.na(trait.data.3$rootDiam.mm)) # 324 missing

# remove individuals where 0 in control or drought

cover.2 = subset(cover.data, cover.data$mean.ctrl.cover > 0)
cover.3 = subset(cover.2, cover.2$mean.drt.cover > 0)

# merge trait data and cover data

all.data = merge(cover.3, trait.data.3, by="Taxon") # 706 data points with site

#### Test for correlation of ALL data ####
cor.test(all.data$leafN.mg.g, all.data$height.m)
cor.test(all.data$leafN.mg.g, all.data$rootN.mg.g) # correlated r = 0.52
cor.test(all.data$leafN.mg.g, all.data$SLA_m2.kg) # correlated r = 0.37
cor.test(all.data$leafN.mg.g, all.data$root.depth_m)
cor.test(all.data$leafN.mg.g, all.data$RTD.g.cm3) # correlated r = -0.20
cor.test(all.data$leafN.mg.g, all.data$SRL_m.g) # correlated r = 0.28
cor.test(all.data$leafN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$height.m, all.data$rootN.mg.g)
cor.test(all.data$height.m, all.data$SLA_m2.kg)
cor.test(all.data$height.m, all.data$root.depth_m)
cor.test(all.data$height.m, all.data$RTD.g.cm3) # correlated r = -0.11
cor.test(all.data$height.m, all.data$SRL_m.g) # correlated r = -0.10
cor.test(all.data$height.m, all.data$rootDiam.mm)
cor.test(all.data$rootN.mg.g, all.data$SLA_m2.kg) # correlated r = 0.26
cor.test(all.data$rootN.mg.g, all.data$root.depth_m)
cor.test(all.data$rootN.mg.g, all.data$RTD.g.cm3) # correlated r = -0.13
cor.test(all.data$rootN.mg.g, all.data$SRL_m.g)
cor.test(all.data$rootN.mg.g, all.data$rootDiam.mm)
cor.test(all.data$SLA_m2.kg, all.data$root.depth_m) # correlated r = -0.12
cor.test(all.data$SLA_m2.kg, all.data$RTD.g.cm3) # correlated r = -0.26
cor.test(all.data$SLA_m2.kg, all.data$SRL_m.g) # correlated r = 0.44
cor.test(all.data$SLA_m2.kg, all.data$rootDiam.mm)
cor.test(all.data$root.depth_m, all.data$RTD.g.cm3) # correlated r = 0.13
cor.test(all.data$root.depth_m, all.data$SRL_m.g) # correlated r =-0.12
cor.test(all.data$root.depth_m, all.data$rootDiam.mm) #  correlated r = 0.25
cor.test(all.data$RTD.g.cm3, all.data$SRL_m.g) # correlated r = -0.28
cor.test(all.data$RTD.g.cm3, all.data$rootDiam.mm)
cor.test(all.data$SRL_m.g, all.data$rootDiam.mm)

cor.mat = cor(all.data[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors for ALL data ####

cor.test(all.data$mean.cover.response, all.data$leafN.mg.g) # r = 0.08
cor.test(all.data$mean.cover.response, all.data$height.m) # r = -0.02
cor.test(all.data$mean.cover.response, all.data$rootN.mg.g) # r = -0.01
cor.test(all.data$mean.cover.response, all.data$SLA_m2.kg) # r = 0.04
cor.test(all.data$mean.cover.response, all.data$root.depth_m) # r = 0.05
cor.test(all.data$mean.cover.response, all.data$RTD.g.cm3) # r = 0.09
cor.test(all.data$mean.cover.response, all.data$SRL_m.g) # r = -0.06
cor.test(all.data$mean.cover.response, all.data$rootDiam.mm) # r = 0.07

#### Data set with Woody removed ####

no.trees = subset(all.data, !all.data$functional_group == "WOODY") # 618 data points

#### test for correlation of ALL data WITHOUT WOODY ####
cor.test(no.trees$leafN.mg.g, no.trees$height.m)
cor.test(no.trees$leafN.mg.g, no.trees$rootN.mg.g) # correlated r = 0.52
cor.test(no.trees$leafN.mg.g, no.trees$SLA_m2.kg) # correlated r = 0.36
cor.test(no.trees$leafN.mg.g, no.trees$root.depth_m)
cor.test(no.trees$leafN.mg.g, no.trees$RTD.g.cm3) # correlated r = -0.21
cor.test(no.trees$leafN.mg.g, no.trees$SRL_m.g) # correlated r = 0.25
cor.test(no.trees$leafN.mg.g, no.trees$rootDiam.mm)
cor.test(no.trees$height.m, no.trees$rootN.mg.g)
cor.test(no.trees$height.m, no.trees$SLA_m2.kg) # correlated r = -0.14
cor.test(no.trees$height.m, no.trees$root.depth_m)
cor.test(no.trees$height.m, no.trees$RTD.g.cm3) # correlated r = -0.15
cor.test(no.trees$height.m, no.trees$SRL_m.g)
cor.test(no.trees$height.m, no.trees$rootDiam.mm)
cor.test(no.trees$rootN.mg.g, no.trees$SLA_m2.kg) # correlated r = 0.25
cor.test(no.trees$rootN.mg.g, no.trees$root.depth_m)
cor.test(no.trees$rootN.mg.g, no.trees$RTD.g.cm3) # correlated r = -0.15
cor.test(no.trees$rootN.mg.g, no.trees$SRL_m.g)
cor.test(no.trees$rootN.mg.g, no.trees$rootDiam.mm)
cor.test(no.trees$SLA_m2.kg, no.trees$root.depth_m) # correlated r = -0.13
cor.test(no.trees$SLA_m2.kg, no.trees$RTD.g.cm3) # correlated r = -0.24
cor.test(no.trees$SLA_m2.kg, no.trees$SRL_m.g) # correlated r = 0.42
cor.test(no.trees$SLA_m2.kg, no.trees$rootDiam.mm)
cor.test(no.trees$root.depth_m, no.trees$RTD.g.cm3) # correlated r = 0.16
cor.test(no.trees$root.depth_m, no.trees$SRL_m.g) # correlated r =-0.13
cor.test(no.trees$root.depth_m, no.trees$rootDiam.mm) #  correlated r = 0.24
cor.test(no.trees$RTD.g.cm3, no.trees$SRL_m.g) # correlated r = -0.23
cor.test(no.trees$RTD.g.cm3, no.trees$rootDiam.mm)
cor.test(no.trees$SRL_m.g, no.trees$rootDiam.mm)

cor.mat = cor(no.trees[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors of ALL data WITHOUT WOODY ####

cor.test(no.trees$mean.cover.response, no.trees$leafN.mg.g) # r = 0.12
cor.test(no.trees$mean.cover.response, no.trees$height.m) # r = -0.06
cor.test(no.trees$mean.cover.response, no.trees$rootN.mg.g) # r = -0.004
cor.test(no.trees$mean.cover.response, no.trees$SLA_m2.kg) # r = 0.05
cor.test(no.trees$mean.cover.response, no.trees$root.depth_m) # r = 0.06
cor.test(no.trees$mean.cover.response, no.trees$RTD.g.cm3) # r = 0.08
cor.test(no.trees$mean.cover.response, no.trees$SRL_m.g) # r = -0.04
cor.test(no.trees$mean.cover.response, no.trees$rootDiam.mm) # r = 0.06

#### Data set with TREE removed ####

trees = subset(all.data, !all.data$local_lifeform == "TREE") # 697 data points

#### test for correlation WITHOUT TREES ####

cor.test(trees$leafN.mg.g, trees$height.m)
cor.test(trees$leafN.mg.g, trees$rootN.mg.g) # correlated r = 0.53
cor.test(trees$leafN.mg.g, trees$SLA_m2.kg) # correlated r = 0.36
cor.test(trees$leafN.mg.g, trees$root.depth_m)
cor.test(trees$leafN.mg.g, trees$RTD.g.cm3) # correlated r = -0.21
cor.test(trees$leafN.mg.g, trees$SRL_m.g) # correlated r = 0.27
cor.test(trees$leafN.mg.g, trees$rootDiam.mm)
cor.test(trees$height.m, trees$rootN.mg.g)
cor.test(trees$height.m, trees$SLA_m2.kg)
cor.test(trees$height.m, trees$root.depth_m)
cor.test(trees$height.m, trees$RTD.g.cm3) # correlated r = -0.11
cor.test(trees$height.m, trees$SRL_m.g) # correlated r = -0.10
cor.test(trees$height.m, trees$rootDiam.mm)
cor.test(trees$rootN.mg.g, trees$SLA_m2.kg) # correlated r = 0.26
cor.test(trees$rootN.mg.g, trees$root.depth_m)
cor.test(trees$rootN.mg.g, trees$RTD.g.cm3) # correlated r = -0.12
cor.test(trees$rootN.mg.g, trees$SRL_m.g)
cor.test(trees$rootN.mg.g, trees$rootDiam.mm)
cor.test(trees$SLA_m2.kg, trees$root.depth_m) # correlated r = -0.12
cor.test(trees$SLA_m2.kg, trees$RTD.g.cm3) # correlated r = -0.26
cor.test(trees$SLA_m2.kg, trees$SRL_m.g) # correlated r = 0.44
cor.test(trees$SLA_m2.kg, trees$rootDiam.mm)
cor.test(trees$root.depth_m, trees$RTD.g.cm3) # correlated r = 0.13
cor.test(trees$root.depth_m, trees$SRL_m.g) # correlated r =-0.11
cor.test(trees$root.depth_m, trees$rootDiam.mm) #  correlated r = 0.25
cor.test(trees$RTD.g.cm3, trees$SRL_m.g) # correlated r = -0.28
cor.test(trees$RTD.g.cm3, trees$rootDiam.mm)
cor.test(trees$SRL_m.g, trees$rootDiam.mm)


cor.mat = cor(trees[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors WITHOUT TREES ####

cor.test(trees$mean.cover.response, trees$leafN.mg.g) # r = 0.09
cor.test(trees$mean.cover.response, trees$height.m) # r = -0.02
cor.test(trees$mean.cover.response, trees$rootN.mg.g) # r = -0.02
cor.test(trees$mean.cover.response, trees$SLA_m2.kg) # r = 0.04
cor.test(trees$mean.cover.response, trees$root.depth_m) # r = 0.05
cor.test(trees$mean.cover.response, trees$RTD.g.cm3) # r = 0.08
cor.test(trees$mean.cover.response, trees$SRL_m.g) # r = -0.06
cor.test(trees$mean.cover.response, trees$rootDiam.mm) # r = 0.07

#### data set split by lifespan ####

# need to read out data and fix the lifespan to particular sites since 
# some species have different lifespan at different sites

# write.csv(all.data, file="./Formatted.Data/all.data.response.0.csv")

all.data.ls = read.csv("./Formatted.Data/all.data.response.0.lifespan.csv", row.names = 1)
table(all.data.ls$local_lifespan)

#### Annual Data ####

annual.data = subset(all.data.ls, all.data.ls$local_lifespan == "ANNUAL") # 127 data points

#### test for correlation with annual data ####
cor.test(annual.data$leafN.mg.g, annual.data$height.m)
cor.test(annual.data$leafN.mg.g, annual.data$rootN.mg.g) # correlated r = 0.43
cor.test(annual.data$leafN.mg.g, annual.data$SLA_m2.kg) # correlated r = 0.28
cor.test(annual.data$leafN.mg.g, annual.data$root.depth_m)
cor.test(annual.data$leafN.mg.g, annual.data$RTD.g.cm3)
cor.test(annual.data$leafN.mg.g, annual.data$SRL_m.g)
cor.test(annual.data$leafN.mg.g, annual.data$rootDiam.mm)
cor.test(annual.data$height.m, annual.data$rootN.mg.g)
cor.test(annual.data$height.m, annual.data$SLA_m2.kg)
cor.test(annual.data$height.m, annual.data$root.depth_m) # correlated r = 0.24
cor.test(annual.data$height.m, annual.data$RTD.g.cm3)
cor.test(annual.data$height.m, annual.data$SRL_m.g)
cor.test(annual.data$height.m, annual.data$rootDiam.mm)
cor.test(annual.data$rootN.mg.g, annual.data$SLA_m2.kg)
cor.test(annual.data$rootN.mg.g, annual.data$root.depth_m)
cor.test(annual.data$rootN.mg.g, annual.data$RTD.g.cm3)
cor.test(annual.data$rootN.mg.g, annual.data$SRL_m.g)
cor.test(annual.data$rootN.mg.g, annual.data$rootDiam.mm)
cor.test(annual.data$SLA_m2.kg, annual.data$root.depth_m)
cor.test(annual.data$SLA_m2.kg, annual.data$RTD.g.cm3)
cor.test(annual.data$SLA_m2.kg, annual.data$SRL_m.g) # correlated r = 0.26
cor.test(annual.data$SLA_m2.kg, annual.data$rootDiam.mm)  # correlated r = 0.27
cor.test(annual.data$root.depth_m, annual.data$RTD.g.cm3)
cor.test(annual.data$root.depth_m, annual.data$SRL_m.g) # correlated r =-0.24
cor.test(annual.data$root.depth_m, annual.data$rootDiam.mm)
cor.test(annual.data$RTD.g.cm3, annual.data$SRL_m.g) # correlated r = -0.30
cor.test(annual.data$RTD.g.cm3, annual.data$rootDiam.mm)
cor.test(annual.data$SRL_m.g, annual.data$rootDiam.mm)

cor.mat = cor(annual.data[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(annual.data$mean.cover.response, annual.data$leafN.mg.g) # r = 0.25
cor.test(annual.data$mean.cover.response, annual.data$height.m) # r = 0.05
cor.test(annual.data$mean.cover.response, annual.data$rootN.mg.g) # r = 0.19
cor.test(annual.data$mean.cover.response, annual.data$SLA_m2.kg) # r = 0.06
cor.test(annual.data$mean.cover.response, annual.data$root.depth_m) # r = 0.21
cor.test(annual.data$mean.cover.response, annual.data$RTD.g.cm3) # r = 0.18
cor.test(annual.data$mean.cover.response, annual.data$SRL_m.g) # r = -0.22
cor.test(annual.data$mean.cover.response, annual.data$rootDiam.mm) # r = 0.10

#### Perennial Data ####

perennial.data = subset(all.data.ls, all.data.ls$local_lifespan == "PERENNIAL") # 563 data points

#### test for correlation with PERENNIAL data ####

cor.test(perennial.data$leafN.mg.g, perennial.data$height.m)
cor.test(perennial.data$leafN.mg.g, perennial.data$rootN.mg.g) # correlated r = 0.51
cor.test(perennial.data$leafN.mg.g, perennial.data$SLA_m2.kg) # correlated r = 0.37
cor.test(perennial.data$leafN.mg.g, perennial.data$root.depth_m)
cor.test(perennial.data$leafN.mg.g, perennial.data$RTD.g.cm3) # correlated r = -0.22
cor.test(perennial.data$leafN.mg.g, perennial.data$SRL_m.g) # correlated r = 0.30
cor.test(perennial.data$leafN.mg.g, perennial.data$rootDiam.mm)
cor.test(perennial.data$height.m, perennial.data$rootN.mg.g)
cor.test(perennial.data$height.m, perennial.data$SLA_m2.kg)
cor.test(perennial.data$height.m, perennial.data$root.depth_m)
cor.test(perennial.data$height.m, perennial.data$RTD.g.cm3)
cor.test(perennial.data$height.m, perennial.data$SRL_m.g) # correlated r = -0.14
cor.test(perennial.data$height.m, perennial.data$rootDiam.mm)
cor.test(perennial.data$rootN.mg.g, perennial.data$SLA_m2.kg) # correlated r = 0.28
cor.test(perennial.data$rootN.mg.g, perennial.data$root.depth_m)
cor.test(perennial.data$rootN.mg.g, perennial.data$RTD.g.cm3)
cor.test(perennial.data$rootN.mg.g, perennial.data$SRL_m.g)
cor.test(perennial.data$rootN.mg.g, perennial.data$rootDiam.mm)
cor.test(perennial.data$SLA_m2.kg, perennial.data$root.depth_m) # correlated r = -0.12
cor.test(perennial.data$SLA_m2.kg, perennial.data$RTD.g.cm3) # correlated r = -0.30
cor.test(perennial.data$SLA_m2.kg, perennial.data$SRL_m.g) # correlated r = 0.47
cor.test(perennial.data$SLA_m2.kg, perennial.data$rootDiam.mm) # correlated r = -0.15
cor.test(perennial.data$root.depth_m, perennial.data$RTD.g.cm3) # correlated r = 0.13
cor.test(perennial.data$root.depth_m, perennial.data$SRL_m.g)
cor.test(perennial.data$root.depth_m, perennial.data$rootDiam.mm) #  correlated r = 0.31
cor.test(perennial.data$RTD.g.cm3, perennial.data$SRL_m.g) # correlated r = -0.31
cor.test(perennial.data$RTD.g.cm3, perennial.data$rootDiam.mm)
cor.test(perennial.data$SRL_m.g, perennial.data$rootDiam.mm)

cor.mat = cor(perennial.data[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors for PERENNIALS ####

cor.test(perennial.data$mean.cover.response, perennial.data$leafN.mg.g) # r = 0.05
cor.test(perennial.data$mean.cover.response, perennial.data$height.m) # r = -0.03
cor.test(perennial.data$mean.cover.response, perennial.data$rootN.mg.g) # r = -0.04
cor.test(perennial.data$mean.cover.response, perennial.data$SLA_m2.kg) # r = 0.02
cor.test(perennial.data$mean.cover.response, perennial.data$root.depth_m) # r = 0.03
cor.test(perennial.data$mean.cover.response, perennial.data$RTD.g.cm3) # r = 0.09
cor.test(perennial.data$mean.cover.response, perennial.data$SRL_m.g) # r = -0.02
cor.test(perennial.data$mean.cover.response, perennial.data$rootDiam.mm) # r = 0.07

#### Perennial data without woody functional group ####

perennial.tree = subset(perennial.data, !perennial.data$functional_group == "WOODY") # 478 data points

#### test for correlation with PERENNIAL data WITHOUT WOODY ####
cor.test(perennial.tree$leafN.mg.g, perennial.tree$height.m)
cor.test(perennial.tree$leafN.mg.g, perennial.tree$rootN.mg.g) # correlated r = 0.52
cor.test(perennial.tree$leafN.mg.g, perennial.tree$SLA_m2.kg) # correlated r = 0.39
cor.test(perennial.tree$leafN.mg.g, perennial.tree$root.depth_m)
cor.test(perennial.tree$leafN.mg.g, perennial.tree$RTD.g.cm3) # correlated r = -0.25
cor.test(perennial.tree$leafN.mg.g, perennial.tree$SRL_m.g) # correlated r = 0.29
cor.test(perennial.tree$leafN.mg.g, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$height.m, perennial.tree$rootN.mg.g)
cor.test(perennial.tree$height.m, perennial.tree$SLA_m2.kg) # correlated r = -0.16
cor.test(perennial.tree$height.m, perennial.tree$root.depth_m)
cor.test(perennial.tree$height.m, perennial.tree$RTD.g.cm3) # correlated r = -0.13
cor.test(perennial.tree$height.m, perennial.tree$SRL_m.g) # correlated r = -0.12
cor.test(perennial.tree$height.m, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$SLA_m2.kg) # correlated r = 0.27
cor.test(perennial.tree$rootN.mg.g, perennial.tree$root.depth_m)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$RTD.g.cm3)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$SRL_m.g)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$rootDiam.mm) # correlated r = 0.14
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$root.depth_m) # correlated r = -0.14
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$RTD.g.cm3) # correlated r = -0.29
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$SRL_m.g) # correlated r = 0.46
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$rootDiam.mm) # correlated r = -0.12
cor.test(perennial.tree$root.depth_m, perennial.tree$RTD.g.cm3) # correlated r = 0.18
cor.test(perennial.tree$root.depth_m, perennial.tree$SRL_m.g) # correlated r =-0.12
cor.test(perennial.tree$root.depth_m, perennial.tree$rootDiam.mm) #  correlated r = 0.29
cor.test(perennial.tree$RTD.g.cm3, perennial.tree$SRL_m.g) # correlated r = -0.27
cor.test(perennial.tree$RTD.g.cm3, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$SRL_m.g, perennial.tree$rootDiam.mm)

cor.mat = cor(perennial.tree[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.tree$mean.cover.response, perennial.tree$leafN.mg.g) # r = 0.08
cor.test(perennial.tree$mean.cover.response, perennial.tree$height.m) # r = -0.09
cor.test(perennial.tree$mean.cover.response, perennial.tree$rootN.mg.g) # r = -0.03
cor.test(perennial.tree$mean.cover.response, perennial.tree$SLA_m2.kg) # r = 0.04
cor.test(perennial.tree$mean.cover.response, perennial.tree$root.depth_m) # r = 0.04
cor.test(perennial.tree$mean.cover.response, perennial.tree$RTD.g.cm3) # r = 0.09
cor.test(perennial.tree$mean.cover.response, perennial.tree$SRL_m.g) # r = 0.003
cor.test(perennial.tree$mean.cover.response, perennial.tree$rootDiam.mm) # r = 0.06

#### PERENNIAL without TREE ####

perennial.no.tree = subset(perennial.data, !perennial.data$local_lifeform == "TREE") # 554 data points

#### test for correlation of PERENNIAL data without TREE ####
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$height.m)
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$rootN.mg.g) # correlated r = 0.51
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$SLA_m2.kg) # correlated r = 0.37
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$RTD.g.cm3) # correlated r = -0.23
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$SRL_m.g) # correlated r = 0.30
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$height.m, perennial.no.tree$rootN.mg.g)
cor.test(perennial.no.tree$height.m, perennial.no.tree$SLA_m2.kg)
cor.test(perennial.no.tree$height.m, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$height.m, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$height.m, perennial.no.tree$SRL_m.g) # correlated r = -0.14
cor.test(perennial.no.tree$height.m, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$SLA_m2.kg) # correlated r = 0.28
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$SRL_m.g)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$root.depth_m) # correlated r = -0.12
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$RTD.g.cm3) # correlated r = -0.30
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$SRL_m.g) # correlated r = 0.46
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$rootDiam.mm) # correlated r = -0.14
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$RTD.g.cm3) # correlated r = 0.13
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$SRL_m.g)
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$rootDiam.mm) #  correlated r = 0.30
cor.test(perennial.no.tree$RTD.g.cm3, perennial.no.tree$SRL_m.g) # correlated r = -0.31
cor.test(perennial.no.tree$RTD.g.cm3, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$SRL_m.g, perennial.no.tree$rootDiam.mm)

cor.mat = cor(perennial.no.tree[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$leafN.mg.g) # r = 0.05
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$height.m) # r = -0.02
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$rootN.mg.g) # r = -0.04
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$SLA_m2.kg) # r = 0.03
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$root.depth_m) # r = 0.03
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$RTD.g.cm3) # r = 0.08
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$SRL_m.g) # r = -0.02
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$rootDiam.mm) # r = 0.07

#### GRASS ####

grass = subset(all.data, all.data$functional_group == "GRASS") # 223 data points

#### test for correlation with GRASS ####
cor.test(grass$leafN.mg.g, grass$height.m)
cor.test(grass$leafN.mg.g, grass$rootN.mg.g) # correlated r = 0.21
cor.test(grass$leafN.mg.g, grass$SLA_m2.kg) # correlated r = 0.43
cor.test(grass$leafN.mg.g, grass$root.depth_m)
cor.test(grass$leafN.mg.g, grass$RTD.g.cm3) # correlated r = -0.31
cor.test(grass$leafN.mg.g, grass$SRL_m.g) # correlated r = 0.61
cor.test(grass$leafN.mg.g, grass$rootDiam.mm)  # correlated r = -0.31
cor.test(grass$height.m, grass$rootN.mg.g)  # correlated r = 0.18
cor.test(grass$height.m, grass$SLA_m2.kg)
cor.test(grass$height.m, grass$root.depth_m)
cor.test(grass$height.m, grass$RTD.g.cm3)
cor.test(grass$height.m, grass$SRL_m.g) # correlated r = 0.16
cor.test(grass$height.m, grass$rootDiam.mm)  # correlated r = 0.16
cor.test(grass$rootN.mg.g, grass$SLA_m2.kg)
cor.test(grass$rootN.mg.g, grass$root.depth_m)
cor.test(grass$rootN.mg.g, grass$RTD.g.cm3)
cor.test(grass$rootN.mg.g, grass$SRL_m.g)
cor.test(grass$rootN.mg.g, grass$rootDiam.mm)
cor.test(grass$SLA_m2.kg, grass$root.depth_m) # correlated r = -0.23
cor.test(grass$SLA_m2.kg, grass$RTD.g.cm3) # correlated r = -0.46
cor.test(grass$SLA_m2.kg, grass$SRL_m.g) # correlated r = 0.42
cor.test(grass$SLA_m2.kg, grass$rootDiam.mm)  # correlated r = -0.22
cor.test(grass$root.depth_m, grass$RTD.g.cm3) # correlated r = 0.23
cor.test(grass$root.depth_m, grass$SRL_m.g) # correlated r =-0.17
cor.test(grass$root.depth_m, grass$rootDiam.mm) #  correlated r = 0.33
cor.test(grass$RTD.g.cm3, grass$SRL_m.g) # correlated r = -0.40
cor.test(grass$RTD.g.cm3, grass$rootDiam.mm)
cor.test(grass$SRL_m.g, grass$rootDiam.mm)  # correlated r = -0.43

cor.mat = cor(grass[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(grass$mean.cover.response, grass$leafN.mg.g) # r = 0.06
cor.test(grass$mean.cover.response, grass$height.m) # r = -0.11
cor.test(grass$mean.cover.response, grass$rootN.mg.g) # r = -0.02
cor.test(grass$mean.cover.response, grass$SLA_m2.kg) # r = 0.01
cor.test(grass$mean.cover.response, grass$root.depth_m) # r = 0.12
cor.test(grass$mean.cover.response, grass$RTD.g.cm3) # r = 0.19
cor.test(grass$mean.cover.response, grass$SRL_m.g) # r = -0.09
cor.test(grass$mean.cover.response, grass$rootDiam.mm) # r = 0.12

#### FORB ####

forb = subset(all.data, all.data$functional_group == "FORB") # 317 data points

#### test for correlation with FORB ####
cor.test(forb$leafN.mg.g, forb$height.m)
cor.test(forb$leafN.mg.g, forb$rootN.mg.g)
cor.test(forb$leafN.mg.g, forb$SLA_m2.kg) # correlated r = 0.27
cor.test(forb$leafN.mg.g, forb$root.depth_m) # correlated r = 0.17
cor.test(forb$leafN.mg.g, forb$RTD.g.cm3)
cor.test(forb$leafN.mg.g, forb$SRL_m.g) # correlated r = 0.19
cor.test(forb$leafN.mg.g, forb$rootDiam.mm)
cor.test(forb$height.m, forb$rootN.mg.g) # correlated r = 0.22
cor.test(forb$height.m, forb$SLA_m2.kg) # correlated r = -0.17
cor.test(forb$height.m, forb$root.depth_m)
cor.test(forb$height.m, forb$RTD.g.cm3) # correlated r = -0.31
cor.test(forb$height.m, forb$SRL_m.g) # correlated r = -0.24
cor.test(forb$height.m, forb$rootDiam.mm)
cor.test(forb$rootN.mg.g, forb$SLA_m2.kg) # correlated r = 0.23
cor.test(forb$rootN.mg.g, forb$root.depth_m)
cor.test(forb$rootN.mg.g, forb$RTD.g.cm3)
cor.test(forb$rootN.mg.g, forb$SRL_m.g) # correlated r = 0.30
cor.test(forb$rootN.mg.g, forb$rootDiam.mm)
cor.test(forb$SLA_m2.kg, forb$root.depth_m)
cor.test(forb$SLA_m2.kg, forb$RTD.g.cm3)
cor.test(forb$SLA_m2.kg, forb$SRL_m.g) # correlated r = 0.40
cor.test(forb$SLA_m2.kg, forb$rootDiam.mm)
cor.test(forb$root.depth_m, forb$RTD.g.cm3)
cor.test(forb$root.depth_m, forb$SRL_m.g)
cor.test(forb$root.depth_m, forb$rootDiam.mm) #  correlated r = 0.21
cor.test(forb$RTD.g.cm3, forb$SRL_m.g)
cor.test(forb$RTD.g.cm3, forb$rootDiam.mm)
cor.test(forb$SRL_m.g, forb$rootDiam.mm)

cor.mat = cor(forb[,c(11:18)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(forb$mean.cover.response, forb$leafN.mg.g) # r = 0.19
cor.test(forb$mean.cover.response, forb$height.m) # r = -0.03
cor.test(forb$mean.cover.response, forb$rootN.mg.g) # r = -0.04
cor.test(forb$mean.cover.response, forb$SLA_m2.kg) # r = 0.09
cor.test(forb$mean.cover.response, forb$root.depth_m) # r = 0.04
cor.test(forb$mean.cover.response, forb$RTD.g.cm3) # r = -0.02
cor.test(forb$mean.cover.response, forb$SRL_m.g) # r = -0.03
cor.test(forb$mean.cover.response, forb$rootDiam.mm) # r = 0.02

#### ANNUAL GRASSES ####

annual.grass = subset(annual.data, annual.data$functional_group == "GRASS") # 39 data points

#### ANNUAL FORB ####
annual.forb = subset(annual.data, annual.data$functional_group == "FORB") # 78 data points

#### PERENNIAL GRASS ####

perennial.grass = subset(perennial.data, perennial.data$functional_group == "GRASS") # 182 data points

#### PERENNIAL FORB ####

perennial.forb = subset(perennial.data, perennial.data$functional_group == "FORB") # 231 data points

#### write out the files ####

write.csv(all.data, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/all.data.2.csv")
write.csv(annual.data, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/annual.data.csv")
write.csv(annual.forb, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/annual.forb.csv")
write.csv(annual.grass, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/annual.grass.csv")
write.csv(forb, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/forb.csv")
write.csv(grass, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/grass.csv")
write.csv(no.trees, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/no.trees.csv")
write.csv(perennial.data, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.data.csv")
write.csv(perennial.forb, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.forb.csv")
write.csv(perennial.grass, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.grass.csv")
write.csv(perennial.no.tree, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.no.tree.csv")
write.csv(perennial.tree, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.tree.csv")
write.csv(trees, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/trees.csv")






#### Generating data files of Drought yr 1 versus Drought yr 2 ####
# read in cover year 2

cover.yr.2 = read.csv("./Formatted.Data/cover.trt.y2.csv")

# read in cover data
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")

# merge cover year 1 with cover year 2

years.cover = merge(cover.data, cover.yr.2, by = c("site_code", "Taxon"))
years.cover = years.cover[,c(1,2,8,17)]
colnames(years.cover)[3] = "mean.drt.cover.yr1"
colnames(years.cover)[4] = "mean.drt.cover.yr2"

# remove 0 if wasn't in plot in year 1 and year 2

years.cover.2 = subset(years.cover, !(years.cover$mean.drt.cover.yr1 == 0))
years.cover.3 = subset(years.cover.2, !(years.cover.2$mean.drt.cover.yr2 == 0)) # 774 data points

years.cover.3$diff.drt.trt = years.cover.3$mean.drt.cover.yr2 - years.cover.3$mean.drt.cover.yr1

# merge new cover data with traits

trait.data.new = read.csv("./Formatted.Data/trait.species.trt.yr1.outlier.2.csv")

# subset traits so they must have SLA
trait.data.2 = trait.data.new[,c(1,7,8,10,12,14,15,18,20,27:29)]
trait.data.3 = subset(trait.data.2, trait.data.2$SLA_m2.kg > 0 ) # 646 data points, 

all.data.year2 = merge(years.cover.3, trait.data.3, by="Taxon") # 563 data points

#### test for correlation with all data drt.v.drt ####
cor.test(all.data.year2$leafN.mg.g, all.data.year2$height.m)
cor.test(all.data.year2$leafN.mg.g, all.data.year2$rootN.mg.g) # correlated r = 0.50
cor.test(all.data.year2$leafN.mg.g, all.data.year2$SLA_m2.kg) # correlated r = 0.33
cor.test(all.data.year2$leafN.mg.g, all.data.year2$root.depth_m)
cor.test(all.data.year2$leafN.mg.g, all.data.year2$RTD.g.cm3) # correlated r = -0.19
cor.test(all.data.year2$leafN.mg.g, all.data.year2$SRL_m.g) # correlated r = 0.25
cor.test(all.data.year2$leafN.mg.g, all.data.year2$rootDiam.mm)
cor.test(all.data.year2$height.m, all.data.year2$rootN.mg.g)
cor.test(all.data.year2$height.m, all.data.year2$SLA_m2.kg)
cor.test(all.data.year2$height.m, all.data.year2$root.depth_m)
cor.test(all.data.year2$height.m, all.data.year2$RTD.g.cm3)
cor.test(all.data.year2$height.m, all.data.year2$SRL_m.g)
cor.test(all.data.year2$height.m, all.data.year2$rootDiam.mm) # correlated r = 0.20
cor.test(all.data.year2$rootN.mg.g, all.data.year2$SLA_m2.kg) # correlated r = 0.24
cor.test(all.data.year2$rootN.mg.g, all.data.year2$root.depth_m)
cor.test(all.data.year2$rootN.mg.g, all.data.year2$RTD.g.cm3) # correlated r = -0.14
cor.test(all.data.year2$rootN.mg.g, all.data.year2$SRL_m.g)
cor.test(all.data.year2$rootN.mg.g, all.data.year2$rootDiam.mm)
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$root.depth_m) # correlated r = -0.10
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$RTD.g.cm3) # correlated r = -0.25
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$SRL_m.g) # correlated r = 0.45
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$rootDiam.mm) # correlated r = -0.16
cor.test(all.data.year2$root.depth_m, all.data.year2$RTD.g.cm3) # correlated r = 0.12
cor.test(all.data.year2$root.depth_m, all.data.year2$SRL_m.g)
cor.test(all.data.year2$root.depth_m, all.data.year2$rootDiam.mm) #  correlated r = 0.20
cor.test(all.data.year2$RTD.g.cm3, all.data.year2$SRL_m.g) # correlated r = -0.30
cor.test(all.data.year2$RTD.g.cm3, all.data.year2$rootDiam.mm)
cor.test(all.data.year2$SRL_m.g, all.data.year2$rootDiam.mm)

cor.mat = cor(all.data.year2[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(all.data.year2$diff.drt.trt, all.data.year2$leafN.mg.g) # r = -0.08
cor.test(all.data.year2$diff.drt.trt, all.data.year2$height.m) # r = 0.006
cor.test(all.data.year2$diff.drt.trt, all.data.year2$rootN.mg.g) # r = -0.05
cor.test(all.data.year2$diff.drt.trt, all.data.year2$SLA_m2.kg) # r = -0.06
cor.test(all.data.year2$diff.drt.trt, all.data.year2$root.depth_m) # r = -0.03
cor.test(all.data.year2$diff.drt.trt, all.data.year2$RTD.g.cm3) # r = -0.003
cor.test(all.data.year2$diff.drt.trt, all.data.year2$SRL_m.g) # r = 0.009
cor.test(all.data.year2$diff.drt.trt, all.data.year2$rootDiam.mm) # r = -0.01

#### data set with WOODY removed ####

no.trees = subset(all.data.year2, !all.data.year2$functional_group == "WOODY") # 485 data points

#### test for correlation with all data without woody drt.v.drt ####
cor.test(no.trees$leafN.mg.g, no.trees$height.m) # correlated r = -0.12
cor.test(no.trees$leafN.mg.g, no.trees$rootN.mg.g) # correlated r = 0.50
cor.test(no.trees$leafN.mg.g, no.trees$SLA_m2.kg) # correlated r = 0.33
cor.test(no.trees$leafN.mg.g, no.trees$root.depth_m)
cor.test(no.trees$leafN.mg.g, no.trees$RTD.g.cm3) # correlated r = -0.18
cor.test(no.trees$leafN.mg.g, no.trees$SRL_m.g) # correlated r = 0.22
cor.test(no.trees$leafN.mg.g, no.trees$rootDiam.mm)
cor.test(no.trees$height.m, no.trees$rootN.mg.g)
cor.test(no.trees$height.m, no.trees$SLA_m2.kg) # correlated r = -0.17
cor.test(no.trees$height.m, no.trees$root.depth_m) # correlated r = 0.11
cor.test(no.trees$height.m, no.trees$RTD.g.cm3) # correlated r = -0.12
cor.test(no.trees$height.m, no.trees$SRL_m.g) # correlated r = -0.12
cor.test(no.trees$height.m, no.trees$rootDiam.mm)
cor.test(no.trees$rootN.mg.g, no.trees$SLA_m2.kg) # correlated r = 0.23
cor.test(no.trees$rootN.mg.g, no.trees$root.depth_m)
cor.test(no.trees$rootN.mg.g, no.trees$RTD.g.cm3) # correlated r = -0.17
cor.test(no.trees$rootN.mg.g, no.trees$SRL_m.g)
cor.test(no.trees$rootN.mg.g, no.trees$rootDiam.mm)
cor.test(no.trees$SLA_m2.kg, no.trees$root.depth_m) # correlated r = -0.10
cor.test(no.trees$SLA_m2.kg, no.trees$RTD.g.cm3) # correlated r = -0.21
cor.test(no.trees$SLA_m2.kg, no.trees$SRL_m.g) # correlated r = 0.42
cor.test(no.trees$SLA_m2.kg, no.trees$rootDiam.mm) # correlated r = -0.14
cor.test(no.trees$root.depth_m, no.trees$RTD.g.cm3) # correlated r = 0.13
cor.test(no.trees$root.depth_m, no.trees$SRL_m.g)
cor.test(no.trees$root.depth_m, no.trees$rootDiam.mm) #  correlated r = 0.17
cor.test(no.trees$RTD.g.cm3, no.trees$SRL_m.g) # correlated r = -0.24
cor.test(no.trees$RTD.g.cm3, no.trees$rootDiam.mm)
cor.test(no.trees$SRL_m.g, no.trees$rootDiam.mm)

cor.mat = cor(no.trees[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(no.trees$diff.drt.trt, no.trees$leafN.mg.g) # r = -0.07
cor.test(no.trees$diff.drt.trt, no.trees$height.m) # r = 0.0002
cor.test(no.trees$diff.drt.trt, no.trees$rootN.mg.g) # r = -0.005
cor.test(no.trees$diff.drt.trt, no.trees$SLA_m2.kg) # r = -0.06
cor.test(no.trees$diff.drt.trt, no.trees$root.depth_m) # r = -0.02
cor.test(no.trees$diff.drt.trt, no.trees$RTD.g.cm3) # r = -0.013
cor.test(no.trees$diff.drt.trt, no.trees$SRL_m.g) # r = 0.02
cor.test(no.trees$diff.drt.trt, no.trees$rootDiam.mm) # r = -0.005

#### data set with TREE removed ####

trees = subset(all.data.year2, !all.data.year2$local_lifeform == "TREE") # 559 data points

#### test for correlation with all data without TREE lifeform drt.v.drt ####

cor.test(trees$leafN.mg.g, trees$height.m)
cor.test(trees$leafN.mg.g, trees$rootN.mg.g) # correlated r = 0.50
cor.test(trees$leafN.mg.g, trees$SLA_m2.kg) # correlated r = 0.33
cor.test(trees$leafN.mg.g, trees$root.depth_m)
cor.test(trees$leafN.mg.g, trees$RTD.g.cm3) # correlated r = -0.19
cor.test(trees$leafN.mg.g, trees$SRL_m.g) # correlated r = 0.25
cor.test(trees$leafN.mg.g, trees$rootDiam.mm)
cor.test(trees$height.m, trees$rootN.mg.g)
cor.test(trees$height.m, trees$SLA_m2.kg)
cor.test(trees$height.m, trees$root.depth_m)
cor.test(trees$height.m, trees$RTD.g.cm3)
cor.test(trees$height.m, trees$SRL_m.g) # correlated r = -0.13
cor.test(trees$height.m, trees$rootDiam.mm)
cor.test(trees$rootN.mg.g, trees$SLA_m2.kg) # correlated r = 0.24
cor.test(trees$rootN.mg.g, trees$root.depth_m)
cor.test(trees$rootN.mg.g, trees$RTD.g.cm3) # correlated r = -0.13
cor.test(trees$rootN.mg.g, trees$SRL_m.g)
cor.test(trees$rootN.mg.g, trees$rootDiam.mm)
cor.test(trees$SLA_m2.kg, trees$root.depth_m) # correlated r = -0.10
cor.test(trees$SLA_m2.kg, trees$RTD.g.cm3) # correlated r = -0.25
cor.test(trees$SLA_m2.kg, trees$SRL_m.g) # correlated r = 0.44
cor.test(trees$SLA_m2.kg, trees$rootDiam.mm) # correlated r = -0.16
cor.test(trees$root.depth_m, trees$RTD.g.cm3) # correlated r = 0.12
cor.test(trees$root.depth_m, trees$SRL_m.g)
cor.test(trees$root.depth_m, trees$rootDiam.mm) #  correlated r = 0.20
cor.test(trees$RTD.g.cm3, trees$SRL_m.g) # correlated r = -0.30
cor.test(trees$RTD.g.cm3, trees$rootDiam.mm)
cor.test(trees$SRL_m.g, trees$rootDiam.mm)

cor.mat = cor(trees[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(trees$diff.drt.trt, trees$leafN.mg.g) # r = -0.08
cor.test(trees$diff.drt.trt, trees$height.m) # r = 0.01
cor.test(trees$diff.drt.trt, trees$rootN.mg.g) # r = -0.05
cor.test(trees$diff.drt.trt, trees$SLA_m2.kg) # r = -0.06
cor.test(trees$diff.drt.trt, trees$root.depth_m) # r = -0.03
cor.test(trees$diff.drt.trt, trees$RTD.g.cm3) # r = -0.004
cor.test(trees$diff.drt.trt, trees$SRL_m.g) # r = 0.01
cor.test(trees$diff.drt.trt, trees$rootDiam.mm) # r = -0.01

#### data set split by lifespan ####

# need to read out data and fix the lifespan to particular sites since some species have different lifespan at different sites

# write.csv(all.data.year2, file="./Formatted.Data/all.data.year2.response.0.csv")

all.data.year2.ls = read.csv("./Formatted.Data/all.data.year2.response.0.lifespan.csv", row.names = 1)
table(all.data.year2.ls$local_lifespan)

#### ANNUAL drt.v.drt ####

annual.data = subset(all.data.year2.ls, all.data.year2.ls$local_lifespan == "ANNUAL") # 83 data points

#### test for correlation with all data ####
cor.test(annual.data$leafN.mg.g, annual.data$height.m)
cor.test(annual.data$leafN.mg.g, annual.data$rootN.mg.g) # correlated r = 0.60
cor.test(annual.data$leafN.mg.g, annual.data$SLA_m2.kg) # correlated r = 0.28
cor.test(annual.data$leafN.mg.g, annual.data$root.depth_m)
cor.test(annual.data$leafN.mg.g, annual.data$RTD.g.cm3)
cor.test(annual.data$leafN.mg.g, annual.data$SRL_m.g)
cor.test(annual.data$leafN.mg.g, annual.data$rootDiam.mm)
cor.test(annual.data$height.m, annual.data$rootN.mg.g)
cor.test(annual.data$height.m, annual.data$SLA_m2.kg)  # correlated r = 0.23
cor.test(annual.data$height.m, annual.data$root.depth_m)  # correlated r = 0.35
cor.test(annual.data$height.m, annual.data$RTD.g.cm3)
cor.test(annual.data$height.m, annual.data$SRL_m.g) # correlated r = 0.32
cor.test(annual.data$height.m, annual.data$rootDiam.mm)
cor.test(annual.data$rootN.mg.g, annual.data$SLA_m2.kg)
cor.test(annual.data$rootN.mg.g, annual.data$root.depth_m)
cor.test(annual.data$rootN.mg.g, annual.data$RTD.g.cm3)
cor.test(annual.data$rootN.mg.g, annual.data$SRL_m.g)
cor.test(annual.data$rootN.mg.g, annual.data$rootDiam.mm)
cor.test(annual.data$SLA_m2.kg, annual.data$root.depth_m)
cor.test(annual.data$SLA_m2.kg, annual.data$RTD.g.cm3) # correlated r = -0.33
cor.test(annual.data$SLA_m2.kg, annual.data$SRL_m.g) # correlated r = 0.29
cor.test(annual.data$SLA_m2.kg, annual.data$rootDiam.mm)
cor.test(annual.data$root.depth_m, annual.data$RTD.g.cm3)
cor.test(annual.data$root.depth_m, annual.data$SRL_m.g)
cor.test(annual.data$root.depth_m, annual.data$rootDiam.mm)
cor.test(annual.data$RTD.g.cm3, annual.data$SRL_m.g) # correlated r = -0.33
cor.test(annual.data$RTD.g.cm3, annual.data$rootDiam.mm)
cor.test(annual.data$SRL_m.g, annual.data$rootDiam.mm)

cor.mat = cor(annual.data[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(annual.data$diff.drt.trt, annual.data$leafN.mg.g) # r = -0.35
cor.test(annual.data$diff.drt.trt, annual.data$height.m) # r = -0.07
cor.test(annual.data$diff.drt.trt, annual.data$rootN.mg.g) # r = -0.22
cor.test(annual.data$diff.drt.trt, annual.data$SLA_m2.kg) # r = -0.18
cor.test(annual.data$diff.drt.trt, annual.data$root.depth_m) # r = -0.22
cor.test(annual.data$diff.drt.trt, annual.data$RTD.g.cm3) # r = 0.02
cor.test(annual.data$diff.drt.trt, annual.data$SRL_m.g) # r = 0.15
cor.test(annual.data$diff.drt.trt, annual.data$rootDiam.mm) # r = 0.05

#### PERENNIAL drt.v.drt ####

perennial.data = subset(all.data.year2.ls, all.data.year2.ls$local_lifespan == "PERENNIAL") # 463 data points

#### test for correlation with all data ####
cor.test(perennial.data$leafN.mg.g, perennial.data$height.m)
cor.test(perennial.data$leafN.mg.g, perennial.data$rootN.mg.g) # correlated r = 0.47
cor.test(perennial.data$leafN.mg.g, perennial.data$SLA_m2.kg) # correlated r = 0.33
cor.test(perennial.data$leafN.mg.g, perennial.data$root.depth_m)
cor.test(perennial.data$leafN.mg.g, perennial.data$RTD.g.cm3) # correlated r = -0.21
cor.test(perennial.data$leafN.mg.g, perennial.data$SRL_m.g) # correlated r = 0.27
cor.test(perennial.data$leafN.mg.g, perennial.data$rootDiam.mm)
cor.test(perennial.data$height.m, perennial.data$rootN.mg.g)
cor.test(perennial.data$height.m, perennial.data$SLA_m2.kg)
cor.test(perennial.data$height.m, perennial.data$root.depth_m)
cor.test(perennial.data$height.m, perennial.data$RTD.g.cm3)
cor.test(perennial.data$height.m, perennial.data$SRL_m.g)
cor.test(perennial.data$height.m, perennial.data$rootDiam.mm) # correlated r = 0.23
cor.test(perennial.data$rootN.mg.g, perennial.data$SLA_m2.kg) # correlated r = 0.24
cor.test(perennial.data$rootN.mg.g, perennial.data$root.depth_m)
cor.test(perennial.data$rootN.mg.g, perennial.data$RTD.g.cm3) # correlated r = -0.13
cor.test(perennial.data$rootN.mg.g, perennial.data$SRL_m.g)
cor.test(perennial.data$rootN.mg.g, perennial.data$rootDiam.mm)
cor.test(perennial.data$SLA_m2.kg, perennial.data$root.depth_m)
cor.test(perennial.data$SLA_m2.kg, perennial.data$RTD.g.cm3) # correlated r = -0.26
cor.test(perennial.data$SLA_m2.kg, perennial.data$SRL_m.g) # correlated r = 0.44
cor.test(perennial.data$SLA_m2.kg, perennial.data$rootDiam.mm)
cor.test(perennial.data$root.depth_m, perennial.data$RTD.g.cm3) # correlated r = -0.29
cor.test(perennial.data$root.depth_m, perennial.data$SRL_m.g)
cor.test(perennial.data$root.depth_m, perennial.data$rootDiam.mm) #  correlated r = 0.24
cor.test(perennial.data$RTD.g.cm3, perennial.data$SRL_m.g) # correlated r = -0.34
cor.test(perennial.data$RTD.g.cm3, perennial.data$rootDiam.mm)
cor.test(perennial.data$SRL_m.g, perennial.data$rootDiam.mm)

cor.mat = cor(perennial.data[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.data$diff.drt.trt, perennial.data$leafN.mg.g) # r = -0.04
cor.test(perennial.data$diff.drt.trt, perennial.data$height.m) # r = 0.02
cor.test(perennial.data$diff.drt.trt, perennial.data$rootN.mg.g) # r = -0.04
cor.test(perennial.data$diff.drt.trt, perennial.data$SLA_m2.kg) # r = -0.06
cor.test(perennial.data$diff.drt.trt, perennial.data$root.depth_m) # r = -0.004
cor.test(perennial.data$diff.drt.trt, perennial.data$RTD.g.cm3) # r = -0.05
cor.test(perennial.data$diff.drt.trt, perennial.data$SRL_m.g) # r = -0.03
cor.test(perennial.data$diff.drt.trt, perennial.data$rootDiam.mm) # r = -0.002

#### Perennial data without woody functional group drt.v.drt ####

perennial.tree = subset(perennial.data, !perennial.data$functional_group == "WOODY") # 388 data points

#### test for correlation with all data ####
cor.test(perennial.tree$leafN.mg.g, perennial.tree$height.m) # correlated r = -0.14
cor.test(perennial.tree$leafN.mg.g, perennial.tree$rootN.mg.g) # correlated r = 0.47
cor.test(perennial.tree$leafN.mg.g, perennial.tree$SLA_m2.kg) # correlated r = 0.34
cor.test(perennial.tree$leafN.mg.g, perennial.tree$root.depth_m)
cor.test(perennial.tree$leafN.mg.g, perennial.tree$RTD.g.cm3) # correlated r = -0.21
cor.test(perennial.tree$leafN.mg.g, perennial.tree$SRL_m.g) # correlated r = 0.25
cor.test(perennial.tree$leafN.mg.g, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$height.m, perennial.tree$rootN.mg.g)
cor.test(perennial.tree$height.m, perennial.tree$SLA_m2.kg) # correlated r = -0.20
cor.test(perennial.tree$height.m, perennial.tree$root.depth_m)
cor.test(perennial.tree$height.m, perennial.tree$RTD.g.cm3)
cor.test(perennial.tree$height.m, perennial.tree$SRL_m.g) # correlated r = -0.16
cor.test(perennial.tree$height.m, perennial.tree$rootDiam.mm) # correlated r = 0.14
cor.test(perennial.tree$rootN.mg.g, perennial.tree$SLA_m2.kg) # correlated r = 0.23
cor.test(perennial.tree$rootN.mg.g, perennial.tree$root.depth_m)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$RTD.g.cm3)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$SRL_m.g)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$root.depth_m)
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$RTD.g.cm3) # correlated r = -0.26
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$SRL_m.g) # correlated r = 0.42
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$rootDiam.mm) # correlated r = -0.22
cor.test(perennial.tree$root.depth_m, perennial.tree$RTD.g.cm3)
cor.test(perennial.tree$root.depth_m, perennial.tree$SRL_m.g)
cor.test(perennial.tree$root.depth_m, perennial.tree$rootDiam.mm) #  correlated r = 0.22
cor.test(perennial.tree$RTD.g.cm3, perennial.tree$SRL_m.g) # correlated r = -0.28
cor.test(perennial.tree$RTD.g.cm3, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$SRL_m.g, perennial.tree$rootDiam.mm)

cor.mat = cor(perennial.tree[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.tree$diff.drt.trt, perennial.tree$leafN.mg.g) # r = -0.02
cor.test(perennial.tree$diff.drt.trt, perennial.tree$height.m) # r = 0.04
cor.test(perennial.tree$diff.drt.trt, perennial.tree$rootN.mg.g) # r = 0.01
cor.test(perennial.tree$diff.drt.trt, perennial.tree$SLA_m2.kg) # r = -0.06
cor.test(perennial.tree$diff.drt.trt, perennial.tree$root.depth_m) # r = 0.006
cor.test(perennial.tree$diff.drt.trt, perennial.tree$RTD.g.cm3) # r = -0.08
cor.test(perennial.tree$diff.drt.trt, perennial.tree$SRL_m.g) # r = -0.02
cor.test(perennial.tree$diff.drt.trt, perennial.tree$rootDiam.mm) # r = 0.003

#### Perennial data without TREE lifeform drt.v.drt ####

perennial.no.tree = subset(perennial.data, !perennial.data$local_lifeform == "TREE") # 459 data points

#### test for correlation with all data ####
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$height.m)
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$rootN.mg.g) # correlated r = 0.47
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$SLA_m2.kg) # correlated r = 0.33
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$RTD.g.cm3) # correlated r = -0.21
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$SRL_m.g) # correlated r = 0.27
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$height.m, perennial.no.tree$rootN.mg.g)
cor.test(perennial.no.tree$height.m, perennial.no.tree$SLA_m2.kg)
cor.test(perennial.no.tree$height.m, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$height.m, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$height.m, perennial.no.tree$SRL_m.g) # correlated r = -0.16
cor.test(perennial.no.tree$height.m, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$SLA_m2.kg) # correlated r = 0.24
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$SRL_m.g)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$RTD.g.cm3) # correlated r = -0.29
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$SRL_m.g) # correlated r = 0.44
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$rootDiam.mm) # correlated r = -0.22
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$SRL_m.g)
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$rootDiam.mm) #  correlated r = 0.24
cor.test(perennial.no.tree$RTD.g.cm3, perennial.no.tree$SRL_m.g) # correlated r = -0.34
cor.test(perennial.no.tree$RTD.g.cm3, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$SRL_m.g, perennial.no.tree$rootDiam.mm)

cor.mat = cor(perennial.no.tree[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$leafN.mg.g) # r = -0.04
cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$height.m) # r = 0.03
cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$rootN.mg.g) # r = -0.04
cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$SLA_m2.kg) # r = -0.06
cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$root.depth_m) # r = -0.004
cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$RTD.g.cm3) # r = -0.05
cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$SRL_m.g) # r = -0.03
cor.test(perennial.no.tree$diff.drt.trt, perennial.no.tree$rootDiam.mm) # r = -0.003

#### GRASS drt.v.drt ####

grass = subset(all.data.year2, all.data.year2$functional_group == "GRASS") # 175 data points

#### test for correlation with all data ####
cor.test(grass$leafN.mg.g, grass$height.m)
cor.test(grass$leafN.mg.g, grass$rootN.mg.g) # correlated r = 0.27
cor.test(grass$leafN.mg.g, grass$SLA_m2.kg) # correlated r = 0.42
cor.test(grass$leafN.mg.g, grass$root.depth_m)
cor.test(grass$leafN.mg.g, grass$RTD.g.cm3) # correlated r = -0.28
cor.test(grass$leafN.mg.g, grass$SRL_m.g) # correlated r = 0.55
cor.test(grass$leafN.mg.g, grass$rootDiam.mm) # correlated r = -0.33
cor.test(grass$height.m, grass$rootN.mg.g) # correlated r = 0.22
cor.test(grass$height.m, grass$SLA_m2.kg)
cor.test(grass$height.m, grass$root.depth_m) # correlated r = 0.19
cor.test(grass$height.m, grass$RTD.g.cm3)
cor.test(grass$height.m, grass$SRL_m.g)
cor.test(grass$height.m, grass$rootDiam.mm) # correlated r = 0.30
cor.test(grass$rootN.mg.g, grass$SLA_m2.kg)
cor.test(grass$rootN.mg.g, grass$root.depth_m)
cor.test(grass$rootN.mg.g, grass$RTD.g.cm3)
cor.test(grass$rootN.mg.g, grass$SRL_m.g)
cor.test(grass$rootN.mg.g, grass$rootDiam.mm)
cor.test(grass$SLA_m2.kg, grass$root.depth_m) # correlated r = -0.19
cor.test(grass$SLA_m2.kg, grass$RTD.g.cm3) # correlated r = -0.52
cor.test(grass$SLA_m2.kg, grass$SRL_m.g) # correlated r = 0.45
cor.test(grass$SLA_m2.kg, grass$rootDiam.mm) # correlated r = -0.35
cor.test(grass$root.depth_m, grass$RTD.g.cm3)
cor.test(grass$root.depth_m, grass$SRL_m.g)
cor.test(grass$root.depth_m, grass$rootDiam.mm) #  correlated r = 0.24
cor.test(grass$RTD.g.cm3, grass$SRL_m.g) # correlated r = -0.36
cor.test(grass$RTD.g.cm3, grass$rootDiam.mm)
cor.test(grass$SRL_m.g, grass$rootDiam.mm) # correlated r = -0.38

cor.mat = cor(grass[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(grass$diff.drt.trt, grass$leafN.mg.g) # r = -0.11
cor.test(grass$diff.drt.trt, grass$height.m) # r = 0.10
cor.test(grass$diff.drt.trt, grass$rootN.mg.g) # r = -0.05
cor.test(grass$diff.drt.trt, grass$SLA_m2.kg) # r = -0.09
cor.test(grass$diff.drt.trt, grass$root.depth_m) # r = 0.10
cor.test(grass$diff.drt.trt, grass$RTD.g.cm3) # r = -0.10
cor.test(grass$diff.drt.trt, grass$SRL_m.g) # r = 0.01
cor.test(grass$diff.drt.trt, grass$rootDiam.mm) # r = -0.03
#### FORB drt.v.drt ####

forb = subset(all.data.year2, all.data.year2$functional_group == "FORB") # 252 data points

#### test for correlation with all data ####
cor.test(forb$leafN.mg.g, forb$height.m)
cor.test(forb$leafN.mg.g, forb$rootN.mg.g)
cor.test(forb$leafN.mg.g, forb$SLA_m2.kg) # correlated r = 0.24
cor.test(forb$leafN.mg.g, forb$root.depth_m) # correlated r = 0.18
cor.test(forb$leafN.mg.g, forb$RTD.g.cm3)
cor.test(forb$leafN.mg.g, forb$SRL_m.g) # correlated r = 0.22
cor.test(forb$leafN.mg.g, forb$rootDiam.mm)
cor.test(forb$height.m, forb$rootN.mg.g)
cor.test(forb$height.m, forb$SLA_m2.kg) # correlated r = -0.22
cor.test(forb$height.m, forb$root.depth_m)
cor.test(forb$height.m, forb$RTD.g.cm3) # correlated r = -0.33
cor.test(forb$height.m, forb$SRL_m.g) # correlated r = -0.26
cor.test(forb$height.m, forb$rootDiam.mm)
cor.test(forb$rootN.mg.g, forb$SLA_m2.kg)
cor.test(forb$rootN.mg.g, forb$root.depth_m)
cor.test(forb$rootN.mg.g, forb$RTD.g.cm3)
cor.test(forb$rootN.mg.g, forb$SRL_m.g) # correlated r = 0.31
cor.test(forb$rootN.mg.g, forb$rootDiam.mm)
cor.test(forb$SLA_m2.kg, forb$root.depth_m)
cor.test(forb$SLA_m2.kg, forb$RTD.g.cm3)
cor.test(forb$SLA_m2.kg, forb$SRL_m.g) # correlated r = 0.38
cor.test(forb$SLA_m2.kg, forb$rootDiam.mm)
cor.test(forb$root.depth_m, forb$RTD.g.cm3)
cor.test(forb$root.depth_m, forb$SRL_m.g)
cor.test(forb$root.depth_m, forb$rootDiam.mm) #  correlated r = 0.20
cor.test(forb$RTD.g.cm3, forb$SRL_m.g)
cor.test(forb$RTD.g.cm3, forb$rootDiam.mm)
cor.test(forb$SRL_m.g, forb$rootDiam.mm)

cor.mat = cor(forb[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(forb$diff.drt.trt, forb$leafN.mg.g) # r = -0.11
cor.test(forb$diff.drt.trt, forb$height.m) # r = -0.08
cor.test(forb$diff.drt.trt, forb$rootN.mg.g) # r = -0.11
cor.test(forb$diff.drt.trt, forb$SLA_m2.kg) # r = -0.06
cor.test(forb$diff.drt.trt, forb$root.depth_m) # r = -0.09
cor.test(forb$diff.drt.trt, forb$RTD.g.cm3) # r = 0.11
cor.test(forb$diff.drt.trt, forb$SRL_m.g) # r = 0.05
cor.test(forb$diff.drt.trt, forb$rootDiam.mm) # r = -0.01

#### annual GRASS drt.v.drt ####

annual.grass = subset(annual.data, annual.data$functional_group == "GRASS") # 31 data points

#### annual forb drt.v.drt ####

annual.forb = subset(annual.data, annual.data$functional_group == "FORB") # 46 data points

#### perennial GRASS drt.v.drt ####

perennial.grass = subset(perennial.data, perennial.data$functional_group == "GRASS") # 143 data points

#### perennial forb drt.v.drt ####

perennial.forb = subset(perennial.data, perennial.data$functional_group == "FORB") # 195 data points


#### write out the files ####

write.csv(all.data.year2, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/all.data.year2.csv")
write.csv(no.trees, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/no.trees.csv")
write.csv(trees, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/trees.csv")
write.csv(annual.data, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/annual.data.csv")
write.csv(annual.forb, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/annual.forb.csv")
write.csv(annual.grass, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/annual.grass.csv")
write.csv(forb, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/forb.csv")
write.csv(grass, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/grass.csv")
write.csv(perennial.data, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.data.csv")
write.csv(perennial.forb, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.forb.csv")
write.csv(perennial.grass, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.grass.csv")
write.csv(perennial.tree, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.tree.csv")
write.csv(perennial.no.tree, file="./Formatted.Data/Drt.yr1.v.drt.yr2.data/perennial.no.tree.csv")









#### Generating data files of Control yr 2 versus Drought yr 2 ####

# read in cover year 2
cover.yr.2 = read.csv("./Formatted.Data/cover.trt.y2.csv")

# remove individuals where 0 in control or drought

cover.2 = subset(cover.yr.2, cover.yr.2$mean.ctrl.cover > 0)
cover.3 = subset(cover.2, cover.2$mean.drt.cover > 0)

# merge new cover data with traits

trait.data.new = read.csv("./Formatted.Data/trait.species.trt.yr1.outlier.2.csv")

# subset traits so they must have SLA
trait.data.2 = trait.data.new[,c(1,7,8,10,12,14,15,18,20,27:29)]
trait.data.3 = subset(trait.data.2, trait.data.2$SLA_m2.kg > 0 ) # 645 data points, 

all.data.year2 = merge(cover.3, trait.data.3, by="Taxon") # 579 data points


#### test for correlation with all data ctr.v.drt.2 ####
cor.test(all.data.year2$leafN.mg.g, all.data.year2$height.m)
cor.test(all.data.year2$leafN.mg.g, all.data.year2$rootN.mg.g) # correlated r = 0.45
cor.test(all.data.year2$leafN.mg.g, all.data.year2$SLA_m2.kg) # correlated r = 0.36
cor.test(all.data.year2$leafN.mg.g, all.data.year2$root.depth_m)
cor.test(all.data.year2$leafN.mg.g, all.data.year2$RTD.g.cm3) # correlated r = -0.16
cor.test(all.data.year2$leafN.mg.g, all.data.year2$SRL_m.g) # correlated r = 0.22
cor.test(all.data.year2$leafN.mg.g, all.data.year2$rootDiam.mm)
cor.test(all.data.year2$height.m, all.data.year2$rootN.mg.g)
cor.test(all.data.year2$height.m, all.data.year2$SLA_m2.kg) # correlated r = -0.09
cor.test(all.data.year2$height.m, all.data.year2$root.depth_m)
cor.test(all.data.year2$height.m, all.data.year2$RTD.g.cm3)
cor.test(all.data.year2$height.m, all.data.year2$SRL_m.g) # correlated r = -0.15
cor.test(all.data.year2$height.m, all.data.year2$rootDiam.mm)
cor.test(all.data.year2$rootN.mg.g, all.data.year2$SLA_m2.kg) # correlated r = 0.18
cor.test(all.data.year2$rootN.mg.g, all.data.year2$root.depth_m)
cor.test(all.data.year2$rootN.mg.g, all.data.year2$RTD.g.cm3)
cor.test(all.data.year2$rootN.mg.g, all.data.year2$SRL_m.g)
cor.test(all.data.year2$rootN.mg.g, all.data.year2$rootDiam.mm)
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$root.depth_m)
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$RTD.g.cm3) # correlated r = -0.24
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$SRL_m.g) # correlated r = 0.42
cor.test(all.data.year2$SLA_m2.kg, all.data.year2$rootDiam.mm) # correlated r = -0.15
cor.test(all.data.year2$root.depth_m, all.data.year2$RTD.g.cm3)
cor.test(all.data.year2$root.depth_m, all.data.year2$SRL_m.g)
cor.test(all.data.year2$root.depth_m, all.data.year2$rootDiam.mm) #  correlated r = 0.18
cor.test(all.data.year2$RTD.g.cm3, all.data.year2$SRL_m.g) # correlated r = -0.29
cor.test(all.data.year2$RTD.g.cm3, all.data.year2$rootDiam.mm)
cor.test(all.data.year2$SRL_m.g, all.data.year2$rootDiam.mm)

cor.mat = cor(all.data.year2[,c(12:19)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(all.data.year2$mean.cover.response, all.data.year2$leafN.mg.g) # r = -0.007
cor.test(all.data.year2$mean.cover.response, all.data.year2$height.m) # r = 0.001
cor.test(all.data.year2$mean.cover.response, all.data.year2$rootN.mg.g) # r = -0.03
cor.test(all.data.year2$mean.cover.response, all.data.year2$SLA_m2.kg) # r = -0.00007
cor.test(all.data.year2$mean.cover.response, all.data.year2$root.depth_m) # r = 0.02
cor.test(all.data.year2$mean.cover.response, all.data.year2$RTD.g.cm3) # r = -0.001
cor.test(all.data.year2$mean.cover.response, all.data.year2$SRL_m.g) # r = -0.02
cor.test(all.data.year2$mean.cover.response, all.data.year2$rootDiam.mm) # r = 0.04


# all outliers removed manually from trait.species.trt.yr1.final.new and made new file trait.species.trt.yr1.outlier
hist(all.data.year2$leafN.mg.g)
boxplot(all.data.year2$leafN.mg.g)
mean = mean(all.data.year2$leafN.mg.g, na.rm = TRUE)
std = sd(all.data.year2$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
all.data.year2$leafN.mg.g[which(all.data.year2$leafN.mg.g <Tmin | all.data.year2$leafN.mg.g > Tmax)]
# removed leafN 58.30000 and 67.58333


#### data set with WOODY removed ####

no.trees = subset(all.data.year2, !all.data.year2$functional_group == "WOODY") # 507 data points

#### test for correlation with all data without woody ctr.v.drt.2  ####
cor.test(no.trees$leafN.mg.g, no.trees$height.m) # correlated r = -0.11
cor.test(no.trees$leafN.mg.g, no.trees$rootN.mg.g) # correlated r = 0.45
cor.test(no.trees$leafN.mg.g, no.trees$SLA_m2.kg) # correlated r = 0.36
cor.test(no.trees$leafN.mg.g, no.trees$root.depth_m) # correlated r = 0.11
cor.test(no.trees$leafN.mg.g, no.trees$RTD.g.cm3) # correlated r = -0.16
cor.test(no.trees$leafN.mg.g, no.trees$SRL_m.g) # correlated r = 0.19
cor.test(no.trees$leafN.mg.g, no.trees$rootDiam.mm)
cor.test(no.trees$height.m, no.trees$rootN.mg.g)
cor.test(no.trees$height.m, no.trees$SLA_m2.kg) # correlated r = -0.16
cor.test(no.trees$height.m, no.trees$root.depth_m) # correlated r = 0.11
cor.test(no.trees$height.m, no.trees$RTD.g.cm3) # correlated r = -0.15
cor.test(no.trees$height.m, no.trees$SRL_m.g) # correlated r = -0.13
cor.test(no.trees$height.m, no.trees$rootDiam.mm)
cor.test(no.trees$rootN.mg.g, no.trees$SLA_m2.kg) # correlated r = 0.16
cor.test(no.trees$rootN.mg.g, no.trees$root.depth_m)
cor.test(no.trees$rootN.mg.g, no.trees$RTD.g.cm3)
cor.test(no.trees$rootN.mg.g, no.trees$SRL_m.g)
cor.test(no.trees$rootN.mg.g, no.trees$rootDiam.mm)
cor.test(no.trees$SLA_m2.kg, no.trees$root.depth_m)
cor.test(no.trees$SLA_m2.kg, no.trees$RTD.g.cm3) # correlated r = -0.22
cor.test(no.trees$SLA_m2.kg, no.trees$SRL_m.g) # correlated r = 0.38
cor.test(no.trees$SLA_m2.kg, no.trees$rootDiam.mm) # correlated r = -0.12
cor.test(no.trees$root.depth_m, no.trees$RTD.g.cm3) # correlated r = 0.12
cor.test(no.trees$root.depth_m, no.trees$SRL_m.g)
cor.test(no.trees$root.depth_m, no.trees$rootDiam.mm) #  correlated r = 0.16
cor.test(no.trees$RTD.g.cm3, no.trees$SRL_m.g) # correlated r = -0.23
cor.test(no.trees$RTD.g.cm3, no.trees$rootDiam.mm)
cor.test(no.trees$SRL_m.g, no.trees$rootDiam.mm)

cor.mat = cor(no.trees[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(no.trees$mean.cover.response, no.trees$leafN.mg.g) # r = 0.002
cor.test(no.trees$mean.cover.response, no.trees$height.m) # r = -0.02
cor.test(no.trees$mean.cover.response, no.trees$rootN.mg.g) # r = -0.03
cor.test(no.trees$mean.cover.response, no.trees$SLA_m2.kg) # r = 0.009
cor.test(no.trees$mean.cover.response, no.trees$root.depth_m) # r = 0.02
cor.test(no.trees$mean.cover.response, no.trees$RTD.g.cm3) # r = -0.02
cor.test(no.trees$mean.cover.response, no.trees$SRL_m.g) # r = -0.002
cor.test(no.trees$mean.cover.response, no.trees$rootDiam.mm) # r = 0.04

#### data set with TREE removed ####

trees = subset(all.data.year2, !all.data.year2$local_lifeform == "TREE") # 576 data points

#### test for correlation with all data without TREE lifeform ctr.v.drt.2  ####

cor.test(trees$leafN.mg.g, trees$height.m)
cor.test(trees$leafN.mg.g, trees$rootN.mg.g) # correlated r = 0.45
cor.test(trees$leafN.mg.g, trees$SLA_m2.kg) # correlated r = 0.35
cor.test(trees$leafN.mg.g, trees$root.depth_m)
cor.test(trees$leafN.mg.g, trees$RTD.g.cm3) # correlated r = -0.16
cor.test(trees$leafN.mg.g, trees$SRL_m.g) # correlated r = 0.21
cor.test(trees$leafN.mg.g, trees$rootDiam.mm)
cor.test(trees$height.m, trees$rootN.mg.g)
cor.test(trees$height.m, trees$SLA_m2.kg) # correlated r = -0.09
cor.test(trees$height.m, trees$root.depth_m)
cor.test(trees$height.m, trees$RTD.g.cm3)
cor.test(trees$height.m, trees$SRL_m.g) # correlated r = -0.15
cor.test(trees$height.m, trees$rootDiam.mm)
cor.test(trees$rootN.mg.g, trees$SLA_m2.kg) # correlated r = 0.18
cor.test(trees$rootN.mg.g, trees$root.depth_m)
cor.test(trees$rootN.mg.g, trees$RTD.g.cm3)
cor.test(trees$rootN.mg.g, trees$SRL_m.g)
cor.test(trees$rootN.mg.g, trees$rootDiam.mm)
cor.test(trees$SLA_m2.kg, trees$root.depth_m)
cor.test(trees$SLA_m2.kg, trees$RTD.g.cm3) # correlated r = -0.24
cor.test(trees$SLA_m2.kg, trees$SRL_m.g) # correlated r = 0.41
cor.test(trees$SLA_m2.kg, trees$rootDiam.mm) # correlated r = -0.14
cor.test(trees$root.depth_m, trees$RTD.g.cm3)
cor.test(trees$root.depth_m, trees$SRL_m.g)
cor.test(trees$root.depth_m, trees$rootDiam.mm) #  correlated r = 0.18
cor.test(trees$RTD.g.cm3, trees$SRL_m.g) # correlated r = -0.28
cor.test(trees$RTD.g.cm3, trees$rootDiam.mm)
cor.test(trees$SRL_m.g, trees$rootDiam.mm)

cor.mat = cor(trees[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(trees$mean.cover.response, trees$leafN.mg.g) # r = -0.007
cor.test(trees$mean.cover.response, trees$height.m) # r = 0.001
cor.test(trees$mean.cover.response, trees$rootN.mg.g) # r = -0.03
cor.test(trees$mean.cover.response, trees$SLA_m2.kg) # r = -0.0004
cor.test(trees$mean.cover.response, trees$root.depth_m) # r = 0.02
cor.test(trees$mean.cover.response, trees$RTD.g.cm3) # r = -0.002
cor.test(trees$mean.cover.response, trees$SRL_m.g) # r = -0.02
cor.test(trees$mean.cover.response, trees$rootDiam.mm) # r = 0.04

#### data set split by lifespan ####

# need to read out data and fix the lifespan to particular sites since some species have different lifespan at different sites

# write.csv(all.data.year2, file="./Formatted.Data/all.data.year2.ctr.v.drt.response.0.csv")

all.data.year2.ls = read.csv("./Formatted.Data/all.data.year2.ctr.v.drt.response.0.lifespan.csv", row.names = 1)
table(all.data.year2.ls$local_lifespan)

#### ANNUAL ctr.v.drt.2  ####

annual.data = subset(all.data.year2.ls, all.data.year2.ls$local_lifespan == "ANNUAL") # 96 data points

#### test for correlation with all data ####
cor.test(annual.data$leafN.mg.g, annual.data$height.m)
cor.test(annual.data$leafN.mg.g, annual.data$rootN.mg.g) # correlated r = 0.48
cor.test(annual.data$leafN.mg.g, annual.data$SLA_m2.kg)  # correlated r = 0.28
cor.test(annual.data$leafN.mg.g, annual.data$root.depth_m)
cor.test(annual.data$leafN.mg.g, annual.data$RTD.g.cm3)
cor.test(annual.data$leafN.mg.g, annual.data$SRL_m.g)
cor.test(annual.data$leafN.mg.g, annual.data$rootDiam.mm)  # correlated r = -0.31
cor.test(annual.data$height.m, annual.data$rootN.mg.g)
cor.test(annual.data$height.m, annual.data$SLA_m2.kg)  # correlated r = 0.31
cor.test(annual.data$height.m, annual.data$root.depth_m)  # correlated r = 0.29
cor.test(annual.data$height.m, annual.data$RTD.g.cm3)
cor.test(annual.data$height.m, annual.data$SRL_m.g)
cor.test(annual.data$height.m, annual.data$rootDiam.mm)
cor.test(annual.data$rootN.mg.g, annual.data$SLA_m2.kg)
cor.test(annual.data$rootN.mg.g, annual.data$root.depth_m)
cor.test(annual.data$rootN.mg.g, annual.data$RTD.g.cm3)
cor.test(annual.data$rootN.mg.g, annual.data$SRL_m.g)
cor.test(annual.data$rootN.mg.g, annual.data$rootDiam.mm)
cor.test(annual.data$SLA_m2.kg, annual.data$root.depth_m)
cor.test(annual.data$SLA_m2.kg, annual.data$RTD.g.cm3)
cor.test(annual.data$SLA_m2.kg, annual.data$SRL_m.g)
cor.test(annual.data$SLA_m2.kg, annual.data$rootDiam.mm)
cor.test(annual.data$root.depth_m, annual.data$RTD.g.cm3)
cor.test(annual.data$root.depth_m, annual.data$SRL_m.g)  # correlated r = -0.27
cor.test(annual.data$root.depth_m, annual.data$rootDiam.mm)
cor.test(annual.data$RTD.g.cm3, annual.data$SRL_m.g) 
cor.test(annual.data$RTD.g.cm3, annual.data$rootDiam.mm)
cor.test(annual.data$SRL_m.g, annual.data$rootDiam.mm)

cor.mat = cor(annual.data[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(annual.data$mean.cover.response, annual.data$leafN.mg.g) # r = -0.04
cor.test(annual.data$mean.cover.response, annual.data$height.m) # r = -0.17
cor.test(annual.data$mean.cover.response, annual.data$rootN.mg.g) # r = 0.04
cor.test(annual.data$mean.cover.response, annual.data$SLA_m2.kg) # r = -0.04
cor.test(annual.data$mean.cover.response, annual.data$root.depth_m) # r = 0.11
cor.test(annual.data$mean.cover.response, annual.data$RTD.g.cm3) # r = 0.11
cor.test(annual.data$mean.cover.response, annual.data$SRL_m.g) # r = -0.11
cor.test(annual.data$mean.cover.response, annual.data$rootDiam.mm) # r = 0.14

#### PERENNIAL ctr.v.drt.2  ####

perennial.data = subset(all.data.year2.ls, all.data.year2.ls$local_lifespan == "PERENNIAL") # 465 data points

#### test for correlation with all data ####
cor.test(perennial.data$leafN.mg.g, perennial.data$height.m)
cor.test(perennial.data$leafN.mg.g, perennial.data$rootN.mg.g) # correlated r = 0.43
cor.test(perennial.data$leafN.mg.g, perennial.data$SLA_m2.kg) # correlated r = 0.35
cor.test(perennial.data$leafN.mg.g, perennial.data$root.depth_m)
cor.test(perennial.data$leafN.mg.g, perennial.data$RTD.g.cm3) # correlated r = -0.24
cor.test(perennial.data$leafN.mg.g, perennial.data$SRL_m.g) # correlated r = 0.27
cor.test(perennial.data$leafN.mg.g, perennial.data$rootDiam.mm)
cor.test(perennial.data$height.m, perennial.data$rootN.mg.g)
cor.test(perennial.data$height.m, perennial.data$SLA_m2.kg) # correlated r = -0.10
cor.test(perennial.data$height.m, perennial.data$root.depth_m)
cor.test(perennial.data$height.m, perennial.data$RTD.g.cm3)
cor.test(perennial.data$height.m, perennial.data$SRL_m.g) # correlated r = -0.18
cor.test(perennial.data$height.m, perennial.data$rootDiam.mm)
cor.test(perennial.data$rootN.mg.g, perennial.data$SLA_m2.kg) # correlated r = 0.19
cor.test(perennial.data$rootN.mg.g, perennial.data$root.depth_m)
cor.test(perennial.data$rootN.mg.g, perennial.data$RTD.g.cm3)
cor.test(perennial.data$rootN.mg.g, perennial.data$SRL_m.g)
cor.test(perennial.data$rootN.mg.g, perennial.data$rootDiam.mm)
cor.test(perennial.data$SLA_m2.kg, perennial.data$root.depth_m)
cor.test(perennial.data$SLA_m2.kg, perennial.data$RTD.g.cm3) # correlated r = -0.31
cor.test(perennial.data$SLA_m2.kg, perennial.data$SRL_m.g) # correlated r = 0.46
cor.test(perennial.data$SLA_m2.kg, perennial.data$rootDiam.mm) # correlated r = -0.18
cor.test(perennial.data$root.depth_m, perennial.data$RTD.g.cm3)
cor.test(perennial.data$root.depth_m, perennial.data$SRL_m.g)
cor.test(perennial.data$root.depth_m, perennial.data$rootDiam.mm) #  correlated r = 0.24
cor.test(perennial.data$RTD.g.cm3, perennial.data$SRL_m.g) # correlated r = -0.34
cor.test(perennial.data$RTD.g.cm3, perennial.data$rootDiam.mm)
cor.test(perennial.data$SRL_m.g, perennial.data$rootDiam.mm)

cor.mat = cor(perennial.data[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.data$mean.cover.response, perennial.data$leafN.mg.g) # r = 0.01
cor.test(perennial.data$mean.cover.response, perennial.data$height.m) # r = 0.02
cor.test(perennial.data$mean.cover.response, perennial.data$rootN.mg.g) # r = -0.03
cor.test(perennial.data$mean.cover.response, perennial.data$SLA_m2.kg) # r = 0.01
cor.test(perennial.data$mean.cover.response, perennial.data$root.depth_m) # r = -0.009
cor.test(perennial.data$mean.cover.response, perennial.data$RTD.g.cm3) # r = -0.05
cor.test(perennial.data$mean.cover.response, perennial.data$SRL_m.g) # r = 0.002
cor.test(perennial.data$mean.cover.response, perennial.data$rootDiam.mm) # r = 0.04

#### Perennial data without woody functional group ctr.v.drt.2  ####

perennial.tree = subset(perennial.data, !perennial.data$functional_group == "WOODY") # 395 data points

#### test for correlation with all data ####
cor.test(perennial.tree$leafN.mg.g, perennial.tree$height.m) # correlated r = -0.14
cor.test(perennial.tree$leafN.mg.g, perennial.tree$rootN.mg.g) # correlated r = 0.43
cor.test(perennial.tree$leafN.mg.g, perennial.tree$SLA_m2.kg) # correlated r = 0.37
cor.test(perennial.tree$leafN.mg.g, perennial.tree$root.depth_m)
cor.test(perennial.tree$leafN.mg.g, perennial.tree$RTD.g.cm3) # correlated r = -0.26
cor.test(perennial.tree$leafN.mg.g, perennial.tree$SRL_m.g) # correlated r = 0.26
cor.test(perennial.tree$leafN.mg.g, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$height.m, perennial.tree$rootN.mg.g)
cor.test(perennial.tree$height.m, perennial.tree$SLA_m2.kg) # correlated r = -0.21
cor.test(perennial.tree$height.m, perennial.tree$root.depth_m)
cor.test(perennial.tree$height.m, perennial.tree$RTD.g.cm3)
cor.test(perennial.tree$height.m, perennial.tree$SRL_m.g) # correlated r = -0.17
cor.test(perennial.tree$height.m, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$SLA_m2.kg) # correlated r = 0.17
cor.test(perennial.tree$rootN.mg.g, perennial.tree$root.depth_m)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$RTD.g.cm3)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$SRL_m.g)
cor.test(perennial.tree$rootN.mg.g, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$root.depth_m)
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$RTD.g.cm3) # correlated r = -0.29
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$SRL_m.g) # correlated r = 0.44
cor.test(perennial.tree$SLA_m2.kg, perennial.tree$rootDiam.mm) # correlated r = -0.15
cor.test(perennial.tree$root.depth_m, perennial.tree$RTD.g.cm3) # correlated r = 0.14
cor.test(perennial.tree$root.depth_m, perennial.tree$SRL_m.g)
cor.test(perennial.tree$root.depth_m, perennial.tree$rootDiam.mm) #  correlated r = 0.21
cor.test(perennial.tree$RTD.g.cm3, perennial.tree$SRL_m.g) # correlated r = -0.28
cor.test(perennial.tree$RTD.g.cm3, perennial.tree$rootDiam.mm)
cor.test(perennial.tree$SRL_m.g, perennial.tree$rootDiam.mm)

cor.mat = cor(perennial.tree[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.tree$mean.cover.response, perennial.tree$leafN.mg.g) # r = 0.03
cor.test(perennial.tree$mean.cover.response, perennial.tree$height.m) # r = -0.001
cor.test(perennial.tree$mean.cover.response, perennial.tree$rootN.mg.g) # r = -0.03
cor.test(perennial.tree$mean.cover.response, perennial.tree$SLA_m2.kg) # r = 0.02
cor.test(perennial.tree$mean.cover.response, perennial.tree$root.depth_m) # r = -0.02
cor.test(perennial.tree$mean.cover.response, perennial.tree$RTD.g.cm3) # r = -0.09
cor.test(perennial.tree$mean.cover.response, perennial.tree$SRL_m.g) # r = 0.03
cor.test(perennial.tree$mean.cover.response, perennial.tree$rootDiam.mm) # r = 0.03

#### Perennial data without TREE lifeform ctr.v.drt.2  ####

perennial.no.tree = subset(perennial.data, !perennial.data$local_lifeform == "TREE") # 462 data points

#### test for correlation with all data ####
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$height.m)
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$rootN.mg.g) # correlated r = 0.43
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$SLA_m2.kg) # correlated r = 0.34
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$RTD.g.cm3) # correlated r = -0.24
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$SRL_m.g) # correlated r = 0.27
cor.test(perennial.no.tree$leafN.mg.g, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$height.m, perennial.no.tree$rootN.mg.g)
cor.test(perennial.no.tree$height.m, perennial.no.tree$SLA_m2.kg) # correlated r = -0.10
cor.test(perennial.no.tree$height.m, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$height.m, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$height.m, perennial.no.tree$SRL_m.g) # correlated r = -0.18
cor.test(perennial.no.tree$height.m, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$SLA_m2.kg) # correlated r = 0.19
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$SRL_m.g)
cor.test(perennial.no.tree$rootN.mg.g, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$root.depth_m)
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$RTD.g.cm3) # correlated r = -0.31
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$SRL_m.g) # correlated r = 0.45
cor.test(perennial.no.tree$SLA_m2.kg, perennial.no.tree$rootDiam.mm) # correlated r = -0.16
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$RTD.g.cm3)
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$SRL_m.g)
cor.test(perennial.no.tree$root.depth_m, perennial.no.tree$rootDiam.mm) #  correlated r = 0.23
cor.test(perennial.no.tree$RTD.g.cm3, perennial.no.tree$SRL_m.g) # correlated r = -0.33
cor.test(perennial.no.tree$RTD.g.cm3, perennial.no.tree$rootDiam.mm)
cor.test(perennial.no.tree$SRL_m.g, perennial.no.tree$rootDiam.mm)

cor.mat = cor(perennial.no.tree[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$leafN.mg.g) # r = 0.01
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$height.m) # r = 0.02
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$rootN.mg.g) # r = -0.03
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$SLA_m2.kg) # r = 0.01
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$root.depth_m) # r = -0.008
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$RTD.g.cm3) # r = -0.05
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$SRL_m.g) # r = 0.002
cor.test(perennial.no.tree$mean.cover.response, perennial.no.tree$rootDiam.mm) # r = 0.04

#### GRASS ctr.v.drt.2  ####

grass = subset(all.data.year2, all.data.year2$functional_group == "GRASS") # 193 data points

#### test for correlation with all data ####
cor.test(grass$leafN.mg.g, grass$height.m)
cor.test(grass$leafN.mg.g, grass$rootN.mg.g) # correlated r = 0.19
cor.test(grass$leafN.mg.g, grass$SLA_m2.kg) # correlated r = 0.40
cor.test(grass$leafN.mg.g, grass$root.depth_m)
cor.test(grass$leafN.mg.g, grass$RTD.g.cm3) # correlated r = -0.32
cor.test(grass$leafN.mg.g, grass$SRL_m.g) # correlated r = 0.58
cor.test(grass$leafN.mg.g, grass$rootDiam.mm) # correlated r = -0.35
cor.test(grass$height.m, grass$rootN.mg.g) # correlated r = 0.18
cor.test(grass$height.m, grass$SLA_m2.kg)
cor.test(grass$height.m, grass$root.depth_m)
cor.test(grass$height.m, grass$RTD.g.cm3)
cor.test(grass$height.m, grass$SRL_m.g)
cor.test(grass$height.m, grass$rootDiam.mm) # correlated r = 0.23
cor.test(grass$rootN.mg.g, grass$SLA_m2.kg)
cor.test(grass$rootN.mg.g, grass$root.depth_m)
cor.test(grass$rootN.mg.g, grass$RTD.g.cm3)
cor.test(grass$rootN.mg.g, grass$SRL_m.g)
cor.test(grass$rootN.mg.g, grass$rootDiam.mm)
cor.test(grass$SLA_m2.kg, grass$root.depth_m) # correlated r = -0.22
cor.test(grass$SLA_m2.kg, grass$RTD.g.cm3) # correlated r = -0.53
cor.test(grass$SLA_m2.kg, grass$SRL_m.g) # correlated r = 0.44
cor.test(grass$SLA_m2.kg, grass$rootDiam.mm) # correlated r = -0.29
cor.test(grass$root.depth_m, grass$RTD.g.cm3) # correlated r = 0.20
cor.test(grass$root.depth_m, grass$SRL_m.g) # correlated r = -0.19
cor.test(grass$root.depth_m, grass$rootDiam.mm) #  correlated r = 0.25
cor.test(grass$RTD.g.cm3, grass$SRL_m.g) # correlated r = -0.35
cor.test(grass$RTD.g.cm3, grass$rootDiam.mm)
cor.test(grass$SRL_m.g, grass$rootDiam.mm) # correlated r = -0.40

cor.mat = cor(grass[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(grass$mean.cover.response, grass$leafN.mg.g) # r = -0.03
cor.test(grass$mean.cover.response, grass$height.m) # r = 0.01
cor.test(grass$mean.cover.response, grass$rootN.mg.g) # r = -0.10
cor.test(grass$mean.cover.response, grass$SLA_m2.kg) # r = -0.04
cor.test(grass$mean.cover.response, grass$root.depth_m) # r = -0.01
cor.test(grass$mean.cover.response, grass$RTD.g.cm3) # r = -0.10
cor.test(grass$mean.cover.response, grass$SRL_m.g) # r = -0.05
cor.test(grass$mean.cover.response, grass$rootDiam.mm) # r = 0.06

#### FORB ctr.v.drt.2  ####

forb = subset(all.data.year2, all.data.year2$functional_group == "FORB") # 254 data points

#### test for correlation with all data ####
cor.test(forb$leafN.mg.g, forb$height.m)
cor.test(forb$leafN.mg.g, forb$rootN.mg.g)
cor.test(forb$leafN.mg.g, forb$SLA_m2.kg) # correlated r = 0.28
cor.test(forb$leafN.mg.g, forb$root.depth_m) # correlated r = 0.20
cor.test(forb$leafN.mg.g, forb$RTD.g.cm3)
cor.test(forb$leafN.mg.g, forb$SRL_m.g)
cor.test(forb$leafN.mg.g, forb$rootDiam.mm)
cor.test(forb$height.m, forb$rootN.mg.g)
cor.test(forb$height.m, forb$SLA_m2.kg) # correlated r = -0.20
cor.test(forb$height.m, forb$root.depth_m)
cor.test(forb$height.m, forb$RTD.g.cm3) # correlated r = -0.38
cor.test(forb$height.m, forb$SRL_m.g) # correlated r = -0.27
cor.test(forb$height.m, forb$rootDiam.mm)
cor.test(forb$rootN.mg.g, forb$SLA_m2.kg) # correlated r = 0.21
cor.test(forb$rootN.mg.g, forb$root.depth_m)
cor.test(forb$rootN.mg.g, forb$RTD.g.cm3)
cor.test(forb$rootN.mg.g, forb$SRL_m.g) # correlated r = 0.28
cor.test(forb$rootN.mg.g, forb$rootDiam.mm)
cor.test(forb$SLA_m2.kg, forb$root.depth_m)
cor.test(forb$SLA_m2.kg, forb$RTD.g.cm3)
cor.test(forb$SLA_m2.kg, forb$SRL_m.g) # correlated r = 0.32
cor.test(forb$SLA_m2.kg, forb$rootDiam.mm)
cor.test(forb$root.depth_m, forb$RTD.g.cm3)
cor.test(forb$root.depth_m, forb$SRL_m.g)
cor.test(forb$root.depth_m, forb$rootDiam.mm)
cor.test(forb$RTD.g.cm3, forb$SRL_m.g)
cor.test(forb$RTD.g.cm3, forb$rootDiam.mm)
cor.test(forb$SRL_m.g, forb$rootDiam.mm)

cor.mat = cor(forb[,c(6:13)],use = "pairwise") 
corrplot(cor.mat, method="number")

#### correlation between response and predictors ####

cor.test(forb$mean.cover.response, forb$leafN.mg.g) # r = 0.08
cor.test(forb$mean.cover.response, forb$height.m) # r = -0.02
cor.test(forb$mean.cover.response, forb$rootN.mg.g) # r = 0.04
cor.test(forb$mean.cover.response, forb$SLA_m2.kg) # r = 0.07
cor.test(forb$mean.cover.response, forb$root.depth_m) # r = 0.05
cor.test(forb$mean.cover.response, forb$RTD.g.cm3) # r = 0.11
cor.test(forb$mean.cover.response, forb$SRL_m.g) # r = 0.02
cor.test(forb$mean.cover.response, forb$rootDiam.mm) # r = 0.02

#### annual GRASS ctr.v.drt.2  ####

annual.grass = subset(annual.data, annual.data$functional_group == "GRASS") # 34 data points

#### annual forb ctr.v.drt.2  ####

annual.forb = subset(annual.data, annual.data$functional_group == "FORB") # 54 data points

#### perennial GRASS ctr.v.drt.2  ####

perennial.grass = subset(perennial.data, perennial.data$functional_group == "GRASS") # 157 data points

#### perennial forb ctr.v.drt.2  ####

perennial.forb = subset(perennial.data, perennial.data$functional_group == "FORB") # 189 data points


#### write out the files ####

write.csv(all.data.year2, file="./Formatted.Data/Ctrl.v.drt.yr2.data/all.data.year2.csv")
write.csv(no.trees, file="./Formatted.Data/Ctrl.v.drt.yr2.data/no.trees.csv")
write.csv(trees, file="./Formatted.Data/Ctrl.v.drt.yr2.data/trees.csv")
write.csv(annual.data, file="./Formatted.Data/Ctrl.v.drt.yr2.data/annual.data.csv")
write.csv(annual.forb, file="./Formatted.Data/Ctrl.v.drt.yr2.data/annual.forb.csv")
write.csv(annual.grass, file="./Formatted.Data/Ctrl.v.drt.yr2.data/annual.grass.csv")
write.csv(forb, file="./Formatted.Data/Ctrl.v.drt.yr2.data/forb.csv")
write.csv(grass, file="./Formatted.Data/Ctrl.v.drt.yr2.data/grass.csv")
write.csv(perennial.data, file="./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.data.csv")
write.csv(perennial.forb, file="./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.forb.csv")
write.csv(perennial.grass, file="./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.grass.csv")
write.csv(perennial.tree, file="./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.tree.csv")
write.csv(perennial.no.tree, file="./Formatted.Data/Ctrl.v.drt.yr2.data/perennial.no.tree.csv")















