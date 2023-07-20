# Script to prep trait data for analyses

library(corrplot)

# read in full data frame 
trait.data = read.csv("./Formatted.Data/trait.species.trt.yr1.final.new.csv")

# subset traits to only those of interest
trait.data.2 = trait.data[,c(1,7,8,10,12,14,15,18,20,27:29)]

# subset out woody species and trees from local_lifeform

trait.data.3 = subset(trait.data.2, !trait.data.2$functional_group == "WOODY")
trait.data.3 = subset(trait.data.3, !trait.data.3$local_lifeform == "TREE")

# subset for species with SLA
trait.data.4 = subset(trait.data.3, trait.data.3$SLA_m2.kg > 0 ) # 561 species

write.csv(trait.data.4, file = "./Formatted.Data/traits.no.woody.SLA.csv")

#### checking for outliers ####
# all outliers removed manually from traits.no.woody.SLA.csv and 
# made new file traits.no.woody.SLA.no.outlier.csv

hist(trait.data.4$leafN.mg.g)
boxplot(trait.data.4$leafN.mg.g)
mean = mean(trait.data.4$leafN.mg.g, na.rm = TRUE)
std = sd(trait.data.4$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$leafN.mg.g[which(trait.data.4$leafN.mg.g <Tmin | trait.data.4$leafN.mg.g > Tmax)])
# removed leafN 58.30000, 67.58333

hist(trait.data.4$height.m)
boxplot(trait.data.4$height.m)
mean = mean(trait.data.4$height.m, na.rm = TRUE)
std = sd(trait.data.4$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$height.m[which(trait.data.4$height.m <Tmin | trait.data.4$height.m > Tmax)])
# removed height 1.996670 - 4.572000

hist(trait.data.4$rootN.mg.g)
boxplot(trait.data.4$rootN.mg.g)
mean = mean(trait.data.4$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data.4$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootN.mg.g[which(trait.data.4$rootN.mg.g <Tmin | trait.data.4$rootN.mg.g > Tmax)])
# remove rootN 33.34667 - 39.26274

hist(trait.data.4$SLA_m2.kg)
boxplot(trait.data.4$SLA_m2.kg)
mean = mean(trait.data.4$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data.4$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SLA_m2.kg[which(trait.data.4$SLA_m2.kg <Tmin | trait.data.4$SLA_m2.kg > Tmax)])
# remove SLA 57.70000 - 98.20000

hist(trait.data.4$root.depth_m)
boxplot(trait.data.4$root.depth_m)
mean = mean(trait.data.4$root.depth_m, na.rm = TRUE)
std = sd(trait.data.4$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$root.depth_m[which(trait.data.4$root.depth_m <Tmin | trait.data.4$root.depth_m > Tmax)])
# remove depth 2.426850 - 3.300000

hist(trait.data.4$RTD.g.cm3)
boxplot(trait.data.4$RTD.g.cm3)
mean = mean(trait.data.4$RTD.g.cm3, na.rm = TRUE)
std = sd(trait.data.4$RTD.g.cm3, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$RTD.g.cm3[which(trait.data.4$RTD.g.cm3 <Tmin | trait.data.4$RTD.g.cm3 > Tmax)])
# remove RTD 0.7450000 - 0.9167417

hist(trait.data.4$SRL_m.g)
boxplot(trait.data.4$SRL_m.g)
mean = mean(trait.data.4$SRL_m.g, na.rm = TRUE)
std = sd(trait.data.4$SRL_m.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SRL_m.g[which(trait.data.4$SRL_m.g <Tmin | trait.data.4$SRL_m.g > Tmax)])
# remove SRL 601.5858 - 929.9185

hist(trait.data.4$rootDiam.mm)
boxplot(trait.data.4$rootDiam.mm)
mean = mean(trait.data.4$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data.4$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootDiam.mm[which(trait.data.4$rootDiam.mm <Tmin | trait.data.4$rootDiam.mm > Tmax)])
# remove diam 1.330000 - 2.015308

#### read in new trait data without outliers and cover data ####

trait.data.new = read.csv("./Formatted.Data/traits.no.woody.SLA.no.outlier.csv")
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")

# number of missing species per trait
table(is.na(trait.data.new$leafN.mg.g)) # 182 missing, 32%
table(is.na(trait.data.new$height.m)) # 72 missing, 13%
table(is.na(trait.data.new$rootN.mg.g)) # 377 missing, 67%
table(is.na(trait.data.new$root.depth_m)) # 247 missing, 44%
table(is.na(trait.data.new$RTD.g.cm3)) # 334 missing, 60%
table(is.na(trait.data.new$SRL_m.g)) # 286 missing, 51%
table(is.na(trait.data.new$rootDiam.mm)) # 279 missing, 50%

# remove individuals where 0 in control or drought

cover.2 = subset(cover.data, cover.data$mean.ctrl.cover > 0)
cover.3 = subset(cover.2, cover.2$mean.drt.cover > 0)
# 682 species with 1008 cover points

# merge trait data and cover data

all.data = merge(cover.3, trait.data.new, by="Taxon") 
# 386 species from 632 data points from 76 sites

# test for correlations among traits
cor.mat.all = cor(all.data[,c(12:19)],use = "pairwise") 
corrplot(cor.mat.all, method="number")

#### data set split by lifespan ####
# need to read out data and fix the lifespan to particular sites since 
# some species have different lifespan at different sites

# write.csv(all.data, file="./Formatted.Data/all.data.response.0.csv")

all.data.ls = read.csv("./Formatted.Data/all.data.response.0.lifespan.csv", row.names = 1)
table(all.data.ls$local_lifespan)
# 490 perennial, 127 annual

annual.data = subset(all.data.ls, all.data.ls$local_lifespan == "ANNUAL") # 127 data points
# 127 data points from 81 species

cor.mat.annual = cor(annual.data[,c(12:19)],use = "pairwise")
corrplot(cor.mat.annual, method="number")

perennial.data = subset(all.data.ls, all.data.ls$local_lifespan == "PERENNIAL") 
# 490 data points of 296 species

cor.mat.perennial = cor(perennial.data[,c(12:19)],use = "pairwise")
corrplot(cor.mat.perennial, method="number")

grass = subset(all.data, all.data$functional_group == "GRASS") 
# 230 data points of 129 species

cor.mat.grass = cor(grass[,c(12:19)],use = "pairwise")
corrplot(cor.mat.grass, method="number")

forb = subset(all.data, all.data$functional_group == "FORB")
# 324 data points of 204 species

cor.mat.forb = cor(forb[,c(12:19)],use = "pairwise")
corrplot(cor.mat.forb, method="number")

#### write out the files ####

#write.csv(all.data, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/all.data.csv")
#write.csv(annual.data, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/annual.data.csv")
#write.csv(forb, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/forb.csv")
#write.csv(grass, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/grass.csv")
#write.csv(perennial.data, file = "./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.data.csv")











