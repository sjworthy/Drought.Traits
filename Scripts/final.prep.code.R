# Script to prep trait data for analyses

library(corrplot)

# read in full data frame 
trait.data = read.csv("./New.dfs/trait.species.trt.yr1.final.new.csv")

# subset traits to only those of interest
trait.data.2 = trait.data[,c(1,7,8,10,12,14,15,18,20,27:29)]

# subset out woody species and trees from local_lifeform

trait.data.3 = subset(trait.data.2, !trait.data.2$functional_group == "WOODY")
trait.data.3 = subset(trait.data.3, !trait.data.3$local_lifeform == "TREE")

# subset for species with SLA
trait.data.4 = subset(trait.data.3, trait.data.3$SLA_m2.kg > 0 ) # 558 taxon

#write.csv(trait.data.4, file = "./New.dfs/traits.no.woody.SLA.csv")

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
# removed leafN 58.30000 67.58333
# percent removed 
table(is.na(trait.data.4$leafN.mg.g))
558-185
2/373*100 #0.54%

hist(trait.data.4$height.m)
boxplot(trait.data.4$height.m)
mean = mean(trait.data.4$height.m, na.rm = TRUE)
std = sd(trait.data.4$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$height.m[which(trait.data.4$height.m <Tmin | trait.data.4$height.m > Tmax)])
# removed height 1.996670 - 4.572000
# percent removed 
table(is.na(trait.data.4$height.m))
558-66
8/492*100 #1.63%

hist(trait.data.4$rootN.mg.g)
boxplot(trait.data.4$rootN.mg.g)
mean = mean(trait.data.4$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data.4$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootN.mg.g[which(trait.data.4$rootN.mg.g <Tmin | trait.data.4$rootN.mg.g > Tmax)])
# remove rootN 33.34667 - 39.26274
# percent removed 
table(is.na(trait.data.4$rootN.mg.g))
558-371
3/187*100 #1.60%

hist(trait.data.4$SLA_m2.kg)
boxplot(trait.data.4$SLA_m2.kg)
mean = mean(trait.data.4$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data.4$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SLA_m2.kg[which(trait.data.4$SLA_m2.kg <Tmin | trait.data.4$SLA_m2.kg > Tmax)])
# remove SLA 53.50000 - 98.20000
# percent removed 
table(is.na(trait.data.4$SLA_m2.kg))
12/558*100 #2.15%

hist(trait.data.4$root.depth_m)
boxplot(trait.data.4$root.depth_m)
mean = mean(trait.data.4$root.depth_m, na.rm = TRUE)
std = sd(trait.data.4$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$root.depth_m[which(trait.data.4$root.depth_m <Tmin | trait.data.4$root.depth_m > Tmax)])
# remove depth 2.426850 - 3.300000
# percent removed 
table(is.na(trait.data.4$root.depth_m))
558-245
9/313*100 #2.88%

hist(trait.data.4$RTD.g.cm3)
boxplot(trait.data.4$RTD.g.cm3)
mean = mean(trait.data.4$RTD.g.cm3, na.rm = TRUE)
std = sd(trait.data.4$RTD.g.cm3, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$RTD.g.cm3[which(trait.data.4$RTD.g.cm3 <Tmin | trait.data.4$RTD.g.cm3 > Tmax)])
# remove RTD 0.7450000 - 0.9167417
# percent removed 
table(is.na(trait.data.4$RTD.g.cm3))
558-331
3/227*100 #1.32%

hist(trait.data.4$SRL_m.g)
boxplot(trait.data.4$SRL_m.g)
mean = mean(trait.data.4$SRL_m.g, na.rm = TRUE)
std = sd(trait.data.4$SRL_m.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SRL_m.g[which(trait.data.4$SRL_m.g <Tmin | trait.data.4$SRL_m.g > Tmax)])
# remove SRL 601.5858 - 929.9185
# percent removed 
table(is.na(trait.data.4$SRL_m.g))
558-281
6/277*100 #2.17%

hist(trait.data.4$rootDiam.mm)
boxplot(trait.data.4$rootDiam.mm)
mean = mean(trait.data.4$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data.4$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootDiam.mm[which(trait.data.4$rootDiam.mm <Tmin | trait.data.4$rootDiam.mm > Tmax)])
# remove diam 1.330000 - 2.015308
# percent removed 
table(is.na(trait.data.4$rootDiam.mm))
558-274
8/284*100 #2.81%

#### read in new trait data without outliers and cover data ####

trait.data.new = read.csv("./New.dfs/traits.no.woody.SLA.no.outlier.csv")
cover.data = read.csv("./New.dfs/cover.response.trt.y1.csv")

# subset trait data for only species in cover data

subset.traits = subset(trait.data.new, trait.data.new$Taxon %in% cover.data$Taxon)

# merge trait data and cover data

all.data = merge(cover.data, subset.traits, by="Taxon") 
# 383 species from 628 data points from 76 sites

# write.csv(all.data, file = "./New.dfs/all.data.response.0.csv")

#### data set split by lifespan ####
# need to read out data and fix the lifespan to particular sites since 
# some species have different lifespan at different sites

all.data.ls = read.csv("./New.dfs/all.data.response.0.lifespan.csv", row.names = 1)
table(all.data.ls$local_lifespan)
# 487 perennial, 126 annual

annual.data = subset(all.data.ls, all.data.ls$local_lifespan == "ANNUAL") # 126 data points
# 126 data points from 80 species

perennial.data = subset(all.data.ls, all.data.ls$local_lifespan == "PERENNIAL") 
# 487 data points of 294 species

grass = subset(all.data.ls, all.data.ls$functional_group == "GRASS") 
# 230 data points of 129 species

forb = subset(all.data.ls, all.data.ls$functional_group == "FORB")
# 320 data points of 201 species

#### write out the files ####

#write.csv(all.data, file = "./New.dfs/all.data.csv")
#write.csv(annual.data, file = "./New.dfs/annual.data.csv")
#write.csv(forb, file = "./New.dfs/forb.csv")
#write.csv(grass, file = "./New.dfs/grass.csv")
#write.csv(perennial.data, file = "./New.dfs/perennial.data.csv")











