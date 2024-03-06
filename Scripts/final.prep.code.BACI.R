# Script to prep trait data for analyses

library(tidyverse)
library(corrplot)

# read in full data frame 
trait.data = read.csv("./New.dfs/trait.species.BACI.final.csv")

# subset traits to only those of interest
trait.data.2 = trait.data[,c(1:9,14:16)]

# subset out woody species and trees from local_lifeform

trait.data.3 = subset(trait.data.2, !trait.data.2$functional_group == "WOODY")
trait.data.3 = subset(trait.data.3, !trait.data.3$local_lifeform == "TREE")

# BACI analysis:
# 244 new species, 29 removed because TREE or WOODY, 1 FERN removed, 4 CACTUS removed
# these were removed prior to adding species and traits to the end of trait.data
# 210 total new species added (if traits available)

# subset for species with SLA
trait.data.4 = subset(trait.data.3, trait.data.3$SLA_m2.kg > 0 ) # 692 taxon

#write.csv(trait.data.4, file = "./New.dfs/BACI.traits.no.woody.SLA.csv")

#### checking for outliers ####
# all outliers removed manually from traits.no.woody.SLA.csv and 
# made new file BACI.traits.no.woody.SLA.no.outlier.csv

hist(trait.data.4$leafN.mg.g)
boxplot(trait.data.4$leafN.mg.g)
mean = mean(trait.data.4$leafN.mg.g, na.rm = TRUE)
std = sd(trait.data.4$leafN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$leafN.mg.g[which(trait.data.4$leafN.mg.g <Tmin | trait.data.4$leafN.mg.g > Tmax)])
# removed leafN 54.82000 58.30000 67.58333
# percent removed 
table(is.na(trait.data.4$leafN.mg.g))
692-207
3/485*100 #0.62%

hist(trait.data.4$height.m)
boxplot(trait.data.4$height.m)
mean = mean(trait.data.4$height.m, na.rm = TRUE)
std = sd(trait.data.4$height.m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$height.m[which(trait.data.4$height.m <Tmin | trait.data.4$height.m > Tmax)])
# removed height 2.600709 - 10
# percent removed 
table(is.na(trait.data.4$height.m))
692-62
7/630*100 #1.11%

hist(trait.data.4$rootN.mg.g)
boxplot(trait.data.4$rootN.mg.g)
mean = mean(trait.data.4$rootN.mg.g, na.rm = TRUE)
std = sd(trait.data.4$rootN.mg.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootN.mg.g[which(trait.data.4$rootN.mg.g <Tmin | trait.data.4$rootN.mg.g > Tmax)])
# remove rootN 34.39000 34.84397 39.26274
# percent removed 
table(is.na(trait.data.4$rootN.mg.g))
692-490
3/202*100 #1.49%

hist(trait.data.4$SLA_m2.kg)
boxplot(trait.data.4$SLA_m2.kg)
mean = mean(trait.data.4$SLA_m2.kg, na.rm = TRUE)
std = sd(trait.data.4$SLA_m2.kg, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SLA_m2.kg[which(trait.data.4$SLA_m2.kg <Tmin | trait.data.4$SLA_m2.kg > Tmax)])
# remove SLA 52.45601 - 98.20000
# percent removed 
table(is.na(trait.data.4$SLA_m2.kg))
13/692*100 #1.88%

hist(trait.data.4$root.depth_m)
boxplot(trait.data.4$root.depth_m)
mean = mean(trait.data.4$root.depth_m, na.rm = TRUE)
std = sd(trait.data.4$root.depth_m, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$root.depth_m[which(trait.data.4$root.depth_m <Tmin | trait.data.4$root.depth_m > Tmax)])
# remove depth 2.500000 - 3.828571
# percent removed 
table(is.na(trait.data.4$root.depth_m))
692-353
9/339*100 #2.65%

hist(trait.data.4$RTD.g.cm3)
boxplot(trait.data.4$RTD.g.cm3)
mean = mean(trait.data.4$RTD.g.cm3, na.rm = TRUE)
std = sd(trait.data.4$RTD.g.cm3, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$RTD.g.cm3[which(trait.data.4$RTD.g.cm3 <Tmin | trait.data.4$RTD.g.cm3 > Tmax)])
# remove RTD 0.7839500 - 0.9167417
# percent removed 
table(is.na(trait.data.4$RTD.g.cm3))
692-458
3/234*100 #1.28%

hist(trait.data.4$SRL_m.g)
boxplot(trait.data.4$SRL_m.g)
mean = mean(trait.data.4$SRL_m.g, na.rm = TRUE)
std = sd(trait.data.4$SRL_m.g, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$SRL_m.g[which(trait.data.4$SRL_m.g <Tmin | trait.data.4$SRL_m.g > Tmax)])
# remove SRL 652.6000 - 929.9185
# percent removed 
table(is.na(trait.data.4$SRL_m.g))
692-380
7/312*100 #2.24%

hist(trait.data.4$rootDiam.mm)
boxplot(trait.data.4$rootDiam.mm)
mean = mean(trait.data.4$rootDiam.mm, na.rm = TRUE)
std = sd(trait.data.4$rootDiam.mm, na.rm = TRUE)
Tmin = mean-(3*std)
Tmax = mean+(3*std)
sort(trait.data.4$rootDiam.mm[which(trait.data.4$rootDiam.mm <Tmin | trait.data.4$rootDiam.mm > Tmax)])
# remove diam 1.242171 - 2.015308
# percent removed 
table(is.na(trait.data.4$rootDiam.mm))
692-387
9/305*100 #2.95%

#### read in new trait data without outliers and cover data ####

trait.data.new = read.csv("./New.dfs/BACI.traits.no.woody.SLA.no.outlier.csv", row.names = 1)
cover.data = read.csv("./New.dfs/BACI.data.final.csv", row.names = 1)

# merge trait data and cover data

all.data = merge(cover.data, trait.data.new, by="Taxon") 
# 692 species from 1234 data points from 83 sites

# write.csv(all.data, file = "./New.dfs/BACI.all.data.csv")

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

table(all.data.ls$functional_group)
# 320 forb, 24 graminoid, 230 grass, 54 legume

grass = subset(all.data.ls, all.data.ls$functional_group == "GRASS") 
# 230 data points of 129 species

table(grass$local_lifespan)
# 41 annuals, 187 perennials

forb = subset(all.data.ls, all.data.ls$functional_group == "FORB")
# 320 data points of 201 species

# grasses and graminoids
graminoid = subset(all.data.ls, all.data.ls$functional_group %in% c("GRASS","GRAMINOID"))

# grass.annuals
grass.annual = subset(grass, grass$local_lifespan == "ANNUAL")
# 41 data points

# grass.perennial
grass.perennial = subset(grass, grass$local_lifespan == "PERENNIAL")
# 187 data points

graminoid.perennial = subset(graminoid, graminoid$local_lifespan == "PERENNIAL")
# 210 data points 
# graminoid annual is same as grass.annual since all graminoids were perennials

#### merge trait data and BACI cover data ####

trait.data.new = read.csv("./New.dfs/traits.no.woody.SLA.no.outlier.csv", row.names = 1)
cover.data = read.csv("./New.dfs/BACI.data.csv", row.names = 1)

# subset trait data for only species in cover data

subset.traits = subset(trait.data.new, trait.data.new$Taxon %in% cover.data$Taxon)

# merge trait data and cover data

all.data = merge(cover.data, subset.traits, by="Taxon") 
# 275 species from 407 data points from 65 sites

# write.csv(all.data, file = "./New.dfs/BACI.all.data.response.0.csv")

# split by lifeform and functional group

all.data.ls = read.csv("./New.dfs/BACI.all.data.response.0.lifespan.csv", row.names = 1)
table(all.data.ls$local_lifespan)
# 321 perennial, 79 annual

annual.data = subset(all.data.ls, all.data.ls$local_lifespan == "ANNUAL")
# 79 data points from 56 species

perennial.data = subset(all.data.ls, all.data.ls$local_lifespan == "PERENNIAL") 
# 321 data points of 212 species

table(all.data.ls$functional_group)
# 201 forb, 15 graminoid, 160 grass, 31 legume

grass = subset(all.data.ls, all.data.ls$functional_group == "GRASS") 
# 160 data points of 105 species

table(grass$local_lifespan)
# 31 annuals, 127 perennials

forb = subset(all.data.ls, all.data.ls$functional_group == "FORB")
# 201 data points of 135 species

# grass.perennial
grass.perennial = subset(grass, grass$local_lifespan == "PERENNIAL")
# 127 data points of 82 species

#### write out the files ####

#write.csv(all.data, file = "./New.dfs/all.data.csv")
#write.csv(annual.data, file = "./New.dfs/annual.data.csv")
#write.csv(forb, file = "./New.dfs/forb.csv")
#write.csv(grass, file = "./New.dfs/grass.csv")
#write.csv(perennial.data, file = "./New.dfs/perennial.data.csv")
#write.csv(graminoid, file = "./New.dfs/graminoid.data.csv")
#write.csv(grass.annual, file = "./New.dfs/grass.annual.data.csv")
#write.csv(grass.perennial, file = "./New.dfs/grass.perennial.data.csv")
#write.csv(graminoid.perennial, file = "./New.dfs/graminoid.perennial.data.csv")

#write.csv(all.data, file = "./New.dfs/all.data.BACI.csv")
#write.csv(annual.data, file = "./New.dfs/annual.data.BACI.csv")
#write.csv(forb, file = "./New.dfs/forb.BACI.csv")
#write.csv(grass, file = "./New.dfs/grass.BACI.csv")
#write.csv(perennial.data, file = "./New.dfs/perennial.data.BACI.csv")
#write.csv(grass.perennial, file = "./New.dfs/grass.perennial.data.BACI.csv")














