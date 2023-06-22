# Plotting raw data for Control vs. Drought year 1

library(ggplot2)
library(cowplot)

#### read in all the data frames needs for the analyses ####

# all data without functional group WOODY
no.trees = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/no.trees.csv", row.names = 1) # 618 data points
# all annual data
annual.data = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/annual.data.csv", row.names = 1) # 127 data points
annual.data = subset(annual.data, !annual.data$functional_group == "WOODY") # 125 data points
# perennial data without functional group WOODY
perennial.tree = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/perennial.tree.csv", row.names = 1) # 478 data points
# grasses
grass = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/grass.csv", row.names = 1) # 223 data points
# forbs
forb = read.csv("./Formatted.Data/Ctrl.v.drt.yr1.data/forb.csv", row.names = 1) # 317 data points

#### no.trees ####

leafN = ggplot(no.trees, aes(y = mean.cover.response, x = leafN.mg.g)) +
  geom_point()+
  theme_classic()
height = ggplot(no.trees, aes(y = mean.cover.response, x = height.m)) +
  geom_point()+
  theme_classic()
rootN = ggplot(no.trees, aes(y = mean.cover.response, x = rootN.mg.g)) +
  geom_point()+
  theme_classic()
SLA = ggplot(no.trees, aes(y = mean.cover.response, x = SLA_m2.kg)) +
  geom_point()+
  theme_classic()
depth = ggplot(no.trees, aes(y = mean.cover.response, x = root.depth_m)) +
  geom_point()+
  theme_classic()
RTD = ggplot(no.trees, aes(y = mean.cover.response, x = RTD.g.cm3)) +
  geom_point()+
  theme_classic()
SRL = ggplot(no.trees, aes(y = mean.cover.response, x = SRL_m.g)) +
  geom_point()+
  theme_classic()
diam = ggplot(no.trees, aes(y = mean.cover.response, x = rootDiam.mm)) +
  geom_point()+
  theme_classic()

all.data.plots = plot_grid(leafN,height,rootN,SLA,depth,RTD,SRL,diam)
ggsave("./Plots/ctrl.v.drt.yr1.no.woody.pdf", height = 10, width = 12)
# height plant is 2 individuals of AMBROSIA TENUIFOLIA
# PHRAGMITES AUSTRALIS has the most negative cover response (-27.5)


