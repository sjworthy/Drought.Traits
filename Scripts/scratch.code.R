# Scratch Code

#### Drought Severity Index ####
# need average drought severity index for each site

yr1.ppt=read.csv("./Raw.Data/cover_ppt_2023-05-10.csv")

# ppt.1 column is the amount of ppt that fell in the plot in the 365 days prior
# to biomass collection
# map column is the mean annual precipitation at the site

drt.trt=read.csv("./Raw.Data/Site_Elev-Disturb.csv")
drt.trt.2 = drt.trt[,c(2,25)]

# drought_trt column has the targeted % water removal for the site

# formula: 
# 1. precip.drt = amount ppt received*% targeted removal (as decimal)
# 2. (precip.drt - MAP)/MAP

# merge drt.trt with yr1.ppt

all.drt.data = merge(yr1.ppt, drt.trt.2, by = c("site_code"))

# turn drought_trt into decimal
all.drt.data$drought_trt_dec = as.numeric(all.drt.data$drought_trt)/100

all.drt.data$precip.drt = all.drt.data$ppt.1*all.drt.data$drought_trt_dec
all.drt.data$drt.sev.index = (all.drt.data$precip.dr - all.drt.data$map)/all.drt.data$map

# slim dataset to match our cover data

all.drt.data.2 = subset(all.drt.data, all.drt.data$n_treat_years == 1)

# get the mean drt.sev.index for each site

site.drt.sev.index = all.drt.data.2 %>%
  group_by(site_code) %>%
  reframe(mean.drt.sev.index = mean(drt.sev.index, na.rm = TRUE))

# write.csv(site.drt.sev.index, "./Formatted.Data/site.drt.dev.index.csv")

