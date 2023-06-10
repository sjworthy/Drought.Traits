# models of year 2 - year 1 control cover

# read in cover year 2
cover.yr.2 = read.csv("./Formatted.Data/cover.trt.y2.csv")

# read in cover data
cover.data = read.csv("./Formatted.Data/cover.response.trt.y1.csv")

# merge cover year 1 with cover year 2

years.cover = merge(cover.data, cover.yr.2, by = c("site_code", "Taxon"))
years.cover = years.cover[,c(1,2,6,15)]
colnames(years.cover)[3] = "mean.ctr.cover.yr1"
colnames(years.cover)[4] = "mean.ctr.cover.yr2"

# remove 0 if wasn't in plot in year 1 and year 2

years.cover.2 = subset(years.cover, !(years.cover$mean.ctr.cover.yr1 == 0))
years.cover.3 = subset(years.cover.2, !(years.cover.2$mean.ctr.cover.yr2 == 0))

years.cover.3$diff.ctr.trt = years.cover.3$mean.ctr.cover.yr2 - years.cover.3$mean.ctr.cover.yr1

# merge new cover data with traits

trait.data.new = read.csv("./Formatted.Data/trait.species.trt.yr1.outlier.2.csv")

# subset traits so they must have SLA
trait.data.2 = trait.data.new[,c(1,7,8,10,12,14,15,18,20,27:29)]
trait.data.3 = subset(trait.data.2, trait.data.2$SLA_m2.kg > 0 ) # 645 data points, 

all.data.year2 = merge(years.cover.3, trait.data.3, by="Taxon") # 619 data points







