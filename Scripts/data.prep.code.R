library(Taxonstand)
library(tidyverse)
library(stringr)

#### Prepping final species list ####

# read in raw IDE cover data
data=read.csv("./Raw.Data/IDE_cover_2023-01-02.csv")

# get taxon list
taxon.df=as.data.frame(table(data$Taxon))
# 1904 species, but some are sp. so will need to be excluded
# 35977 rows

# excluded species with sp.

data.2=subset(data, !(data$Taxon %in% c("ADESMIA SP.(eea.br)","AGOSERIS SP.(OREAC.US)","AGROSTIS sp.","AGROSTIS SP.",
                                "AGROSTIS SP.(SCRUZH.US)","ALLIUM SP.", "ALLIUM SP.(BFL.US)","ALOPECURUS  SP.(LCSOUTH.CL)",
                                "ALOPECURUS SP.(LCNORTH.CL)","ALOPECURUS SP.(PURDUE.US)","AMBROSIA SP.(PURDUE.US)",
                                "ANTHOXANTHUM  SP.(LCSOUTH.CL)","ARISTIDA SP.","ASTER SP.(BFL.US)","ASTER SP.1(MILPARINKA.AU)",
                                "ASTER SP.7(QUILPIE.AU)","ASTRAGALUS  SP.(LCSOUTH.CL)","ASTRAGALUS SP.(LCNORTH.CL)",
                                "ATRIPLEX SP","ATRIPLEX SP.(MILPARINKA.AU)","AUSTROSTIPA SP.(JILPANGER.AU)","AVENA SP.",
                                "AXONOPUS SP.(BIVENSARM.US)","BACCHARIS SP.(chacra.ar)","BOTHRIOCHLOA SP.(NAPOSTA.AR)",
                                "BRACHYCOME SP.","BROMUS SP.","BRYOPHYTE","BRYOPHYTE SP.(CHACRA.AR)","BRYOPHYTE SP.(KERNB.CA)",
                                "BRYOPHYTE SP.(KERNNU.CA)","BRYOPHYTE SP.(LYGRAINT.NO)","BRYOPHYTE SP.(lygraold.no)",
                                "BRYOPHYTE SP.1(LYGRAINT.NO)","BRYOPHYTE SP.1(lygraold.no)","BRYOPHYTE SP.2(LYGRAINT.NO)",
                                "BRYOPHYTE SP.2(lygraold.no)","BRYOPHYTE SP.3(lygraold.no)","BULBOSTYLIS SP.(UKULINGADRT.ZA)",
                                "CAREX SP.","CAREX SP.(BIDDULPH.CA)","CAREX SP.(BIVENSARM.US)","CAREX SP.(KERNB.CA)",
                                "CAREX SP.(KERNNU.CA)","CAREX SP.(KINSELLA.CA)","CAREX SP.(LYGRAINT.NO)","CAREX SP.(lygraold.no)",
                                "CAREX SP.(MATADOR.CA)","CAREX SP.(MATTHEIS.CA)","CAREX SP.(OREAC.US)","CAREX SP.1(LYGRAINT.NO)",
                                "CASTILLEJA SP.(OKLAH.US)","CENTAUREA sp.","CENTAUREA SP.","CENTAUREA SP.(ayora.es)",
                                "CERASTIUM  SP.(LCSOUTH.CL)","CERASTIUM SP.(LCNORTH.CL)","CLADONIA SP.(LYGRAINT.NO)",
                                "CLADONIA SP.1(LYGRAINT.NO)","CONVOLVULUS SP.(JILPANGER.AU)","CRASSULA SP.(JILPANGER.AU)",
                                "CREPIS SP.(OREAA.US)","CREPIS SP.(OREAC.US)","CRYPTANTHA SP.(JRNCHI.US)","CYPERUS SP.",
                                "CYPERUS SP.(BIVENSARM.US)","DANTHONIA SP.(eea.br)","DELPHINIUM SP.(OREAC.US)","DESCURAINIA SP.(SPVDRT.AR)",
                                "DIANTHUS SP.(qdtsouth.cl)","DICHANTHELIUM SP.(SLP.US)","DICHANTHELIUM SP.2(SLP.US)","DICRANUM SP.(LYGRAINT.NO)",
                                "DICRANUM SP.(lygraold.no)","DIOSCOREA  S2P.(LCSOUTH.CL)","DIOSCOREA  SP1.(LCSOUTH.CL)","DIOSCOREA SP1.(LCNORTH.CL)",
                                "DIOSCOREA SP1.(QDTNORTH.CL)","DIOSCOREA SP1.(qdtsouth.cl)","DIOSCOREA SP2.(qdtsouth.cl)","ELEOCHARIS SP.(OKLAH.US)",
                                "ERAGROSTIS SP.(OKLAH.US)","ERIGERON SP.(KERNB.CA)","ERIGERON SP.(KERNNU.CA)","ERIGERON SP.(MATADOR.CA)",
                                "ERIOGONUM SP.(NNSS.US)","EUONYMUS SP.","FORB sp.","FORB SP.(KERNB.CA)","FORB SP.(KERNNU.CA)","GAILLARDIA SP.(bfl.us)",
                                "GAMOCHAETA SP.(CHACRA.AR)","GAMOCHAETA SP.(NAPOSTA.AR)","GAMOCHAETA SP.(SPVDRT.AR)","GERANIUM  SP.(LCSOUTH.CL)",
                                "GERANIUM SP.(JILPANGER.AU)","GERANIUM SP.(LCNORTH.CL)","GERANIUM SP.(MILPARINKA.AU)","GILIA SP.(NNSS.US)",
                                "GILIA SP.(OREAA.US)","GILIA SP.(OREAC.US)","GUTIERREZIA SP.(BFL.US)","HELIANTHEMUM  SP.(LCSOUTH.CL)",
                                "HELIANTHEMUM SP.(LCNORTH.CL)","HELIANTHUS SP.(BROOKDALE.CA)","HIBBERTIA SP.(JILPANGER.AU)","HIERACIUM SP.(ELVADRT.EE)",
                                "HYMENOPAPPUS SP. (cdpt_drt.us)","HYPOCHAERIS SP.","HYPOCHAERIS SP.(SCRUZH.US)","HYPOCHAERIS SP.(SCRUZM.US)",
                                "IPOMOEA SP.(BIVENSARM.US)","JUNCUS SP.","JUNCUS SP.(BIVENSARM.US)","JUNCUS SP.(LYGRAINT.NO)","JUNCUS SP.(SLP.US)",
                                "JUNELLIA SP.(CHACRA.AR)","LATHYRUS SP.(BROOKDALE.CA)","LEONTODON SP.","LEPIDIUM SP.(NNSS.US)","LESQUERELLA SP.",
                                "LEUCOCORYNE SP1.(QDTNORTH.CL)","LEUCOCORYNE SP2.(QDTNORTH.CL)","LICHEN ","LICHEN","LICHEN SP.(LYGRAINT.NO)","LICHEN SP.(lygraold.no)",
                                "LICHEN SP.1(lygraold.no)","LINUM SP.(BFL.US)","LITHOPHRAGMA SP.(OREAA.US)","LITHOPHRAGMA SP.(OREAC.US)","LOLIUM SP.(CHACRA.AR)",
                                "LOMATIUM SP.(OREAA.US)","LOMATIUM SP.(OREAC.US)","LOTUS","LOTUS SP.","LUZULA SP.","LUZULA SP.(FALLS.AU)","MADIA SP.(OREAC.US)",
                                "MALVA SP.","MALVA SP.(SCRUZL.US)","MELILOTUS SP.(ESW.CA)","MINURIA SP.(credoj.au)","MINURIA SP.(CREDOM.AU)","MONTIA SP.(OREAA.US)",
                                "MONTIA SP.(OREAC.US)","OLEARIA SP.","OPHIOGLOSSUM SP.(JILPANGER.AU)","OXALIS  SP.(LCSOUTH.CL)","OXALIS SP.(LCNORTH.CL)",
                                "OXALIS SP.(QDTNORTH.CL)","OXALIS SP.(qdtsouth.cl)","PAEPALANTHUS SP.GUARIBAS.BR","PANICUM SP.(SLP.US)","PASPALUM SP.",
                                "PASPALUM SP.(SLP.US)","PECTOCARYA SP.","PLANTAGO SP.","PLANTAGO SP.(OKLAH.US)","PLANTAGO SP.(SLP.US)",
                                "POA SP.(LYGRAINT.NO)","POLYGALA SP.","POLYGALA SP.GUARIBAS.BR","POLYGONUM SP.(OREAA.US)","POLYGONUM SP.(OREAC.US)",
                                "POLYPOGON SP.","POTENTILLA SP.(KERNNU.CA)","PRUNELLA SP.","PRUNUS SP.","PTILOTUS SP.","PTILOTUS SP.(MILPARINKA.AU)",
                                "RANUNCULUS sp.","RANUNCULUS SP.","ROSA SP.","RUBUS SP.(BIVENSARM.US)","RUBUS SP.(OKLAH.US)","RUMEX SP.",
                                "RYTIDOSPERMA SP.(JILPANGER.AU)","SCLEROLAENA SP.","SCLEROLAENA SP.1(MILPARINKA.AU)","SCLEROLAENA SP.2(MILPARINKA.AU)",
                                "SCUTELLARIA SP.(SLP.US)","SENECIO SP.(OREAC.US)","SIDA SP.","SIDA SP.(CREDOJ.AU)","SIDA SP.(CREDOM.AU)",
                                "SILENE SP.(KERNNU.CA)","SMILAX SP.(BIVENSARM.US)","SOLIDAGO SP.(NAPOSTA.AR)","SOLIDAGO SP.(OKLAH.US)",
                                "SONCHUS  SP1.(LCSOUTH.CL)","SONCHUS  SP2.(LCSOUTH.CL)","SPHAERALCEA SP.","SPOROBOLUS SP.","Sporobulus sp.",
                                "STELLARIA  SP.(LCSOUTH.CL)","STELLARIA SP.(LCNORTH.CL)","TALINUM SP.(CERRILLOS.AR)","TARAXACUM SP.",
                                "TARAXACUM SP.(ELVADRT.EE)","THELYMITRA SP.(JILPANGER.AU)","TOXICOSCORDION SP.","TRIFOLIUM SP.","TRIFOLIUM SP.(BROOKDALE.CA)",
                                "TRIFOLIUM sp.(jilpanger.au)","TRIFOLIUM SP.(PURDUE.US)","TRIFOLIUM SP.(SCRUZH.US)","ULMUS SP.(OKLAH.US)","UNKNOWN ",
                                "UNKNOWN  SP.(LCSOUTH.CL)","UNKNOWN  SP2.(LCSOUTH.CL)","UNKNOWN A(COWIDRT.CA)","UNKNOWN AMARYLLIDACEAE SP.(HARD.US)",
                                "UNKNOWN ASTERACEAE ","UNKNOWN ASTERACEAE  SP4.(LCSOUTH.CL)","UNKNOWN ASTERACEAE SP.(CREDOM.AU)","UNKNOWN ASTERACEAE SP.(HARD.US)",
                                "UNKNOWN ASTERACEAE SP.(OREAC.US)","UNKNOWN ASTERACEAE SP.2(MILPARINKA.AU)","UNKNOWN ASTERACEAE SP.2(OREAC.US)","UNKNOWN ASTERACEAE SP1.(qdtsouth.cl)",
                                "UNKNOWN ASTERACEAE SP2.(qdtsouth.cl)","UNKNOWN ASTERACEAE SP3.(LCNORTH.CL)","UNKNOWN ASTERACEAE SP3.(QDTNORTH.CL)",
                                "UNKNOWN ASTERACEAE SP3.(qdtsouth.cl)","UNKNOWN ASTERACEAE SP4.(LCNORTH.CL)","UNKNOWN BRASSICACEAE SP.(HARD.US)",
                                "UNKNOWN D(COWIDRT.CA)","UNKNOWN FABACEAE SP.(HARD.US)","UNKNOWN FORB(BIDDULPH.CA)","UNKNOWN G(COWIDRT.CA)","UNKNOWN GRASS",
                                "UNKNOWN GRASS ","UNKNOWN GRASS SP.","UNKNOWN GRASS(COWIDRT.CA)","UNKNOWN H(COWIDRT.CA)","UNKNOWN MINT(COWIDRT.CA)",
                                "UNKNOWN ONAGRACEAE SP.(QDTNORTH.CL)","UNKNOWN POACEAE SP.(HARD.US)","UNKNOWN POACEAE SP.(JILPANGER.AU)","UNKNOWN POACEAE SP.(MATADOR.CA)",
                                "UNKNOWN POACEAE SP.(MILPARINKA.AU)","UNKNOWN POACEAE SP.(QUILPIE.AU)","UNKNOWN POACEAE SP.(SPVDRT.AR)",
                                "UNKNOWN POACEAE SP.1(QUILPIE.AU)","UNKNOWN POACEAE SP1.(qdtsouth.cl)","UNKNOWN POACEAE SP2.(QDTNORTH.CL)",
                                "UNKNOWN POACEAE SP3.(qdtsouth.cl)","UNKNOWN POLEMONIACEAE SP.(OREAA.US)","UNKNOWN POLEMONIACEAE SP.(OREAC.US)",
                                "UNKNOWN SCLEROLAENA ","UNKNOWN SP.","UNKNOWN SP.(CERRILLOS.AR)","UNKNOWN SP.(GUARIBAS.BR)","UNKNOWN SP.(HARD.US)",
                                "UNKNOWN sp.(jilpanger.au)","UNKNOWN SP.(JRNCHI.US)","UNKNOWN SP.(KONZADRT.US)","UNKNOWN SP.(LCNORTH.CL)",
                                "UNKNOWN SP.(LYGRAINT.NO)","UNKNOWN SP.(MATADOR.CA)","UNKNOWN SP.(NNSS.US)","UNKNOWN SP.(OKLAH.US)","SELAGINELLA DENSA",
                                "UNKNOWN SP.(OREAA.US)","UNKNOWN SP.(OREAC.US)","UNKNOWN SP.(SCRUZH.US)","UNKNOWN SP.(SLP.US)","UNKNOWN SP.1(OKLAH.US)",
                                "UNKNOWN SP.11(OKLAH.US)","UNKNOWN SP.2(MILPARINKA.AU)","UNKNOWN SP.2(OKLAH.US)","UNKNOWN SP.4(OKLAH.US)","UNKNOWN SP.4(SPVDRT.AR)",
                                "UNKNOWN SP.5(SPVDRT.AR)","UNKNOWN SP.6(OKLAH.US)","UNKNOWN SP.6(SPVDRT.AR)","UNKNOWN SP.7(OKLAH.US)","UNKNOWN SP.9(OKLAH.US)",
                                "UNKNOWN SP.9(SPVDRT.AR)","UNKNOWN SP.D20(OKLAH.US)","UNKNOWN SP.D22(OKLAH.US)","UNKNOWN SP.D25(OKLAH.US)","UNKNOWN SP.D29(OKLAH.US)",
                                "UNKNOWN SP3.(LCNORTH.CL)","UNKNOWN VIOLACEAE  SP.(LCSOUTH.CL)","UNKNOWN WEED(COWIDRT.CA)","UTRICULARIA SP.(GUARIBAS.BR)",
                                "VERNONIA SP.(PURDUE.US)","VERONICA SP.(OREAA.US)","VICIA SP.","VIOLA sp.","VIOLA SP.","VIOLA SP.(ayora.es)","VIOLA SP.(OREAC.US)",
                                "WAHLENBERGIA SP.","XYRIS SP.","XYRIS SP.GUARIBAS.BR","ACALYPHA SP.(SLP.US)","ATRIPLEX SP.","LASERPITIUM SP.")))

# get taxon list
taxon.2=as.data.frame(table(data.2$Taxon))
# 1605 species, but some are sp. so will need to be excluded
# 33973 rows

# write.csv(taxon.2, "IDE.species.check.csv")

#### verify species names ####
# https://cran.r-project.org/web/packages/Taxonstand/Taxonstand.pdf
taxon.check=Taxonstand::TPL(taxon.2$Var1)
# write.csv(taxon.check, file="taxon.check.csv")


#### calculate cover change ####
# split by drought and control

drought=subset(data.2, data.2$trt == "Drought")
# remove species with negative n_treat_days

drought.2=subset(drought, drought$n_treat_days > 0)
# only keep cover data for after first year of experiment

control=subset(data.2, data.2$trt == "Control")

control.2=subset(control, control$n_treat_days > 0)
# only keep cover data for after first year of experiment

# group by site, Taxon and number of years of treatment
drought.3 = drought.2 %>%
  group_by(site_code, Taxon, n_treat_years) %>%
  summarise(sum.drt.cover = sum(max_cover),
            mean.drt.cover = mean(max_cover)) 

# group by site, Taxon and number of years of treatment
control.3 = control.2 %>%
  group_by(site_code, Taxon, n_treat_years) %>%
  summarise(sum.ctrl.cover = sum(max_cover),
            mean.ctrl.cover = mean(max_cover)) 


# merge data frames together
# match up site code, taxon and treatment years, but put NA if not a match

cover.data=merge(control.3, drought.3, by=c("site_code","Taxon","n_treat_years"),all = TRUE)
# make NA = 0 

cover.data.2=cover.data

cover.data.2[is.na(cover.data.2)] = 0

# calculate cover change

cover.data.2$sum.cover.response = cover.data.2$sum.drt.cover-cover.data.2$sum.ctrl.cover
cover.data.2$mean.cover.response = cover.data.2$mean.drt.cover-cover.data.2$mean.ctrl.cover

# write.csv(cover.data.2, file = "cover.response.all.csv")

# subset for only treatment year = 1

cover.trt.y1=subset(cover.data.2, cover.data.2$n_treat_years == 1)

# write.csv(cover.trt.y1, file = "cover.response.trt.y1.csv")

# get final species list for traits

trait.species=as.data.frame(unique(cover.trt.y1$Taxon))

# write.csv(trait.species, file="trait.species.trt.yr1.csv")

# subset for other treatment years

cover.response.all=read.csv("./Formatted.Data/cover.response.all.csv")

cover.trt.y2=subset(cover.response.all, cover.response.all$n_treat_years == 2)
cover.trt.y3=subset(cover.response.all, cover.response.all$n_treat_years == 3)
cover.trt.y4=subset(cover.response.all, cover.response.all$n_treat_years == 4)
cover.trt.y5=subset(cover.response.all, cover.response.all$n_treat_years == 5)
cover.trt.y6=subset(cover.response.all, cover.response.all$n_treat_years == 6)
cover.trt.y7=subset(cover.response.all, cover.response.all$n_treat_years == 7)
cover.trt.y8=subset(cover.response.all, cover.response.all$n_treat_years == 8)

# get trait species for other year

trait.species.trt.2=as.data.frame(unique(cover.trt.y2$Taxon))
trait.species.trt.3=as.data.frame(unique(cover.trt.y3$Taxon))
trait.species.trt.4=as.data.frame(unique(cover.trt.y4$Taxon))
trait.species.trt.5=as.data.frame(unique(cover.trt.y5$Taxon))
trait.species.trt.6=as.data.frame(unique(cover.trt.y6$Taxon))
trait.species.trt.7=as.data.frame(unique(cover.trt.y7$Taxon))
trait.species.trt.8=as.data.frame(unique(cover.trt.y8$Taxon))

# get number of trt years for each site

trt.yrs.site = as.data.frame(unclass(table(cover.response.all$site_code, cover.response.all$n_treat_years)))

write.csv(trt.yrs.site, file="./Formatted.Data/trt.years.per.site.csv")

# species metadata
data.3=data.2[,c(1,4,12,13,16,17,18)]

data.4=unique(data.3)

write.csv(data.4, file="./Formatted.Data/species.metadata.csv")

#### Austraits ####

install.packages("remotes")
remotes::install_github("traitecoevo/austraits", dependencies = TRUE, upgrade = "ask")
library(austraits) 
austraits <- load_austraits(version = "4.1.0", path = "austraits")

setwd("/Users/sjworthy/Documents/Courses/ECL271/austraits")
austraits.2=readRDS("austraits-4.1.0.rds")

taxon.df=as.data.frame(unique(austraits.2$traits$taxon_name))

# read in species names
setwd("/Users/sjworthy/Documents/GitHub/Drought.Traits/Formatted.Data")

cover.response.all=read.csv("cover.response.all.csv")
cover.species=as.data.frame(unique(cover.response.all$Taxon))

cover.species.list=as.vector(cover.species$`unique(cover.response.all$Taxon)`)
cover.species.list.2=str_to_sentence(cover.species.list)

austraits.subset <- extract_taxa(austraits.2, taxon_name = cover.species.list.2)
austraits.traits <- extract_trait(austraits.subset, c("leaf_area","leaf_C_per_dry_mass",
                                                      "leaf_dry_matter_content","leaf_N_per_dry_mass",
                                                      "leaf_thickness","plant_height","root_C_per_dry_mass",
                                                      "root_diameter","root_N_per_dry_mass","root_shoot_ratio",
                                                      "root_specific_root_length","seed_dry_mass",
                                                      "leaf_mass_per_area"))
austraits.subset.traits=austraits.traits$traits
write.csv(austraits.subset.traits, file="AusTraits.subset.traits.csv")



