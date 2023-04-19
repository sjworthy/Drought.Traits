library(stringr)
library(dplyr)
library(tidyverse)

try.data = read.csv("/Users/sjworthy/Documents/Courses/ECL271/TRYTRAIT.csv")

# remove rows that don't include trait data - most include metadata we don't need
try.data.2 <- subset(try.data, try.data$TraitName != "")

# get species list

yr.1.species = read.csv("./Formatted.Data/trait.species.trt.yr1.csv")

sp.list=str_to_sentence(yr.1.species$Taxon)

# subset try data with year 1 drought species

try.data.3=subset(try.data.2, try.data.2$AccSpeciesName %in% sp.list)

#### Leaf Area ####
# Condensing to a single Leaf Area variable

# we chose to combine measurements with and without petiole that used leaflet measurements rather than compound leaf measurements
leaf_area_keep <- c('Leaf area (in case of compound leaves: leaflet, petiole excluded)',
                    'Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)',
                    "Leaf area (in case of compound leaves: leaflet, petiole included)")

leaf_area_toss <- c("Leaf area (in case of compound leaves: leaf, petiole excluded)", 
                    "Leaf area (in case of compound leaves: leaf, petiole included)", 
                    "Leaf area per plant")

# reassign TraitName to "Leaf area"
try.data.3$TraitName[try.data.3$TraitName %in% leaf_area_keep] <- "Leaf area"

# remove other non-usable "leaf area" columns
try.data.3 <- subset(try.data.3, !try.data.3$TraitName %in% leaf_area_toss)

#### Specific Leaf Area ####
# we chose to combine all SLA measurements regardless of whether petiole was included or excluded
SLA_keep <- c("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", 
              "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", 
              "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included")

try.data.3$TraitName[try.data.3$TraitName %in% SLA_keep] <- "Specific leaf area"

# remove entries that correspond to maximum or minimum trait values measured
#keeping everything that is either a single obs or a mean/median/best estimate
measurements_keep <- c("Best estimate", "Mean", "Median", "Single", "Site specific mean", "Species mean")

try.data.3 <- subset(try.data.3, try.data.3$ValueKindName %in% measurements_keep)

#### only keep traits we need ####

keep.traits = c("Coarse root nitrogen (N) content per coarse root dry mass",
                "Fine root carbon (C) content per fine root dry mass",
                "Fine root length per fine root dry mass (specific fine root length, SRL)",
                "Fine root nitrogen (N) content per fine root dry mass",
                "Leaf area","Leaf carbon (C) content per leaf dry mass",
                "Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)",
                "Leaf nitrogen (N) content per leaf dry mass","Leaf thickness",
                "Plant height vegetative","Root carbon (C) content per root dry mass",
                "Root diameter","Root nitrogen (N) content per root dry mass",
                "Root rooting depth","Seed dry mass","Specific leaf area")

try.data.4 = subset(try.data.3, try.data.3$TraitName %in% keep.traits)

### summarize ####
# summarize means by study
summary_data <- try.data.4 %>% group_by(DatasetID, AccSpeciesName, TraitName) %>% 
  summarise(mean_value=mean(StdValue, na.rm = T),
            .groups = 'drop')

# summarize means across study
summary <- summary_data %>% group_by(AccSpeciesName, TraitName) %>% 
  summarise(mean_value=mean(mean_value, na.rm = T),
            .groups = 'drop')

summary_wide <- pivot_wider(data = summary, names_from = TraitName, values_from = mean_value)

write.csv(summary_wide, file="./Formatted.Data/try.traits.yr1.redo.csv")


