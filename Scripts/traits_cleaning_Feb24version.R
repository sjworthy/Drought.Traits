# Dylan, Regina, Ben, Jenny's ECL 271 Group
# taking TRY output with cleaned species list
# Cleaning traits and then summarizing
# Averaged within study and then between studies

library(dplyr)
library(tidyverse)

data <- read.csv("~/Downloads/TRY_droughtnetsp.csv")

# remove rows that don't include trait data - most include metadata we don't need
data <- subset(data, data$TraitName != "")

# Specific Leaf Area
# Condensing to a single Leaf Area variable

# we chose to combine measurements with and without petiole that used leaflet measurements rather than compound leaf measurements
leaf_area_keep <- c('Leaf area (in case of compound leaves: leaflet, petiole excluded)',
                    'Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)',
                    "Leaf area (in case of compound leaves: leaflet, petiole included)")

leaf_area_toss <- c("Leaf area (in case of compound leaves: leaf, petiole excluded)", 
                    "Leaf area (in case of compound leaves: leaf, petiole included)", 
                    "Leaf area per plant")

# reassign TraitName to "Leaf area"
data$TraitName[data$TraitName %in% leaf_area_keep] <- "Leaf area"

# remove other non-usable "leaf area" columns
data <- subset(data, !data$TraitName %in% leaf_area_toss)

# Specific Leaf Area
# we chose to combine all SLA measurements regardless of whether petiole was included or excluded
SLA_keep <- c("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", 
              "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", 
              "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included")

data$TraitName[data$TraitName %in% SLA_keep] <- "Specific leaf area"

# remove entries that correspond to maximum or minimum trait values measured
#keeping everything that is either a single obs or a mean/median/best estimate
measurements_keep <- c("Best estimate", "Mean", "Median", "Single", "Site specific mean", "Species mean")

data <- subset(data, data$ValueKindName %in% measurements_keep)


# Removing traits that don't have good coverage
# Trying a cutoff at 10% of species with data for that trait

# how many unique species are represented among each trait?
#unique <- as.data.frame(table(unique(data[c("AccSpeciesName","TraitName")])$TraitName))

#traits_keep <- subset(unique$Var1, unique$Freq > 0.1*length(unique(data$AccSpeciesName)))

# remove traits that are not well represented
#data <- subset(data, data$TraitName %in% traits_keep)

# dropping drought tolerance trait because it's categorical
data <- subset(data, data$TraitName != "Species tolerance to drought")

# we now have 12 traits remaining
# only rooting depth remains of root traits

# summarize means by study
summary_data <- data %>% group_by(DatasetID, AccSpeciesName, TraitName) %>% 
  summarise(mean_value=mean(StdValue, na.rm = T),
            .groups = 'drop')

# summarize means across study
summary <- summary_data %>% group_by(AccSpeciesName, TraitName) %>% 
  summarise(mean_value=mean(mean_value, na.rm = T),
            .groups = 'drop')

# number of species remaining - 1084
length(unique(summary$AccSpeciesName))

# distribution of number of traits across species 
table(table(summary$AccSpeciesName))

# distribution of number of species across traits
trait_coverage <- as.data.frame(table(summary$TraitName))

trait_coverage$Freq <- trait_coverage$Freq/length(unique(summary$AccSpeciesName))

View(trait_coverage)

# pivot the dataframe to have rows = species, columns = traits
summary_wide <- pivot_wider(data = summary, names_from = TraitName, values_from = mean_value)

write.csv(summary_wide, "traits_species_summarized.csv")
