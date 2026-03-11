################################################################################
# Lichen Trait analysis

# Author: Jonas Dotterweich
# Date: 2026-03-06
# Purpose: Validation/support of lichen grouping and selection for utilizing as OGF indicators
################################################################################

# Load packages
library(tidyverse)
library(readxl)
library(dplyr)





# Load data

traits <- read_excel("Lichens/traits/lichen_traits_spruce.xlsx")

# Explore data
names(traits)

# Make column names easier to type
colnames(traits) <- c(
  "species",
  "macrolichen",
  "thallus_type",
  "photobiont",
  "metabolites",
  "veget_rep",    # Was "Vegetative", now "vegetative reproduction"
  "conidia",
  "fruit_body_type",
  "fruit_body_area",
  "spore_septation",
  "spore_pigmentation",
  "spore_shape",
  "spore_volume"
)

# Check it worked
colnames(traits)

# For each important column, see unique values

# Thallus type
unique(traits$thallus_type)
table(traits$thallus_type)

# Photobiont type
unique(traits$photobiont)
table(traits$photobiont)

# Reproduction mode
unique(traits$veget_rep)
table(traits$veget_rep)

# Spore volume (numeric)
summary(traits$spore_volume)
hist(traits$spore_volume, breaks=30, main="Spore Volume Distribution")


# Count missing values per column
colSums(is.na(traits))

# Which species have missing spore volume?
traits$species[is.na(traits$spore_volume)]


# Your 30 indicator species from all groups

my_indicators <- c(
  # Group 1: Parmelia
  "Parmelia saxatilis",
  "Parmelia sulcata",
  
  # Group 2: Mycoblastus
  "Mycoblastus sanguinarius",
  "Mycoblastus affinis",
  
  # Group 3: Ochrolechia
  "Ochrolechia alboflavescens",
  "Ochrolechia microstictoides",
  "Ochrolechia androgyna",
  "Ochrolechia mahluensis",
  
  # Group 4: Calicioids
  "Calicium viride",
  "Calicium glaucellum",
  "Calicium trabinellum",
  "Chaenotheca sphaerocephala",
  "Chaenotheca ferruginea",
  "Chaenotheca chrysocephala",
  "Chaenotheca xyloxena",
  "Chaenotheca trichialis",
  "Chaenotheca brunneola",
  "Chaenotheca furfuracea",
  "Chaenotheca stemonea",
  "Chaenothecopsis pusilla",
  "Mycocalicium subtile",
  
  # Group 5: Xylographa
  "Xylographa vitiligo",
  "Xylographa parallela",
  
  # Group 6: Elite Rare (some overlap with above)
  # Already included above
  
  # Group 7: Core OGF
  "Lecanora subintricata",
  "Micarea globulosella"
)

# Remove duplicates
my_indicators <- unique(my_indicators)

# How many indicator species?
length(my_indicators)

# Filter trait data to only your indicators
trait_indicators <- traits %>%
  filter(species %in% my_indicators)


nrow(trait_indicators)

# Check which indicators are missing from trait data
missing_indicators <- setdiff(my_indicators, trait_indicators$species)

missing_indicators

# Add group assignments to your filtered data
trait_indicators <- trait_indicators %>%
  mutate(
    group = case_when(
      species %in% c("Parmelia saxatilis", "Parmelia sulcata") ~ "Parmelia",
      species %in% c("Mycoblastus sanguinarius", "Mycoblastus affinis") ~ "Mycoblastus",
      species %in% c("Ochrolechia alboflavescens", "Ochrolechia microstictoides", 
                     "Ochrolechia androgyna", "Ochrolechia mahluensis") ~ "Ochrolechia",
      species %in% c("Calicium viride", "Calicium glaucellum", "Calicium trabinellum",
                     "Chaenotheca sphaerocephala", "Chaenotheca ferruginea",
                     "Chaenotheca chrysocephala", "Chaenotheca xyloxena",
                     "Chaenotheca trichialis", "Chaenotheca brunneola",
                     "Chaenotheca furfuracea", "Chaenotheca stemonea",
                     "Chaenothecopsis pusilla", "Mycocalicium subtile") ~ "Calicioids",
      species %in% c("Xylographa vitiligo", "Xylographa parallela") ~ "Xylographa",
      species %in% c("Lecanora subintricata", "Micarea globulosella") ~ "Other_OGF",
      TRUE ~ "Ungrouped"
    )
  )

traits <- traits %>%
  mutate(
    indicator_status = ifelse(species %in% my_indicators, 
                              "OGF_indicator", 
                              "Non_indicator")
  )

# Check it worked
table(traits$indicator_status)


# Simple view of all your indicators with key traits
trait_indicators %>%
  select(species, group, thallus_type, photobiont, 
         veget_rep, spore_volume) %>%
  arrange(group, species) %>%
  print(n = 30) 


## simlifying the vegetative reproduction column to just "sexual" vs "other" (which includes asexual and both)
#for both traits and trait indicators

traits <- traits %>%
  mutate(
    reproduction_simple = case_when(
      veget_rep == "absent" ~ "Sexual_only",
      veget_rep %in% c("soredia", "isidia") ~ "Asexual",
      TRUE ~ "Other"
    )
  )

trait_indicators <- trait_indicators %>%
  mutate(
    reproduction_simple = case_when(
      veget_rep == "absent" ~ "Sexual_only",
      veget_rep %in% c("soredia", "isidia") ~ "Asexual",
      TRUE ~ "Other"
    )
  )


#Indicators vs Non-Indicators - Simple Comparison

# === REPRODUCTION MODE ===
cat("=== REPRODUCTION MODE ===\n")
table(traits$reproduction_simple, traits$indicator_status)

# Percentages
round(prop.table(table(traits$reproduction_simple, traits$indicator_status), 
                 margin = 2) * 100, 1)

# === PHOTOBIONT ===
cat("\n=== PHOTOBIONT TYPE ===\n")
traits_lich <- filter(traits, photobiont != "absent")
table(traits_lich$photobiont, traits_lich$indicator_status)

round(prop.table(table(traits_lich$photobiont, traits_lich$indicator_status), 
                 margin = 2) * 100, 1)

# === SPORE VOLUME ===
cat("\n=== SPORE VOLUME ===\n")
traits %>%
  group_by(indicator_status) %>%
  summarize(
    n = n(),
    mean_spore = round(mean(spore_volume), 1),
    median_spore = round(median(spore_volume), 1)
  )




### Puts out a table with the grouping of traits by indicator group
# Summary by group - 
trait_indicators %>%
  group_by(group) %>%
  summarize(
    n_species = n(),
    
    # Reproduction - CORRECTED
    n_sexual = sum(veget_rep == "absent"),
    pct_sexual = round(mean(veget_rep == "absent") * 100, 0),
    
    # List what reproduction types present
    repro_types = paste(unique(veget_rep), collapse = ", "),
    
    # Photobiont - SIMPLIFIED (just show what you have)
    photobiont_types = paste(unique(photobiont), collapse = ", "),
    
    # Thallus - most common type
    main_thallus = names(sort(table(thallus_type), decreasing = TRUE))[1],
    
    # Spores
    mean_spore = round(mean(spore_volume), 1),
    median_spore = round(median(spore_volume), 1)
  ) %>%
  arrange(desc(pct_sexual))


### with these chuncks the single groups can be investigated: # CALICIOIDS
cat("=== CALICIOIDS ===\n")
trait_indicators %>%
  filter(group == "Calicioids") %>%
  select(species, veget_rep, photobiont, spore_volume) %>%
  print(n = 20)

# XYLOGRAPHA
cat("\n=== XYLOGRAPHA ===\n")
trait_indicators %>%
  filter(group == "Xylographa") %>%
  select(species, veget_rep, photobiont, spore_volume)

# PARMELIA
cat("\n=== PARMELIA ===\n")
trait_indicators %>%
  filter(group == "Parmelia") %>%
  select(species, veget_rep, photobiont, spore_volume)

# MYCOBLASTUS
cat("\n=== MYCOBLASTUS ===\n")
trait_indicators %>%
  filter(group == "Mycoblastus") %>%
  select(species, veget_rep, photobiont, spore_volume)

# OCHROLECHIA
cat("\n=== OCHROLECHIA ===\n")
trait_indicators %>%
  filter(group == "Ochrolechia") %>%
  select(species, veget_rep, photobiont, spore_volume)
