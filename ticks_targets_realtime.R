# Download and clean tick data 
Sys.setenv("NEONSTORE_HOME" = "/home/rstudio/data/neonstore")
#Sys.setenv("NEONSTORE_DB" = "/home/rstudio/data/neonstore")
#Sys.setenv("NEONSTORE_DB")

### Script to create tidy dataset of NEON tick abundances 
### To be used for the RCN Tick Forecasting Challenge 
# Tick Data
# link data from lab and field
# filter out poor quality data
# only included training data 
# through end of 2020
# A. amercanum nymphs
# standardize nymphs across forested plots
# export a .csv

### Load packages
# devtools::install_github("NEONScience/NEON-utilities/neonUtilities")
# library(neonUtilities) # for downloading data (use GitHub version)

# library(usethis)

library(tidyverse) # for data wrangling and piping (dplyr probably ok)
library(lubridate) # for finding year from dates
library(stringr) # for searching within character strings 
library(here) # for working within subdirectories
library(parallel) # for using more than one core in download
library(uuid) # for unique IDs
library(MMWRweek) # for converting from date to MMWR week

# select target species and life stage
target.species <- "Amblyomma americanum" # NEON species name
target.lifestage <- "Nymph"

sites.df <- read_csv("Ticks_NEON_Field_Site_Metadata_20210928.csv")
target.sites <- sites.df %>% pull(field_site_id)

if(!"neonstore" %in% installed.packages()){
  library(remotes)
  remotes::install_github("cboettig/neonstore", ref = "patch/api-updates")
}

library(neonstore)

efi_server <- TRUE

# get data from neon
#product <- "DP1.10093.001"
#neon_download(product = product,
#              site = target.sites)

tick.field.raw <- neon_read("tck_fielddata-basic", keep_filename = TRUE)
tick.taxon.raw <- neon_read("tck_taxonomyProcessed-basic")

# there are lots of reasons why sampling didn't occur (logistics, too wet, too cold, etc.)
# so, keep records when sampling occurred
tick.field <- tick.field.raw %>% 
  filter(totalSampledArea > 0) %>% 
  mutate(time = floor_date(collectDate, unit = "day")) %>% 
  unite(namedLocation, time, col = "occasionID", sep = "_")

# combine adults into single category and make wide to get zero counts
tick.taxon.wide <- tick.taxon.raw %>% 
  filter(sampleCondition == "OK") %>% # remove taxonomy samples with quality issues
  mutate(sexOrAge = if_else(sexOrAge == "Female" | sexOrAge == "Male", 
                            "Adult",     # convert to Adult
                            sexOrAge),
         time = floor_date(collectDate, unit = "day")) %>% 
  unite(namedLocation, time, col = "occasionID", sep = "_") %>% 
  pivot_wider(id_cols = occasionID, # make wide by species and life stage
              names_from = c(acceptedTaxonID, sexOrAge),
              values_from = individualCount, 
              names_sep = "_",
              values_fn = {sum}, # duplicates occur because of Adults that where F/M - add them 
              values_fill = 0)

# join taxonomy and field data
tick.joined <- left_join(tick.taxon.wide, tick.field, by = "occasionID") %>% 
  select(-NA_NA, -geodeticDatum, -samplingImpractical, -targetTaxaPresent,
         -adultCount, -nymphCount, -larvaCount, -samplingProtocolVersion, 
         -measuredBy, -sampleCode, -biophysicalCriteria, -plotType)

# all the species column names
spp.cols <- tick.joined %>% 
  select(contains("Larva"), contains("Nymph"), contains("Adult")) %>% 
  colnames()

# get matching taxon ids
taxon.ids <- tick.taxon.raw %>%
  filter(!is.na(acceptedTaxonID)) %>% 
  select(acceptedTaxonID, scientificName, taxonRank) %>% 
  distinct() 

# make longer
tick.long <- tick.joined %>% 
  pivot_longer(cols = all_of(spp.cols), 
               names_to = "taxonAge",
               values_to = "processedCount",
               values_drop_na = TRUE) %>% 
  separate(col = taxonAge, into = c("acceptedTaxonID", "lifeStage"), sep = "_")

# add taxon ids
tick.long <- left_join(tick.long, taxon.ids, by = "acceptedTaxonID") 

# standardize the data and subset to targets
tick.standard <- tick.long %>% 
  filter(siteID %in% target.sites, # sites we want
         lifeStage == target.lifestage, # life stage we want
         scientificName == target.species, # species we want
         grepl("Forest", nlcdClass)) %>%  # forest plots
  mutate(date = floor_date(collectDate, unit = "day"),
         date = ymd(date),
         year = year(date),
         iso_week = ISOweek::ISOweek(collectDate),
         time = ISOweek::ISOweek2date(paste0(iso_week, "-1"))) %>% 
  select(time, processedCount, totalSampledArea, siteID) %>%
  mutate(totalSampledArea = as.numeric(totalSampledArea)) %>% 
  group_by(siteID, time) %>%
  summarise(totalCount = sum(processedCount), # all counts in a week
            totalArea = sum(totalSampledArea),# total area surveyed in a week
            amblyomma_americanum = totalCount / totalArea * 1600) %>% # scale to the size of a plot
  mutate(iso_week = ISOweek::ISOweek(time)) %>% 
  arrange(siteID, time) %>% 
  filter() %>% 
  select(time, iso_week, siteID, amblyomma_americanum)


tick.targets <- tick.standard %>% 
  #filter(time < challenge.time) |> 
  rename(site_id = siteID,
         observed = amblyomma_americanum) |> 
  mutate(variable = "amblyomma_americanum") |> 
  select(time, site_id, variable, observed, iso_week)
  

# write targets to csv
write_csv(tick.targets,
          file = "ticks-targets.csv.gz")

aws.s3::put_object(file = "ticks-targets.csv.gz", 
                   object = "ticks/ticks-targets.csv.gz",
                   bucket = "neon4cast-targets")

