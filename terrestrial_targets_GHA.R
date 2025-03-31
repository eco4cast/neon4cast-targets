#renv::restore()
print(paste0("Running Creating Terrestrial Targets at ", Sys.time()))

#Sys.setenv("NEONSTORE_HOME" = "/home/rstudio/data/neonstore")

Sys.setenv("NEONSTORE_HOME" = "/home/rstudio/data/neonstore_temp")

readRenviron("~/.Renviron") # compatible with littler

non_store_dir <- "/home/rstudio/data/neon_flux_data"
use_5day_data <- TRUE

library(neonUtilities)
library(neonstore)
library(tidyverse)
library(lubridate)
library(contentid)
library(ncdf4)
library(arrow)

sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv", show_col_types = FALSE) |> 
  dplyr::filter(terrestrial == 1)

site_names <- sites$field_site_id

s3 <- arrow::s3_bucket("bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D",
                       endpoint_override = "sdsc.osn.xsede.org",
                       access_key = Sys.getenv("OSN_KEY"),
                       secret_key = Sys.getenv("OSN_SECRET"))


s3_neon_portal <- arrow::s3_bucket("bio230014-bucket01/flux_staging/neon_portal",
                                   endpoint_override = "sdsc.osn.xsede.org",
                                   access_key = Sys.getenv("OSN_KEY"),
                                   secret_key = Sys.getenv("OSN_SECRET"))

s3_current_month <- arrow::s3_bucket("bio230014-bucket01/flux_staging/current_month",
                                     endpoint_override = "sdsc.osn.xsede.org",
                                     access_key = Sys.getenv("OSN_KEY"),
                                     secret_key = Sys.getenv("OSN_SECRET"))


existing_targets_daily <- arrow::read_csv_arrow(s3$path("terrestrial_daily-targets-test.csv.gz"))

start_date <- as_date(max(targets_daily$datetime) - months(3))

# Terrestrial

neonstore::neon_download(product = "DP4.00200.001", type = "basic", start_date = start_date)
neonstore::neon_store(product = "DP4.00200.001")
#neon_disconnect()

#DP4.00200.001 & DP1.00094.001
#neon_store(product = "DP4.00200.001") 
flux_data <- neon_table(table = "nsae-basic", site = site_names, lazy = TRUE) %>% 
  mutate(timeBgn = as_datetime(timeBgn),
         timeEnd = as_datetime(timeEnd)) |> 
  select(timeBgn, timeEnd, data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux, qfqm.fluxCo2.turb.qfFinl) |> 
  collect()

flux_data |> 
  mutate(year = year(timeBgn),
         month = month(timeBgn)) |> 
  arrow::write_dataset(s3_neon_portal, partitioning = c("siteID", "year", "month"))

#Get the current unpublished flux data (5-day latency)

if(use_5day_data){
  #files <- readr::read_csv("https://storage.googleapis.com/neon-sae-files/ods/dataproducts/sae_file_url_unpublished.csv", show_col_types = FALSE)
  files <- readr::read_csv("https://storage.googleapis.com/neon-sae-files/ods/sae_files_unpublished/sae_file_url_unpublished.csv", show_col_types = FALSE)
  #Convert old S3 links to GCS
  files$url <- base::gsub("https://s3.data.neonscience.org", "https://storage.googleapis.com", files$url)
  files <- files %>% 
    filter(site %in% site_names) %>% 
    mutate(file_name = basename(url)) |> 
    filter(date > max(lubridate::as_date(flux_data$timeBgn)))
  
  fs::dir_create(file.path(non_store_dir,"current_month"), recurse = TRUE)
  fs::dir_create(file.path(non_store_dir,"current_month_parquet"), recurse = TRUE)
  
  fn_parquet_s3 <- s3_current_month$ls()
  
  for(i in 1:nrow(files)){
    destfile <- file.path(non_store_dir,"current_month",files$file_name[i])
    if(!paste0(files$file_name[i],".parquet") %in% fn_parquet_s3){
      download.file(files$url[i], destfile = destfile)
    }
  }
  
  
  fn_parquet_s3 <- s3_current_month$ls()
  
  #Checking for files that have been updated
  if(length(fn_parquet_s3) > 0){
    d <- tibble(file_full = fn_parquet_s3,
                file = tools::file_path_sans_ext(tools::file_path_sans_ext(fn_parquet_s3)),
                date = lubridate::as_datetime(stringr::str_split(fn_parquet_s3, "\\.", 12, simplify = TRUE)[, 10])) |> 
      arrange(file, date) |> 
      group_by(file) |> 
      mutate(max_date = max(date),
             delete = ifelse(date != max(date), 1, 0)) |> 
      ungroup()
    
    for(i in 1:nrow(d)){
      if(d$delete[i] == 1){
        s3_current_month$DeleteFile(d$file[i])
      }
    }
    #remove files that are no longer in the unpublished s3 bucket because they are now in NEON portal
    
    for(i in 1:length(fn_parquet_s3)){
      if(!tools::file_path_sans_ext(basename(fn_parquet_s3[i])) %in% files$file_name){
        message(paste0("removing: ", basename(fn_parquet_s3[i])))
        s3_current_month$DeleteFile(fn_parquet_s3[i])
      }
    }
  }
  
  fn_parquet_s3 <- s3_current_month$ls()
  fn <- list.files(file.path(non_store_dir,"current_month"), full.names = TRUE)
  
  if(length(fn) > 0){
    message(paste0("reading in ", length(fn), " non-NEON portal files"))
    purrr::walk(1:length(fn), function(i, fn){
      message(paste0(i, " of ", length(fn), " reading file ",fn[i]))
      df <- neonstore::neon_read(files = fn[i]) |> 
        select(timeBgn, timeEnd, data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux, qfqm.fluxCo2.turb.qfFinl) 
      arrow::write_parquet(x = df, s3_current_month$path(paste0(basename(fn[i]),".parquet")))
      unlink(fn[i])
    },
    fn = fn)
  }
  
  flux_data_curr <- arrow::open_dataset(s3_current_month) |> 
    select(timeBgn, timeEnd, data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux, qfqm.fluxCo2.turb.qfFinl) |> 
    collect()
  
  #Combined published and unpublished
  
  flux_data <- bind_rows(flux_data, flux_data_curr)
}

#flux_data |> 
#  group_by(siteID) |> 
#  summarize(min = min(timeBgn),
#            max = max(timeBgn), .groups = "drop") |> 
#  arrange(min) |> 
#  print(n = 50) 

flux_data <- flux_data %>% 
  mutate(time = as_datetime(timeBgn))

co2_data <- flux_data %>% 
  filter(qfqm.fluxCo2.turb.qfFinl == 0 & data.fluxCo2.turb.flux > -50 & data.fluxCo2.turb.flux < 50) %>% 
  select(time,data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux) %>% 
  rename(nee = data.fluxCo2.turb.flux,
         le = data.fluxH2o.turb.flux,
         site_id = siteID) |> 
  mutate(nee = ifelse(site_id == "OSBS" & year(time) < 2019, NA, nee),
         nee = ifelse(site_id == "SRER" & year(time) < 2019, NA, nee),
         nee = ifelse(site_id == "BARR" & year(time) < 2019, NA, nee),
         le = ifelse(site_id == "OSBS" & year(time) < 2019, NA, le),
         le = ifelse(site_id == "SRER" & year(time) < 2019, NA, le),
         le = ifelse(site_id == "BARR" & year(time) < 2019, NA, le)) |> 
  pivot_longer(-c("time","site_id"), names_to = "variable", values_to = "observation")

flux_target_30m <- co2_data #left_join(full_time, co2_data, by = c("time", "site_id", "variable"))

valid_dates <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  filter(!is.na(observation)) %>%
  group_by(date, site_id, variable) %>% 
  summarise(count = n(), .groups = "drop")

flux_target_daily <- flux_target_30m %>% 
  mutate(date = as_date(time)) %>% 
  group_by(date, site_id, variable) %>% 
  summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |> 
  left_join(valid_dates, by = c("date","site_id", "variable")) |> 
  mutate(observation = ifelse(count > 24, observation, NA),
         observation = ifelse(is.nan(observation), NA, observation)) %>% 
  rename(time = date) %>% 
  select(-count) |> 
  mutate(observation = ifelse(variable == "nee", (observation * 12 / 1000000) * (60 * 60 * 24), observation))

#flux_target_daily %>% 
#  filter(year(time) > 2021) %>% 
#  ggplot(aes(x = time, y = observation)) + 
#  geom_point() +
#  facet_grid(variable~site_id, scale = "free")


readRenviron("~/.Renviron") # compatible with littler
# Write 30 minute data

#flux_target_30m <- flux_target_30m |>
#  select(time, site_id, variable, observation, type) |> 
#  rename(datetime = time) |> 
#  mutate(datetime = lubridate::as_datetime(datetime),
#         duration = "PT30M",
#         project_id = "neon4cast") |>
#  select(project_id, site_id, datetime, duration, variable, observation, type) |> 
#  arrange(variable, datetime) |> 
#  na.omit()


#combined_30_min <- bind_rows(flux_target_30m, targets_30min) |> 
#  arrange(project_id, site_id, datetime, duration, variable, type) |> 
#  group_by(project_id, site_id, datetime, duration, variable) |> 
#  slice(1) |> 
#  ungroup()

#s3 <- arrow::s3_bucket("bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=PT30M",
#                       endpoint_override = "sdsc.osn.xsede.org",
#                       access_key = Sys.getenv("OSN_KEY"),
#                       secret_key = Sys.getenv("OSN_SECRET"))

#arrow::write_csv_arrow(combined_30_min, sink = s3$path("terrestrial_30min-targets.csv.gz"))

#Write daily data

flux_target_daily <- flux_target_daily |> 
  select(time, site_id, variable, observation) |> 
  rename(datetime = time) |> 
  mutate(datetime = lubridate::as_datetime(datetime),
         duration = "P1D",
         project_id = "neon4cast") |>
  select(project_id, site_id, datetime, duration, variable, observation) |> 
  arrange(variable, datetime) |> 
  na.omit()


min_year <- year(min(as_date(flux_target_daily$datetime)))
min_month <- month(min(as_date(flux_target_daily$datetime)))

old_fluxes <- existing_targets_daily |> 
  filter(datetime < min(as_date(flux_target_daily$datetime)))
         
combined_daily <- bind_rows(flux_target_daily, old_fluxes) |> 
  arrange(project_id, site_id, datetime, duration, variable)

s3 <- arrow::s3_bucket("bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D",
                       endpoint_override = "sdsc.osn.xsede.org",
                       access_key = Sys.getenv("OSN_KEY"),
                       secret_key = Sys.getenv("OSN_SECRET"))


arrow::write_csv_arrow(combined_daily, sink = s3$path("terrestrial_daily-targets-test.csv.gz"))

message(paste0("Completed Terrestrial Target at ", Sys.time()))
