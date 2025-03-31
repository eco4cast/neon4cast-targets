#renv::restore()
print(paste0("Running Creating Terrestrial Targets at ", Sys.time()))

#Sys.setenv("NEONSTORE_HOME" = "/home/rstudio/data/neonstore")

Sys.setenv("NEONSTORE_HOME" = "/home/rstudio/data/neonstore_temp")

pecan_flux_uncertainty <- "../pecan/modules/uncertainty/R/flux_uncertainty.R"
readRenviron("~/.Renviron") # compatible with littler

non_store_dir <- "/home/rstudio/data/neon_flux_data"
use_5day_data <- TRUE
add_TERN <- FALSE

library(neonUtilities)
library(neonstore)
library(tidyverse)
library(lubridate)
library(contentid)
library(ncdf4)
library(arrow)

sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(terrestrial == 1)

site_names <- sites$field_site_id

# Terrestrial

neonstore::neon_download(product = "DP4.00200.001", type = "basic")
neonstore::neon_store(product = "DP4.00200.001")
#neon_disconnect()

#DP4.00200.001 & DP1.00094.001
#neon_store(product = "DP4.00200.001") 
flux_data <- neon_table(table = "nsae-basic", site = site_names, lazy = TRUE) %>% 
  mutate(timeBgn = as_datetime(timeBgn),
         timeEnd = as_datetime(timeEnd)) |> 
  select(timeBgn, timeEnd, data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux, qfqm.fluxCo2.turb.qfFinl) |> 
  mutate(year = year(timeBgn),
         month = month(timeBgn)) |> 
  collect()

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
  
  for(i in 1:nrow(files)){
    destfile <- file.path(non_store_dir,"current_month",files$file_name[i])
    parquet_file <- file.path(non_store_dir,"current_month_parquet",paste0(files$file_name[i],".parquet"))
    if(!(file.exists(parquet_file))){
      download.file(files$url[i], destfile = destfile)
    }
  }
  
  fn_parquet <- list.files(file.path(non_store_dir,"current_month_parquet"))
  
  #Checking for files that have been updated
  if(length(fn_parquet) > 0){
    d <- tibble(file_full = fn_parquet,
                file = tools::file_path_sans_ext(tools::file_path_sans_ext(fn_parquet)),
                date = lubridate::as_datetime(stringr::str_split(fn_parquet, "\\.", 12, simplify = TRUE)[, 10])) |> 
      arrange(file, date) |> 
      group_by(file) |> 
      mutate(max_date = max(date),
             delete = ifelse(date != max(date), 1, 0)) |> 
      ungroup()
    
    for(i in 1:nrow(d)){
      if(d$delete[i] == 1){
        unlink(file.path(non_store_dir,"current_month_parquet", d$file_full[i]))
      }
    }
    #remove files that are no longer in the unpublished s3 bucket because they are now in NEON portal
    
    for(i in 1:length(fn_parquet)){
      if(!tools::file_path_sans_ext(basename(fn_parquet[i])) %in% files$file_name){
        message(paste0("removing: ", basename(fn_parquet[i])))
        unlink(file.path(non_store_dir,"current_month_parquet",fn_parquet[i]))
      }
    }
  }
  
  fn <- list.files(file.path(non_store_dir,"current_month"), full.names = TRUE)
  #fn_parquet <- file.path(non_store_dir,"current_month_parquet", paste0(basename(fn[i]),".parquet"))
  
  if(length(fn) > 0){
    message(paste0("reading in ", length(fn), " non-NEON portal files"))
    
    #future::plan("future::multisession", workers = 4)
    
    #new_files <- fn[which(!(basename(fn) %in% tools::file_path_sans_ext(basename(fn_parquet))))]
    purrr::walk(1:length(fn), function(i, fn){
      message(paste0(i, " of ", length(fn), " reading file ",fn[i]))
      df <- neonstore::neon_read(files = fn[i]) |> 
        select(timeBgn, timeEnd, data.fluxCo2.turb.flux, siteID, data.fluxH2o.turb.flux, qfqm.fluxCo2.turb.qfFinl)
      arrow::write_parquet(x = df, file.path(non_store_dir,"current_month_parquet", paste0(basename(fn[i]),".parquet")))
      unlink(fn[i])
    },
    fn = fn)
  }
  
  #s3 <- arrow::SubTreeFileSystem$create(file.path(non_store_dir,"current_month_parquet"))
  
  flux_data_curr <- arrow::open_dataset(file.path(non_store_dir,"current_month_parquet")) |> 
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

#co2_data %>% 
#  filter(variable == "nee") |> 
#  ggplot(aes(x = time, y = observation)) +
#  geom_point() +
#  facet_wrap(~site_id)

#earliest <- min(as_datetime(c(co2_data$time)), na.rm = TRUE)
#latest <- max(as_datetime(c(co2_data$time)), na.rm = TRUE)


#full_time_vector <- seq(min(c(co2_data$time), na.rm = TRUE), 
#                        max(c(co2_data$time), na.rm = TRUE), 
#                        by = "30 min")

#full_time <- NULL
#for(i in 1:length(site_names)){
#  df_nee <- tibble(time = full_time_vector,
#                   site_id = rep(site_names[i], length(full_time_vector)),
#                   variable = "nee")
#  df_le <- tibble(time = full_time_vector,
#                  site_id = rep(site_names[i], length(full_time_vector)),
#                  variable = "le")
#  full_time <- bind_rows(full_time, df_nee, df_le)
#  
#}

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

flux_target_30m <- flux_target_30m |> 
  select(time, site_id, variable, observation) |> 
  rename(datetime = time)

flux_target_daily <- flux_target_daily |> 
  select(time, site_id, variable, observation) |> 
  rename(datetime = time)


if(add_TERN){
  
  sites <- read_csv("tern_field_site_metadata.csv", show_col_types = FALSE) |> 
    filter(!is.na(data_url))
  
  read_tern_site <- function(i, sites){
    nc <- ncdf4::nc_open(sites$data_url[i])
    
    co2 <- ncdf4::ncvar_get(nc, "Fco2")
    le <- ncdf4::ncvar_get(nc, "Fe")
    time <- ncdf4::ncvar_get(nc, "time")
    
    time <- lubridate::as_datetime("1800-01-01 00:00:00.0") + lubridate::seconds(time * 86400)
    
    co2 <- tibble::tibble(datetime = time,
                          site_id = sites$`EFI ID`[i],
                          variable = "nee",
                          observation = co2) |> 
      mutate(observation = ifelse(observation < -100, NA, observation))
    
    le <- tibble::tibble(datetime = time,
                         site_id = sites$`EFI ID`[i],
                         variable = "le",
                         observation = le) |> 
      mutate(observation = ifelse(observation < -100, NA, observation))
    
    bind_rows(co2, le)
  }
  
  tern_flux_target_30m <- map_dfr(1:nrow(sites), read_tern_site, sites)
  
  full_time <- NULL
  for(i in 1:nrow(sites)){
    
    site_targets <- tern_flux_target_30m |> filter(site_id == sites$`EFI ID`[i])
    
    full_time_vector <- seq(min(c(site_targets$datetime), na.rm = TRUE), 
                            max(c(site_targets$datetime), na.rm = TRUE), 
                            by = "30 min")
    
    df_nee <- tibble(datetime = full_time_vector,
                     site_id = rep(sites$`EFI ID`[i], length(full_time_vector)),
                     variable = "nee")
    df_le <- tibble(datetime = full_time_vector,
                    site_id = rep(sites$`EFI ID`[i], length(full_time_vector)),
                    variable = "le")
    full_time <- bind_rows(full_time, df_nee, df_le)
    
  }
  
  tern_flux_target_30m <- left_join(full_time, tern_flux_target_30m, by = c("datetime", "site_id", "variable"))
  
  valid_dates <- tern_flux_target_30m %>% 
    mutate(date = as_date(datetime)) %>% 
    filter(!is.na(observation)) %>%
    group_by(date, site_id, variable) %>% 
    summarise(count = n(), .groups = "drop")
  
  
  tern_flux_target_daily <- tern_flux_target_30m %>% 
    mutate(date = as_date(datetime)) %>% 
    group_by(date, site_id, variable) %>% 
    summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |> 
    left_join(valid_dates, by = c("date","site_id", "variable")) |> 
    mutate(observation = ifelse(count > 24, observation, NA),
           observation = ifelse(is.nan(observation), NA, observation)) %>% 
    rename(datetime = date) %>% 
    select(-count) |> 
    mutate(observation = ifelse(variable == "nee", (observation * 12 / 1000000) * (60 * 60 * 24), observation))
  
  tern_flux_target_daily %>% 
    filter(year(datetime) > 2021) %>% 
    ggplot(aes(x = datetime, y = observation)) + 
    geom_point() +
    facet_grid(variable~site_id, scale = "free")
  
  ggplot(tern_flux_target_daily, aes(x = datetime, y = observation)) + 
    geom_point() + 
    facet_wrap(~variable, scale = "free")
  
  
  flux_target_30m <- bind_rows(flux_target_30m, tern_flux_target_30m)
  
  flux_target_daily <- bind_rows(flux_target_daily, tern_flux_target_daily)
  
}

#write_csv(flux_target_30m, "terrestrial_30min-targets.csv.gz")
#write_csv(flux_target_daily, "terrestrial_daily-targets.csv.gz")

readRenviron("~/.Renviron") # compatible with littler
# Write 30 minute data

s3 <- arrow::s3_bucket("neon4cast-targets/terrestrial_30min",
                       endpoint_override = "data.ecoforecast.org",
                       access_key = Sys.getenv("AWS_ACCESS_KEY"),
                       secret_key = Sys.getenv("AWS_SECRET_ACCESS_KEY"))

flux_target_30m |> 
  arrange(variable, datetime) |> 
  arrow::write_csv_arrow(sink = s3$path("terrestrial_30min-targets.csv.gz"))

flux_target_30m2 <- flux_target_30m |>
  mutate(datetime = lubridate::as_datetime(datetime),
         duration = "PT30M",
         project_id = "neon4cast") |>
  select(project_id, site_id, datetime, duration, variable, observation) |> 
  arrange(variable, datetime) |> 
  na.omit()

s3 <- arrow::s3_bucket("bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=PT30M",
                       endpoint_override = "sdsc.osn.xsede.org",
                       access_key = Sys.getenv("OSN_KEY"),
                       secret_key = Sys.getenv("OSN_SECRET"))

arrow::write_csv_arrow(flux_target_30m2, sink = s3$path("terrestrial_30min-targets.csv.gz"))

#Write daily data

s3 <- arrow::s3_bucket("neon4cast-targets/terrestrial_daily",
                       endpoint_override = "data.ecoforecast.org",
                       access_key = Sys.getenv("AWS_ACCESS_KEY"),
                       secret_key = Sys.getenv("AWS_SECRET_ACCESS_KEY"))

flux_target_daily |> 
  arrange(variable, datetime) |> 
  arrow::write_csv_arrow(sink = s3$path("terrestrial_daily-targets.csv.gz"))

flux_target_daily2 <- flux_target_daily |>
  mutate(datetime = lubridate::as_datetime(datetime),
         duration = "P1D",
         project_id = "neon4cast") |>
  select(project_id, site_id, datetime, duration, variable, observation) |> 
  arrange(variable, datetime) |> 
  na.omit()

s3 <- arrow::s3_bucket("bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D",
                       endpoint_override = "sdsc.osn.xsede.org",
                       access_key = Sys.getenv("OSN_KEY"),
                       secret_key = Sys.getenv("OSN_SECRET"))


arrow::write_csv_arrow(flux_target_daily2, sink = s3$path("terrestrial_daily-targets.csv.gz"))

message(paste0("Completed Terrestrial Target at ", Sys.time()))
