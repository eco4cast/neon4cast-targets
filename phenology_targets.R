#remotes::install_deps()

library(tidyverse)


source("R/downloadPhenoCam.R")
source("R/calculatePhenoCamUncertainty.R")

sites <- readr::read_csv("NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(phenology == 1)

allData <- data.frame(matrix(nrow = 0, ncol = 5))

message(paste0("Downloading and generating phenology targets ", Sys.time()))

for(i in 1:nrow(sites)){
  siteName <- sites$phenocam_code[i]
  site_roi <- sites$phenocam_roi[i]
  message(siteName)
  ##URL for daily summary statistics
  URL_gcc90 <- paste('https://phenocam.nau.edu/data/archive/',siteName,"/ROI/",siteName,"_",site_roi,"_1day.csv",sep="") 
  ##URL for individual image metrics
  URL_individual <- paste('https://phenocam.nau.edu/data/archive/',siteName,"/ROI/",siteName,"_",site_roi,"_roistats.csv",sep="") 
  
  phenoData <- download.phenocam(URL = URL_gcc90)
  dates <- unique(phenoData$date)
  phenoData_individual <- download.phenocam(URL=URL_individual,skipNum = 17)
  ##Calculates standard deviations on daily gcc90 values
  gcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates)    
  rcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates,target="rcc")
  
  subPhenoData <- phenoData %>% 
    mutate(site_id = stringr::str_sub(siteName, 10, 13), 
           time = date) %>% 
    select(time, site_id, gcc_90, rcc_90) |> 
    pivot_longer(-c("time", "site_id"), names_to = "variable", values_to = "observed") |> 
    mutate(sd = ifelse(variable == "gcc_90", gcc_sd, rcc_sd))

  allData <- rbind(allData,subPhenoData)
  
}

full_time <- seq(min(allData$time),max(allData$time), by = "1 day")
combined <- NULL

for(i in 1:nrow(sites)){
  
  full_time_curr1 <- tibble(time = full_time,
                           site_id = rep(sites$field_site_id[i],length(full_time)),
                           variable = "gcc_90")
  
  full_time_curr2 <- tibble(time = full_time,
                            site_id = rep(sites$field_site_id[i],length(full_time)),
                            variable = "rcc_90")
  
  combined <- bind_rows(combined, full_time_curr1, full_time_curr2)
}



allData2 <- left_join(combined, allData, by = c("time", "site_id", "variable"))

readr::write_csv(allData2, "phenology-targets.csv.gz")

aws.s3::put_object(file = "phenology-targets.csv.gz", 
                   object = "phenology/phenology-targets.csv.gz",
                   bucket = "neon4cast-targets")

unlink("phenology-targets.csv.gz")
