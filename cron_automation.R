#remotes::install_github("rqthomas/cronR")
#remotes::install_deps()
library(cronR)

home_dir <-  path.expand("~")
log_dir <- path.expand("~/log/cron")

targets_repo <- "neon4cast-targets"

## Phenocam Download and Target Generation

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, targets_repo, "phenology_targets.R"),
                           rscript_log = file.path(log_dir, "phenology-targets.log"),
                           log_append = FALSE,
                           workdir = file.path(home_dir, targets_repo),
                           trailing_arg = "curl -fsS -m 10 --retry 5 -o /dev/null https://hc-ping.com/f5d48d96-bb41-4c21-b028-930fa2b01c5a")
#cronR::cron_add(command = cmd,  frequency = '0 */2 * * *', id = 'phenocam_download')
cronR::cron_add(command = cmd,  frequency = 'daily', at = '2PM', id = 'phenocam-targets')

## Aquatics Targets

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, targets_repo,"aquatics_targets.R"),
                           rscript_log = file.path(log_dir, "aquatics-target.log"),
                           log_append = FALSE,
                           workdir = file.path(home_dir, targets_repo),
                           cmd = "/usr/local/bin/r", # use litter, more robust on CLI
                           trailing_arg = "curl -fsS -m 10 --retry 5 -o /dev/null https://hc-ping.com/1267b13e-8980-4ddf-8aaa-21aa7e15081c")
cronR::cron_add(command = cmd, frequency = 'daily', at = "7AM", id = 'aquatics-targets')

## Beetles

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, targets_repo,"beetles_targets.R"),
                           rscript_log = file.path(log_dir, "beetles-targets.log"),
                           log_append = FALSE,
                           workdir = file.path(home_dir, targets_repo),
                           trailing_arg = "curl -fsS -m 10 --retry 5 -o /dev/null https://hc-ping.com/ed35da4e-01d3-4750-ae5a-ad2f5dfa6e99")
cronR::cron_add(command = cmd, frequency = "0 10 * * SUN", id = 'beetles-targets')

## Terrestrial targets

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, targets_repo,"terrestrial_targets.R"),
                           rscript_log = file.path(log_dir, "terrestrial-targets.log"),
                           log_append = FALSE,
                           cmd = "/usr/local/bin/r", # use litter, more robust on CLI
                           workdir = file.path(home_dir, targets_repo),
                           trailing_arg = "curl -fsS -m 10 --retry 5 -o /dev/null https://hc-ping.com/c1fb635f-95f8-4ba2-a348-98924548106c")
cronR::cron_add(command = cmd, frequency = 'daily', at = "9AM", id = 'terrestrial-targets')

## Ticks

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, targets_repo,"ticks_targets.R"),
                           rscript_log = file.path(log_dir, "ticks-targets.log"),
                           log_append = FALSE,
                           workdir = file.path(home_dir, targets_repo),
                           trailing_arg = "curl -fsS -m 10 --retry 5 -o /dev/null https://hc-ping.com/09c7ab10-eb4e-40ef-a029-7a4addc3295b")
cronR::cron_add(command = cmd, frequency = "0 11 * * SUN", id = 'ticks-targets')


cronR::cron_ls()


