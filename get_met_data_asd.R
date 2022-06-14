# Script to download met, soil and water level data from the NOBV database and 
# then aggregate to 30 min timesteps.

library(tidyverse)
library(lubridate)
library(tictoc)

source('src/db_get_data.R')
source('src/db_config.R')
source('src/db_process_data.R')

tz = 'UTC'
Sys.setenv(TZ=tz)

# Helper functions for processing EC data
make_dir <- function(path) {
  # Makes a directory if it doesn't exist
  if(dir.exists(path) == FALSE) dir.create(path, recursive = TRUE)
}

site_id = 'asd_mp_ec03'
site = 'ASD'
plot_mt1 <- 'RF'
plot_spg <- 'MP'
plot_wlg <- 'MP'

# Number of months to download
start_date = as_date('2022-03-01')
end_date   = as_date('2022-05-01')

# Dir to save csv files
data_dir <- file.path('~/Data/ECProcessing/met_data', site_id)

dates = seq.Date(start_date, end_date-1, 1)
#------------------------------------------------------------------------------#
# Download raw data
#------------------------------------------------------------------------------#
download_met_data()

#------------------------------------------------------------------------------#
# Process data
#------------------------------------------------------------------------------#
process_met_data(overwrite = TRUE)
