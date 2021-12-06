library(tidyverse)
library(lubridate)
library(dbplyr)

tz = 'UTC'
Sys.setenv(tz=tz)

#------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------#
get_fluxnet <- function(x) {
  src = file.path(eddypro_dir, year(x), format(x, '%m'))
  f = list.files(src, 'fluxnet')
  if(length(f) == 1) {
    df = read_csv(file.path(src, f), locale=locale(tz=tz)) %>% 
      rename(datetime = TIMESTAMP_END) %>% 
      select(datetime, any_of(default_fluxnet_variables)) %>% 
      mutate(datetime = as_datetime(as.character(datetime),
                                    format = '%Y%m%d%H%M',
                                    tz = 'UTC')) %>% 
      filter(datetime >= start_date, datetime < end_date)
  } else {
    df = tibble()
  }
  return(df)
}

get_meteo <- function(x) {
  src = file.path(db_data_dir, year(x), format(x, '%m'))
  
  mt1 <- read_csv(file.path(src, paste0(strftime(x, '%Y-%m_MT1.csv'))), locale=locale(tz=tz))
  spg <- read_csv(file.path(src, paste0(strftime(x, '%Y-%m_SPG.csv'))), locale=locale(tz=tz))
  wlg <- read_csv(file.path(src, paste0(strftime(x, '%Y-%m_WLG.csv'))), locale=locale(tz=tz))
  
  # WLG doesn't always exist
  if(nrow(wlg) > 0) {
    df <- reduce(list(mt1,spg,wlg), left_join)
  } else {
    df <- left_join(mt1, spg)
  }
  
  # Filter dates
  df <- df %>% 
    rename(datetime = dt) %>% 
    filter(datetime >= start_date, datetime < end_date)
  
  return(df)
}

default_fluxnet_variables <- 
  c('H', 'LE', 'ET', 'FC', 'FCH4', 'SC_SINGLE', 'WS', 'WS_MAX', 'WD', 'WD_SIGMA', 'USTAR', 'MO_LENGTH',
    'T_SONIC', 'TA_EP', 'PA_EP', 'RH_EP', 'VPD_EP', 'TDEW', 'U_SIGMA', 'V_SIGMA', 'W_SIGMA',
    'FC_SSITC_TEST', 'FCH4_SSITC_TEST', 'BADM_LOCATION_LAT', 'BADM_LOCATION_LONG','BADM_INST_HEIGHT_SA', 
    'BADM_HEIGHTC', 'DISPLACEMENT_HEIGHT', 'ROUGHNESS_LENGTH', 'CUSTOM_RSSI_77_MEAN', 'CUSTOM_CO2_SIGNAL_STRENGTH_7500_MEAN')
#------------------------------------------------------------------------------#

# Dates to extract and write to DB
start_date <- as_date('2021-08-01') # inclusive
end_date <- as_date('2021-11-21') # exclusive
# Path to eddy pro output (fluxnet files)
eddypro_dir = '~/Data/Assendelft/processed'
# Path to data extracted from NOBV database
db_data_dir = '~/Data/Assendelft/Data/processed'

# Path/name of DB to create/write to
db = "data/Assendelft.db"

# Check if the DB exists
db_exists = file.exists(db)
# Create DB or connect to an existing one
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db)

if(db_exists) {
  # Update the start date if the db exists
  if('data' %in% DBI::dbListTables(con)) {
    start_date <- tbl(con, 'data') %>%
      select(datetime) %>% 
      filter(datetime == max(datetime)) %>% 
      collect() %>% 
      mutate(datetime = as_datetime(datetime, tz='UTC')) %>% 
      pull(datetime)
    start_date <- as_date(start_date) + 1
  }
}

# Function expects data to be arranged yyyy/mm underneath the provided path
months = seq(start_date, end_date, 1) %>% 
  round_date('month') %>% 
  unique()

# Bind list. Replace -9999 with NA and parse datetimes
fluxnet <- lapply(months, get_fluxnet) %>% 
  bind_rows() %>% 
  mutate(across(where(is.numeric), ~na_if(., -9999)))

# Read in monthly dfs of meteo, soil and water level data
db_data = lapply(months, get_meteo) %>% 
  bind_rows() %>% 
  distinct()

# Join tables together and pivot to long format
df = full_join(fluxnet, db_data) %>%
  pivot_longer(-datetime)

id_name = df %>%
  select(name) %>%
  distinct() %>%
  mutate(id_name = 1L:n()) %>%
 select(id_name, name)

df = left_join(df, id_name) %>%
  select(datetime, id_name, value) %>%
  mutate(datetime = as.integer(datetime))

if(!db_exists) {
  # Create the table
  DBI::dbCreateTable(con, 'id_name', id_name)
  DBI::dbCreateTable(con, 'data', df)
}

# Write/append data to DB
DBI::dbAppendTable(con, 'id_name', id_name)
DBI::dbAppendTable(con, 'data', df)

# Close connection
DBI::dbDisconnect(con)
