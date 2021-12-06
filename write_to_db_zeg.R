library(tidyverse)
library(lubridate)
library(dbplyr)

#------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------#
get_fluxnet <- function(x) {
  src = file.path(eddypro_dir, year(x), format(x, '%m'))
  f = list.files(src, 'fluxnet')
  if(length(f) == 1) {
    df = read_csv(file.path(src, f)) %>% 
      rename(datetime = TIMESTAMP_END) %>% 
      select(datetime, any_of(default_fluxnet_variables)) %>% 
      mutate(datetime = as_datetime(as.character(datetime),
                                    format = '%Y%m%d%H%M',
                                    tz = 'UTC')) %>% 
      filter(datetime >= start_dt, datetime < end_dt)
  } else {
    df = tibble()
  }
  return(df)
}

get_meteo <- function(x) {
  src = file.path(db_data_dir, year(x), format(x, '%m'))
  
  mt1 <- read_csv(file.path(src, paste0(strftime(x, '%Y-%m_MT1.csv'))))
  spg <- read_csv(file.path(src, paste0(strftime(x, '%Y-%m_SPG.csv'))))
  # WLG doesn't always exist
  wlg_file = file.path(src, paste0(strftime(x, '%Y-%m_WLG.csv')))
  if(file.exists(wlg_file)) {
    wlg <- read_csv(wlg_file)
    df <- reduce(list(mt1,spg,wlg), left_join)
  } else {
    df <- left_join(mt1, spg)
  }
  
  # Filter dates
  df <- df %>% 
    rename(datetime = dt) %>% 
    filter(datetime >= start_dt, datetime < end_dt)
  
  return(df)
}

default_fluxnet_variables <- 
  c('H', 'LE', 'ET', 'FC', 'FCH4', 'SC_SINGLE', 'WS', 'WS_MAX', 'WD', 'WD_SIGMA', 'USTAR', 'MO_LENGTH',
    'T_SONIC', 'TA_EP', 'PA_EP', 'RH_EP', 'VPD_EP', 'TDEW', 'U_SIGMA', 'V_SIGMA', 'W_SIGMA',
    'FC_SSITC_TEST', 'FCH4_SSITC_TEST', 'BADM_LOCATION_LAT', 'BADM_LOCATION_LONG','BADM_INST_HEIGHT_SA', 
    'BADM_HEIGHTC', 'DISPLACEMENT_HEIGHT', 'ROUGHNESS_LENGTH', 'CUSTOM_RSSI_77_MEAN', 'CUSTOM_CO2_SIGNAL_STRENGTH_7500_MEAN')
#------------------------------------------------------------------------------#

# Dates to extract and write to DB
start_dt <- as_datetime('2020-05-01T00:00:00')
end_dt   <- as_datetime('2021-12-01T00:00:00')

# Path to eddy pro output (fluxnet files)
eddypro_dir = '~/Data/Zegveld/processed'
# Path to data extracted from NOBV database
db_data_dir = '~/Data/Zegveld/Data/processed'

# Path/name of DB to create/write to
db = "data/Zegveld.db"

# Check if the DB exists
db_exists = file.exists(db)
# Create DB or connect to an existing one
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db)

if(db_exists) {
  # Update the start date if the db exists
  if('data' %in% DBI::dbListTables(con)) {
    start_dt <- tbl(con, 'data') %>%
      select(datetime) %>% 
      filter(datetime == max(datetime)) %>% 
      collect() %>% 
      distinct() %>% 
      mutate(datetime = as_datetime(datetime, tz='UTC')) %>% 
      pull(datetime)
    start_dt <- start_dt + 1
  }
}

# Function expects data to be arranged yyyy/mm underneath the provided path
if(strftime(end_dt, '%d%H%M%S') == '01000000') {
  interval = interval(floor_date(start_dt, 'month')+1, end_dt-1) 
} else {
  interval = interval(floor_date(start_dt, 'month')+1, end_dt) 
}
months = floor_date(start_dt, 'month') %m+% months(0:(interval %/% months(1)))

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
