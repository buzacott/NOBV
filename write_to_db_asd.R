library(tidyverse)
library(lubridate)
library(dbplyr)

source('src/db_helper.R')

site_name = 'Assendelft'
site_id = 'asd_mp_ec03'
#------------------------------------------------------------------------------#
# Dates to extract and write to DB
start_dt <- as_datetime('2021-08-01T00:00:01')
end_dt   <- as_datetime('2022-05-01T00:00:00')

# Path to eddy pro output (fluxnet files)
eddypro_dir = file.path('~/Data/ECProcessing/processed', site_id)
# Path to data extracted from NOBV database
db_data_dir = file.path('~/Data/ECProcessing/met_data', site_id, 'processed')

# Path/name of DB to create/write to
db = file.path('data', paste0(site_id, '.db'))

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
if(strftime(end_dt, '%d%H%M') == '010000') {
  interval = interval(floor_date(start_dt, 'month')+1, end_dt-1) 
} else {
  interval = interval(floor_date(start_dt, 'month')+1, end_dt) 
}
months = floor_date(start_dt, 'month') %m+% months(0:(interval %/% months(1)))

# Bind list
fluxnet <- lapply(months, get_fluxnet) %>% 
  bind_rows()

# Read in monthly dfs of meteo, soil and water level data
db_data = lapply(months, get_meteo) %>% 
  bind_rows()

df = tibble(datetime = seq(min(c(fluxnet$datetime, db_data$datetime)),
                           max(c(fluxnet$datetime, db_data$datetime)),
                           1800)) %>% 
  left_join(fluxnet) %>% 
  left_join(db_data) %>% 
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

