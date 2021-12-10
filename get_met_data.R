# Script to download met, soil and water level data from the NOBV database and 
# then aggregate to 30 min timesteps.

library(tidyverse)
library(lubridate)

source('src/database_util.R')
source('src/db_config.R')

# Helper functions for processing EC data
make_dir <- function(path) {
  # Makes a directory if it doesn't exist
  if(dir.exists(path) == FALSE) dir.create(path, recursive = TRUE)
}

tz = 'UTC'
Sys.setenv(TZ=tz)

site = 'ZEG'

# Number of months to download
start_date = as_date('2021-04-01')
end_date   = as_date('2021-05-01')

interval = interval(start_date, end_date)
months = interval %/% months(1)

# + 1 / - 1 to get obs recorded within the month
start_months <- as_datetime((floor_date(start_date, 'month') %m+% months(0:months))) + 1
end_months <- start_months %m+% months(1) - 1

# Dir to save csv files
data_dir <- '~/Data/Zegveld/Data'

for(i in seq_along(start_months)) {
  # Make dir to save files
  dst <- file.path(data_dir,
                    'raw',
                    year(start_months[i]),
                    strftime(start_months[i], '%m'))
  make_dir(dst)
  
  # Meteo data
  mt1_file = file.path(dst, paste0(strftime(start_months[i], '%Y-%m_MT1.csv')))
  if(!file.exists(mt1_file)) {
    mt1 = get_db_data(site = site,
                      plot = 'RF',
                      sensor_group = 'MT1',
                      start_dt = start_months[i],
                      end_dt = end_months[i])
    mt1 %>% 
      pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
      write_csv(mt1_file)
  }

  # Soil data
  spg_file = file.path(dst, paste0(strftime(start_months[i], '%Y-%m_SPG.csv')))
  if(!file.exists(spg_file)) {
    spg = get_db_data(site = site,
                      plot = 'RF',
                      sensor_group = 'SPG',
                      start_dt = start_months[i],
                      end_dt = end_months[i])
    spg %>% 
      pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
      write_csv(spg_file)
  }
  
  # Water level data
  wlg_file = file.path(dst, paste0(strftime(start_months[i], '%Y-%m_WLG.csv')))
  if(!file.exists(wlg_file)) {
    wlg = get_db_data(site = site,
                      plot = 'MP16',
                      sensor_group = 'WLG',
                      start_dt = start_months[i],
                      end_dt = end_months[i])
    wlg %>% 
      pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
      write_csv(wlg_file)
  }
}

# Take 30 min averages
mt1 <- vector('list', length(start_months))
spg <- vector('list', length(start_months))
wlg <- vector('list', length(start_months))
for(i in seq_along(start_months)) {
  src <- file.path(data_dir,
                   'raw',
                   year(start_months[i]),
                   strftime(start_months[i], '%m'))
  
  mt1[[i]] = read_csv(file.path(src, paste0(strftime(start_months[i], '%Y-%m_MT1.csv'))))
  spg[[i]] = read_csv(file.path(src, paste0(strftime(start_months[i], '%Y-%m_SPG.csv'))))
  wlg[[i]] = read_csv(file.path(src, paste0(strftime(start_months[i], '%Y-%m_WLG.csv')))) %>% 
    mutate(dt = as_datetime(dt, tz='UTC'))
} 


# Take 30 min average
mt1 <- lapply(mt1, function(df) {
  df %>% 
    mutate(dt = ceiling_date(dt, '30 min')) %>%
    group_by(dt) %>%
    # Aggregate variables
    summarise(
      # Better averaging of wind direction, use wind vectors
      u_east  = mean(MT1_WINS_1_H_200_Avg * sin(MT1_WIND_1_H_200_Avg / 180*pi), na.rm=TRUE),
      u_north = mean(MT1_WINS_1_H_200_Avg * cos(MT1_WIND_1_H_200_Avg / 180*pi), na.rm=TRUE),
      # Take mean of variables except wind dir and rain
      across(-c(MT1_RAIN_1_H_040_Tot, MT1_WIND_1_H_200_Avg), mean, na.rm=TRUE),
      # Sum rain
      MT1_RAIN_1_H_040_Tot = sum(MT1_RAIN_1_H_040_Tot, na.rm=TRUE),
    ) %>%
    mutate(MT1_WIND_1_H_200_Avg = (atan2(u_east, u_north) + pi) * 180 / pi) %>%
    select(-c(u_east, u_north))
})


# SPG data already half-hourly

# WLG data is hourly, do linear interp to half hourly
# Remove empty months
wlg <- wlg[ sapply(wlg, function(x) nrow(x)>1) ]

wlg <- bind_rows(wlg) %>% 
  right_join(tibble(dt = seq(min(.$dt), max(.$dt), 1800))) %>% 
  arrange(dt) %>% 
  mutate(across(where(is.numeric), zoo::na.approx, na.rm=FALSE),
         month = floor_date(dt-1, 'month')) %>% 
  group_split(month, .keep=FALSE)


lapply(mt1, function(df) {
  dt = df$dt[1]
  dst <- file.path(data_dir,
                   'processed',
                   year(dt),
                   strftime(dt, '%m'))
  make_dir(dst)
  
  write_csv(df,  file.path(dst, paste0(strftime(dt, '%Y-%m_MT1.csv'))))
})

lapply(spg, function(df) {
  dt = df$dt[1]
  dst <- file.path(data_dir,
                   'processed',
                   year(dt),
                   strftime(dt, '%m'))
  make_dir(dst)
  
  write_csv(df,  file.path(dst, paste0(strftime(dt, '%Y-%m_SPG.csv'))))
})

lapply(wlg, function(df) {
  dt = df$dt[1]
  dst <- file.path(data_dir,
                   'processed',
                   year(dt),
                   strftime(dt, '%m'))
  make_dir(dst)
  
  write_csv(df,  file.path(dst, paste0(strftime(dt, '%Y-%m_WLG.csv'))))
})