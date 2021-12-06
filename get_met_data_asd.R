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

site = 'ASD'
start_dt = lubridate::as_datetime('2021-08-01T00:00:01', tz=tz)
end_dt   = lubridate::as_datetime('2021-11-20T00:00:00', tz=tz)

# Number of months to download
interval = interval(start_dt, end_dt)
months = interval(start_dt, end_dt) %/% months(1)

start_months <- round_date(start_dt, 'month') %m+% months(0:months)
end_months <- start_months %m+% months(1) - 1
start_months[1] <- start_dt
end_months[length(end_months)] <- end_dt

# Dir to save csv files
data_dir <- '~/Data/Assendelft/Data'

for(i in seq_along(start_months)) {
  # Make dir to save files
  dst <- file.path(data_dir,
                   'raw',
                   year(start_months[i]),
                   strftime(start_months[i], '%m'))
  make_dir(dst)
  
  # Meteo data
  mt1 = get_db_data(site = site,
                    plot = 'RF',
                    sensor_group = 'MT1',
                    start_dt = start_months[i],
                    end_dt = end_months[i])
  # Soil data
  spg = get_db_data(site = site,
                    plot = 'MP',
                    sensor_group = 'SPG',
                    start_dt = start_months[i],
                    end_dt = end_months[i])
  # Water level data
  wlg = get_db_data(site = site,
                    plot = 'MP',
                    sensor_group = 'WLG',
                    start_dt = start_months[i],
                    end_dt = end_months[i])
  
  # Write data
  mt1 %>% 
    pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
    write_csv(file.path(dst, paste0(strftime(start_months[i], '%Y-%m_MT1.csv'))))
  
  spg %>% 
    pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
    write_csv(file.path(dst, paste0(strftime(start_months[i], '%Y-%m_SPG.csv'))))
  
  wlg %>% 
    pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
    write_csv(file.path(dst, paste0(strftime(start_months[i], '%Y-%m_WLG.csv'))))
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
  
  mt1[[i]] = read_csv(file.path(src, paste0(strftime(start_months[i], '%Y-%m_MT1.csv'))), guess_max = 100)
  spg[[i]] = read_csv(file.path(src, paste0(strftime(start_months[i], '%Y-%m_SPG.csv'))), guess_max = 100)
  wlg[[i]] = read_csv(file.path(src, paste0(strftime(start_months[i], '%Y-%m_WLG.csv'))), guess_max = 100) %>% 
    mutate(dt = as_datetime(dt, tz='UTC'))
} 
# Take 30 m average
mt1 <- bind_rows(mt1) %>%
  # Group by 30 mins
  mutate(dt = ceiling_date(dt, '30 min')) %>%
  group_by(dt) %>%
  # Aggregate variables
  summarise(
    # Better averaging of wind direction, use wind vectors
    u_east  = mean(MT1_WINS_1_H_200_Avg * sin(MT1_WIND_1_H_200_Avg / 180*pi), na.rm=TRUE),
    u_north = mean(MT1_WINS_1_H_200_Avg * cos(MT1_WIND_1_H_200_Avg / 180*pi), na.rm=TRUE),
    # Sum rain
    MT1_RAIN_1_H_040_Tot = sum(MT1_RAIN_1_H_040_Tot, na.rm=TRUE),
    # Take mean of variables except wind dir and rain
    across(-c(MT1_RAIN_1_H_040_Tot, MT1_WIND_1_H_200_Avg), mean, na.rm=TRUE),
  ) %>%
  mutate(MT1_WIND_1_H_200_Avg = (360 + (atan2(u_east, u_north) * 180/pi)) %% 360) %>%
  select(-c(u_east, u_north))


spg <- bind_rows(spg) %>%
  # Group by 30 mins
  mutate(dt = ceiling_date(dt, '30 min')) %>%
  group_by(dt) %>%
  # Aggregate variables
  summarise(across(.fns = mean, na.rm=TRUE))

# WLG data is hourly, do linear interp to half hourly
wlg <- bind_rows(wlg) %>% 
  right_join(tibble(dt = sort(c(.$dt, .$dt + 60*30)))) %>%
  arrange(dt) %>% 
  mutate(across(where(is.numeric), zoo::na.approx, na.rm=FALSE))

for(i in seq_along(start_months)) {
  dst <- file.path(data_dir,
                   'processed',
                   year(start_months[i]),
                   strftime(start_months[i], '%m'))
  make_dir(dst)
  
  # Write data
  mt1 %>%
    filter(floor_date(dt, 'month') == floor_date(start_months[i], 'month')) %>% 
    write_csv(file.path(dst, paste0(strftime(start_months[i], '%Y-%m_MT1.csv'))))
  
  spg %>%
    filter(floor_date(dt, 'month') == floor_date(start_months[i], 'month')) %>% 
    write_csv(file.path(dst, paste0(strftime(start_months[i], '%Y-%m_SPG.csv'))))
  
  wlg %>% 
    filter(floor_date(dt, 'month') == floor_date(start_months[i], 'month')) %>% 
    write_csv(file.path(dst, paste0(strftime(start_months[i], '%Y-%m_WLG.csv'))))
}
