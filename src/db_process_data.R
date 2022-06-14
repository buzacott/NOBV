download_met_data <- function() {
  for(i in seq_along(dates)) {
    # Make dir to save files
    dst <- file.path(data_dir,
                     'raw',
                     format(dates[i], '%Y'),
                     format(dates[i], '%m'))
    make_dir(dst)
    
    # Meteo data
    mt1_file = file.path(dst, paste0(strftime(dates[i], '%Y-%m-%d_MT1.csv')))
    if(!file.exists(mt1_file)) {
      mt1 = get_db_data(site = site,
                        plot = plot_mt1,
                        sensor_group = 'MT1',
                        start_dt = as_datetime(dates[i])+1,
                        end_dt = as_datetime(dates[i]+1))
      if(nrow(mt1) > 0) {
        mt1 %>% 
          pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
          arrange(dt) %>% 
          write_csv(mt1_file)
      }
    }
    
    # Soil data
    spg_file = file.path(dst, paste0(strftime(dates[i], '%Y-%m-%d_SPG.csv')))
    if(!file.exists(spg_file)) {
      spg = get_db_data(site = site,
                        plot = plot_spg,
                        sensor_group = 'SPG',
                        start_dt = as_datetime(dates[i])+1,
                        end_dt = as_datetime(dates[i]+1))
      if(nrow(spg) > 0) {
        spg %>% 
          pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
          arrange(dt) %>% 
          write_csv(spg_file)
      }
    }
    
    # Water level data
    wlg_file = file.path(dst, paste0(strftime(dates[i], '%Y-%m-%d_WLG.csv')))
    if(!file.exists(wlg_file)) {
      wlg = get_db_data(site = site,
                        plot = plot_wlg,
                        sensor_group = 'WLG',
                        start_dt = as_datetime(dates[i])+1,
                        end_dt = as_datetime(dates[i]+1))
      
      if(nrow(wlg) > 0) {
        wlg %>% 
          pivot_wider(id_cols=dt, names_from=name, values_from=value) %>% 
          arrange(dt) %>% 
          write_csv(wlg_file)
      }
    }
  }
}

process_met_data <- function(overwrite = FALSE) {
  # Averages the met data to 30m blocks
  # Take 30 min averages
  yearmon <- unique(floor_date(dates, 'month'))
  
  mt1 <- vector('list', length(yearmon))
  spg <- vector('list', length(yearmon))
  wlg <- vector('list', length(yearmon))
  
  for(i in seq_along(yearmon)) {
    dst <- file.path(data_dir,
                     'processed',
                     format(yearmon[i], '%Y'),
                     format(yearmon[i], '%m'))
    
    if(dir.exists(dst) == FALSE | overwrite == TRUE) {
      src <- file.path(data_dir,
                       'raw',
                       year(yearmon[i]),
                       strftime(yearmon[i], '%m'))
      
      mt1_files <- list.files(src, pattern='MT1', recursive = TRUE)
      spg_files <- list.files(src, pattern='SPG', recursive = TRUE)
      wlg_files <- list.files(src, pattern='WLG', recursive = TRUE)
      
      mt1[[i]] = lapply(mt1_files, function(f) {
        read_csv(file.path(src, f),
                 col_types = cols(.default = "d", dt = "T"))
      }) %>% 
        bind_rows()
      spg[[i]] = lapply(spg_files, function(f) {
        read_csv(file.path(src, f),
                 col_types = cols(.default = "d", dt = "T"))
      }) %>% 
        bind_rows()
      wlg[[i]] = lapply(wlg_files, function(f) {
        read_csv(file.path(src, f),
                 col_types = cols(.default = "d", dt = "T"))
      }) %>% 
        bind_rows()
    }
  }
  
  # Take 30 min average
  if(overwrite == TRUE) {
    mt1 <- lapply(mt1, function(df) {
      df %>% 
        mutate(dt = ceiling_date(dt, '30 min')) %>%
        group_by(dt) %>%
        # Aggregate variables
        summarise(
          # Better averaging of wind direction, use wind vectors
          u_east  = -mean((MT1_WINS_1_H_200_Avg+0.001) * sin(MT1_WIND_1_H_200_Avg * pi/180), na.rm=TRUE),
          u_north = -mean((MT1_WINS_1_H_200_Avg+0.001) * cos(MT1_WIND_1_H_200_Avg * pi/180), na.rm=TRUE),
          # Take mean of variables except wind dir and rain
          across(-c(MT1_RAIN_1_H_040_Tot, MT1_WIND_1_H_200_Avg), mean, na.rm=TRUE),
          # Arithmetic mean wind
          MT1_WINS_1_H_200_Avg_am = mean(MT1_WINS_1_H_200_Avg, na.rm=TRUE),
          # Sum rain
          MT1_RAIN_1_H_040_Tot = sum(MT1_RAIN_1_H_040_Tot, na.rm=TRUE),
        ) %>%
        mutate(MT1_WIND_1_H_200_Avg = (atan2(u_east, u_north) + pi) * 180 / pi) %>%
        select(-c(u_east, u_north))
    })
    
    lapply(mt1, function(df) {
      dt = df$dt[1]
      dst <- file.path(data_dir,
                       'processed',
                       year(dt),
                       strftime(dt, '%m'))
      make_dir(dst)
      
      write_csv(df,  file.path(dst, paste0(strftime(dt, '%Y-%m_MT1.csv'))))
    })
    
    

    
    # WLG data is hourly, do linear interp to half hourly
    # Remove empty months
    wlg <- wlg[ sapply(wlg, function(x) nrow(x)>1) ]
    if(length(wlg) > 0) {
      wlg <- bind_rows(wlg) %>% 
        right_join(tibble(dt = seq(min(.$dt), max(.$dt), 1800))) %>% 
        arrange(dt) %>% 
        mutate(across(where(is.numeric), zoo::na.approx, na.rm=FALSE),
               month = floor_date(dt-1, 'month')) %>% 
        group_split(month, .keep=FALSE)
      
      lapply(wlg, function(df) {
        dt = df$dt[1]
        dst <- file.path(data_dir,
                         'processed',
                         year(dt),
                         strftime(dt, '%m'))
        make_dir(dst)
        
        write_csv(df,  file.path(dst, paste0(strftime(dt, '%Y-%m_WLG.csv'))))
      })
    }
    
    # SPG data already half-hourly
    lapply(spg, function(df) {
      dt = df$dt[1]
      dst <- file.path(data_dir,
                       'processed',
                       year(dt),
                       strftime(dt, '%m'))
      make_dir(dst)
      
      write_csv(df,  file.path(dst, paste0(strftime(dt, '%Y-%m_SPG.csv'))))
    })
    

  }
}
