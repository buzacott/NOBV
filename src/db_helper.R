#------------------------------------------------------------------------------#
# Functions for writing to db
#------------------------------------------------------------------------------#
get_fluxnet <- function(x) {
  src = file.path(eddypro_dir, year(x), format(x, '%m'))
  f = list.files(src, 'fluxnet')
  if(length(f) == 1) {
    df = read_csv(file.path(src, f)) %>% 
      rename(datetime = TIMESTAMP_END) %>% 
      select(datetime, any_of(default_fluxnet_variables),
             CO2, CH4, AIR_MV, BADM_INST_HEIGHT_SA) %>% 
      mutate(datetime = as_datetime(as.character(datetime),
                                    format = '%Y%m%d%H%M',
                                    tz = 'UTC')) %>% 
      filter(datetime >= start_dt, datetime <= end_dt) %>% 
      # NaN values
      mutate(across(where(is.numeric), ~na_if(., -9999))) %>% 
      # Recalculate storage fluxes due to daily processing (missing midnight
      # storage). Method from EddyPro code
      mutate(SC_SINGLE_2 = c(NaN, diff(CO2)) / AIR_MV *  BADM_INST_HEIGHT_SA / 1800,
             SCH4_SINGLE_2 = c(NaN, diff(CH4)) / AIR_MV *  BADM_INST_HEIGHT_SA / 1800) %>% 
      select(-c(CO2, CH4, AIR_MV, BADM_INST_HEIGHT_SA))
  } else {
    df = tibble()
  }
  return(df)
}

get_meteo <- function(x) {
  src = file.path(db_data_dir, year(x), format(x, '%m'))
  
  mt1_file <- file.path(src, paste0(strftime(x, '%Y-%m_MT1.csv')))
  if(file.exists(mt1_file)) {
    mt1 <- read_csv(mt1_file)
  } else {
    mt1 <- tibble(dt = POSIXct(tz = 'UTC'))
  }
  spg_file <- file.path(src, paste0(strftime(x, '%Y-%m_SPG.csv')))
  if(file.exists(spg_file)) {
    spg <- read_csv(spg_file)
  } else {
    spg <- tibble(dt = POSIXct(tz = 'UTC'))
  }
  wlg_file = file.path(src, paste0(strftime(x, '%Y-%m_WLG.csv')))
  if(file.exists(wlg_file)) {
    wlg <- read_csv(wlg_file)
  } else {
    wlg <- tibble(dt = POSIXct(tz = 'UTC'))
  }
  
  df <- reduce(list(mt1, spg, wlg), full_join) %>% 
    rename(datetime = dt) %>% 
    filter(datetime >= start_dt, datetime <= end_dt)
  
  return(df)
}

default_fluxnet_variables <- 
  c('H', 'LE', 'ET', 'FC', 'FCH4', 'SC_SINGLE', 'WS', 'WS_MAX', 'WD',
    'WD_SIGMA', 'USTAR', 'MO_LENGTH', 'T_SONIC', 'TA_EP', 'PA_EP', 'RH_EP',
    'VPD_EP', 'TDEW', 'U_SIGMA', 'V_SIGMA', 'W_SIGMA', 'FC_SSITC_TEST',
    'FCH4_SSITC_TEST', 'BADM_LOCATION_LAT', 'BADM_LOCATION_LONG',
    'BADM_INST_HEIGHT_SA',  'BADM_HEIGHTC', 'DISPLACEMENT_HEIGHT',
    'ROUGHNESS_LENGTH', 'ROUGHNESS_LENGTH', 'FETCH_MAX', 'FETCH_OFFSET',
    paste0('FETCH_', c(10, 30, 50, 70, 80, 90)), 'CUSTOM_RSSI_77_MEAN',
    'CUSTOM_CO2_SIGNAL_STRENGTH_7500_MEAN', 'ALBEDO_1_1_1', 'LW_IN_1_1_1',
    'LW_OUT_1_1_1', 'PTEMP_1_1_1', 'RH_1_1_1', 'RH_1_2_1', 'RH_1_3_1',
    'SW_IN_1_1_1', 'SW_OUT_1_1_1', 'TA_1_1_1', 'TA_1_2_1', 'TA_1_3_1',
    'TS_1_1_1', 'TS_1_2_1', 'VIN_1_1_1', 'WL_1_1_1', 'TW_1_1_1',
    'U_UNROT', 'V_UNROT', 'W_UNROT', 'DRYAIR_MV', 'AIR_MV',
    'CO2', 'CH4', 'CO2_MIXING_RATIO', 'CH4_MIXING_RATIO', 
    'U_VM97_TEST', 'V_VM97_TEST', 'W_VM97_TEST', 'CO2_VM97_TEST', 
    'CH4_VM97_TEST')