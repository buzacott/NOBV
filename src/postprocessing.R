#------------------------------------------------------------------------------#
# Scripts to postprocess NOBV eddy covariance data
# R Code written by Alexander Buzacott, with adaptations from hidmu's 
# Python code
#------------------------------------------------------------------------------#

library(dplyr)
library(dbplyr)
library(readr)
library(tidyr)
library(DBI)
library(ggplot2)
library(raster, exclude = 'select')
library(lubridate)
library(randomForest)
library(sf)
library(doParallel)
library(REddyProc)
# Import Python libraries via reticulate for neural network gapfilling
# Comment out or ignore if you don't want to setup reticulate or use neural
# network gapfilling.
library(reticulate)
#conda_create(envname = 'RStudio',
#             packages = c('python', 'scikit-learn', 'pandas', 'numpy'))
reticulate::use_condaenv(condaenv = 'RStudio')
sklearn <- reticulate::import('sklearn')
model_selection <- reticulate::import('sklearn.model_selection')
preproc <- reticulate::import("sklearn.preprocessing")
pipeline <- reticulate::import("sklearn.pipeline")
neural_network <- reticulate::import("sklearn.neural_network")
pd <- reticulate::import('pandas')

Sys.setenv(TZ='UTC')

set.seed(123)

# Import functions
source('src/postprocessing_functions.R')
# Import Kljun footprint model 2015
source('src/calc_footprint_FFP_climatology.R')

#------------------------------------------------------------------------------#
# Postprocessing class
#------------------------------------------------------------------------------#
Postprocessing <- setRefClass(
  "Postprocessing",
  fields = list(
    site_id = 'character',
    default_variables = "character",
    db_location = "character",
    data = "ANY",
    EProc = "ANY",
    ch4_rf = "ANY",
    latlon = "numeric",
    ustar_threshold = "numeric"
  ),
  methods = list(
    initialize = function(site_id, latlon)
    {
      site_id <<- site_id
      latlon <<- latlon
      default_variables <<- c(
        'FC','SC_SINGLE', 'FCH4', 'USTAR', 'DISPLACEMENT_HEIGHT', 
        'ROUGHNESS_LENGTH', 'MT1_ATMP_1_H_200_Avg', 'WS', 'WD', 'MO_LENGTH',
        'V_SIGMA', 'MT1_SWIN_1_H_180','FC_SSITC_TEST', 'FCH4_SSITC_TEST',
        'MT1_PAR_1_H_180', 'MT1_WIND_1_H_200_Avg', 'MT1_WINS_1_H_200_Avg',
        'MT1_NETL_1_H_180', 'MT1_NETS_1_H_180', 'SPG_STMP_1_D_005_Avg',
        'MT1_RHUM_1_H_200_Avg', 'CUSTOM_RSSI_77_MEAN',
        'CUSTOM_CO2_SIGNAL_STRENGTH_7500_MEAN')
      
      db_location <<- file.path('data', 'sample.db')
    },
    load_data = function(start_dt, end_dt, db = NULL, variables = NULL)
    {
      # Check to see if there is a given list of variables, otherwise use default
      # variables
      if(is.null(variables)) {
        variables <- default_variables
      }
      if(is.null(db)) {
        db <- db_location
      }

      con <- DBI::dbConnect(RSQLite::SQLite(), db)
      tryCatch({
        id_name <- tbl(con, 'id_name')
        df <- tbl(con, 'data') %>% 
          left_join(id_name) %>% 
          filter(name %in% !!variables,
                 datetime >= !!as.integer(start_dt),
                 datetime <= !!as.integer(end_dt)) %>% 
          select(datetime, name, value) %>% 
          collect() %>% 
          # Pivot to wide format
          pivot_wider(id_col=datetime, names_from=name, values_from=value)
      },
      error = function(e) message(e),
      finally = DBI::dbDisconnect(con))
      
      col_name_change <- c(
        'ATMP' = 'MT1_ATMP_1_H_200_Avg',
        'SWIN' = 'MT1_SWIN_1_H_180',
        'FC_flag' = 'FC_SSITC_TEST',
        'FCH4_flag' = 'FCH4_SSITC_TEST',
        'PAR_met' = 'MT1_PAR_1_H_180' ,
        'WD_met' = 'MT1_WIND_1_H_200_Avg',
        'WS_met' = 'MT1_WINS_1_H_200_Avg',
        'NETL_met' = 'MT1_NETL_1_H_180',
        'NETS_met' = 'MT1_NETS_1_H_180',
        'STMP' = 'SPG_STMP_1_D_005_Avg',
        'RELH' = 'MT1_RHUM_1_H_200_Avg',
        'CO2_SS' = 'CUSTOM_CO2_SIGNAL_STRENGTH_7500_MEAN',
        'CH4_SS' = 'CUSTOM_RSSI_77_MEAN'
      )
      
      # Rename columns
      df <- df %>% 
        rename(any_of(col_name_change))
      
      # Some adjustments
      df <- df %>% 
        # Add NEE, night time and month columns
        mutate(datetime = as_datetime(datetime, tz='UTC'),
               date = floor_date(datetime, 'day'), # temporary for sr/ss
               NEE = FC + SC_SINGLE,
               FCH4 = FCH4 / 1000, # nmol/s/m2 -> µmol/s/m2
               year = year(datetime - 1),
               month = month(datetime - 1),
               season = quarter(datetime - 1),
               month_sin = sin(pi * (month - 1) / 12),
               day_sin = cos(2 * pi * (yday(datetime)-172)/365.25)) %>%  
        arrange(datetime)
      
      # Use astronomical day/night to define night, r_g thresholds used later
      crds <- matrix(rev(latlon), nrow=1)
      crds <- sp::SpatialPoints(crds,
                                proj4string=CRS("+proj=longlat +datum=WGS84"))
      
      sr_ss_df <- tibble(date = unique(df$date)) %>% 
        mutate(sunrise = maptools::sunriset(crds = crds,
                                            dateTime = date,
                                            direction = c('sunrise'),
                                            POSIXct.out = TRUE)[,2],
               sunset = maptools::sunriset(crds = crds,
                                           dateTime = date,
                                           direction = c('sunset'),
                                           POSIXct.out = TRUE)[,2])
      
      df <- df %>% 
        left_join(sr_ss_df) %>% 
        mutate(night = datetime < sunrise | datetime > sunset) %>% 
        select(-c(date, sunrise, sunset))
      
      # Assign to field
      data <<- df
    },
    filter_co2 = function(vm97 = FALSE, CO2_SS_threshold = 80)
    {
      # Filter the CO2 data using EddyPro quality flags
      # NEE is also updated at the end of the function
      
      print('Filtering FC data')
      n_before = sum(!is.na(data$FC))
      print(paste('Data points before filtering:', n_before))
      
      # Quality flag filtering, possible to also add 1 if you want to be strict
      # Flags
      # 0 is good
      # 1 is ok
      # 2 is bad
      # Can customise flagging system in EddyPro
      data <<- data %>% 
        mutate(FC_f = if_else(FC_flag > 1, NaN, FC))
      
      # Apply Vickers and Mahrt 1997 tests
      if(vm97) {
        vm97_colnames <- c('spikes_hf',
                           'amplitude_hf',
                           'drop_out_hf',
                           'absolute_limits_hf',
                           'skewness_kurtosis_hf',
                           'skewness_kurtosis_sf',
                           'discontinuities_hf',
                           'discontinuities_sf')
        # Missing value index
        na_idx <- is.na(data$CO2_VM97_TEST)
        co2_vm97_test <- matrix(nrow=length(na_idx), ncol=8)
        co2_vm97_test[!na_idx,] <- t(sep_digits(data$CO2_VM97_TEST[!na_idx]))[,-1]
        # Select only hardflags for now
        # TODO: add here soft/hard flag choosing
        co2_vm97_test <- co2_vm97_test[,c(-6, -8)]
        data <<- data %>% 
          mutate(FC_VM97_flag = rowSums(co2_vm97_test),
                 FC_f = if_else(FC_VM97_flag > 1, NaN, FC_f))
      }
      
      # Signal strength filter
      data <<- data %>% 
        mutate(FC_f = if_else(CO2_SS < CO2_SS_threshold, NaN, FC_f))
      
      # Update NEE
      data <<- data %>% 
        mutate(NEE_f = if_else(is.na(FC_f), NaN, NEE))
      
      # Print n data points lost
      n_after = sum(!is.na(data$FC_f))
      print(paste('Data points after filtering:', n_after))
      print(paste('Data points lost:', n_before-n_after,
                  paste0('(', (1 - round(n_after/n_before, 2)) * 100, '%)')))
    },
    filter_ch4 = function(vm97 = FALSE, CH4_SS_threshold)
    {
      # Filter methane data for quality flags and outliers
      print('Filtering FCH4 data')
      n_before = sum(!is.na(data$FCH4))
      print(paste('Data points before filtering:', n_before))
      
      # Quality flags filter
      data <<- data %>% 
        mutate(FCH4_f = if_else(FCH4_flag > 1, NaN, FCH4))
      
      # Apply Vickers and Mahrt 1997 tests
      if(vm97) {
        vm97_colnames <- c('spikes_hf',
                           'amplitude_hf',
                           'drop_out_hf',
                           'absolute_limits_hf',
                           'skewness_kurtosis_hf',
                           'skewness_kurtosis_sf',
                           'discontinuities_hf',
                           'discontinuities_sf')
        
        # Missing value index
        na_idx <- is.na(data$CH4_VM97_TEST)
        ch4_vm97_test <- matrix(nrow=length(na_idx), ncol=8)
        ch4_vm97_test[!na_idx,] <- t(sep_digits(data$CH4_VM97_TEST[!na_idx]))[,-1]
        # Select only hardflags for now
        # TODO: add here soft/hard flag choosing
        ch4_vm97_test <- ch4_vm97_test[,c(-6, -8)]
        
        data <<- data %>% 
          mutate(FCH4_VM97_flag = rowSums(ch4_vm97_test),
                 FCH4_f = if_else(FCH4_VM97_flag > 0, NaN, FCH4_f))
        
        # Signal strength filter
        data <<- data %>% 
          mutate(FCH4_f = if_else(CH4_SS < CH4_SS_threshold, NaN, FCH4_f))
      }
      
      n_after = sum(!is.na(data$FCH4_f))
      print(paste('Data points after filtering:', n_after))
      print(paste('Data points lost:', n_before-n_after,
                  paste0('(', (1 - round(n_after/n_before, 2)) * 100, '%)')))
    },
    filter_rain = function() {
      # Removes timesteps before, and after rain
      print('Filtering timesteps with, before and after rain')
      n_fc_before = sum(!is.na(data$FC_f))
      n_fch4_before = sum(!is.na(data$FCH4_f))
      print(paste('n FC before:', n_fc_before))
      print(paste('n FCH4 before:', n_fch4_before))
      
      n_fc <- sum(is.na(data$FC_f))
      n_fch4 <- sum(is.na(data$FCH4_f))
        
      rain_idx <- which(data$MT1_RAIN_1_H_040_Tot > 0)
      n <- length(data$MT1_RAIN_1_H_040_Tot)
      rain_idx <- unique(c(rain_idx, rain_idx-1, rain_idx+1))
      rain_idx <- 1:n %in% rain_idx
      
      data <<- data %>% 
        mutate(FC_f = if_else(rain_idx, NaN, FC_f),
               NEE_f = if_else(rain_idx, NaN, NEE_f),
               FCH4_f = if_else(rain_idx, NaN, FCH4_f))
      
      n_fc_after = sum(!is.na(data$FC_f))
      n_fch4_after = sum(!is.na(data$FCH4_f))
      print(paste('n FC after', n_fc_after, '(% lost:',
                  paste0((1 - round(n_fc_after/n_fc_before, 2)) * 100, '%)')))
      print(paste('n FCH4 after', n_fch4_after, '(% lost:',
                  paste0((1 - round(n_fch4_after/n_fch4_before, 2)) * 100, '%)')))
    },
    filter_outlier_co2 = function(lower = 0.01, upper = 0.99) {
      # Quantile filtering of NEE outliers if desired
      data <<- data %>% 
        mutate(across(c(FC_f, NEE_f),
                      ~if_else(.x <= quantile(.x, lower, na.rm=TRUE) |
                                 .x >= quantile(.x, upper, na.rm=TRUE),
                               NaN, .x)))
    },
    filter_outlier_ch4 = function(lower = 0.01, upper = 0.99) {
      # Quantile filtering of CH4 outliers if desired
      data <<- data %>% 
        mutate(FCH4_f = if_else(FCH4_f <= quantile(FCH4_f, lower, na.rm=TRUE) |
                                  FCH4_f >= quantile(FCH4_f, upper, na.rm=TRUE),
                                NaN, FCH4_f))
    },
    filter_ustar = function(min_ustar_threshold = 0.1,
                            force_min_threshold = FALSE,
                            manual_threshold = NULL) {
      #------------------------------------------------------------------------#
      # uStar filtering using REddyProc
      #------------------------------------------------------------------------#
      # Create data for REddyProc
      eproc_data <- data %>% 
        select(DateTime = datetime,
               NEE = NEE_f,
               Rg =  SWIN,
               Tair = ATMP,
               rH = RELH,
               Ustar = USTAR,
               VPD) %>% 
        # Add VPD
        mutate(Rg = if_else(Rg < 0, 0, Rg)) %>% 
        as.data.frame()
      
      # Create REddyProc class for ustar filtering
      EProc <<- sEddyProc$new(
        site_id,
        eproc_data,
        c('NEE', 'Rg', 'Tair', 'rH', 'Ustar', 'VPD')
      )

      EProc$sSetLocationInfo(LatDeg = latlon[1],
                             LongDeg = latlon[2],
                             TimeZoneHour = 1)
      
      # Ustar filtering, use MPT from REddyProc
      # Estimate threshold with 100 bootstrap iterations
      if(is.null(manual_threshold)) {
        EProc$sEstimateUstarScenarios(nSample=100L, isVerbose=FALSE)
        
        # Use single ustar threshold
        ustar_threshold <<- EProc$sUSTAR$uStar[1]
        message(paste('MPT u* threshold is', ustar_threshold))
        if(ustar_threshold < min_ustar_threshold) {
          message(paste('Ustar threshold is below the minimum, defaulting
                      to minimum threshold of', min_ustar_threshold))
          ustar_threshold <<- min_ustar_threshold
        }
      } else {
        suppressWarnings({
          EProc$sUSTAR <<- EProc$sEstUstarThold()
          
        })
        EProc$sUSTAR$uStar[1] <<- manual_threshold
        ustar_threshold <<- manual_threshold
      }
      
      # Force threshold
      if(force_min_threshold) {
        # Use single ustar threshold
        EProc$sUSTAR$uStar[1] <<- min_ustar_threshold
        ustar_threshold <<- min_ustar_threshold
      }
      
      data <<- data %>% 
        mutate(across(c(FC_f, NEE_f, FCH4_f),
                      ~if_else(USTAR < ustar_threshold, NaN, .x)))
      
    },
    flux_partitioning = function(t_ref = 15,
                                 t_0 = -46.02,
                                 rg_threshold = 20,
                                 annual = TRUE,
                                 use_gf_vars = FALSE,
                                 plot = FALSE) {
      #------------------------------------------------------------------------#
      # Respiration
      #------------------------------------------------------------------------#
      # Method follows Reichstein et al. 2005
      
      # Subset variables and add date
      df <- data %>% 
        select(year, datetime, night, ATMP, SWIN, NEE=NEE_f) %>% 
        mutate(date = floor_date(datetime, unit='day'))
      
      if(use_gf_vars) {
        df <- data %>% 
          select(year, datetime, night, ATMP = ATMP_gf, SWIN = SWIN_gf, NEE=NEE_f) %>% 
          mutate(date = floor_date(datetime, unit='day'))
      }
      
      # Calculate Reco per year)
      if((length(unique(df$date)) > 365) & annual) {
        # First see how many days of data are in each year
        # If there are too few it will often fail to find E0
        ndays <- df %>%
          select(year, datetime) %>% 
          mutate(datetime = floor_date(datetime, 'day')) %>% 
          distinct() %>% 
          group_by(year) %>% 
          tally()
        
        # Assign year groups
        ndays$grp <- 1
        nyears <- nrow(ndays)
        for(i in seq_len(nyears)) {
          if(i > 1) {
            if(ndays$n[i] > 182) {
              ndays$grp[i:nyears] <- ndays$grp[i:nyears] + 1
            }
          }
        }
        
        df <- left_join(df, ndays %>% select(-n))
        
        # Calculate respiration in groups
        suppressWarnings({
          res <- df %>% 
            group_by(grp) %>% 
            nest() %>% 
            mutate(data = purrr::map(data, function(x) {
              fit_respiration_curve(x,
                                    t_ref = t_ref,
                                    t_0 = t_0,
                                    rg_threshold = rg_threshold)
            })) %>% 
            unnest(data) %>% 
            ungroup() %>% 
            select(-grp)
        })
      } else {
        # Try and fit once
        res <- fit_respiration_curve(df,
                                     t_ref = t_ref,
                                     t_0 = t_0,
                                     rg_threshold = rg_threshold)
      }
      
      data <<- left_join(data, res)
      
      #------------------------------------------------------------------------#
      # GPP
      #------------------------------------------------------------------------#
      # Select data for fitting
      if(use_gf_vars) {
        gpp_data = data %>% 
          select(datetime, year, month, season, SWIN=SWIN_gf, NEE_f, Reco)
      } else {
        gpp_data = data %>% 
          select(datetime, year, month, season, SWIN=SWIN, NEE_f, Reco)
      }
      gpp_data = gpp_data %>% 
        mutate(month_group = ceiling(month / 2)) %>% 
        mutate(GPP = NEE_f - Reco,
               SWIN = if_else(SWIN > 0, SWIN, 0))
      
      variables = gpp_data %>% 
        drop_na() %>% 
        group_by(year, season) %>% 
        # group_by(year, month_group) %>% 
        nest() %>% 
        mutate(variables = purrr::map(data, function(x) {
          variables = fit_gpp_curve(x)
          return(tibble(beta1 = variables[[1]],
                        beta2 = variables[[2]]))
        })) %>% 
        select(-data) %>% 
        unnest(variables)
      
      gpp_data <- gpp_data %>% 
        left_join(variables) %>% 
        mutate(across(c(beta1, beta2), ~replace_na(.x, mean(.x))),
               GPP = gpp(SWIN, beta1, beta2))
        
      
      data <<- data %>% 
        left_join(gpp_data %>% select(datetime, beta1, beta2, GPP))
      
      if(plot) {
        data %>% 
          select(datetime, NEE_f, GPP, Reco) %>% 
          pivot_longer(-datetime) %>% 
          ggplot(aes(datetime, value, col=name)) +
          geom_line() +
          labs(y = expression(value~'('*µmol.s^-1*.m^-2*')'))
      }
    },
    flux_partitioning_ll = function(t_ref = 15,
                                    t_0 = -46.02,
                                    rg_threshold_resp = 20,
                                    rg_threshold_gpp = 4,
                                    annual = TRUE,
                                    use_gf_vars = FALSE,
                                    plot = FALSE) {
      #------------------------------------------------------------------------#
      # Method follows Lasslop et al. 2010
      #------------------------------------------------------------------------#
      data <<- data %>% 
        select(-any_of(c('GPP_ll', 'Reco_ll', 'NEE_ll')))
      
      # Subset variables and add date
      df <- data %>% 
        select(year, datetime, night, ATMP, SWIN, NEE=NEE_f) %>% 
        mutate(date = floor_date(datetime, unit='day'))
      
      if(use_gf_vars) {
        df <- data %>% 
          select(year, datetime, night, ATMP=ATMP_gf, SWIN=SWIN_gf, NEE=NEE_f) %>% 
          mutate(date = floor_date(datetime, unit='day'))
      }
      
      if((length(unique(df$date)) > 365) & annual) {
        # First see how many days of data are in each year
        # If there are too few it will often fail to find E0
        ndays <- df %>%
          select(year, datetime) %>% 
          mutate(datetime = floor_date(datetime, 'day')) %>% 
          distinct() %>% 
          group_by(year) %>% 
          tally()
        
        # Assign year groups
        ndays$grp <- 1
        nyears <- nrow(ndays)
        for(i in seq_len(nyears)) {
          if(i > 1) {
            if(ndays$n[i] > 182) {
              ndays$grp[i:nyears] <- ndays$grp[i:nyears] + 1
            }
          }
        }
        
        df <- left_join(df, ndays %>% select(-n))
      
        suppressWarnings({
          gpp_ll <- df %>% 
            group_by(grp) %>% 
            nest() %>% 
            mutate(data = purrr::map(data, function(x) {
              gpp_lasslop_2010(x,
                               t_ref = t_ref,
                               t_0 = t_0,
                               rg_threshold_resp = 20,
                               rg_threshold_gpp = 4)
            })) %>% 
            unnest(data) %>% 
            ungroup() %>% 
            select(-grp)
        })
      } else {
        gpp_ll <- gpp_lasslop_2010(df,
                                   t_ref = t_ref,
                                   t_0 = t_0,
                                   rg_threshold_resp = 20,
                                   rg_threshold_gpp = 4)
      }
      
      gpp_ll <- gpp_ll %>% 
        select(datetime,
               alpha_ll = alpha,
               beta_ll = beta,
               r_ref_ll = r_ref,
               e_0_ll = e_0,
               GPP_ll,
               Reco_ll,
               NEE_ll)
      
      data <<- left_join(data, gpp_ll)
      
      if(plot) {
        data %>% 
          select(datetime, NEE_ll, GPP_ll, Reco_ll) %>% 
          pivot_longer(-datetime) %>% 
          ggplot(aes(datetime, value, col=name)) +
          geom_line() +
          labs(y = expression(value~'('*µmol.s^-1*.m^-2*')'))
      }
    },
    flux_partitioning_eproc = function(plot = FALSE) 
    {
      #------------------------------------------------------------------------#
      # Using REddyProc to gapfill NEE
      #------------------------------------------------------------------------#
      EProc$sMDSGapFillAfterUstar('NEE', FillAll=TRUE, uStarTh = ustar_threshold)
      
      # Gapfill temperature and VPD data for flux partitioning
      EProc$sMDSGapFill('Tair', FillAll = FALSE,  minNWarnRunLength = NA)
      EProc$sMDSGapFill('VPD', FillAll = FALSE,  minNWarnRunLength = NA)
      # Fill longer gaps still present in VPD_f
      EProc$sFillVPDFromDew()
      # Partition fluxes
      EProc$sMRFluxPartition(suffix = 'uStar')
      #EProc$sMRFluxPartitionUStarScens()
      eproc_results <- EProc$sExportResults()
      
      data <<- data %>% 
        mutate(NEE_eproc = eproc_results$NEE_uStar_f,
               GPP_eproc = eproc_results$GPP_uStar_f,
               Reco_eproc = eproc_results$Reco_uStar)
    },
    gapfilling_nlr = function(plot = FALSE)
    {
      # Non-linear regression approach to gapfilling
      data <<- data %>% 
        mutate(nlr_NEE = Reco + GPP,
               NEE_GF_NLR = if_else(is.na(NEE_f), nlr_NEE, NEE_f))
    },
    gapfilling_nn = function(training_variables,
                             plot = FALSE)
    {
      # Gapfill NEE using a neural network
      if('nn_NEE' %in% colnames(data)) {
        data <<- data %>% 
          select(-nn_NEE)
      }
      
      if(missing(training_variables)) {
        training_variables <- c('ATMP', 'SWIN', 'VPD', 'day_sin')
      }
      
      training_data = data %>%
        select(datetime, NEE_f, all_of(training_variables)) %>% 
        drop_na()
      
      x = training_data %>%
        select(all_of(training_variables))
      y = training_data %>%
        pull(NEE_f)
      
      train_test = model_selection$train_test_split(x, y)
      
      pipe = sklearn$pipeline$Pipeline(list(
        tuple('scaler', preproc$StandardScaler()),
        tuple('nn', neural_network$MLPRegressor(
          hidden_layer_sizes = tuple(5L,5L,5L,5L,5L),
          activation='logistic',
          learning_rate_init = 0.05,
          max_iter = 1000L))
      ))
      
      pipe$fit(train_test[[1]], train_test[[3]])
      
      score = pipe$score(train_test[[2]], train_test[[4]])
      print(paste('The accuracy is', round(score, 2)))
      
      # Simulate NEE
      pred_data = data %>% 
        select(datetime, all_of(training_variables)) %>% 
        drop_na() 
      pred_data$nn_NEE = pipe$predict(pred_data %>% select(-datetime))
      
      # Merge to data
      data <<- data %>% 
        left_join(pred_data %>% select(datetime, nn_NEE)) %>% 
        mutate(NEE_GF_NN = if_else(is.na(NEE_f), nn_NEE, NEE_f))
      if(plot == TRUE) {
        ggplot(data, aes(datetime, NEE_GF_NN)) +
          geom_line()
      }
    },
    gapfill_ch4_nn = function(training_variables,
                              plot = FALSE)
    {
      #------------------------------------------------------------------------#
      # Gapfill methane using a neural network
      #------------------------------------------------------------------------#
      if('FCH4_NN' %in% colnames(data)) {
        data <<- data %>% 
          select(-FCH4_NN)
      }
      
      if(missing(training_variables)) {
        training_variables = c('ATMP', 'SWIN', 'PAR_met', 'day_sin', 'WS', 'WD')
      }
      
      training_data = data %>%
        select(datetime, FCH4_f, all_of(training_variables)) %>%
        drop_na()

      x = training_data %>%
        select(all_of(training_variables))
      y = training_data %>%
        pull(FCH4_f)
      
      train_test = model_selection$train_test_split(x, y)
      
      # TODO: evaluate NN design
      pipe = sklearn$pipeline$Pipeline(list(
        tuple('scaler', preproc$StandardScaler()),
        tuple('nn', neural_network$MLPRegressor(
          hidden_layer_sizes = tuple(10L, 10L, 10L),
          activation='logistic',
          solver = 'lbfgs',
          learning_rate_init = 0.01,
          max_iter = 15000L))
      ))
      
      pipe$fit(train_test[[1]], train_test[[3]])
      
      score = pipe$score(train_test[[2]], train_test[[4]])
      print(paste('The accuracy is', round(score, 2)))
      
      # Simulate FCH4
      pred_data = data %>% 
        select(datetime, all_of(training_variables)) %>% 
        drop_na() 
      pred_data$FCH4_NN = pipe$predict(pred_data %>% select(-datetime))
      
      # Merge to data
      data <<- left_join(data,
                         pred_data %>% select(datetime, FCH4_NN)) %>% 
        mutate(FCH4_GF_NN = if_else(is.na(FCH4_f), FCH4_NN, FCH4_f))
      
      if(plot == TRUE) {
        ggplot(data, aes(datetime, FCH4_GF_NN)) +
          geom_line()
      }
    },
    gapfill_ch4_rf = function(training_variables,
                              plot = FALSE)
    {
      #------------------------------------------------------------------------#
      # Gapfill methane using Random Forest
      #------------------------------------------------------------------------#
      if('FCH4_RF' %in% colnames(data)) {
        data <<- data %>% 
          select(-FCH4_RF)
      }
      
      if(missing(training_variables)) {
        # Default training variables
        training_variables <- c('ATMP', 'SWIN', 'WS', 'STMP', 'day_sin')
      }
      
      rf_data = data %>%
        select(datetime, FCH4_f, all_of(training_variables)) %>%
        drop_na()
      
      # Training/validation data, take 70/30 split
      training_data = rf_data %>% 
        slice_sample(prop = 0.7)
      test_data = rf_data %>% 
        filter(!datetime %in% training_data$datetime)
      
      # Fit RandomForest
      # Irvin et al. 2021 has a good overview of parameters to use
      # https://doi.org/10.1016/j.agrformet.2021.108528
      rf = randomForest(FCH4_f ~ .,
                        data=training_data %>% select(-datetime),
                        importance = TRUE)
      
      ch4_rf <<- rf
      
      if(plot) {
        # Plot variable importance
        varImpPlot(rf)
      }
      
      # Gap filling
      predict_data = data %>%
        select(datetime, all_of(training_variables)) %>%
        drop_na() %>%
        mutate(FCH4_RF = predict(rf, newdata=.))
      
      # Join back and fill in final missing values (missing predictor data)
      # with linear interp
      data <<- data %>% 
        left_join(predict_data %>% select(datetime, FCH4_RF)) %>% 
        mutate(# FCH4_RF = zoo::na.approx(FCH4_RF, na.rm=FALSE),
               FCH4_GF_RF = if_else(is.na(FCH4_f), FCH4_RF, FCH4_f),
               RF_src = ifelse(datetime %in% training_data$datetime, 'train',
                               ifelse(datetime %in% test_data$datetime, 'test',
                                      ifelse(is.na(FCH4_RF), 'GF-lm', 'GF-rf'))),
               RF_src = factor(RF_src,
                               levels = c('train', 'test', 'GF-rf', 'GF-lm')))
      
      if(plot) {
        ggplot(data) +
          geom_line(aes(datetime, FCH4_f, col='obs')) + 
          geom_line(aes(datetime, FCH4_GF_RF, col='GF RF')) 
      }
    },
    upload_data = function(db, variables, overwrite=FALSE) {
      if(missing(db)) {
        db <- db_location
      }
      if(missing(variables)) {
        variables = c('NEE_f', 'GPP', 'Reco', 'NEE_GF_NLR', 'NEE_GF_NN',
                      'FCH4_GF_RF', 'Contribution')
      }
      
      con <- dbConnect(RSQLite::SQLite(), db, append=FALSE)
      
      # Check if new IDs need to be generated
      id_name <- tbl(con, 'id_name') %>% 
        collect()
      
      id_append = variables[!variables %in% id_name$name]
      if(length(id_append) > 0) {
        id_append = tibble(id_name = (1L:length(id_append)) + max(id_name$id_name),
                           name = id_append)
        DBI::dbAppendTable(con, 'id_name', id_append)
        id_name <- tbl(con, 'id_name') %>% 
          collect()
      }
      
      # Upload data to postprocessing table
      pp_data <- data %>% 
        select(datetime, any_of(variables)) %>%
        pivot_longer(-datetime) %>%
        left_join(id_name) %>%
        mutate(datetime = as.integer(datetime)) %>% 
        select(-name)
      
      # Check if table exists
      table_list = DBI::dbListTables(con)
      
      if(!'postprocessing' %in% table_list) {
        message('Postprocessing table created in db')
        DBI::dbCreateTable(con, 'postprocessing', pp_data)
        DBI::dbAppendTable(con, 'postprocessing', pp_data)
      } else if(overwrite) {
        DBI::dbWriteTable(con, 'postprocessing', pp_data, overwrite=TRUE)
      } else {
        # This doesn't work so well at the moment, just overwrite for now
        pp_data_db <- tbl(con, 'postprocessing') %>% 
          collect()
        
        pp_data_append <- anti_join(pp_data, pp_data_db) %>% 
          filter(!is.nan(value))
        
        if(nrow(pp_data_append) > 0) {
          DBI::dbAppendTable(con, 'postprocessing', pp_data_append)
        }
      }
      
      DBI::dbDisconnect(con)
      message('Postprocessing data uploaded')
    },
    load_pp_data = function(db, start_dt, end_dt, variables)
    {
      if(missing(db)) {
        db <- db_location
      }
      if(missing(variables)) {
        variables = c('NEE_f', 'GPP', 'Reco', 'NEE_GF_NLR', 'NEE_GF_NN',
                      'FCH4_f', 'FCH4_GF_RF', 'Contribution')
      }
      
      con <- DBI::dbConnect(RSQLite::SQLite(), db)
      tryCatch({
        id_name <- tbl(con, 'id_name')
        df <- tbl(con, 'postprocessing') %>% 
          left_join(id_name) %>% 
          filter(name %in% !!variables,
                 datetime >= !!as.integer(start_dt),
                 datetime <= !!as.integer(end_dt)) %>% 
          select(datetime, name, value) %>% 
          collect() %>% 
          # Pivot to wide format
          pivot_wider(id_col=datetime, names_from=name, values_from=value) %>% 
          mutate(datetime = as_datetime(datetime)) %>% 
          arrange(datetime)
        
        #if('NEE' %in% variables) data <<- data %>% select(-NEE)
        #if('FCH4' %in% variables) data <<- data %>% select(-FCH4)
        
        # Assign to field
        data <<- left_join(data, df)
        message('postprocessing data loaded')
      },
      error = function(e) message(e),
      finally = DBI::dbDisconnect(con))
      
    }
  )
)
