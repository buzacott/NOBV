#------------------------------------------------------------------------------#
# Scripts to postprocess NOBV eddy covariance data
# Code written by Alexander Buzacott, with large parts adapted from hidmu's 
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
# Import Python libraries via reticulate for neural network gapfilling
# Comment out or ignore if you don't want to setup reticulate or use neural
# network gapfilling.
library(reticulate)
reticulate::use_condaenv(condaenv = 'RStudio')
sklearn <- reticulate::import('sklearn')
model_selection <- reticulate::import('sklearn.model_selection')
preproc <- reticulate::import("sklearn.preprocessing")
pipeline <- reticulate::import("sklearn.pipeline")
neural_network <- reticulate::import("sklearn.neural_network")
pd <- reticulate::import('pandas')

#------------------------------------------------------------------------------#
# Function definitions
#------------------------------------------------------------------------------#
respiration <- function(x, b, r_ref) {
  # Loyd and Taylor (1994) Respiration temperature dependence
  r_ref * b^((x - 10) / 10 )
}

gpp <- function(x, b1, b2) {
  # GPP light dependence curve
  (b1 * x) / ( x + b2)
}

fit_respiration_curve <- function(data) {
  # Fit the repsiration to air temperature
  # Needs NEE and ATMP in dataframe
  # Calculate the reference respiration at 10 degrees
  # see Boone et al. 1998 for equation
  
  if(sum(!is.na(data$NEE)) < 100) {
    return(c(NaN, NaN))
  }
  
  data = data %>% 
    drop_na()
  
  # Calculate reference respiration
  r_ref = data %>% 
    filter(ATMP < 11 & ATMP > 9) %>% 
    summarise(NEE = mean(NEE, na.rm=TRUE)) %>% 
    pull(NEE)
  # TODO: improve this
  
  
  if(is.na(r_ref)) r_ref = 2.0
  
  # Get initial estimates of params to avoid fitting errors
  # https://stats.stackexchange.com/questions/183653/getting-the-right-starting-values-for-an-nls-model-in-r
  mod <- suppressWarnings(
    lm(I(log(NEE)-r_ref) ~ 0 + ATMP, data = data)
  )
  
  mod <- nls(
    NEE ~ r_ref * b^( (ATMP-10) / 10 ),
    data = data,
    start = list(b = exp(mod$coefficients[1]),
                 r_ref = r_ref),
    lower = c(-Inf, r_ref-0.1),
    upper = c(Inf,  r_ref+0.1),
    algorithm = 'port'
  )
  
  fitted_pars <- mod$m$getPars()
  
  return(fitted_pars)
}

fit_gpp_curve <- function(data) {
  if(sum(!is.na(data$GPP)) < 100) {
    var = c(NaN, NaN)
  } else {
    mod <- nls(
      GPP ~ b1 * SWIN / (SWIN + b2),
      data = data,
      start = list(b1 = -1, b2 = 10),
      control = list(maxiter=100)
    )
    var = mod$m$getPars()
  }
  return(var)
}

# Import Kljun model
source('src/calc_footprint_FFP_climatology.R')
#------------------------------------------------------------------------------#
# Postprocessing class
#------------------------------------------------------------------------------#
Postprocessing <- setRefClass(
  "Postprocessing",
  fields = list(
    default_variables = "character",
    db_location = "character",
    data = "ANY",
    ch4_rf = "ANY"
  ),
  methods = list(
    initialize = function()
    {
      default_variables <<- c(
        'FC','SC_SINGLE', 'FCH4', 'USTAR', 'DISPLACEMENT_HEIGHT',
        'MT1_ATMP_1_H_200_Avg', 'WS', 'WD', 'MO_LENGTH', 'V_SIGMA',
        'MT1_SWIN_1_H_180','FC_SSITC_TEST', 'FCH4_SSITC_TEST',
        'MT1_PAR_1_H_180', 'MT1_WIND_1_H_200_Avg', 'MT1_WINS_1_H_200_Avg',
        'MT1_NETL_1_H_180', 'MT1_NETS_1_H_180', 'SPG_STMP_1_D_005_Avg')
      
      db_location <<- file.path('data', 'sample.db')
    },
    load_data = function(start_dt, end_dt, db = NULL, variable_list = NULL)
    {
      # Check to see if there is a given list of variables, otherwise use default
      # variables
      if(is.null(variable_list)) {
        variable_list <- default_variables
      }
      if(is.null(db)) {
        db <- db_location
      }

      con <- DBI::dbConnect(RSQLite::SQLite(), db)
      tryCatch({
        id_name <- tbl(con, 'id_name')
        df <- tbl(con, 'data') %>% 
          left_join(id_name) %>% 
          filter(name %in% !!variable_list,
                 datetime >= !!as.integer(start_dt),
                 datetime <= !!as.integer(end_dt)) %>% 
          select(datetime, name, value) %>% 
          collect() %>% 
          # Pivot to wide format
          pivot_wider(id_col=datetime, names_from=name, values_from=value) %>% 
          # Add NEE, night time and month columns
          mutate(datetime = as_datetime(datetime, tz='UTC'),
                 NEE = FC + SC_SINGLE,
                 night = MT1_SWIN_1_H_180 < 10,
                 month = month(datetime),
                 FCH4 = FCH4 / 1000) %>% # nmol/s/m2 -> µmol/s/m2
          arrange(datetime)
        
        # Rename columns
        df <- df %>% 
          rename('ATMP' = 'MT1_ATMP_1_H_200_Avg',
                 'SWIN' = 'MT1_SWIN_1_H_180',
                 'FC_flag' = 'FC_SSITC_TEST',
                 'FCH4_flag' = 'FCH4_SSITC_TEST',
                 'PAR_met' = 'MT1_PAR_1_H_180' ,
                 'WD_met' = 'MT1_WIND_1_H_200_Avg',
                 'WS_met' = 'MT1_WINS_1_H_200_Avg',
                 'NETL_met' = 'MT1_NETL_1_H_180',
                 'NETS_met' = 'MT1_NETS_1_H_180',
                 'STMP' = 'SPG_STMP_1_D_005_Avg')
        
        # Assign to field
        data <<- df
      },
      error = function(e) message(e),
      finally = DBI::dbDisconnect(con))
    },
    filter_co2 = function(ustar_threshold = 0.10)
    {
      # Filter the co2 data for outliers, friction velocity and quality flags
      
      variables <- c('FC', 'NEE')
      # Ustar filtering
      # The ustar filtering can be improved to account of seasons / literature
      # 0.1 was chosen by Hidde as it removed most bad data, but some stuff still gets through
      data <<- data %>% 
        mutate(across(all_of(variables), ~if_else(USTAR < ustar_threshold, NaN, .x)))
      
      # Quality flag filtering, possible to also add 1 if you want to be strict
      # Flags
      # 0 is good
      # 1 is ok
      # 2 is bad
      # Can customise flagging system in EddyPro
      data <<- data %>% 
        mutate(across(all_of(variables), ~if_else(FC_flag > 1, NaN, .x)))
      
      # Filter outliers
      # Note: To be replace with robust outlier detection
      data <<- data %>% 
        mutate(across(all_of(variables), ~if_else(FC > quantile(FC, 0.97, na.rm=TRUE), NaN, .x)),
               across(all_of(variables), ~if_else(FC < quantile(FC, 0.01, na.rm=TRUE), NaN, .x)))
    },
    filter_ch4 = function()
    {
      # Filter methane data for quality flags and outliers
      
      # Quality flags filter
      data <<- data %>% 
        mutate(FCH4 = if_else(FCH4_flag > 1, NaN, FCH4))
      
      # Outlier detection
      # Note replace with robust detection algorithm
      data <<- data %>% 
        mutate(FCH4 = if_else(FCH4 > quantile(FCH4, 0.99, na.rm=TRUE), NaN, FCH4),
               FCH4 = if_else(FCH4 < quantile(FCH4, 0.01, na.rm=TRUE), NaN, FCH4))
    },
    flux_partitioning = function(plot = FALSE) 
    {
      #Partition fluxes
      
      data <<- data %>% 
        mutate(season = quarter(datetime))
      
      # Respiration
      # Select night data with only respiration
      night_flux = data %>%
        filter(night == TRUE & NEE > 0) %>%
        select(datetime, NEE, ATMP, month, season) %>%
        drop_na()
      
      variables = night_flux %>% 
        group_by(season) %>% 
        nest() %>% 
        mutate(variables = purrr::map(data, function(x) {
          variables = fit_respiration_curve(x)
          return(tibble(r_param = variables[[1]],
                        r_ref = variables[[2]]))
        })) %>% 
        select(-data) %>% 
        unnest(variables)
      
      # Fill up empty season with averages and calculate respiration for whole
      # dataset
      data <<- data %>% 
        left_join(variables) %>% 
        mutate(across(c(r_ref, r_param), ~replace_na(.x, mean(.x))),
               Respiration = respiration(ATMP, r_param, r_ref))
      
      # GPP
      # select data for fitting
      gpp_data = data %>% 
        select(datetime, SWIN, NEE, Respiration, month, season) %>% 
        drop_na() %>% 
        mutate(GPP = NEE - Respiration)
      
      variables = gpp_data %>% 
        group_by(season) %>% 
        nest() %>% 
        mutate(variables = purrr::map(data, function(x) {
          variables = fit_gpp_curve(x)
          return(tibble(beta1 = variables[[1]],
                        beta2 = variables[[2]]))
        })) %>% 
        select(-data) %>% 
        unnest(variables)
      
      data <<- data %>% 
        left_join(variables) %>% 
        mutate(across(c(beta1, beta2), ~replace_na(.x, mean(.x))),
               GPP = gpp(SWIN, beta1, beta2))
      
      if(plot) {
        data %>% 
          select(datetime, NEE, GPP, Respiration) %>% 
          pivot_longer(-datetime) %>% 
          ggplot(aes(datetime, value, col=name)) +
          geom_line() +
          labs(y = expression(value~'('*µmol.s^-1*.m^-2*')'))
      }
    },
    gapfilling_nlr = function(plot = FALSE)
    {
      data <<- data %>% 
        mutate(nlr_NEE = Respiration + GPP,
               FC_GF_NLR = if_else(is.na(NEE), nlr_NEE, NEE))
    },
    gapfilling_nn = function(plot = FALSE)
    {
      data$month_sin <<- sin(pi * (data$month - 1) / 12)
      
      training_data = data %>%
        select(datetime, NEE, ATMP, SWIN, month_sin) %>%
        drop_na()
      
      x = training_data %>%
        select(ATMP, SWIN, month_sin)
      y = training_data %>%
        pull(NEE)
      
      # set.seed(123)
      # 
      # train_idx = sample.int(nrow(x)*0.75)
      # 
      # x_train = x[train_idx,]
      # y_train = y[train_idx]
      # x_test = x[-train_idx,]
      # y_test = y[-train_idx]
      # 
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
        select(datetime, ATMP, SWIN, month_sin) %>% 
        drop_na() 
      pred_data$FC_GF_NN = pipe$predict(pred_data %>% select(-datetime))
      
      # Merge to data
      data <<- left_join(data,
                         pred_data %>% select(datetime, FC_GF_NN)) %>% 
        mutate(FC_GF_NN = if_else(is.na(NEE), FC_GF_NN, NEE))
      
      if(plot == TRUE) {
        ggplot(data, aes(datetime, FC_GF_NN)) +
          geom_line()
      }
    },
    gapfill_ch4_nn = function(plot = FALSE)
    {
      
      training_variables = c('ATMP', 'SWIN','month_sin', #'FC_GF_NN',
                            'WD_met', 'WS_met', 'PAR_met', 'NETL_met', 'NETS_met')
      
      training_data = data %>%
        select(datetime, FCH4, all_of(training_variables)) %>%
        drop_na()

      x = training_data %>%
        select(all_of(training_variables))
      y = training_data %>%
        pull(FCH4)
      
      train_test = model_selection$train_test_split(x, y)
      
      pipe = sklearn$pipeline$Pipeline(list(
        tuple('scaler', preproc$StandardScaler()),
        tuple('nn', neural_network$MLPRegressor(
          hidden_layer_sizes = tuple(100L, 100L, 100L),
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
      pred_data$FCH4_GF_NN = pipe$predict(pred_data %>% select(-datetime))
      
      # Merge to data
      data <<- left_join(data,
                         pred_data %>% select(datetime, FCH4_GF_NN)) %>% 
        mutate(FCH4_GF_NN = if_else(is.na(FCH4), FCH4_GF_NN, FCH4))
      
      if(plot == TRUE) {
        ggplot(data, aes(datetime, FCH4_GF_NN)) +
          geom_line()
      }
    },
    gapfill_ch4_rf = function(plot = FALSE)
    {
      # Gapfill methane using Random Forest
      # Based on Irvin et al. 2021
      # https://doi.org/10.1016/j.agrformet.2021.108528
      
      data$month_sin <<- sin(pi * (data$month - 1) / 12)
      
      training_variables <- c('ATMP', 'SWIN', 'PAR_met', 'WS_met', #NEE,
                              'NETL_met', 'NETS_met', 'STMP', 'month_sin')
      
      rf_data = data %>%
        select(datetime, FCH4, all_of(training_variables)) %>%
        drop_na()
      
      # Training/validation data, take 70/30 split
      training_data = rf_data %>% 
        slice_sample(prop = 0.7)
      test_data = rf_data %>% 
        filter(!datetime %in% training_data$datetime)
      
      # Fit RandomForest
      # TODO compare RF params to Irvin et al.
      rf = randomForest(FCH4 ~ .,
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
        mutate(RF_src = ifelse(datetime %in% training_data$datetime, 'train',
                               ifelse(datetime %in% test_data$datetime, 'test',
                                      ifelse(is.na(FCH4_RF), 'GF-lm', 'GF-rf'))),
               RF_src = factor(RF_src,
                               levels = c('train', 'test', 'GF-rf', 'GF-lm'))) %>% 
        mutate(FCH4_RF = zoo::na.approx(FCH4_RF, na.rm=FALSE),
               FCH4_GF_RF = if_else(is.na(FCH4), FCH4_RF, FCH4))
      
      if(plot) {
        ggplot(data) +
          geom_line(aes(datetime, FCH4, col='obs')) + 
          geom_line(aes(datetime, FCH4_GF_RF, col='GF RF')) 
      }
    },
    upload_data = function(db, variables, overwrite=FALSE) {
      if(missing(db)) {
        db <- db_location
      }
      if(missing(variables)) {
        variables = c('NEE', 'GPP', 'Respiration', 'FC_GF_NLR', 'FC_GF_NN',
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
        variables = c('NEE', 'GPP', 'Respiration', 'FC_GF_NLR', 'FC_GF_NN',
                      'FCH4_GF_RF', 'Contribution')
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
        
        # Assign to field
        data <<- left_join(data, df)
        message('postprocessing data loaded')
      },
      error = function(e) message(e),
      finally = DBI::dbDisconnect(con))
      
    }
  )
)

Contribution <- setRefClass(
  'Contribution',
  fields = list(
    cdata = "ANY",
    xy = "numeric",
    shapefile = "ANY",
    epsg = "ANY",
    kljun_params = "ANY",
    footprints = "ANY"
  ),
  methods = list(
    initialize = function(data,
                          xy,
                          epsg = 28992,
                          shapefile,
                          kljun_params = list(
                            domain = c(-120, 120,-120, 120),
                            nx = 240,
                            ny = 240
                          ))
    {
      # Footprint contribution class initialisation
      # Inputs
      #   data <df>: dataframe/tibble from the postprocesing class
      #   xy <numeric>: projected coordinates of the EC tower
      #   epsg <int/character>: EPSG code for xy and for use with the shapefile
      #                         and rasters
      #   shapefile <character>: path to the shapefile to calculate the percent
      #                          contribution within an AOI
      #   kljun_params <list>: Kljun footprint parameters
      
      # Select variables that are needed for footprint
      contribution_variables = c('datetime', 'USTAR', 'DISPLACEMENT_HEIGHT',
                                 'WS', 'WD', 'MO_LENGTH', 'V_SIGMA')
      cdata <<- data %>% 
        select(all_of(contribution_variables))
    
      # Assign variables to class
      epsg <<- epsg
      xy <<- xy
      xy_sf <- st_point(xy) %>% 
        st_sfc(crs=epsg)
      lonlat <- st_transform(xy_sf, crs=4326) %>% 
        st_coordinates()
      shapefile <<- read_sf(shapefile) %>% 
        st_transform(epsg)
      
      kljun_params <<- kljun_params
      
      #Calculate parameters that are needed for contribution calculations
      latitude = lonlat[2] * (pi / 180)        # Degrees -> radians
      angular_velocity = 7.2921159 * 10^-5     # rad/s
      coriolis_parameter = 2 * angular_velocity * sin(latitude)
      # coefficient for calculating planetary boundary layer height
      c_n = 0.3  
      # planetary boundary layer height
      cdata$PBLH <<- c_n * cdata$USTAR / coriolis_parameter 
    },
    run_kljun = function(row, verbosity=FALSE)
    {
      # Run the Kljun footprint model. Code has been sourced from the kljun
      # website: https://footprint.kljun.net
      # Inputs
      #   row <vec>: row of data / timeseries slice from the cdata df
      #   verbosity <bool>: flag to suppress function output
      # Returns
      #   <list> output of FPP climatology function
      # TODO: the verbosity flag won't work on windows
      if(!verbosity) sink('/dev/null')
      ffp_clim = calc_footprint_FFP_climatology(
        zm = (3 - row[['DISPLACEMENT_HEIGHT']]),
        z0 = NaN,
        rs = NaN,
        umean = row[['WS']],
        h = row[['PBLH']],
        ol = row[['MO_LENGTH']],
        sigmav = row[['V_SIGMA']],
        ustar = row[['USTAR']],
        wind_dir = row[['WD']],
        domain = kljun_params$domain,
        nx = kljun_params$nx,
        ny = kljun_params$ny
      )
      if(!verbosity) sink()
      return(ffp_clim)
    },
    make_raster = function(row)
    {
      # Internal function to create a raster from the Kljun footprint matrix
      # Inputs
      #   row <vec>: row of data / timeseries slice from the cdata df
      # Returns
      #   <raster> footprint raster
      m <- matrix(NaN, kljun_params$ny+1, kljun_params$nx+1)
      if( !any(is.na(row)) ) {
        ffp_clim = run_kljun(row)
        if(ffp_clim$flag_err != 1) {
          m <- ffp_clim$fclim_2d
        }
      }
      suppressWarnings({
        # Raster expects left to right, top to bottom
        r <- raster(m[nrow(m):1,],
                    xmn = xy[1] + kljun_params$domain[1] + 0.5,
                    xmx = xy[1] + kljun_params$domain[2] + 0.5,
                    ymn = xy[2] + kljun_params$domain[3] + 0.5,
                    ymx = xy[2] + kljun_params$domain[4] + 0.5,
                    crs = epsg)
      })
      return(r)
    },
    calculate_footprints = function(parallel = NULL)
    {
      # Get a timeseries of EC footprints
      # Inputs:
      #   parallel <int>: number of cores to use. NULL to not use parallel
      #                   processing
      # Returns:
      #   <raster brick>
      message('Calculating footprints')
      tic()
      if(is.null(parallel)) {
        
        footprints <<- apply(cdata %>% select(-datetime),
                        MARGIN = 1,
                        simplify = FALSE,
                        make_raster)
      } else {
        registerDoParallel()
        cl = makeCluster(parallel)
        footprints <<- foreach(row=iterators::iter(cdata, by='row')) %dopar% {
          make_raster(row)
        }
        stopCluster(cl) 
      }
      toc()
      footprints <<- brick(footprints)
      names(footprints) <<- format_ISO8601(cdata$datetime)
    },
    calculate_contribution = function(keep_footprints = TRUE, parallel = FALSE)
    {
      # Calculate the percent contribution of flux coming from a defined area.
      # The shapefile, EPSG code, and Kljun model parameters are specified
      # during class creation
      # Inputs
      #   parallel <int>: number of cores to use. NULL to not use parallel
      #                   processing
      # Returns
      #   adds a vector of percent contribution to cdata
      
      tic()
      if(keep_footprints) {
        calculate_footprints(parallel)
        
        message('Calculating contribution percentage')
        r_sum <- cellStats(footprints, 'sum')
        suppressWarnings({
          ee_sum <- exactextractr::exact_extract(footprints,
                                                 shapefile,
                                                 fun='sum')
        })
        contribution <- unname( t(ee_sum)[,1] / r_sum)
      } else {
        if(parallel) {
          registerDoParallel()
          cl = makeCluster(4)
          tic()
          contribution <- foreach(row=iterators::iter(cdata, by='row'),
                                  .combine = c) %dopar% {
                                    r <- make_raster(row)
                                    r_sum = cellStats(r, sum)
                                    ext_sum = exactextractr::exact_extract(r, shapefile, fun='sum')
                                    cont = ext_sum/r_sum
                                  }
        } else {
          contribution <- apply(cdata %>% select(-datetime), 1, function(row) {
            r <- make_raster(row)
            r_sum = cellStats(r, sum)
            ext_sum = exactextractr::exact_extract(r, shapefile, fun='sum')
            cont = ext_sum/r_sum
            return(cont)
          })
        }
      }
      toc()
      cdata$Contribution <<- contribution
    }
  )
)
