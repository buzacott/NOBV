#------------------------------------------------------------------------------#
# Function definitions
#------------------------------------------------------------------------------#
respiration <- function(r_10, q_10, theta) {
  # Loyd and Taylor (1994) Respiration temperature dependence
  r_10 * q_10^((theta - 10) / 10)
}

respiration_lloyd_taylor <- function(
  r_ref, 
  e_0,
  t,
  t_ref = 15, # Reference temperature of 15C
  t_0 = -46.02 
) {
  # Soil respiration function
  # Inputs:
  #   r_10: Respiration rate at reference temperature
  #   e_0: Activation energy parameter 
  #   t: soil temperature (K)
  #   t_ref: reference temperature (K) (10C)
  #   t_0: temperature as fitted by LloydTaylor (1994)
  # Returns:
  #   r_eco: value of respiration
  r_eco <- r_ref * exp(e_0 * (1/(t_ref - t_0) - 1/(t - t_0)))
  
  return(r_eco)
}

fit_respiration_curve_dep <- function(data) {
  # Fit the respiration to air temperature
  # Needs NEE and ATMP in dataframe
  # Calculate the reference respiration at 10 degrees
  # see Boone et al. 1998 for equation
  
  if(sum(!is.na(data$NEE)) < 100) {
    return(c(NaN, NaN))
  }
  
  data = data %>% 
    drop_na()
  
  # Calculate reference (R10) respiration
  r_10 = data %>% 
    filter(ATMP > 9 & ATMP < 11) %>% 
    summarise(NEE = mean(NEE, na.rm=TRUE)) %>% 
    pull(NEE)
  
  if(is.na(r_10)) r_10 = 2.0
  
  # Get initial estimates of params to avoid fitting errors
  # https://stats.stackexchange.com/questions/183653/getting-the-right-starting-values-for-an-nls-model-in-r
  mod <- suppressWarnings(
    lm(I(log(NEE)-log(r_10)) ~ I((ATMP-10)/10) + 0, data=data)
  )
  
  mod <- nls(
    NEE ~ r_10 * q_10^( (ATMP-10) / 10 ),
    data = data,
    start = list(r_10 = r_10,
                 q_10 = exp(mod$coefficients[1])),
    lower = c(r_10-0.1, -Inf),
    upper = c(r_10+0.1, Inf),
    algorithm = 'port'
  )
  
  fitted_pars <- mod$m$getPars()
  
  return(fitted_pars)
}

fit_respiration_curve <- function(df,
                                  t_ref = 15,
                                  t_0 = -46.02,
                                  rg_threshold = 10) {
  # Estimate the LLoyd Taylor respiration in the postprocessing class
  # Input
  #   df <tbl>: dataframe from postprocessing
  #   t_ref <num>: reference temperature to use (e.g. 10 or 15C)
  #   t_0 <num>: temperature as fitted by Lloyd and Taylor (1994)
  #   rg_threshold <num>: radiation threshold to determine nighttime
  
  # Create respiration df and filter to night time conditions
  r_df <- df %>% 
    # Apply filters
    filter(night == TRUE, SWIN < rg_threshold) %>% 
    filter(!is.na(NEE),
           NEE > 0)
  
  # Find the limits to start the sliding window
  dt_min <- min(r_df$datetime)
  dt_max <- max(r_df$datetime)
  
  # Initial timestep
  dt_i <- floor_date(dt_min, 'day') + 7*86400
  # List to store e_0 results
  e_0_results <- list()
  # Counter to store the e_0 results
  i <- 1
  while(dt_i < dt_max) {
    # Collect the current day 15 window
    r_df_sub <- r_df %>% 
      filter(datetime >= dt_i - 7 * 86400,
             datetime < dt_i + 8 * 86400) %>% 
      # Temperature range must be more than 5 degrees
      filter(diff(range(ATMP)) > 5)
    
    # Check that there are more than 6 data points in the window
    if(nrow(r_df_sub) > 6) {
      # Do the first fit to find e_0
      mod <- nls(
        formula = NEE ~ respiration_lloyd_taylor(r_ref,
                                                 e_0,
                                                 ATMP,
                                                 t_ref = t_ref),
        trace = FALSE,
        data = r_df_sub,
        start = list(r_ref = 2, e_0 = 200) 
        # Check why I use these starting parameters
        # Sometimes in the mean of the window is used for e_0
      )
      
      # Repeat fit with trimmed dataset (5-95% residuals)
      r_df_sub <- r_df_sub %>%
        filter(residuals(mod) > quantile(residuals(mod), 0.05) |
                 residuals(mod) < quantile(residuals(mod), 0.95))
      
      # Fit again
      mod <- nls(
        formula = NEE ~ respiration_lloyd_taylor(r_ref,
                                                 e_0,
                                                 ATMP,
                                                 t_ref = t_ref),
        trace = FALSE,
        data = r_df_sub,
        start = list(r_ref = 2, e_0 = 200)
      )
      
      # Get the summary of the regression to store the results
      s_mod <- summary(mod)
      
      e_0_results[[i]] <- tibble(
        dt_i = dt_i,
        dt_start = dt_i - 7 * 86400,
        dt_end = dt_i + 8 * 86400,
        r_ref = coef(s_mod)['r_ref', 'Estimate'],
        r_ref_std = coef(s_mod)['r_ref', 'Std. Error'],
        e_0 = coef(s_mod)['e_0', 'Estimate'],
        e_0_std = coef(s_mod)['e_0', 'Std. Error']
      )
      # Old results return if addiational data is wanted
      # e_0_results[[i]] <- r_df_sub %>% 
      #   summarise(dt_start = dt_i,
      #             dt_end = dt_i + 15 * 86400,
      #             dt_min = min(datetime),
      #             dt_max = max(datetime),
      #             ATMP = mean(ATMP),
      #             NEE = mean(NEE)) %>% 
      #   mutate(r_ref = coef(s_mod)['r_ref', 'Estimate'],
      #          r_ref_std = coef(s_mod)['r_ref', 'Std. Error'],
      #          e_0 = coef(s_mod)['e_0', 'Estimate'],
      #          e_0_std = coef(s_mod)['e_0', 'Std. Error'])
      i <- i + 1
    }
    # dt_i <- dt_i + 5 * 86400 # Reichstein use sliding 5 day windows here
    dt_i <- dt_i + 15 * 86400 # Wutzler et al. 2018 use consec 15 day windows
  }
  
  # Combine the results and arrange the data
  e_0_results <- bind_rows(e_0_results)
  
  e_0_est <- e_0_results %>% 
    # Filter to where std error is less than 50% of the estimate
    # and e_0 is within an accepted range (30-450) 
    filter(e_0_std < 0.5 * e_0,
           e_0 > 30, 
           e_0 < 450) %>% 
    arrange(e_0_std)
  
  # Take the mean of the 3 best e_0 parameters
  e_0_est <- mean(e_0_est$e_0[1:3], na.rm=TRUE)
  df$e_0 <- e_0_est
  
  # Estimate r_ref
  # 7 day windows shifted consecutively by 4 days
  dt_i <- floor_date(dt_min, 'day') + 3 * 86400
  r_ref_results <- list()
  i <- 1
  while(dt_i < dt_max) {
    r_df_sub <- r_df %>% 
      filter(datetime >= dt_i - 3 * 86400,
             datetime < dt_i + 4 * 86400)
    
    if(nrow(r_df_sub) > 6) {
      mod <- nls(
        formula = NEE ~ respiration_lloyd_taylor(r_ref,
                                                 e_0=e_0_est,
                                                 t = ATMP,
                                                 t_ref = t_ref),
        trace = FALSE,
        data = r_df_sub,
        start = list(r_ref = 2) # mean(r_df_sub$NEE, na.rm=TRUE))
      )
      
      s_mod <- summary(mod)
      
      r_ref_results[[i]] <- r_df_sub %>% 
        summarise(dt_i = dt_i,
                  dt_start = dt_i - 3 * 86400,
                  dt_end = dt_i + 4 * 86400,
                  dt_min = min(datetime),
                  dt_max = max(datetime),
                  ATMP = mean(ATMP),
                  NEE = mean(NEE)) %>% 
        mutate(r_ref = coef(s_mod)['r_ref', 'Estimate'],
               r_ref_std = coef(s_mod)['r_ref', 'Std. Error'])
      i <- i + 1
    }
    dt_i <- dt_i + 4 * 86400
  }
  
  # Linearly interpolate r_ref results
  r_ref_results <- bind_rows(r_ref_results)
  r_ref_interp <- tibble(datetime = seq.POSIXt(min(df$datetime), max(df$datetime), 1800)) %>% 
    left_join(r_ref_results %>% 
                select(datetime = dt_i, # Assign to midpoint
                       r_ref) %>% 
                mutate(datetime = datetime + 2*86400)) 
  # Assign the first result value to the first data point as the linear interpo-
  # lation won't go backwards
  r_ref_interp$r_ref[1] = r_ref_results$r_ref[1]
  r_ref_interp$r_ref[nrow(r_ref_interp)] = r_ref_results$r_ref[nrow(r_ref_results)]
  # Linearly interpolate
  r_ref_interp <- r_ref_interp %>% 
    mutate(r_ref = zoo::na.approx(r_ref, na.rm=FALSE))
  
  # Calculate r_eco using complete parameter sets
  df <- left_join(df, r_ref_interp) %>% 
    mutate(Reco = respiration_lloyd_taylor(r_ref = r_ref,
                                           e_0 = e_0_est,
                                           t = ATMP,
                                           t_ref = t_ref)) %>% 
    select(datetime, r_ref, e_0, Reco)
  
  return(df)
}

gpp <- function(x, b1, b2) {
  # GPP light dependence curve
  (b1 * x) / (x + b2)
}

fit_gpp_curve <- function(data) {
  
  if(sum(!is.na(data$GPP)) < 100) {
    var = c(NaN, NaN)
  } else {
    var <- tryCatch({
      mod <- nls(
        GPP ~ b1 * SWIN / (SWIN + b2),
        data = data,
        start = list(b1 = -1, b2 = 10),
        control = list(maxiter=100)
      )
      mod$m$getPars()
    },
    error = function(e) {
      message(e)
      return(c(NaN, NaN))
    })
  }
  return(var)
}

season = function(dt) {
  DJF = c(12, 1, 2)
  MAM = 3:5
  JJA = 6:8
  SON = 9:11
  
  m = lubridate::month(dt)
  
  s = ifelse(m %in% DJF, 'DJF',
             ifelse(m %in% MAM, 'MAM',
                    ifelse(m %in% JJA, 'JJA',
                           ifelse(m %in% SON, 'SON', NA))))
  s = factor(s, levels = c('DJF', 'MAM', 'JJA', 'SON'))
  return(s)
}

kljun04 <- function(zm, z0, ustar, sigma_w) {
  # Inputs
  #   zm: Measurement height above displacement height (i.e. z-d) [m]
  #   z0: roughness length [m]
  #   ustar: friction velocity [m/s]
  #   sigma_w: standard deviation of vertical velocity fluctuations
  #            (after coordinate rotation)
  # Returns
  #   Peak contribution distance and distances of contribution percentages [m]
  
  # Constants
  Ac = 4.277
  B = 3.418
  Ad = 1.685
  
  # Taken from EddyPro Engine Code, which was provided to LICOR by Kljun
  L = c(0.000000, 0.302000, 0.368000, 0.414000, 0.450000, 0.482000, 0.510000,
        0.536000, 0.560000, 0.579999, 0.601999, 0.621999, 0.639999, 0.657998,
        0.675998, 0.691998, 0.709998, 0.725998, 0.741997, 0.755997, 0.771997,
        0.785997, 0.801997, 0.815996, 0.829996, 0.843996, 0.857996, 0.871996,
        0.885995, 0.899995, 0.911995, 0.925995, 0.939995, 0.953995, 0.965994,
        0.979994, 0.993994, 1.005994, 1.019994, 1.033994, 1.045993, 1.059993,
        1.073993, 1.085993, 1.099993, 1.113993, 1.127992, 1.141992, 1.155992,
        1.169992, 1.183992, 1.197991, 1.211991, 1.225991, 1.239991, 1.253991,
        1.269991, 1.283990, 1.299990, 1.315990, 1.329990, 1.345990, 1.361989,
        1.379989, 1.395989, 1.411989, 1.429989, 1.447988, 1.465988, 1.483988,
        1.501988, 1.521987, 1.539987, 1.559987, 1.581987, 1.601986, 1.623986,
        1.647986, 1.669985, 1.693985, 1.719985, 1.745984, 1.773984, 1.801984,
        1.831983, 1.863983, 1.895983, 1.931982, 1.969982, 2.009982, 2.053984,
        2.101986, 2.153988, 2.211991, 2.279994, 2.355998)
  
  # From Kljun et al. 2004 eq 13-16
  c = Ac*(B - log(z0))
  d = Ad*(B - log(z0))
  
  Xstar = c(1, L[c(11, 31, 51, 71, 91)]) * c - d
  X = Xstar * (z-zd) * (sigma_w / ustar)^-0.8
  X = round(X, 4)
  names(X) = c('x_peak', paste0('x_', c(10, 30, 50, 70, 90), '%'))
  
  return(X)
}

# For Vickers and Mahrt 1997 tests
sep_digits <- Vectorize(function(x) {
  n_digits <- floor(log10(x))
  y <- vector(length = n_digits)
  for(i in 0:n_digits) {
    y[i+1] <- (x %/% 10^(i)) %% 10^1
  }
  return(rev(y))
})

gpp_lasslop_2010 <- function(df,
                             t_ref = 15,
                             t_0 = -46.02,
                             rg_threshold_resp = 20,
                             rg_threshold_gpp = 4) {
  
  # Estimate respiration
  respiration <- fit_respiration_curve(df = df,
                                       t_ref = t_ref,
                                       t_0 = t_0,
                                       rg_threshold = 20)
  
  # GPP model
  rg_threshold = 4
  alpha = 0.01
  e_0 = respiration$e_0[1]
  
  dt_min <- min(df$datetime)
  dt_max <- max(df$datetime)
  
  dt_i <- floor_date(dt_min, 'day')
  opt_results <- list()
  i <- 1
  
  while(dt_i < dt_max) {
    # Day df
    df_d = df %>% 
      # 4 day moving window for daytime
      filter(datetime >= dt_i, datetime <= dt_i + 4 * 86400) %>% 
      filter(!night, SWIN > rg_threshold) %>% 
      drop_na()
    # Night df
    df_n = df %>% 
      # 12 day moving window for nighttime
      filter(datetime >= dt_i, datetime <= dt_i + 12 * 86400) %>% 
      filter(night, SWIN < rg_threshold,
             NEE > 0) %>% 
      drop_na()
    
    if(nrow(df_n) > 6) {
      # First estimate E0 using night data
      try({
        mod <- nls(
          formula = NEE ~ respiration_lloyd_taylor(r_ref,
                                                   e_0,
                                                   ATMP,
                                                   t_ref = t_ref),
          trace = FALSE,
          data = df_n,
          start = list(r_ref = mean(df_n$NEE), e_0 = 100)
        )
        e_0_fit <- coef(mod)[2]
        if(i == 1) {
          e_0 <- if_else(e_0_fit < 50, 50, e_0_fit)
          e_0 <- if_else(e_0_fit > 400, 400, e_0)
        } else {
          if(e_0_fit >= 50 & e_0_fit <= 400) e_0 <- e_0_fit
        }
        
        initial_par = c(alpha = 0.01,
                        beta = unname( abs( quantile(df_d$NEE, 0.03) - quantile(df_d$NEE, 0.97) ) ),
                        r_ref = mean(df_n$NEE))
        
        opt <- optim(initial_par,
                     nee_lrc_optim,
                     method = 'BFGS',
                     NEE = df_d$NEE,
                     SWIN = df_d$SWIN,
                     ATMP = df_d$ATMP,
                     e_0 = e_0,
                     t_ref = t_ref,
                     t_0 = t_0,
                     hessian = TRUE)
        
        std_err <- sqrt(abs(diag(solve(-opt$hessian))))
        
        opt_results[[i]] <- df_d %>% 
          summarise(dt_start = dt_i,
                    dt_end = dt_i + 4 * 86400,
                    dt_min = min(datetime),
                    dt_max = max(datetime),
                    NEE = mean(NEE),
                    ATMP = mean(ATMP),
                    SWIN = mean(SWIN)) %>% 
          mutate(alpha = opt$par[1],
                 beta = opt$par[2],
                 r_ref = opt$par[3],
                 e_0 = e_0,
                 alpha_std = std_err[1],
                 beta_std = std_err[2],
                 r_ref_std = std_err[3])
        
        i <- i + 1
      })
    }
    dt_i <- dt_i + 2 * 86400
  }
  
  # Bind together and filter parameter sets
  opt_results = bind_rows(opt_results) %>% 
    filter(alpha > 0 & alpha <= 0.22,
           beta > 0 & beta < 250,
           beta <= 100 | beta_std < beta,
           r_ref > 0) %>% 
    select(datetime = dt_start,
           alpha, beta, r_ref, e_0) %>% 
    # Assign dt to midpoint
    mutate(datetime = datetime + 2 * 86400)
  
  # Linearly interpolate gpp results
  interp <- tibble(datetime = seq.POSIXt(dt_min, dt_max, 1800)) %>% 
    left_join(opt_results) 
  # No data for first rows, simply backfill first valid result
  interp[1, c('alpha', 'beta', 'r_ref', 'e_0')] = opt_results[1, c('alpha', 'beta', 'r_ref', 'e_0')]
  # Same at the end
  interp[nrow(interp), c('alpha', 'beta', 'r_ref', 'e_0')] = opt_results[nrow(opt_results), c('alpha', 'beta', 'r_ref', 'e_0')]
  interp <- interp %>% 
    mutate(across(-datetime, ~zoo::na.approx(.x, na.rm=FALSE)))
  
  # Finally join the parameter sets to the initial df and apply functions
  df <- df %>% 
    left_join(interp) %>% 
    mutate(GPP_ll = -nee_lrc(alpha, beta, SWIN, 0),
           Reco_ll = respiration_lloyd_taylor(r_ref, e_0, ATMP),
           NEE_ll = -GPP_ll + Reco_ll)
  
  return(df)
}

nee_lrc <- function(alpha,
                    beta,
                    SWIN,
                    gamma) {
  # NEE light response curve
  NEE <- - (alpha * beta * SWIN) / (alpha * SWIN + beta) + gamma
  return(NEE)
}

nee_lrc_optim <- function(x, NEE, SWIN, ATMP, e_0, t_ref, t_0) {
  alpha <- x[1]
  beta <- x[2]
  r_ref <- x[3]
  
  fit <- nee_lrc(alpha = alpha,
                 beta = beta,
                 SWIN = SWIN,
                 respiration_lloyd_taylor(r_ref, e_0, ATMP, t_ref, t_0))
  
  obj_fun <- sum( (NEE - fit)^2 / sd(NEE) )
  
  return(obj_fun)
}
