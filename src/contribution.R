#------------------------------------------------------------------------------#
# Contribution class script
# Applies the Kljun et al. 2015 FFP climatology to estimate source contribution
# within the footprint
#------------------------------------------------------------------------------#
Contribution <- setRefClass(
  'Contribution',
  fields = list(
    cdata = "ANY",
    tower_height = "numeric",
    xy = "numeric",
    shapefile = "ANY",
    epsg = "ANY",
    kljun_params = "ANY",
    footprints = "ANY",
    mean_footprint = "ANY"
  ),
  methods = list(
    initialize = function(data,
                          xy,
                          epsg = 28992,
                          shapefile,
                          tower_height = 3,
                          kljun_params = list(
                            domain = c(-120, 120,-120, 120),
                            nx = 241,
                            ny = 241
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
      #   tower_height <character>: height of the EC tower/anemometer
      #   kljun_params <list>: Kljun footprint parameters
      
      # Select variables that are needed for footprint
      contribution_variables = c('datetime', 'USTAR', 'DISPLACEMENT_HEIGHT',
                                 'ROUGHNESS_LENGTH', 'WS', 'WD', 'MO_LENGTH',
                                 'V_SIGMA', 'PBLH')
      cdata <<- data %>% 
        select(any_of(contribution_variables))
      
      # Assign variables to class
      epsg <<- epsg
      xy <<- xy
      tower_height <<- tower_height
      xy_sf <- st_point(xy) %>% 
        st_sfc(crs=epsg)
      lonlat <- st_transform(xy_sf, crs=4326) %>% 
        st_coordinates()
      shapefile <<- read_sf(shapefile) %>% 
        st_transform(epsg)
      
      kljun_params <<- kljun_params
      
      # If the PBL height is not supplied then estimate it
      # Kljun et al 2015 appendix B talks about the approaches and caveats
      # Try and use ceilometer or other observed data where available
      if(!'PBLH' %in% colnames(cdata)) {
        # Hanna and Chang 1993 approach
        latitude = lonlat[2] * (pi / 180)        # Degrees -> radians
        angular_velocity = 7.2921159 * 10^-5     # rad/s
        coriolis_parameter = 2 * angular_velocity * sin(latitude)
        # coefficient for calculating planetary boundary layer height
        c_n = 0.3
        # planetary boundary layer height
        cdata$PBLH <<- c_n * cdata$USTAR / coriolis_parameter
        
        # Niuewstadt 1981
        # # angular velocity (rad/s)
        # omega = 7.2921150 * 10^-5 
        # # lat (rad)
        # phi = lonlat[2] * (pi / 180)
        # # Coriolis parameter
        # f <- 2 * omega * sin(phi)
        # # PBLH
        # cdata$PBLH <<- cdata$MO_LENGTH/3.8 * ( -1 + ( 1 + 2.28*(cdata$USTAR/(f*cdata$MO_LENGTH)) )^0.5 )
        # 
      }
    },
    run_kljun = function(row, verbose=FALSE)
    {
      # Run the Kljun footprint model. Code has been sourced from the kljun
      # website: https://footprint.kljun.net
      # Inputs
      #   row <vec>: row of data / timeseries slice from the cdata df
      #   verbosity <bool>: flag to suppress function output
      # Returns
      #   <list> output of FPP climatology function
      # TODO: the verbosity flag won't work on windows
      if(!verbose) sink('/dev/null')
      ffp_clim = calc_footprint_FFP_climatology(
        zm = (tower_height - row[['DISPLACEMENT_HEIGHT']]),
        z0 = row[['ROUGHNESS_LENGTH']],
        # rs = NaN,
        umean = row[['WS']],
        h = row[['PBLH']],
        ol = row[['MO_LENGTH']],
        sigmav = row[['V_SIGMA']],
        ustar = row[['USTAR']],
        wind_dir = row[['WD']],
        domain = kljun_params$domain,
        nx = kljun_params$nx-1,
        ny = kljun_params$ny-1,
        rslayer = 1
      )
      if(!verbose) sink()
      return(ffp_clim)
    },
    make_raster = function(row)
    {
      # Internal function to create a raster from the Kljun footprint matrix
      # Inputs
      #   row <vec>: row of data / timeseries slice from the cdata df
      # Returns
      #   <raster> footprint raster
      m <- matrix(NaN, kljun_params$ny, kljun_params$nx)
      if( !any(is.na(row)) ) {
        ffp_clim = run_kljun(row)
        if(ffp_clim$flag_err != 1) {
          m <- ffp_clim$fclim_2d
        }
      }
      suppressWarnings({
        # Raster expects left to right, top to bottom
        r <- raster::raster(t(m[,ncol(m):1]),
                            xmn = xy[1] + kljun_params$domain[1] - 0.5,
                            xmx = xy[1] + kljun_params$domain[2] + 0.5,
                            ymn = xy[2] + kljun_params$domain[3] - 0.5,
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
        
        fp <- apply(cdata %>% select(-datetime),
                    MARGIN = 1,
                    simplify = FALSE,
                    make_raster)
      } else {
        message('processing in parallel')
        tryCatch({
          registerDoParallel(cores=parallel)
          fp <- foreach(row=iterators::iter(cdata %>% select(-datetime), by='row')) %dopar% 
            {
              make_raster(row)
            }
        },
        error = function(e) {
          stopImplicitCluster()
          message(e)
        },
        finally = {
          stopImplicitCluster()
        })
      }
      toc()
      footprints <<- brick(fp)
      names(footprints) <<- format_ISO8601(cdata$datetime)
    },
    calculate_contribution_dep = function(keep_footprints = FALSE, parallel = NULL)
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
        contribution <- unname( t(ee_sum)[,1] / r_sum )
      } else {
        if(is.null(parallel)) {
          # Not parallel
          contribution <- apply(cdata %>% select(-datetime), 1, function(row) {
            r <- make_raster(row)
            r_sum = cellStats(r, sum)
            ext_sum = exactextractr::exact_extract(r, shapefile, fun='sum')
            cont = ext_sum/r_sum
            return(cont)
          })
        } else {
          # Parallel
          tic()
          tryCatch({
            message('starting parallel processing')
            registerDoParallel(cores=parallel)
            contribution <- foreach(row=iterators::iter(cdata, by='row'),
                                    .combine = c) %dopar% 
              {
                r <- make_raster(row)
                r_sum = raster::cellStats(r, sum)
                ext_sum = exactextractr::exact_extract(r, shapefile, fun='sum')
                cont = ext_sum/r_sum
              }
          },
          error = function(e) {
            stopImplicitCluster()
            message(e)
          },
          finally = {
            stopImplicitCluster()
          })
        }
      }
      toc()
      cdata$Contribution <<- contribution
    },
    calculate_contribution = function(keep_footprints = FALSE,
                                      parallel = NULL)
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
        contribution <- unname( t(ee_sum)[,1] / r_sum )
        r_mean <- mean(footprints[[which(!is.na(contribution))]])
      } else {
        if(is.null(parallel)) {
          # Not parallel
          contribution <- vector('numeric', length=nrow(cdata))
          # Create raster for average footprint
          r_mean <- raster::raster(matrix(0, kljun_params$ny, kljun_params$nx),
                                   xmn = xy[1] + kljun_params$domain[1] - 0.5,
                                   xmx = xy[1] + kljun_params$domain[2] + 0.5,
                                   ymn = xy[2] + kljun_params$domain[3] - 0.5,
                                   ymx = xy[2] + kljun_params$domain[4] + 0.5,
                                   crs = epsg)
          # Counter for valid footprints
          n <- 0
          for(i in seq_len(nrow(cdata))) {
            r <- make_raster(cdata[i,])
            r_sum = cellStats(r, sum)
            if(r_sum > 0) {
              ext_sum = exactextractr::exact_extract(r, shapefile, fun='sum')
              contribution[i] = ext_sum/r_sum
              n <- n + 1
              r_mean <- r_mean + r
            } else {
              contribution[i] <- NaN
            }
          }
          r_mean <- r_mean / n
        } else {
          # Parallel
          tryCatch({
            message('starting parallel processing')
            registerDoParallel(cores=parallel)
            # Split data into chunks
            par_data <- cdata %>% 
              mutate(group = rep(1:parallel, each=ceiling(nrow(.)/parallel), length.out=nrow(.))) %>% 
              group_split(group, .keep=FALSE)
            par_data <- foreach(cdata=par_data) %dopar% 
            {
              contribution <- vector('numeric', length=nrow(cdata))
              # Create raster for average footprint
              r_mean <- raster::raster(matrix(0, kljun_params$ny, kljun_params$nx),
                                       xmn = xy[1] + kljun_params$domain[1] - 0.5,
                                       xmx = xy[1] + kljun_params$domain[2] + 0.5,
                                       ymn = xy[2] + kljun_params$domain[3] - 0.5,
                                       ymx = xy[2] + kljun_params$domain[4] + 0.5,
                                       crs = epsg)
              for(i in seq_len(nrow(cdata))) {
                r <- make_raster(cdata[i,])
                r_sum = cellStats(r, sum)
                if(r_sum > 0) {
                  ext_sum = exactextractr::exact_extract(r, shapefile, fun='sum')
                  contribution[i] = ext_sum/r_sum
                  r_mean <- r_mean + r
                } else {
                  contribution[i] <- NaN
                }
              }
              return(list(contribution = contribution,
                          r_mean = r_mean))
              }
            contribution <- unlist(lapply(par_data, function(l) l$contribution))
            r_mean <- raster::brick(lapply(par_data, function(l) l$r_mean))
            r_mean <- sum(r_mean) / length(contribution[!is.na(contribution)])
          },
          error = function(e) {
            stopImplicitCluster()
            message(e)
          },
          finally = {
            stopImplicitCluster()
          })
        }
      }
      cdata$Contribution <<- contribution
      mean_footprint <<- r_mean
    },
    plot_mean_footprint = function() {
      
      # Make a tbl of the mean footprint raster
      r_tbl <- as_tibble(raster::rasterToPoints(mean_footprint)) %>% 
        arrange(layer) %>% 
        mutate(level = cumsum(layer))
      
      f_sort <- r_tbl$layer
      f_cs <- r_tbl$level
      rs <- seq(0.1, 0.8, 0.1)
      fr <- sapply(rs, function(r) {
        diff <- abs(f_cs - r)
        idx <- which.min(diff)
        return(f_sort[idx])
      })
      
      xy_sf <- st_point(xy) %>% 
        st_sfc(crs=epsg)
      
      ggplot() +
        geom_sf(data=shapefile) +
        geom_contour_filled(data=r_tbl, aes(x, y, z=layer), breaks=c(fr,Inf), alpha=0.5) +
        geom_sf(data=xy_sf, aes(col='EC tower'), size=1) +
        coord_sf(datum=epsg) +
        scale_fill_manual(values=viridis::viridis(length(rs)),
                          labels=rev(rs)*100) +
        guides(fill = guide_legend(reverse = TRUE)) +
        labs(fill='FP dist %')
    }
  )
)

