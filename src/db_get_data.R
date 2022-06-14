#------------------------------------------------------------------------------#
# Retrieve data from the NOBV database
#------------------------------------------------------------------------------#

library(dbplyr)
library(dplyr)
library(stringr)

get_db_data <- function(site,
                        plot,
                        sensor_group,
                        sensor_type,
                        start_dt,
                        end_dt) {
  # Function to retrieve data from the NOBV database
  # Params
  #   site <character>: abbreviated name of the site, e.g. ZEG
  #   plot <character>: abbreviated name of the plot, e.g. RF
  #   sensor_group <character>: name of the sensor group, e.g. MT1
  #   sensor_type <character>: name of the sensor type, e.g. WIND. Leave blank 
  #   to return all sensors
  #   start_dt <datetime>: timestamp when to start extraction
  #   end_dt <datetime>: timestamp when to end extraction
  # Returns
  #   <tbl> timeseries of data and metadata
  
  # Handle case where all sensors in a group want to be returned
  if(missing(sensor_type)) sensor_type = ''
  # Handle case if multiple sensor group/types are provided
  if(length(sensor_group) > 1) sensor_group <- paste(sensor_group, collapse='|')
  if(length(sensor_type) > 1) sensor_type <- paste(sensor_type, collapse='|')

  tryCatch(
    {
      # Make database connection
      con <- DBI::dbConnect(RPostgreSQL::PostgreSQL(),
                            dbname = db_name,
                            host = db_host,
                            port = db_port,
                            user = db_user,
                            password = db_password)
      
      # Get the site id
      site_id = tbl(con, in_schema("cdr", "sites")) %>% 
        filter(shortname == !!paste(site, plot, sep='_')) %>% 
        collect()
      
      # Filter to the sensor we want to get
      sensor_id = tbl(con, in_schema('cdr', 'logvalproviders')) %>% 
        filter(site == !!site_id$id,
               str_detect(name, sensor_group) &
               str_detect(name, sensor_type) &
               !str_detect(name, 'Std|mv')) %>% 
        collect()

      # The units of the sensor
      sensor_units = tbl(con, in_schema('cdr', 'units')) %>% 
        filter(id %in% !!sensor_id$unit) %>% 
        collect()
      
      # Retrieve the actual data
      sensor_data = tbl(con, in_schema('cdr', 'pointdata')) %>% 
        filter(logicid %in% !!sensor_id$id,
               dt >= start_dt,
               dt <= end_dt) %>% 
        select(id = logicid, dt, value) %>% 
        mutate(dt = timezone('UTC', dt)) %>% 
        collect()
      
      # Add the sensor info to the tbl
      sensor_data = left_join(sensor_id %>% select(id, name, unit_id = unit),
                      sensor_units %>% select(unit_id = id, unit = abbreviation)) %>% 
              right_join(sensor_data) %>% 
        select(name, dt, value, unit)
    },
    error = function(e) {
      message(paste("ERROR:", e))
      return(NULL)
    },
    finally = {
      # Safely disconnect
      DBI::dbDisconnect(con)
    }
  )
  return(sensor_data)
}
