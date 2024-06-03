library(tidyverse)
library(stars)

gcm <- "CNRM-CERFACS/CNRM-CM6-1"
mem <- "r1i1p1f2"
var <- "pr"


# Z500

root_dsn <- str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/Amon/{var}/gr/v20180917"/')

# Pressure level vector
dsn = str_glue('{root_dsn}:/plev')
plev_vector <- read_stars(dsn) %>% pull() %>% as.vector()

# Time vector
dsn = str_glue('{root_dsn}:/time')
a <- gdal_utils("info", dsn, quiet = T)
o <- a %>% str_extract("time_origin.*\n") %>% str_extract("....-..-..")
time_vector <- read_stars(dsn) %>% pull() %>% as.vector() %>% as_date(origin = o)

# Data
dsn = str_glue('{root_dsn}:zg.zarr/')
d = read_mdim(dsn, 
              offset = c(0,0,
                         which(plev_vector == 50000)-1,
                         first(which(year(time_vector) == 1940))-1), 
              count = c(NA,NA,
                        1,
                        NA)) %>% 
  adrop()




# PRECIP

root_dsn <- str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/Amon/{var}/gr/v20180917"/')

# Time vector
dsn = str_glue('{root_dsn}:/time')
a <- gdal_utils("info", dsn, quiet = T)
o <- a %>% str_extract("time_origin.*\n") %>% str_extract("....-..-..")
time_vector <- read_stars(dsn) %>% pull() %>% as.vector() %>% as_date(origin = o)

# Data
dsn = str_glue('{root_dsn}:zg.zarr/')
y <- first(which(year(time_vector) == 1940))-1
d = read_mdim(dsn, 
              offset = c(0, 0, y), 
              count = c(NA, NA, NA)) %>% 
  adrop()

d <- 
  d %>% 
  mutate(pr = 
           pr %>% 
           units::set_units(kg/m^2/d) %>% 
           units::set_units(NULL) %>% 
           {. * days_in_month(tail(time_vector, -y))})
 
 


