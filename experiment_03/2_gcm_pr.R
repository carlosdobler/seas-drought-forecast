
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")

gcm <- "CNRM-CERFACS/CNRM-CM6-1"
mems <- str_glue("r{seq(1,30,5)}i1p1f2")
var <- "pr"
prov <- "Amon"
grid <- "gr"


# TIME

# root path
r <- 
  str_glue("gcloud storage ls gs://cmip6/CMIP6/CMIP/{gcm}/historical/{mems[1]}/{prov}/{var}/{grid}/") %>% 
  system(intern = T) %>% 
  str_sub(start = 6, end = -2)

dsn <- 
  str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:{var}.zarr/')

# format time
ret <-
  gdal_utils("mdiminfo", dsn, quiet = TRUE)

time_array <- 
  jsonlite::fromJSON(ret)$arrays$time

time_vector <- 
  str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:/time') %>% 
  read_stars() %>% 
  pull() %>% 
  as.vector() %>% 
  {. * 24 * 60 * 60} %>% # pcict works in seconds
  PCICt::as.PCICt(cal = time_array$attributes$calendar,
                  origin = time_array$attributes$time_origin)

time_i <- 
  first(which(year(time_vector) == 1940))

time_vector_sub <- 
  time_vector %>% 
  tail(-time_i+1) %>% 
  str_sub(end = 8) %>% 
  paste0("01") %>% 
  as_date()


# proxy 
s_proxy <- 
  read_mdim(dsn, 
            count = c(NA,NA,1)) %>% 
  adrop()

cell_x <- fn_get_cell_pos(s_proxy, 1, 25)
cell_y <- fn_get_cell_pos(s_proxy, 2, 0.5)

# cell_x <- fn_get_cell_pos(s_proxy, 1, 14)
# cell_y <- fn_get_cell_pos(s_proxy, 2, -4)
# 
# cell_x <- fn_get_cell_pos(s_proxy, 1, 26)
# cell_y <- fn_get_cell_pos(s_proxy, 2, -2)


tb_pr <- 
  map_dfr(mems, function(mem){
    
    print(str_glue("Processing member {mem}"))
    
    r <- 
      str_glue("gcloud storage ls gs://cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/{prov}/{var}/{grid}/") %>% 
      system(intern = T) %>% 
      str_sub(start = 6, end = -2)
    
    dsn <- 
      str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:{var}.zarr/')
    
    
    
    
    s <- 
      read_mdim(dsn, 
                offset = c(cell_x-1, cell_y-1, time_i-1),
                count = c(1,1,NA))
    
    tb_var <-
      s %>% 
      as_tibble() %>% 
      mutate(pr = pr %>% units::set_units(kg/m^2/d)) %>% 
      units::drop_units() %>% 
      
      # filter(month(time) %in% c(9,10,11)) %>% 
      # group_by(time = year(time)) %>% 
      # summarise(pr = mean(pr)) %>% 
      mutate(time = paste0(str_sub(time, end = 8), "01") %>% as_date(),
             mem = mem %>% as.character()) %>% 
      select(-1,-2)
    
    
    return(tb_var)
    
  })



tb_pr %>% 
  write_csv("experiment_03/data/pr_1loc_gcm.csv")





