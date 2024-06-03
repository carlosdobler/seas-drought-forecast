
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")



pts <- 
  seq(70, -70, -20) %>% 
  map(function(lat) {
    
    m <- matrix(c(-180,lat,180,lat), ncol = 2, byrow = T)
    l <- st_linestring(m) %>% st_sfc(crs = 4326)
    ln <- units::drop_units(st_length(l) * 1e12)
    n <- round(ln/120)
    
    l_proj <- st_transform(l, 3857)
    pt <- st_line_sample(l_proj, n = n)  
    st_transform(pt, 4326) %>% st_cast("POINT")
    
  })

pts <- 
  do.call(c, pts) %>% 
  st_as_sf() %>% 
  mutate(id = row_number())





# 
d <- "/mnt/pers_disk/data_forecast/"
fs::dir_create(d)



# load
# s_proxy <- 
#   d %>% 
#   fs::dir_ls(regexp = "geopotential-500") %>% 
#   first() %>% 
#   read_ncdf(ncsub = cbind(start = c(1, 1, 1),
#                           count = c(NA,NA,1))) %>% 
#   adrop()




vars <- c(`mean-sst` = "sst",
          `geopotential-200` = "z200",
          `geopotential-500` = "z500",
          `geopotential-850` = "z850",
          `specific-humidity-200` = "sh200",
          `specific-humidity-850` = "sh850",
          `uwind-200` = "uwind200",
          `uwind-850` = "uwind850",
          `vwind-200` = "vwind200",
          `vwind-850` = "vwind850"
          )



tb <- 
  imap_dfc(vars, function(var, var_d) {
    
    
    # DOWNLOAD
    "/mnt/bucket_mine/era/monthly/" %>% 
      fs::dir_ls(regexp = var_d) %>% 
      fs::dir_ls() %>% 
      str_subset(str_flatten(1940:2023, "|")) %>% 
      future_walk(function(f) {
        
        f_gs <- str_replace(f, "/mnt/bucket_mine", "gs://clim_data_reg_useast1")
        
        str_glue("gsutil cp {f_gs} {d}") %>% 
          system(ignore.stderr = T, ignore.stdout = T)
        
      })
    
    
    # LOAD DATA
    ff <- 
      d %>% 
      fs::dir_ls()
    
    pts_data <- 
      future_map_dfc(seq(nrow(pts)), function(r) {
        
        cell_x <- fn_get_cell_pos(s_proxy, 1, pts[r,] %>% st_coordinates() %>% .[,1] %>% {.+360})
        cell_y <- fn_get_cell_pos(s_proxy, 2, pts[r,] %>% st_coordinates() %>% .[,2])
        
        col_name <- str_glue("loc{str_pad(r, 2, 'left', '0')}_{var}")
        
        v <- 
          ff %>%  
          map(function(f) {
            
            f %>% 
              read_ncdf(ncsub = cbind(start = c(cell_x, cell_y, 1),
                                      count = c(1,      1,      NA))) %>% 
              setNames(var) %>% 
              suppressMessages() %>% 
              pull() %>%
              as.vector()
            
          }) %>% 
          do.call(c, .) %>% 
          unname() 
        
        if(all(!is.na(v))) {
          
          tibble({{col_name}} := v)
          
        }
        
      })
    
    ff %>% 
      fs::file_delete()
    
    return(pts_data)
    
  })



tibble(time = seq(as_date("1940-01-01"), as_date("2023-12-31"), by = "1 month"),
       tb) %>% 
  
  write_csv("experiment_02/predictors_data.csv")


  


