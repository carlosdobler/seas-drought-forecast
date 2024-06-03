
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")

d <- "/mnt/pers_disk/data_forecast/"
fs::dir_create(d)





"/mnt/bucket_mine/era/monthly/mean-daily-precip/" %>% 
  fs::dir_ls() %>% 
  str_subset(str_flatten(1940:2023, "|")) %>% 
  future_walk(function(f) {
    
    f_gs <- str_replace(f, "/mnt/bucket_mine", "gs://clim_data_reg_useast1")
    
    str_glue("gsutil cp {f_gs} {d}") %>% 
      system(ignore.stderr = T, ignore.stdout = T)
    
  })



s_proxy <- 
  "/mnt/bucket_mine/era/monthly/mean-daily-precip/era5_mon_mean-daily-precip_1940.nc" %>%  
  read_ncdf(ncsub = cbind(start = c(1, 1, 1),
                          count = c(NA,NA,1))) %>% 
  adrop()

cell_x <- fn_get_cell_pos(s_proxy, 1, 25)
cell_y <- fn_get_cell_pos(s_proxy, 2, 0.5)




v <- 
  d %>% 
  fs::dir_ls() %>% 
  map_dfr(function(f) {
    
    f %>% 
      read_ncdf(ncsub = cbind(start = c(cell_x, cell_y, 1),
                              count = c(1,      1,      NA))) %>% 
      suppressMessages() %>%
      as_tibble() %>% 
      mutate(tp = tp %>% units::set_units(mm),
             time = as_date(time))
    
  })


vf <- 
  v %>% 
  mutate(tp = tp * days_in_month(time) %>% round()) %>% 
  group_by(mon = month(time)) %>% 
  mutate(tp_anom = tp - mean(tp),
         tp_anom_std = tp_anom/sd(tp)) %>% 
  ungroup() %>% 
  select(time, tp, tp_anom, tp_anom_std) %>% 
  units::drop_units()


vf %>% 
  write_csv("experiment_02/response_data.csv")

d %>% 
  fs::dir_ls() %>% 
  fs::file_delete()









