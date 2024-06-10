
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")

d <- "/mnt/pers_disk/data_forecast/"


# Download precip

fs::dir_create(d)

"/mnt/bucket_mine/era/monthly/mean-daily-precip/" %>% 
  fs::dir_ls() %>% 
  str_subset(str_flatten(1940:2023, "|")) %>% 
  future_walk(function(f) {
    
    f_gs <- str_replace(f, "/mnt/bucket_mine", "gs://clim_data_reg_useast1")
    
    str_glue("gsutil cp {f_gs} {d}") %>% 
      system(ignore.stderr = T, ignore.stdout = T)
    
  })


# Select 1 location

s_proxy <- 
  "/mnt/bucket_mine/era/monthly/mean-daily-precip/era5_mon_mean-daily-precip_1940.nc" %>%  
  read_ncdf(ncsub = cbind(start = c(1, 1, 1),
                          count = c(NA,NA,1))) %>% 
  adrop()


cell_x <- fn_get_cell_pos(s_proxy, 1, 25)
cell_y <- fn_get_cell_pos(s_proxy, 2, 0.5)

cell_x <- fn_get_cell_pos(s_proxy, 1, 14)
cell_y <- fn_get_cell_pos(s_proxy, 2, -4)

cell_x <- fn_get_cell_pos(s_proxy, 1, 26)
cell_y <- fn_get_cell_pos(s_proxy, 2, -2)



# Load precip data

v <- 
  d %>% 
  fs::dir_ls() %>% 
  map_dfr(function(f) {
    
    f %>% 
      read_ncdf(ncsub = cbind(start = c(cell_x, cell_y, 1),
                              count = c(1,      1,      NA))) %>% 
      suppressMessages() %>%
      as_tibble() %>% 
      mutate(time = as_date(time),
             tp = tp %>% units::set_units(mm))
    
  })


vf <- 
  v %>% 
  units::drop_units() %>% 
  rename(ds = time,
         y = tp) %>% 
  select(ds, y)


# vf %>% 
#   write_csv("experiment_04/response_data.csv")

# d %>% 
#   fs::dir_ls() %>% 
#   fs::file_delete()




library(prophet)

# m <- prophet(vf, changepoint.prior.scale = 1, yearly.seasonality = 20)


vf %>% 
  mutate(is_01 = if_else(month(ds) == 1, 1, 0),
         is_02 = if_else(month(ds) == 2, 1, 0),
         is_03 = if_else(month(ds) == 3, 1, 0),
         is_04 = if_else(month(ds) == 4, 1, 0),
         is_05 = if_else(month(ds) == 5, 1, 0),
         is_06 = if_else(month(ds) == 6, 1, 0),
         is_07 = if_else(month(ds) == 7, 1, 0),
         is_08 = if_else(month(ds) == 8, 1, 0),
         is_09 = if_else(month(ds) == 9, 1, 0),
         is_10 = if_else(month(ds) == 10, 1, 0),
         is_11 = if_else(month(ds) == 11, 1, 0),
         is_12 = if_else(month(ds) == 12, 1, 0)) -> vf2

m <- prophet(yearly.seasonality=F)
m <- add_regressor(m, "is_01")
m <- add_regressor(m, "is_02")
m <- add_regressor(m, "is_03")
m <- add_regressor(m, "is_04")
m <- add_regressor(m, "is_05")
m <- add_regressor(m, "is_06")
m <- add_regressor(m, "is_07")
m <- add_regressor(m, "is_08")
m <- add_regressor(m, "is_09")
m <- add_regressor(m, "is_10")
m <- add_regressor(m, "is_11")
m <- add_regressor(m, "is_12")

m <- fit.prophet(m, vf2)

m <- prophet()
m <- fit.prophet(m, vf)


# future <- make_future_dataframe(m, periods = 12, freq = 'month')
fcst <- predict(m, vf)

fcst %>% 
  ggplot(aes(x = as_date(ds))) +
  
  geom_line(data = vf, aes(y = y), color = "red") +
  geom_line(aes(y = yhat))

plot(m, fcst)
prophet_plot_components(m, fcst)











