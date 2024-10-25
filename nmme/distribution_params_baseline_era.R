
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

source("https://raw.github.com/carlosdobler/spatial-routines/master/general_tools.R")


dir_data <- "/mnt/pers_disk/tmp"


ff <- 
  "gsutil ls gs://clim_data_reg_useast1/era5/monthly_means/total_precipitation/" %>% 
  system(intern = T) %>% 
  str_subset(str_flatten(1991:2020, "|"))


# DOWNLOAD ERA DATA

ff %>% 
  future_walk(\(f){
    
    str_glue("gsutil cp {f} {dir_data}") %>% 
      system(ignore.stderr = T, ignore.stdout = T)
    
  })

ff <- 
  ff %>% 
  fs::path_file() %>% 
  {str_glue("{dir_data}/{.}")}




# CLIMATOLOGY BASED ON DISTRIBUTION

str_pad(seq(12), 2, "left", "0") %>% 
  walk(\(m){
    
    print(str_glue("PROCESSING MONTH {m}"))
    
    # month files
    ff_m <- 
      ff %>% 
      str_subset(str_glue("-{m}-"))
    
    # read data
    s <- 
      ff_m %>% 
      future_map(read_ncdf) %>% 
      suppressMessages()
    
    # concatenate and convert units
    s <- 
      do.call(c, c(s, along = "time")) %>% 
      mutate(tp = tp %>% units::set_units(mm))
    
    # fit distribution per grid cell
    s_gamma <- 
      s %>% 
      st_apply(c(1,2), \(x){
        
        if (all(x == 0)) {
          c(-9999,-9999)
        } else {
          
          x[x == 0] <- 0.001 # avoid errors in gamma fitting
          
          # fit gamma
          p <- MASS::fitdistr(x, "gamma")
          
          # shape and rate
          c(p$estimate[1], p$estimate[2])
          
        }
        
      },
      .fname = "params",
      FUTURE = T) %>% 
      aperm(c(2,3,1))
    
    # format and export
    s_gamma %>% 
      split("params") %>% 
      setNames(c("shape", "rate")) %>% 
      rt_write_nc(str_glue("{dir_data}/era5_total-precipitation_mon_gamma-params_1991-2020_{m}.nc"))
    
  })

# transfer to bucket
dir_data %>% 
  fs::dir_ls() %>% 
  str_subset("gamma-params") %>% 
  walk(\(f){
    
    str_glue("gsutil cp {f} gs://clim_data_reg_useast1/era5/climatologies/") %>% 
      system()
    
  })


# clean up
dir_data %>% 
  fs::dir_ls() %>% 
  str_subset("gamma-params", negate = T) %>% 
  walk(fs::file_delete)







# TEST DISTRIBUTION

# # random grid cell
# s %>% units::drop_units() %>% pull() %>% .[sample(1440,1),sample(721,1),] -> x
# x[x == 0] <- 0.001 # avoid errors
# 
# # fit distribution
# p <- MASS::fitdistr(x, "gamma")
# 
# # observed quantiles
# xq <- quantile(x, seq(0,0.99,0.01)) %>% unname()
# # theoretical quantiles
# pq <- qgamma(seq(0,0.99,0.01), p$estimate[1], p$estimate[2])
# 
# # q-q plot
# tibble(obs = xq, fitted = pq) %>%
#   ggplot(aes(obs, fitted)) +
#   geom_point() +
#   coord_equal() +
#   geom_abline()
# 
# # random sequence based on distr param
# r <- rgamma(length(x), p$estimate[1], p$estimate[2])
# 
# tibble(obs = sort(x), fitted = sort(r)) %>%
#   ggplot(aes(obs, fitted)) +
#   geom_point() +
#   coord_equal() +
#   geom_abline()

