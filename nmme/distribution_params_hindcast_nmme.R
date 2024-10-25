
library(tidyverse)
library(stars)
library(furrr)
library(ncdf4)

options(future.fork.enable = T)
plan(multicore)

source("https://raw.github.com/carlosdobler/spatial-routines/master/general_tools.R")


dir_data <- "/mnt/pers_disk/tmp"


# model <- "gfdl-spear"
# model <- "canesm5"
model <- "gem5p2-nemo"


ff <- 
  "gsutil ls gs://clim_data_reg_useast1/nmme/monthly/{model}/precipitation/" %>%
  str_glue() %>% 
  system(intern = T) %>% 
  str_subset(str_flatten(1991:2020, "|"))


# DOWNLOAD DATA

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
      future_map(\(f){
        
        nc <- nc_open(f)
        dim_names <- nc$var$prec$dim %>% map_chr(pluck, "name")
        ic <- nc$var$prec$dim[[which(dim_names == "S")]]
        s_val <- ic$units %>% str_sub(14) %>% as_date() %>% {. + months(ic$vals)}
        
        s <- read_ncdf(f, proxy = F)
        names(st_dimensions(s)) <- dim_names
        s <-
          s %>%
          st_set_dimensions("S", values = year(s_val)) %>% 
          st_set_dimensions("L", values = seq(s_val, s_val + months(5), by = "1 month") %>% month()) %>% 
          st_set_dimensions("M", values = 1:10)
        
        
      }) %>% 
      suppressMessages()
    
    # concatenate and convert units
    s <- 
      do.call(c, c(s, along = "S"))
    
    # fit distribution per grid cell
    # loop through leads
    
    leads <- st_get_dimension_values(s, "L")
    
    seq_along(leads) %>% 
      walk(\(lead_in){
        
        print(str_glue("   PROCESSING LEAD {lead_in}"))
        
        s_gamma <- 
          s %>% 
          slice(L, lead_in) %>% 
          st_apply(c(1,2), \(x){
            
            xv <- as.vector(x)
            
            if(all(is.na(xv))) {
              c(NA,NA)
            } else if (all(xv == 0)) {
              c(-9999,-9999)
            } else {
              
              xv[xv == 0] <- 0.001 # avoid errors in gamma fitting
              
              # fit gamma
              p <- MASS::fitdistr(xv, "gamma")
              
              # shape and rate
              c(p$estimate[1], p$estimate[2])
              
            }
            
          },
          .fname = "params",
          FUTURE = T) %>% 
          aperm(c(2,3,1))
        
        # format and export
        
        r <- str_glue("nmme_{model}_precipitation_mon_gamma-params_1991-2020_{m}_lead-{lead_in-1}.nc")
        
        s_gamma %>% 
          split("params") %>% 
          setNames(c("shape", "rate")) %>% 
          rt_write_nc(str_glue("{dir_data}/{r}"))
        
        
        str_glue("gsutil mv {dir_data}/{r} gs://clim_data_reg_useast1/nmme/climatologies/{model}/") %>% 
          system()
        
        
      })
    
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

