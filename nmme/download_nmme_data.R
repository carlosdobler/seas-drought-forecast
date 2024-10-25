
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore, workers = 5)

dir_data <- "/mnt/pers_disk/tmp"
fs::dir_create(dir_data)



dates <- seq(as_date("1991-01-01"), as_date("2020-12-01"), by = "1 month")
model <- "nasa-geoss2s"  #"gem5p2-nemo"

future_walk(dates, \(d){
  
  # print(d)
  
  # url <- 
  #   str_glue(
  #     "http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.GFDL-SPEAR/.HINDCAST/.MONTHLY/.prec/L/%280.5%29%285.5%29RANGEEDGES/S/%280000%201%20{format(d, '%b')}%20{year(d)}%29VALUES/M/%281.0%29%2810.0%29RANGEEDGES/data.nc"
  #   )
  
  # url <- 
  #   str_glue(
  #     "http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CanSIPS-IC4/.CanESM5/.HINDCAST/.MONTHLY/.prec/L/%280.5%29%285.5%29RANGEEDGES/S/%280000%201%20{format(d, '%b')}%20{year(d)}%29VALUES/M/%281.0%29%2810.0%29RANGEEDGES/data.nc"
  #   )
  
  # url <- 
  #   str_glue(
  #     "http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CanSIPS-IC4/.GEM5.2-NEMO/.HINDCAST/.MONTHLY/.prec/L/%280.5%29%285.5%29RANGEEDGES/S/%280000%201%20{format(d, '%b')}%20{year(d)}%29VALUES/M/%281.0%29%2810.0%29RANGEEDGES/data.nc"
  #   )
  
  url <- 
    str_glue(
      "http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NASA-GEOSS2S/.HINDCAST/.MONTHLY/.prec/L/%280.5%29%285.5%29RANGEEDGES/S/%280000%201%20{format(d, '%b')}%20{year(d)}%29VALUES/data.nc"
    )
  
  f <-
    str_glue("{dir_data}/nmme_{model}_precipitation_mon_{d}_plus5.nc")
  
  
  download.file(url, f, method = "wget", quiet = T)
  
  str_glue("gsutil mv {f} gs://clim_data_reg_useast1/nmme/{model}/precipitation/{fs::path_file(f)}") %>% 
    system(ignore.stdout = T, ignore.stderr = T)
  
})




