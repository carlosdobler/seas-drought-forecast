
start <- lubridate::as_date("2024-09-01")



# ***************

library(tidyverse)
library(stars)
library(furrr)
library(ncdf4)

options(future.fork.enable = T)
plan(multicore)

source("https://raw.github.com/carlosdobler/spatial-routines/master/general_tools.R")

dir_data <- "/mnt/pers_disk/tmp"
dir_gs_era <- "gs://clim_data_reg_useast1/era5/climatologies"
dir_gs_nmme <- "gs://clim_data_reg_useast1/nmme/climatologies"


dates_fcst <- seq(start, start + months(5), by = "1 month")


# ***********


# download ERA parameters
ff_era <- 
  str_glue("gsutil ls {dir_gs_era}") %>% 
  system(intern = T) %>% 
  str_subset("gamma-params") %>% 
  str_subset("1991-2020") %>% 
  str_subset(str_flatten(str_glue("_{str_pad(month(dates_fcst), 2, 'left', '0')}"), "|"))

# rearrange
ff_era <- 
  str_pad(month(dates_fcst), 2, 'left', '0') %>% 
  map_chr(\(m){
    
    ff_era %>% 
      str_subset(str_glue("_{m}"))
    
  })
  
ff_era %>% 
  future_walk(\(f) {
    
    str_glue("gsutil cp {f} {dir_data}") %>% 
      system(ignore.stdout = T, ignore.stderr = T)
    
  })

ff_era <- 
  ff_era %>% 
  fs::path_file() %>% 
  {str_glue("{dir_data}/{.}")}


# load era

# aoi
ncs_era <- 
  ff_era %>% 
  first() %>% 
  read_ncdf(proxy = T) %>% 
  suppressMessages() %>% 
  rt_from_coord_to_ind(7,-14,34,9)

era_params <- 
  ff_era %>%
  set_names(ff_era %>% str_sub(-5,-4)) %>% 
  map(read_ncdf, ncsub = cbind(start = c(ncs_era$x_start, ncs_era$y_start),
                               count = c(ncs_era$x_count, ncs_era$y_count))) %>% 
  suppressMessages() %>% 
  {do.call(c, c(., along = "L"))} %>% 
  st_set_dimensions("L", values = dates_fcst)
  




# ***********


models <- c("canesm5",
            "gfdl-spear",
            "gem5p2-nemo")

urls <- c(str_glue("http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CanSIPS-IC4/.CanESM5/.FORECAST/.MONTHLY/.prec/L/%280.5%29%285.5%29RANGEEDGES/S/%280000%201%20{format(start, '%b')}%20{year(start)}%29VALUES/M/%281.0%29%2810.0%29RANGEEDGES/data.nc"),
          str_glue("http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.GFDL-SPEAR/.FORECAST/.MONTHLY/.prec/M/%281.0%29%2810.0%29RANGEEDGES/L/%280.5%29%285.5%29RANGEEDGES/S/%280000%201%20{format(start, '%b')}%20{year(start)}%29VALUES/data.nc"),
          str_glue("http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CanSIPS-IC4/.GEM5.2-NEMO/.FORECAST/.MONTHLY/.prec/L/%280.5%29%285.5%29RANGEEDGES/S/%280000%201%20{format(start, '%b')}%20{year(start)}%29VALUES/M/%281.0%29%2810.0%29RANGEEDGES/data.nc"))



# loop through models
rr <- 
  
  map2(models, urls, \(model, url){
    
    print(model)
    
    # model = models[2]
    # url = urls[2]
    
    
    # download forecast
    f <- str_glue("{dir_data}/nmme_{model}_{start}.nc")
    download.file(url, method = "wget", destfile = f, quiet = T)
    
    # read_forecast
    # confirm start date is correct
    nc <- nc_open(f)
    dim_names <- nc$var$prec$dim %>% map_chr(pluck, "name")
    ic <- nc$var$prec$dim[[which(dim_names == "S")]]
    s_val <- ic$units %>% str_sub(14) %>% as_date() %>% {. + months(ic$vals)}
    if(!s_val == start){
      print(str_glue("wrong dates!"))
    }
    
    # read data and format
    fcst <- 
      read_ncdf(f, proxy = F) %>% 
      suppressMessages() %>% 
      suppressWarnings()
    
    names(st_dimensions(fcst)) <- dim_names
    
    fcst <-
      fcst %>%
      st_set_dimensions("L", values = dates_fcst) %>% 
      st_set_dimensions("M", values = 1:10) %>% 
      adrop() %>% 
      st_set_crs(4326)
    
    
    
    # ***********
    
    
    # download model parameters
    
    ff_nmme <- 
      str_glue("gsutil ls {dir_gs_nmme}/{model}") %>% 
      system(intern = T) %>% 
      str_subset("gamma-params") %>% 
      str_subset("1991-2020") %>% 
      str_subset(str_glue("_{str_pad(month(start), 2, 'left', '0')}"))
    
    
    ff_nmme %>% 
      future_walk(\(f){
        
        str_glue("gsutil cp {f} {dir_data}") %>% 
          system(ignore.stdout = T, ignore.stderr = T)
        
      })
    
    ff_nmme <- 
      ff_nmme %>% 
      fs::path_file() %>% 
      {str_glue("{dir_data}/{.}")}
    
    
    
    # ***********
    
    
    # PROCESS
    
    
    # loop through lead dates
    r <- 
      dates_fcst %>% 
      imap(\(d,d_in){
        
        print(d)
        
        # d_in = 1
        # d = dates_fcst[d_in]
        
        # era params for 1 lead month
        era_params_1mon <- 
          era_params %>%
          slice(L, d_in)
        
        # nmme params for 1 lead month
        nmme_params_1mon <- 
          ff_nmme %>% 
          str_subset(str_glue("lead-{d_in-1}")) %>% 
          read_ncdf() %>%
          suppressMessages() %>% 
          st_warp(era_params_1mon)
        
        # forecast for 1 lead month
        fcst_1mon <- 
          fcst %>%
          units::drop_units() %>% 
          slice(L, d_in) %>% 
          st_warp(era_params_1mon)
        
        
        # bias adjustment process
        fcst_biasadj <- 
          map(1:10, \(mem){
            
            # percentile
            nmme_percentile <- 
              fcst_1mon %>%
              slice(M, mem) %>% 
              c(nmme_params_1mon) %>% 
              mutate(percentile = pgamma(prec, shape = shape, rate = rate),
                     percentile = if_else(percentile > 0.999, 0.999, percentile) # prevents errors in qgamma
                     ) %>% 
              select(percentile) #%>% 
              # starsExtra::focal2(w = matrix(1,3,3),
              #                    fun = "mean",
              #                    na.rm = T,
              #                    mask = T)
            
            # forecasted amount
            fcst_precip <- 
              nmme_percentile %>% 
              c(era_params_1mon) %>% 
              mutate(precip = if_else(shape == -9999, 0, qgamma(percentile, shape = shape, rate = rate) %>% round(1))) %>% 
              suppressWarnings() %>% 
              select(precip) #%>% 
              # starsExtra::focal2(w = matrix(1,3,3),
              #                    fun = "mean",
              #                    na.rm = T,
              #                    mask = T)
            
            
            # forecasted percentile
            fcst_percentile <- 
              fcst_precip %>% 
              c(era_params_1mon) %>% 
              mutate(percentile = if_else(shape == -9999, -9999, pgamma(precip, shape = shape, rate = rate))) %>% 
              select(percentile)
            
            
            c(fcst_precip, fcst_percentile) #nmme_percentile)
            
            
          }) %>% 
          
          {do.call(c, c(., along = "M"))}
        
        
        # test
        {
          
          # real_precip <- 
          #   "/mnt/bucket_mine/era5/monthly_means/total_precipitation/era5_total-precipitation_mon_{d}.nc" %>%
          #   str_glue() %>%
          #   read_ncdf() %>% 
          #   suppressMessages() %>% 
          #   st_warp(era_params_1mon) %>%
          #   mutate(tp = tp %>% units::set_units(mm)) %>%
          #   units::drop_units()
          # 
          # nmme_precip <- 
          #   fcst %>%
          #   units::drop_units() %>%
          #   slice(L, d_in) %>%
          #   st_warp(era_params_1mon) %>%
          #   st_apply(c(1,2), mean)
          # 
          # 
          # fcst_biasadj %>% 
          #   select(precip) %>% 
          #   st_apply(c(1,2), mean, .fname = "precip") %>% 
          #   {. - real_precip} %>%
          #   pull() %>%
          #   abs() %>%
          #   quantile(c(seq(0,0.9,by = 0.1), 0.99)) %>%
          #   round(4)
          # 
          # (nmme_precip - real_precip) %>%
          #   pull() %>%
          #   abs() %>%
          #   quantile(c(seq(0,0.9,by = 0.1), 0.99)) %>%
          #   round(4)
          
          
        }
        
        return(fcst_biasadj)
        
      }) %>% 
      
      {do.call(c, c(., along = "L"))} %>% 
      st_set_dimensions("L", values = dates_fcst)
    
    
    c(f, ff_nmme) %>% 
      fs::file_delete()
    
    return(r)
    
  })
  

# merge everything
rr <- 
  do.call(c, c(rr, along = "M")) %>% 
  st_warp(era_params)


# mean forecasted precipitation across members
rr_precip <- 
  rr %>%
  select(precip) %>% 
  st_apply(c(1,2,4), mean, .fname = "precip") %>% 
  mutate(precip = round(precip, 1))

# rr_precip %>%
#   slice(L, 4) %>%
#   as_tibble() %>% 
#   ggplot(aes(longitude, latitude, fill = precip)) +
#   geom_raster() +
#   coord_equal() +
#   colorspace::scale_fill_continuous_sequential(n.breaks = 11,
#                                            na.value = "transparent",
#                                            rev = T,
#                                            trans = "sqrt",
#                                            limits = c(0,10),
#                                            oob = scales::squish)


# climatological mean
era_mean_pr <-
  era_params %>%
  mutate(precip_mean = if_else(shape == -9999, 0, shape * 1/rate))



# percentile ---> what the map will show
rr_percentile <- 
  map(seq(6), \(lead){
    
    rr_precip %>% 
      c(era_mean_pr) %>% 
      slice(L, lead) %>% 
      mutate(percentile = round(pgamma(precip, shape = shape, rate = rate)*100),
             percentile = if_else(precip_mean < 0.15 & abs(precip - precip_mean) < 0.15, 50, percentile)) %>% 
      
      select(percentile) %>%
      starsExtra::focal2(w = matrix(1,3,3),
                         fun = "mean",
                         na.rm = F,
                         mask = T)
    
  }) %>% 
  {do.call(c, c(., along = "L"))} %>% 
  st_set_dimensions("L", values = dates_fcst) %>% 
  st_warp(era_params)


rr_percentile %>%
  slice(L, lead_) %>%
  as_tibble() %>% #filter(percentile < 0.01)
  ggplot(aes(longitude, latitude, fill = percentile)) +
  geom_raster() +
  coord_equal() +
  colorspace::scale_fill_binned_diverging(mid = 50,
                                          n.breaks = 11,
                                          na.value = "transparent",
                                          rev = T,
                                          limits = c(0,100))



# agreement ---> 2nd layer of map
rr_agreement <- 
  map(seq(6), \(lead){
    
    rr %>% 
      select(percentile) %>% 
      slice(L, lead) %>%
      
      st_apply(c(1,2), \(x){
        
        x_terc <- cut(x, seq(0,1, length.out = 4), labels = F, include.lowest = T)
        mode_terc <- seq(3)[which.max(tabulate(match(x_terc, seq(3))))]
        round(unname(sum(x_terc == mode_terc)/length(x_terc)*100))
        
      },
      .fname = "agreement") %>% 
      starsExtra::focal2(w = matrix(1,3,3),
                         fun = "mean",
                         na.rm = F,
                         mask = T)
    
  }) %>% 
  {do.call(c, c(., along = "L"))} %>% 
  st_set_dimensions("L", values = dates_fcst) %>% 
  st_warp(era_params)

# rr_agreement %>%
#   slice(L, 2) %>%
#   as_tibble() %>%
#   ggplot(aes(longitude, latitude, fill = agreement)) +
#   geom_raster() +
#   coord_equal() +
#   colorspace::scale_fill_binned_sequential(palette = "viridis",
#                                            n.breaks = 9,
#                                            na.value = "transparent",
#                                            rev = T,
#                                            limits = c(33,100)
#   )



# anom ---> how much it normally rains and what percent will fall   

rr_percentnorm <-
  era_mean_pr %>% 
  select(precip_mean) %>% 
  c(rr_precip) %>% 
  mutate(percent = if_else(precip_mean < 0.15 & abs(precip - precip_mean) < 0.15, 100, round(precip/precip_mean*100)),
         precip_mean = precip_mean * 30) %>%
  select(precip_mean, percent)

rr_percentnorm <- 
  map(seq(6), \(lead){
    
    c(
      rr_percentnorm %>% 
        select(precip_mean) %>% 
        slice(L, lead) %>% 
        starsExtra::focal2(w = matrix(1,3,3),
                           fun = "mean",
                           na.rm = F,
                           mask = T),
      
      rr_percentnorm %>% 
        select(percent) %>% 
        slice(L, lead) %>% 
        starsExtra::focal2(w = matrix(1,3,3),
                           fun = "mean",
                           na.rm = F,
                           mask = T)
      
    )
    
  }) %>% 
  {do.call(c, c(., along = "L"))} %>% 
  st_set_dimensions("L", values = dates_fcst) %>% 
  st_warp(era_params)



# rr_percentnorm %>%
#   slice(L, lead_) %>%
#   select(percent) %>%
#   as_tibble() %>% #filter(percent == 0)
#   ggplot(aes(longitude, latitude, fill = percent)) +
#   geom_raster() +
#   coord_equal() +
#   colorspace::scale_fill_binned_diverging(mid = 100, 
#                                           n.breaks = 9,
#                                           na.value = "transparent",
#                                           rev = T,
#                                           limits = c(0,150),
#                                           oob = scales::squish
#   )
# 
# rr_percentnorm %>%
#   slice(L, lead_) %>%
#   select(precip_mean) %>%
#   as_tibble() %>%
#   ggplot(aes(longitude, latitude, fill = precip_mean)) +
#   geom_raster() +
#   coord_equal() +
#   colorspace::scale_fill_continuous_sequential("plasma",
#                                                na.value = "transparent",
#                                                rev = F,
#                                                trans = "sqrt"
#                                                # limits = c(0,250),
#                                                # oob = scales::squish
#                                                )




# EXPORT
c(rr_percentile,
  rr_agreement,
  rr_percentnorm)

# ...

# clean up
dir_data %>%
  fs::dir_ls() %>% 
  fs::file_delete()












