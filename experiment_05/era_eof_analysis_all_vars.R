
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

dir_data <- "/mnt/pers_disk_/data"

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")



tb_extents <- 
  bind_rows(
    
    tibble(var = "daily-precip", dom = "drc", xmin = 10, ymin = -14, xmax = 32, ymax = 6),
    
    tibble(var = "sst", dom = "tp", xmin = 125, ymin = -47, xmax = 287, ymax = 30),
    tibble(var = "sst", dom = "io", xmin = 40, ymin = -40, xmax = 110, ymax = 20),
    tibble(var = "sst", dom = "sao", xmin = 315, ymin = -39, xmax = 13, ymax = 13),
    tibble(var = "sst", dom = "nao", xmin = 290, ymin = 24, xmax = 350, ymax = 65),
    tibble(var = "sst", dom = "np", xmin = 135, ymin = 20, xmax = 237, ymax = 60),
    
    # tibble(var = "geopotential-200", dom = "tp", xmin = 125, ymin = -47, xmax = 287, ymax = 30),
    # tibble(var = "geopotential-200", dom = "io", xmin = 40, ymin = -40, xmax = 110, ymax = 20),
    # tibble(var = "geopotential-200", dom = "sao", xmin = 315, ymin = -39, xmax = 13, ymax = 13),
    # tibble(var = "geopotential-200", dom = "nao", xmin = 290, ymin = 24, xmax = 350, ymax = 65),
    # tibble(var = "geopotential-200", dom = "np", xmin = 135, ymin = 20, xmax = 237, ymax = 60),
    # 
    # tibble(var = "geopotential-500", dom = "tp", xmin = 125, ymin = -47, xmax = 287, ymax = 30),
    # tibble(var = "geopotential-500", dom = "io", xmin = 40, ymin = -40, xmax = 110, ymax = 20),
    # tibble(var = "geopotential-500", dom = "sao", xmin = 315, ymin = -39, xmax = 13, ymax = 13),
    # tibble(var = "geopotential-500", dom = "nao", xmin = 290, ymin = 24, xmax = 350, ymax = 65),
    # tibble(var = "geopotential-500", dom = "np", xmin = 135, ymin = 20, xmax = 237, ymax = 60),
    # 
    # tibble(var = "geopotential-850", dom = "tp", xmin = 125, ymin = -47, xmax = 287, ymax = 30),
    # tibble(var = "geopotential-850", dom = "io", xmin = 40, ymin = -40, xmax = 110, ymax = 20),
    # tibble(var = "geopotential-850", dom = "sao", xmin = 315, ymin = -39, xmax = 13, ymax = 13),
    # tibble(var = "geopotential-850", dom = "nao", xmin = 290, ymin = 24, xmax = 350, ymax = 65),
    # tibble(var = "geopotential-850", dom = "np", xmin = 135, ymin = 20, xmax = 237, ymax = 60)#,
    
    )


s_proxy <- 
  "/mnt/bucket_mine/era/monthly/mean-daily-precip/era5_mon_mean-daily-precip_1950.nc" %>% 
  read_ncdf(proxy = T)



map(unique(tb_extents$var), function(var){
  
  # DOWNLOAD FILES
  
  d <-
    "gsutil ls gs://clim_data_reg_useast1/era/monthly" %>%
    system(intern = T) %>%
    str_subset(var)

  "gsutil -m cp {d}* {dir_data}" %>%
    str_glue() %>%
    system()
  
  
  # LOOP THROUGH REGIONS
  tb_regs <- 
    tb_extents %>% 
    filter(var == {{var}})
  
  map(seq(nrow(tb_regs)), function(i){
    
    
    # LOAD DATA
    
    # extent
    xmin <- fn_get_cell_pos(s_proxy, 1, tb_regs$xmin[i])
    ymin <- fn_get_cell_pos(s_proxy, 2, tb_regs$ymax[i])
    xmax <- fn_get_cell_pos(s_proxy, 1, tb_regs$xmax[i])
    ymax <- fn_get_cell_pos(s_proxy, 2, tb_regs$ymin[i])
    
    
    # load
    # if reg does not crosses 360 line
    if (xmax-xmin+1 > 0) {
      
      nc <- cbind(start = c(xmin, 
                            ymin, 
                            1),
                  count = c(xmax-xmin+1,
                            ymax-ymin+1,
                            NA))
      
      s <- 
        dir_data %>% 
        fs::dir_ls() %>% 
        str_subset(str_flatten(seq(1970, 2022), "|")) %>%
        future_map(read_ncdf,
                   ncsub = nc,
                   downsample = c(3,3,0)) %>% 
        suppressMessages() %>% 
        do.call(c, .)
      
    } else {
      
      nc <- cbind(start = c(xmin, 
                            ymin, 
                            1),
                  count = c(NA,
                            ymax-ymin+1,
                            NA))
      
      s1 <- 
        dir_data %>% 
        str_subset(str_flatten(seq(1950, 2022), "|")) %>%
        future_map(read_ncdf,
                   ncsub = nc,
                   downsample = c(3,3,0)) %>% 
        suppressMessages() %>% 
        do.call(c, .)
      
      nc <- cbind(start = c(1, 
                            ymin, 
                            1),
                  count = c(xmax,
                            ymax-ymin+1,
                            NA))
      
      s2 <- 
        dir_data %>% 
        str_subset(str_flatten(seq(1950, 2022), "|")) %>%
        future_map(read_ncdf,
                   ncsub = nc,
                   downsample = c(3,3,0)) %>% 
        suppressMessages() %>% 
        do.call(c, .)
      
      s <- 
        c(s1,s2, along = 1) %>% 
        st_set_crs(4326)
      
    }
    
    
    
    # anomalies
    s_anom <- 
      s %>% 
      st_apply(c(1,2), function(x) {
        
        if(all(is.na(x))) {
          rep(NA, length(x))
        } else {
          
          xx <- 
            x %>% 
            zoo::rollmean(3, fill = NA, align = "right") 
          
          mm <- 
            xx %>% 
            matrix(ncol = 12, byrow = T) %>% 
            apply(2, mean, na.rm = T)
          
          xx-mm
          
          # m <- matrix(x, ncol = 12, byrow = T)
          # mm <- apply(m, 2, mean, na.rm = T)
          # as.vector(t(x-mm))
          
        }
        
      },
      .fname = "time",
      FUTURE = F) %>% 
      aperm(c(2,3,1))
    
    st_dimensions(s_anom) <- st_dimensions(s)
  
    
    # pca
    tb_anom <- 
      s_anom %>% 
      st_set_dimensions("time", 
                        values = 
                          st_get_dimension_values(s_anom, "time") %>% 
                          str_sub(end = 7) %>% 
                          str_replace("-", "_")) %>% 
      as_tibble() %>% 
      na.omit()
    
    tb_anom <- 
      tb_anom %>% 
      mutate(time = as.character(time)) %>% 
      rename("v" = 4)
    
    tb_for_pca <- 
      tb_anom %>%
      group_by(time) %>% 
      mutate(gcell = row_number()) %>% 
      ungroup() %>%
      select(time, v, gcell) %>% 
      pivot_wider(names_from = "gcell", names_prefix = "g_", values_from = v)
    
    
    pca <- 
      tb_for_pca %>% 
      select(-1) %>% 
      
      apply(1, scale, scale = F) %>% # scaled divides by sd
      aperm(c(2,1)) %>%
      
      prcomp(scale. = F)
    
    
    # save res
    
    
    
  })
  
  
  dir_data %>% 
    fs::dir_ls() %>% 
    fs::file_delete()
  
})


