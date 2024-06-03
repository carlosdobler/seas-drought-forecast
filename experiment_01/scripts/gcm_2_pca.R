

library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)


source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")


tb_extents <- 
  bind_rows(
    tibble(var = "tos", prov = "Omon", grid = "gr1", dom = "tp", xmin = 99, ymin = -47, xmax = 287, ymax = 28),
    tibble(var = "tos", prov = "Omon", grid = "gr1", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "tos", prov = "Omon", grid = "gr1", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    
    tibble(var = "zg", lev = 500, prov = "Amon", grid = "gr", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "zg", lev = 500, prov = "Amon", grid = "gr", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    tibble(var = "zg", lev = 500, prov = "Amon", grid = "gr", dom = "caf", xmin = 9, ymin = -19, xmax = 37, ymax = 12),
    
    tibble(var = "ua", lev = 200, prov = "Amon", grid = "gr", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "ua", lev = 200, prov = "Amon", grid = "gr", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    tibble(var = "ua", lev = 200, prov = "Amon", grid = "gr", dom = "caf", xmin = 9, ymin = -19, xmax = 37, ymax = 12),
    
    tibble(var = "va", lev = 200, prov = "Amon", grid = "gr", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "va", lev = 200, prov = "Amon", grid = "gr", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    tibble(var = "va", lev = 200, prov = "Amon", grid = "gr", dom = "caf", xmin = 9, ymin = -19, xmax = 37, ymax = 12)
    
    )


gcm <- "CNRM-CERFACS/CNRM-CM6-1"
mems <- str_glue("r{seq(1,30,5)}i1p1f2")




l <- vector("list", nrow(tb_extents))
names(l) <- 
  tb_extents %>% 
  mutate(lev = if_else(is.na(lev), "", as.character(lev)),
         v = str_glue("{var}{lev}")) %>% 
  pull(v)


for (i in seq(nrow(tb_extents))) {
  
  print(str_glue("Loading var-domain {i} / {nrow(tb_extents)}"))
  
  # root path
  r <- 
    str_glue("gcloud storage ls gs://cmip6/CMIP6/CMIP/{gcm}/historical/{mems[1]}/{tb_extents$prov[i]}/{tb_extents$var[i]}/{tb_extents$grid[i]}/") %>% 
    system(intern = T) %>% 
    str_sub(start = 6, end = -2)
  
  dsn <- 
    str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:{tb_extents$var[i]}.zarr/')
  
  
  # get time
  ret = gdal_utils("mdiminfo", dsn, quiet = TRUE)
  time_array <- jsonlite::fromJSON(ret)$arrays$time
  
  if (tb_extents$grid[i] == "gr1") {
    
    time_vector <- 
      str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:/time') %>% 
      read_stars() %>% 
      pull() %>% 
      as.vector() %>% 
      {(. + 360) * 60 * 60} %>% # pcict works in seconds; 360 for it to be mid month
      PCICt::as.PCICt(cal = time_array$attributes$calendar,
                      origin = time_array$attributes$time_origin)
    
  } else if (tb_extents$grid[i] == "gr") {
    
    time_vector <- 
      str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:/time') %>% 
      read_stars() %>% 
      pull() %>% 
      as.vector() %>% 
      {. * 24 * 60 * 60} %>% # pcict works in seconds
      PCICt::as.PCICt(cal = time_array$attributes$calendar,
                      origin = time_array$attributes$time_origin)
    
  }
  
  time_i <- first(which(year(time_vector) == 1940))
  
  
  
  # get spatial extents
  if (is.na(tb_extents$lev[i])) {
    
    s_proxy <- 
      read_mdim(dsn,
                count = c(NA,NA,1)) %>% 
      adrop()
    
  } else {
    
    s_proxy <- 
      read_mdim(dsn,
                count = c(NA,NA,1,1)) %>% 
      adrop()
    
    # get levels
    plev_vector <- 
      str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:/plev') %>% 
      read_stars() %>% 
      pull() %>% 
      as.vector() %>% 
      {./100}
    
  }
  
  
  off_x <- fn_get_cell_pos(s_proxy, 1, tb_extents$xmin[i])
  count_x <- fn_get_cell_pos(s_proxy, 1, tb_extents$xmax[i]) - off_x
  off_y <- fn_get_cell_pos(s_proxy, 2, tb_extents$ymin[i])
  count_y <- fn_get_cell_pos(s_proxy, 2, tb_extents$ymax[i]) - off_y
  
  
  # load data
  
  ll <- 
    future_map(set_names(mems), function(mem) {
      
      # print(str_glue("   member: {mem}"))
      
      r <- 
        str_glue("gcloud storage ls gs://cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/{tb_extents$prov[i]}/{tb_extents$var[i]}/{tb_extents$grid[i]}/") %>% 
        system(intern = T) %>% 
        str_sub(start = 6, end = -2)
      
      dsn <- 
        str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:{tb_extents$var[i]}.zarr/')
      
      
      if (count_x > 0) {
        
        if (is.na(tb_extents$lev[i])) {
          
          s <- 
            read_mdim(dsn, 
                      offset = c(off_x, off_y, time_i-1), 
                      count = c(count_x, count_y, NA))
          
        } else {
          
          s <-
            read_mdim(dsn, 
                      offset = c(off_x, 
                                 off_y,
                                 which(plev_vector == tb_extents$lev[i])-1,
                                 time_i-1),
                      
                      count = c(count_x, 
                                count_y,
                                1,
                                NA)) %>% 
            adrop()
          
        }
        
      } else {
        
        if (is.na(tb_extents$lev[i])) {
          
          s1 <- 
            read_mdim(dsn, 
                      offset = c(off_x, off_y, time_i-1), 
                      count = c(NA, count_y, NA)) %>% 
            st_set_dimensions(1, 
                              values = st_get_dimension_values(., 1, where = "start") - 360)
          
          s2 <- 
            read_mdim(dsn,
                      offset = c(0, off_y, time_i-1),
                      count = c(fn_get_cell_pos(s_proxy, 1, tb_extents$xmax[i]),
                                count_y,
                                NA))
        
          s <- 
            c(s1,s2, along = 1) %>% 
            st_set_crs(4326)
          
            
        } else {
          
          s1 <- 
            read_mdim(dsn, 
                      offset = c(off_x, 
                                 off_y, 
                                 which(plev_vector == tb_extents$lev[i])-1,
                                 time_i-1),
                      
                      count = c(NA, 
                                count_y, 
                                1,
                                NA)) %>% 
            
            st_set_dimensions(1, 
                              values = st_get_dimension_values(., 1, where = "start") - 360)
          
          s2 <- 
            read_mdim(dsn,
                      offset = c(0, 
                                 off_y, 
                                 which(plev_vector == tb_extents$lev[i])-1,
                                 time_i-1),
                      
                      count = c(fn_get_cell_pos(s_proxy, 1, tb_extents$xmax[i]),
                                count_y,
                                1,
                                NA))
          
          s <- 
            c(s1,s2, along = 1) %>% 
            st_set_crs(4326) %>% 
            adrop()
          
        }
        
      }
      
      
      if (tb_extents$grid[i] == "gr") {
        
        # reg dimension
        s <- 
          s %>%
          st_set_dimensions(2, 
                            values = seq(st_get_dimension_values(s, 2, where = "start") %>% first(),
                                         st_get_dimension_values(s, 2, where = "start") %>% last(),
                                         length.out = dim(s)[2])) %>%
          st_set_crs(4326)
        
      }
      
      return(s)
      
    })
  
  l[[i]] <- ll
  
}





# anomalies

ss_anom <- 
  
  l %>% 
  map(function(ss) {
    
    s <-
      do.call(c, c(ss, along = "time"))

    s_anom <-
      s %>%
      st_apply(c(1,2), function(x) {

        if(all(is.na(x))) {
          rep(NA, length(x))
        } else {

          m <- matrix(x, ncol = 12, byrow = T)
          mm <- apply(m, 2, mean, na.rm = T)
          as.vector(t(x-mm))

        }

      },
      .fname = "time",
      FUTURE = F) %>%
      aperm(c(2,3,1))

    return(s_anom)
    
    # ss %>%
    #   map(function(s) {
    # 
    #     s %>%
    #       st_apply(c(1,2), function(x) {
    # 
    #         if(all(is.na(x))) {
    #           rep(NA, length(x))
    #         } else {
    # 
    #           m <- matrix(x, ncol = 12, byrow = T)
    #           mm <- apply(m, 2, mean, na.rm = T)
    #           as.vector(t(x-mm))
    # 
    #         }
    # 
    #       },
    #       .fname = "time",
    #       FUTURE = F) %>%
    #       aperm(c(2,3,1))
    # 
    # 
    #   }) %>%
    # 
    #   {do.call(c, c(., along = "time"))}
    
  })




# pca

tb_time_mem <- 
  
  tidyr::expand_grid(
    
    mem = 
      mems,
    
    time = 
      l[[1]][[1]] %>% 
      st_get_dimension_values("time") %>% 
      str_sub(end = 7)
    
  )

time_mem_vector <- 
  tb_time_mem %>% 
  mutate(tm = str_glue("{mem}_{time}")) %>% 
  pull(tm)



tbs <- 
  pmap(list(ss_anom, names(ss_anom), tb_extents$dom), function(s, var, dom){
    
    # s = ss_anom[[3]]
    
    tb_anom <- 
      s %>% 
      st_set_dimensions("time", 
                        values = time_mem_vector) %>% 
      as_tibble() %>% 
      na.omit()
    
    tb_anom <- 
      tb_anom %>% 
      mutate(time = as.character(time)) %>% # from factor
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
      
      apply(1, scale, scale = F) %>%
      aperm(c(2,1)) %>%
      
      prcomp(scale. = T, # T
             center = F) # virtually no diff
    
    
    # # plot (spatial)
    # pca$rotation[,"PC3"] %>%
    #   as_tibble() %>%
    #   bind_cols(
    #     s %>%
    #       slice(time, 3) %>%
    #       as_tibble() %>%
    #       na.omit() %>%
    #       select(1:2)
    #   ) %>%
    #   ggplot(aes(lon, lat, fill = value)) +
    #   geom_raster() +
    #   scale_fill_gradient2() +
    #   coord_equal()
    # 
    # 
    # # plot (temporal)
    # pca$x[, "PC1"] %>%
    #   as_tibble() %>%
    #   mutate(time = time_mem_vector) %>%
    #   separate_wider_delim(time, "_", names = c("mem", "time")) %>%
    #   mutate(time = as_date(str_glue("{time}-01"))) %>%
    # 
    #   filter(mem == "r1i1p1f2") %>%
    # 
    #   ggplot(aes(time, value)) +
    #   geom_line() +
    #   # scale_y_reverse() +
    #   scale_x_date(date_breaks = "5 years",
    #                date_labels = "%Y",
    #                ) +
    #   geom_hline(yintercept = 0, linetype = "2222")
    
    
    
    
    # final table
    
    tb <- 
      tibble(time = time_mem_vector) %>% 
      separate_wider_delim(time, "_", names = c("mem", "time")) %>% 
      mutate(time = as_date(str_glue("{time}-01"))) %>% 
      bind_cols(pca$x[,1:5] %>% 
                  as_tibble() %>% 
                  set_names(str_glue("{dom}_{var}_{seq(1,5)}")))
    
    
    return(tb)
    
  })


bind_cols(
  tbs[[1]],
  tbs[-1] %>% 
    map(select, -c(1:2))
) %>% 
  
  write_csv("gcm_tb_eof.csv")


