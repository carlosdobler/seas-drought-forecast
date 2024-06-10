
# EOF analysis of several regions with several 
# variables to obtain the components' time series,
# which will represent the model's predictors.

# This script is (almost) identical to emperiment_01/gcm_2_pca.R
# but includes the EOF analysis of precipitation over the DRC



library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

sf_use_s2(F)



cmip_var_specs <- 
  bind_rows(
    
    tibble(var = "pr",
           prov = "Amon",
           grid = "gr",
           regs = list(list(
             drc = c(xmin = 10, ymin = -14, xmax = 32, ymax = 6)))),
    
    tibble(var = "tos", 
           prov = "Omon", 
           grid = "gr1", 
           regs = list(list(
             tp = c(xmin = 125, ymin = -47, xmax = 287, ymax = 30),
             io = c(xmin = 40, ymin = -40, xmax = 110, ymax = 20),
             sao = c(xmin = 315, ymin = -39, xmax = 13, ymax = 13),
             nao = c(xmin = 290, ymin = 24, xmax = 350, ymax = 65),
             np = c(xmin = 135, ymin = 20, xmax = 237, ymax = 60)))),
    
    tibble(var = "zg",
           lev = list(c(500)),#list(c(200, 500, 850)),
           prov = "Amon",
           grid = "gr",
           regs = list(list(
             tp = c(xmin = 125, ymin = -47, xmax = 287, ymax = 30),
             io = c(xmin = 40, ymin = -40, xmax = 110, ymax = 20),
             sao = c(xmin = 315, ymin = -39, xmax = 13, ymax = 13),
             nao = c(xmin = 290, ymin = 24, xmax = 350, ymax = 65),
             np = c(xmin = 135, ymin = 20, xmax = 237, ymax = 60),
             caf = c(xmin = 4, ymin = -21, xmax = 42, ymax = 15)))),
    # 
    # tibble(var = "hus", 
    #        lev = list(c(200, 850)), 
    #        prov = "Amon", 
    #        grid = "gr",
    #        regs = list(list(
    #          #tp = c(xmin = 125, ymin = -47, xmax = 287, ymax = 30),
    #          io = c(xmin = 40, ymin = -40, xmax = 110, ymax = 20),
    #          sao = c(xmin = 315, ymin = -39, xmax = 13, ymax = 13),
    #          #nao = c(xmin = 290, ymin = 24, xmax = 350, ymax = 65),
    #          #np = c(xmin = 135, ymin = 20, xmax = 237, ymax = 60),
    #          caf = c(xmin = 4, ymin = -21, xmax = 42, ymax = 15)))),
    # 
    # tibble(var = "ua",
    #        lev = list(c(200, 850)),
    #        prov = "Amon",
    #        grid = "gr",
    #        regs = list(list(
    #          #tp = c(xmin = 125, ymin = -47, xmax = 287, ymax = 30),
    #          io = c(xmin = 40, ymin = -40, xmax = 110, ymax = 20),
    #          sao = c(xmin = 315, ymin = -39, xmax = 13, ymax = 13),
    #          #nao = c(xmin = 290, ymin = 24, xmax = 350, ymax = 65),
    #          #np = c(xmin = 135, ymin = 20, xmax = 237, ymax = 60),
    #          caf = c(xmin = 4, ymin = -21, xmax = 42, ymax = 15)))),
    # 
    # tibble(var = "va", 
    #        lev = list(c(200, 850)), 
    #        prov = "Amon", 
    #        grid = "gr",
    #        regs = list(list(
    #          #tp = c(xmin = 125, ymin = -47, xmax = 287, ymax = 30),
    #          io = c(xmin = 40, ymin = -40, xmax = 110, ymax = 20),
    #          sao = c(xmin = 315, ymin = -39, xmax = 13, ymax = 13),
    #          #nao = c(xmin = 290, ymin = 24, xmax = 350, ymax = 65),
    #          #np = c(xmin = 135, ymin = 20, xmax = 237, ymax = 60),
    #          caf = c(xmin = 4, ymin = -21, xmax = 42, ymax = 15))))
  )


gcm <- "CNRM-CERFACS/CNRM-CM6-1"
mems <- str_glue("r{seq(1,30,5)}i1p1f2")



# Get time

tb_for_time_dims <- 
  cmip_var_specs %>%
  group_by(grid) %>% 
  slice_head(n = 1)

time_dims <- 
  tb_for_time_dims %>% 
  pmap(function(var, prov, grid, ...) {
    
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
    
    if (grid == "gr1") {
      
      time_vector <- 
        str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:/time') %>% 
        read_stars() %>% 
        pull() %>% 
        as.vector() %>% 
        {(. + 360) * 60 * 60} %>% # pcict works in seconds; 360 to ensure it to be mid month
        PCICt::as.PCICt(cal = time_array$attributes$calendar,
                        origin = time_array$attributes$time_origin)
      
    } else if (grid == "gr") {
      
      time_vector <- 
        str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:/time') %>% 
        read_stars() %>% 
        pull() %>% 
        as.vector() %>% 
        {. * 24 * 60 * 60} %>% # pcict works in seconds
        PCICt::as.PCICt(cal = time_array$attributes$calendar,
                        origin = time_array$attributes$time_origin)
      
    }
    
    time_i <- 
      first(which(year(time_vector) == 1940))
    
    time_vector_sub <- 
      time_vector %>% 
      tail(-time_i+1) %>% 
      str_sub(end = 8) %>% 
      paste0("01") %>% 
      as_date()
    
    return(list(time_i = time_i,
                time_vector_sub = time_vector_sub))
    
  }) %>% 
  set_names(tb_for_time_dims$grid)


tm_vector <- 
  tidyr::expand_grid(
    mem = 
      mems,
    time = time_dims$gr$time_vector_sub
  ) %>%
  mutate(tm = str_glue("{mem}_{time}")) %>% 
  pull(tm)




# Prepare data

reference_grid <-
  st_bbox(c(xmin = 0, ymin = -90, xmax = 360, ymax = 90), crs = 4326) %>% 
  st_as_stars(dx = 1.5, values = NA_real_)


tb_pca <- 
  
  pmap_dfc(cmip_var_specs, function(var, prov, grid, lev, regs) {
    
    # var = cmip_var_specs$var[3]
    # prov = cmip_var_specs$prov[3]
    # grid = cmip_var_specs$grid[3]
    # lev = cmip_var_specs$lev[3]
    # regs = cmip_var_specs$regs[3]
    
    time_dim <- time_dims %>% pluck(grid)
    
    print(str_glue("Processing var {var}"))
    
    ss <-
      map(mems %>% set_names(), function(mem){
        
        print(str_glue("   member {mem}"))
        
        r <- 
          str_glue("gcloud storage ls gs://cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/{prov}/{var}/{grid}/") %>% 
          system(intern = T) %>% 
          str_sub(start = 6, end = -2)
        
        dsn <- 
          str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:{var}.zarr/')
        
        
        
        if (is.null(unlist(lev))) {
          
          # load zarr file
          s <- 
            read_mdim(dsn, 
                      offset = c(0, 0, time_dim$time_i-1))
          
          
          # irregular lat?
          if(st_get_dimension_values(s, 2) %>% diff() %>% unique() %>% length() %>% {. > 1}) {
            
            s <- 
              s %>% 
              st_set_dimensions(2,
                                values = seq(st_get_dimension_values(., 2, where = "start") %>% first(),
                                             st_get_dimension_values(., 2, where = "start") %>% last(),
                                             length.out = dim(.)[2])) %>%
              st_set_crs(4326)
            
          }
          
          s <- 
            s %>%
            st_warp(reference_grid)
          
          
        } else { # variable has plevels
          
          # load zarr file
          s <- 
            read_mdim(dsn, 
                      offset = c(0, 0, 0, time_dim$time_i-1))
          
          
          # irregular lat?
          if(st_get_dimension_values(s, 2) %>% diff() %>% unique() %>% length() %>% {. > 1}) {
            
            s <- 
              s %>% 
              st_set_dimensions(2,
                                values = seq(st_get_dimension_values(., 2, where = "start") %>% first(),
                                             st_get_dimension_values(., 2, where = "start") %>% last(),
                                             length.out = dim(.)[2])) %>%
              st_set_crs(4326)
            
          }
          
          
          plev_val <- st_get_dimension_values(s, "plev")/100
          
          plev_pos <- which(plev_val %in% unlist(lev))
          
          s <- 
            s %>% 
            slice(plev, plev_pos)
          
          s <- 
            s %>% 
            st_warp(reference_grid)
          
        }
        
        
        
        # crop regions
        s_regs <- 
          regs %>% #unlist(recursive = F) %>% 
          map(function(reg) {
            
            if (reg["xmin"] > 180 & reg["xmax"] < 180){
              
              reg1 <- reg
              reg1["xmax"] <- 359.5
              
              s1 <- 
                st_crop(s, st_bbox(reg1, crs = 4326)) %>% 
                suppressMessages() %>% 
                st_set_dimensions(1, 
                                  values = st_get_dimension_values(., 1, where = "start") - 360)
              
              reg2 <- reg
              reg2["xmin"] <- 0.5
              
              s2 <- 
                st_crop(s, st_bbox(reg2, crs = 4326)) %>% 
                suppressMessages()
              
              c(s1,s2, along = 1) %>%
                st_set_crs(4326)
              
              
            } else {
              
              s %>% 
                st_crop(st_bbox(reg, crs = 4326)) %>% 
                suppressMessages()
              
            }
            
          })
        
        
        # split plev
        if (s_regs[[1]] %>% dim() %>% length() > 3) {
          
          s_regs <- 
            s_regs %>% 
            map(function(s_r) {
              
              plev_val <- st_get_dimension_values(s_r, "plev")/100
              
              seq_len(dim(s_r)[3]) %>% 
                map(function(pp) {
                  
                  s_r %>% 
                    slice(plev, pp) %>% 
                    setNames(str_glue("{var}{plev_val[pp]}"))
                  
                })
              
            }) %>% 
            unlist(recursive = F)
          
        }
        
        
        # rolling means
        # s_regs_roll <-
        #   s_regs %>%
        #   map(function(s_reg){
        # 
        #     s_reg %>%
        #       st_apply(c(1,2), function(x) {
        # 
        #         if(all(is.na(x))) {
        #           rep(NA, length(x))
        #         } else {
        # 
        #           zoo::rollmean(x, 3, fill = NA, align = "right")
        # 
        #         }
        # 
        #       },
        #       .fname = "time") %>%
        #       aperm(c(2,3,1))
        # 
        #   })
        # 
        # return(s_regs_roll)
        
        return(s_regs)
        
        
      })
    
    
    print(str_glue(" Calculating PCAs"))
    
    tb_var <- 
      ss %>% 
      transpose() %>% # members 2nd level; regions 1st level
      imap_dfc(function(ssr, i){
        
        print(str_glue("   {names(ssr[[1]])} {str_extract(i, '[a-z]+')}"))
        
        # anomalies
        s_anom <- 
          ssr %>% 
          {do.call(c, c(., along = 3))} %>% 
          
          st_apply(c(1,2), function(x){
            
            if(all(is.na(x))) {
              rep(NA, length(x))
            } else {
              
              mm <- 
                x %>% 
                matrix(ncol = 12, byrow = T) %>% 
                apply(2, mean, na.rm = T)
              
              x-mm
              
            }
            
          },
          .fname = "time",
          FUTURE = F) %>%
          aperm(c(2,3,1))
        
        
        # pca prep
        
        tb <- 
          s_anom %>%
          st_set_dimensions(3, values = tm_vector) %>% 
          as_tibble() %>% 
          mutate(time = as.character(time)) %>% 
          rename("value" = 4)
        
        tb_for_pca <- 
          tb %>% 
          drop_na() %>%
          group_by(time) %>% 
          mutate(gcell = row_number()) %>%
          ungroup() %>%
          select(time, value, gcell) %>% 
          pivot_wider(names_from = "gcell", names_prefix = "g_", values_from = value) %>%
          select(where(~!any(is.na(.))))
        
        pca <- 
          tb_for_pca %>% 
          select(-1) %>% 
          apply(1, scale, scale = F) %>%
          aperm(c(2,1)) %>%
          prcomp(scale. = F)
        
        {
          # # plot (spatial)
          # pca$rotation[,"PC2"] %>%
          #   as_tibble() %>%
          #   bind_cols(
          #     tb %>%
          #       filter(time == "r1i1p1f2_1941-01-01") %>%
          #       drop_na() %>%
          #       select(1:2)
          #   ) %>%
          #   ggplot(aes(x, y, fill = value)) +
          #   geom_raster() +
          #   scale_fill_gradient2() +
          #   coord_equal()
          # 
          # 
          # # plot (temporal)
          # pca$x[, "PC1"] %>%
          #   as_tibble() %>%
          #   mutate(time = tb_for_pca$time) %>%
          #   separate_wider_delim(time, "_", names = c("mem", "time")) %>%
          #   mutate(time = as_date(time)) %>%
          # 
          #   filter(mem == "r1i1p1f2") %>%
          #   # filter(month(time) == 11) %>% 
          # 
          #   ggplot(aes(time, value)) +
          #   geom_line() +
          #   scale_x_date(date_breaks = "5 years",
          #                date_labels = "%Y",
          #                ) +
          #   scale_y_reverse() +
          #   geom_hline(yintercept = 0, linetype = "2222")
          # 
          # 
          # # # composite (high)
          # pca$x[, "PC1"] %>%
          #   as_tibble() %>%
          #   mutate(time = tm_vector) %>%
          #   slice_max(value, n = 20) %>%
          #   select(-value) %>%
          #   left_join(tb, by = "time") %>%
          #   separate_wider_delim(time, "_", names = c("member", "time")) %>%
          #   mutate(time = as_date(time)) %>% 
          #   ggplot(aes(x, y, fill = value)) +
          #   geom_raster() +
          #   colorspace::scale_fill_continuous_diverging(rev = T, na.value = "transparent") +
          #   facet_wrap(~time, ncol = 5)
          # 
          # # composite (low)
          # pca$x[, "PC1"] %>%
          #   as_tibble() %>%
          #   mutate(time = tm_vector) %>%
          #   slice_min(value, n = 20) %>%
          #   select(-value) %>%
          #   left_join(tb, by = "time") %>%
          #   separate_wider_delim(time, "_", names = c("member", "time")) %>%
          #   mutate(time = as_date(time)) %>% 
          #   ggplot(aes(x, y, fill = value)) +
          #   geom_raster() +
          #   colorspace::scale_fill_continuous_diverging(rev = T, na.value = "transparent") +
          #   facet_wrap(~time, ncol = 5)
          
          }
        
        n_comps <- 3
        
        bind_cols(pca$x[,1:n_comps] %>% 
                    as_tibble() %>% 
                    mutate(across(everything(), \(x) as.vector(scale(x, center = F)))) %>% 
                    set_names(str_glue("{names(ssr[[1]])}_{str_extract(i, '[a-z]+')}_eof{seq(1,n_comps)}")))
        
      })
    
    return(tb_var)
    
  })



# Assemble table

tb_f <- 
  tibble(tm = tm_vector) %>% 
  separate_wider_delim(tm, "_", names = c("mem", "time")) %>%
  mutate(time = as_date(time)) %>% 
  bind_cols(tb_pca)

write_csv(tb_f, "experiment_05/data/gcm_tb_eof.csv")



