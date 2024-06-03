
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)

# source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")




cmip_var_specs <- 
  bind_rows(
    tibble(var = "tos", prov = "Omon", grid = "gr1"),
    tibble(var = "zg", lev = list(c(200, 500, 850)), prov = "Amon", grid = "gr"),
    tibble(var = "hus", lev = list(c(200, b = 850)), prov = "Amon", grid = "gr"),
    tibble(var = "ua", lev = list(c(200, 850)), prov = "Amon", grid = "gr"),
    tibble(var = "va", lev = list(c(200, 850)), prov = "Amon", grid = "gr")
  )



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


gcm <- "CNRM-CERFACS/CNRM-CM6-1"
mems <- str_glue("r{seq(1,30,5)}i1p1f2")



# l <- vector("list", nrow(cmip_var_specs))
# names(l) <- 
#   cmip_var_specs %>% 
#   mutate(lev = if_else(is.na(lev), "", as.character(lev)),
#          v = str_glue("{var}{lev}")) %>% 
#   pull(v)



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




tbs <- 
  map(mems, function(mem){
    
    print(str_glue("Processing member {mem}"))
    
    pmap_dfc(cmip_var_specs, function(var, prov, grid, lev) {
      
      # var = cmip_var_specs$var[2]
      # prov = cmip_var_specs$prov[2]
      # grid = cmip_var_specs$grid[2]
      # lev = cmip_var_specs$lev[2]
      
      time_dim <- time_dims %>% pluck(grid)
      
      
      print(str_glue("   var {var}"))
      
      r <- 
        str_glue("gcloud storage ls gs://cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/{prov}/{var}/{grid}/") %>% 
        system(intern = T) %>% 
        str_sub(start = 6, end = -2)
      
      dsn <- 
        str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:{var}.zarr/')
      
      
      
      if (is.null(unlist(lev))) {
        
        s <- 
          read_mdim(dsn, 
                    offset = c(0, 0, time_dim$time_i-1))
        
        tb_var <-
          s %>% 
          st_extract(pts) %>% 
          as_tibble() %>% 
          mutate(id = rep(pts$id, length(time_dim$time_vector_sub)),
                 id = str_pad(id, 2, "left", "0")) %>% 
          select(-x) %>% 
          pivot_wider(names_from = id, 
                      values_from = var,
                      names_prefix = "loc") %>% 
          select(-time) %>% 
          rename_with(~str_glue("{.x}_{var}")) %>% 
          select(where(~sum(!is.na(.x)) > 0)) %>% 
          units::drop_units()
        
        
      } else {
        
        s <-
          read_mdim(dsn, 
                    offset = c(0,0,0, time_dim$time_i-1)) %>% 
          st_set_dimensions(2, 
                            values = seq(st_get_dimension_values(., 2, where = "start") %>% first(),
                                         st_get_dimension_values(., 2, where = "start") %>% last(),
                                         length.out = dim(.)[2])) %>% 
          st_set_crs(4326)
        
        plev_val <- st_get_dimension_values(s, "plev")/100
        
        plev_pos <- which(plev_val %in% unlist(lev))
        
        tb_var <- 
          map_dfc(plev_pos, function(pp) {
            
            v <- str_glue("{var}{plev_val[pp]}")
            
            s %>% 
              slice(plev, pp) %>% 
              st_extract(pts) %>% 
              as_tibble() %>% 
              mutate(id = rep(pts$id, length(time_dim$time_vector_sub)),
                     id = str_pad(id, 2, "left", "0")) %>% 
              select(-x) %>% 
              pivot_wider(names_from = id, 
                          values_from = var,
                          names_prefix = "loc") %>% 
              select(-time) %>% 
              rename_with(~str_glue("{.x}_{v}")) %>% 
              select(where(~sum(!is.na(.x)) > 0)) %>% 
              units::drop_units()
            
          }) 
        
      }
      
      return(tb_var)
      
    }) %>% 
      
      {tibble(member = mem,
              time = time_dims$gr$time_vector_sub,
              .)}
    
  })



write_rds(tbs, "experiment_02/predictors_data_gcm.rds")





