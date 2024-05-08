
library(tidyverse)
library(stars)
library(furrr)

plan(multisession)

dir_data <- "/mnt/pers_disk/data"

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")


tb_extents <- 
  bind_rows(
    tibble(var = "sst", dom = "tp", xmin = 99, ymin = -47, xmax = 287, ymax = 28),
    tibble(var = "sst", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "sst", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    
    tibble(var = "geopotential-500", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "geopotential-500", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    tibble(var = "geopotential-500", dom = "caf", xmin = 9, ymin = -19, xmax = 37, ymax = 12),
    
    tibble(var = "uwind-200", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "uwind-200", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    tibble(var = "uwind-200", dom = "caf", xmin = 9, ymin = -19, xmax = 37, ymax = 12),
    
    tibble(var = "vwind-200", dom = "io", xmin = 44, ymin = -27, xmax = 104, ymax = 24),
    tibble(var = "vwind-200", dom = "ao", xmin = 336, ymin = -22, xmax = 13, ymax = 5),
    tibble(var = "vwind-200", dom = "caf", xmin = 9, ymin = -19, xmax = 37, ymax = 12))




# DOWNLOAD TO LOCAL DISK ----

# fs::dir_create(dir_data)

# tb_extents %>% 
#   pull(var) %>% 
#   unique() %>%
#   walk(function(v) {
#     
#     d <- 
#       "gsutil ls gs://clim_data_reg_useast1/era/monthly" %>%
#       system(intern = T) %>% 
#       str_subset(v)
#       
#     "gsutil -m cp {d}* {dir_data}" %>%
#       str_glue() %>% 
#       system()
# 
#   })






# PROCESS ----
# loop through domains

s_proxy <- 
  dir_data %>% 
  fs::dir_ls() %>% 
  first() %>% 
  read_ncdf(proxy = T)

# loop
l <- vector("list", nrow(tb_extents))

for (i in seq(nrow(tb_extents))) {
  
  print(str_glue("PROCESSING {tb_extents$dom[i]} {tb_extents$var[i]}"))
  
  
  ## LOAD DATA ----
  
  # extent
  xmin <- fn_get_cell_pos(s_proxy, 1, tb_extents$xmin[i])
  ymin <- fn_get_cell_pos(s_proxy, 2, tb_extents$ymax[i])
  xmax <- fn_get_cell_pos(s_proxy, 1, tb_extents$xmax[i])
  ymax <- fn_get_cell_pos(s_proxy, 2, tb_extents$ymin[i])
  
  
  if (xmax-xmin+1 > 0) {
    
    nc <- cbind(start = c(xmin, 
                          ymin, 
                          1),
                count = c(xmax-xmin+1,
                          ymax-ymin+1,
                          NA))
    
    s <- 
      dir_data %>% 
      fs::dir_ls(regexp = tb_extents$var[i]) %>% 
      str_subset(str_flatten(seq(1950, 2022), "|")) %>%
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
      fs::dir_ls(regexp = tb_extents$var[i]) %>% 
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
      fs::dir_ls(regexp = tb_extents$var[i]) %>% 
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
  
  
  
  
  
  ## PRE-PROCESS DATA ----
  
  # 3-month rolling mean
  
  s_3mon <- 
    s %>% 
    st_apply(c(1,2), function(x) {
      
      zoo::rollmean(x, k = 3, align = "right", na.pad = T)
      
    },
    .fname = "time",
    FUTURE = T) %>% 
    aperm(c(2,3,1))
  
  st_dimensions(s_3mon) <- st_dimensions(s)
  
  
  # no de-trending necessary (no significant trend identified)
  # e.g.
  
  # s_3mon %>%
  #   pull() %>%
  #   .[50,35,] %>%
  #   plot(type = "l")
  
  
  # anomalies (by calendar month)
  s_anom <- 
    s_3mon %>% 
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
  
  st_dimensions(s_anom) <- st_dimensions(s)
  
  
  
  ## PREPARE DATA FOR PCA ----
  
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
    
    apply(1, scale, scale = F) %>%
    aperm(c(2,1)) %>%
    
    prcomp(scale. = T, # T
           center = F) # virtually no diff
  
  
  # # % var explained
  # pca$sdev^2/sum(pca$sdev^2)
  # (summary(pca)$importance[2,])[1:10]
  
  
  # # plot (spatial)
  # pca$rotation[,"PC2"] %>%
  #   as_tibble() %>%
  #   bind_cols(
  #     s_anom %>%
  #       slice(time, 3) %>%
  #       as_tibble() %>%
  #       na.omit() %>%
  #       select(1:2)
  #   ) %>%
  #   ggplot(aes(longitude, latitude, fill = value)) +
  #   geom_raster() +
  #   scale_fill_gradient2() +
  #   coord_equal()
  # 
  # 
  # # plot (temporal)
  # pca$x[, "PC1"] %>%
  #   as_tibble() %>%
  #   mutate(time =
  #            tail(st_get_dimension_values(s, 3), -2) %>%
  #            as_date()) %>%
  #   ggplot(aes(time, value)) +
  #   geom_line() +
  #   # scale_y_reverse() +
  #   scale_x_date(date_breaks = "5 years",
  #                date_labels = "%Y",
  #                # limits = c(as_date("1980-01-01"),
  #                #            NA)
  #                ) +
  #   geom_hline(yintercept = 0, linetype = "2222")
  
  
  
  # # test prediction
  # tb_for_pca %>% 
  #   slice(1) %>% 
  #   select(-1) %>% 
  #   apply(1, scale, scale = F) %>%
  #   aperm(c(2,1)) %>% 
  #   {predict(pca, .)} %>% .[,"PC1"]
  # 
  # pca$x[1,"PC1"]
  
  
  
  tb <- 
    tibble(time = 
             st_get_dimension_values(s_anom, 3) %>% 
             as_date() %>% 
             tail(-2)) %>% # remove first 2 months: rolling mean
    bind_cols(pca$x[,1:3] %>% 
                as_tibble() %>% 
                set_names(str_glue("{tb_extents$dom[i]}_{str_remove(tb_extents$var[i], '-')}_{seq(1,3)}")))
  
  
  l[[i]] <- tb
  
}


bind_cols(
  l[[1]],
  l[-1] %>% 
    map(select, -1)
) %>% 
  
  write_csv("tb_eof_ts.csv")



