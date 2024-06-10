
# Cluster analysis of rainfall anomalies in the DRC.
# Resulting clusters will represent the response variable
# in a classification task model.


library(tidyverse)
library(stars)
library(furrr)

plan(multisession)

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")


land <- 
  "/mnt/bucket_mine/misc_data/ne_110m_land/" %>% 
  read_sf(quiet = T) %>% 
  mutate(a = 1) %>% 
  select(a)




# prepare vars to load data

gcm <- "CNRM-CERFACS/CNRM-CM6-1"
prov <- "Amon"
var <- "pr"
grid <- "gr"

mems <- str_glue("r{seq(1,30,5)}i1p1f2")

root_dsn <- 
  str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/cmip6/CMIP6/CMIP/{gcm}/historical/{mems[1]}/{prov}/{var}/{grid}/v20180917"/')

dsn <- 
  str_glue('{root_dsn}:{var}.zarr/')

ret = gdal_utils("mdiminfo", dsn, quiet = TRUE)
time_array <- jsonlite::fromJSON(ret)$arrays$time

time_vector <- 
  str_glue('{root_dsn}:/time') %>% 
  read_stars() %>% 
  pull() %>% 
  as.vector() %>% 
  {. * 24 * 60 * 60} %>% # pcict works in seconds
  PCICt::as.PCICt(cal = time_array$attributes$calendar,
                  origin = time_array$attributes$time_origin)


time_i <- first(which(year(time_vector) == 1940))

gcm_proxy <- 
  read_mdim(dsn,
            count = c(NA,NA,1)) %>% 
  adrop()


off_x <- fn_get_cell_pos(gcm_proxy, 1, 5)
count_x <- fn_get_cell_pos(gcm_proxy, 1, 37) - off_x
off_y <- fn_get_cell_pos(gcm_proxy, 2, -19)
count_y <- fn_get_cell_pos(gcm_proxy, 2, 8) - off_y


gcm_proxy_ext <- 
  read_mdim(dsn, 
            offset = c(off_x, off_y, 0), 
            count = c(count_x, count_y, 1)) %>% 
  adrop()

new_grid <- 
  st_as_stars(st_bbox(gcm_proxy_ext), dx = 1.4, values = NA_real_, inside = T)

land_r <- 
  land %>% 
  st_rasterize(new_grid) #%>% 
  # st_warp(new_grid)




# load each member

s_all_mem <- 
  map(mems %>% set_names(), function(mem) {
    
    print(mem)
    
    r <- 
      str_glue("gcloud storage ls gs://cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/{prov}/{var}/{grid}/") %>% 
      system(intern = T) %>% 
      str_sub(start = 6, end = -2)
    
    dsn <- 
      str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/{r}"/:{var}.zarr/')
    
    pr_gcm <- 
      read_mdim(dsn, 
                offset = c(off_x, off_y, time_i-1), 
                count = c(count_x, count_y, NA))
      
    
    # regularize dimensions
    pr_gcm <- 
      pr_gcm %>%
      st_set_dimensions(2, values = seq(st_bbox(pr_gcm)[2],
                                        st_bbox(pr_gcm)[4],
                                        length.out = dim(pr_gcm)[2])) %>%
      st_set_crs(4326) %>% 
      st_warp(new_grid)
    
    # change units
    pr_gcm <- 
      pr_gcm %>% 
      mutate(pr = 
               pr %>% 
               units::set_units(kg/m^2/d) %>% 
               units::set_units(NULL))
    
    
    
    # apply land mask
    pr_gcm[is.na(land_r)] <- NA
    
    return(pr_gcm)
    
  })



# rolling window
s_all_mem_seas <- 
  s_all_mem %>% 
  map(function(s) {
    
    r <- 
      s %>% 
      st_apply(c(1,2), function(x) {
        
        if(any(is.na(x))){
          rep(NA, length(x))
        } else {
          
          zoo::rollmean(x, 3, align = "right", fill = NA)
          
        }
        
      },
      .fname = "time") %>% 
      aperm(c(2,3,1))
    
    st_dimensions(r) <- st_dimensions(s)
    return(r)
  })





# APPROACH 1: NORMALIZE EACH MEMBER ON ITS OWN
{ 
# # normalize
# s_all_mem_norm <- 
#   s_all_mem_seas %>% 
#   map(function(s) {
#     
#     r <- 
#       s %>% 
#       st_apply(c(1,2), function(x) {
#         
#         if(all(is.na(x))){
#           rep(NA, length(x))
#         } else {
#           
#           x %>% 
#             matrix(ncol = 12, byrow = T) %>% 
#             apply(2, function(xx) {
#               ecdf(xx)(xx)
#             }) %>% 
#             t() %>% 
#             as.vector()
#           
#           
#           # x_mean <- 
#           #   x %>% 
#           #   matrix(ncol = 12, byrow = T) %>% 
#           #   apply(2, mean, na.rm = T)
#           # 
#           # x_sd <- 
#           #   x %>% 
#           #   matrix(ncol = 12, byrow = T) %>% 
#           #   apply(2, sd, na.rm = T)
#           # 
#           # (x-x_mean)/x_sd
#           
#           
#         }
#         
#       },
#       .fname = "time") %>% 
#       aperm(c(2,3,1))
#     
#     st_dimensions(r) <- st_dimensions(s)
#     return(r)
#   })
# 
# 
# # tb for clustering
# tb <-
#   s_all_mem_norm %>%
#   imap_dfr(function(s, i) {
# 
#     s %>%
#       as_tibble() %>%
#       group_by(time) %>%
#       mutate(time = as_date(time),
#              r = row_number()) %>%
#       ungroup() %>%
#       mutate(mem = i) %>% 
#       drop_na()
# 
#   }) %>% 
#   mutate(time = paste0(mem, "_", time)) %>% 
#   select(-mem)
# 
# tb_for_clustering <-
#   tb %>%
#   select(-x, -y) %>%
#   pivot_wider(names_from = r,
#               names_prefix = "g_",
#               values_from = "pr")
# 
# 
# set.seed(1)
# km <-
#   tb_for_clustering %>%
#   select(-time) %>%
#   kmeans(4, iter.max = 50)
# 
# # plot
# tb_for_clustering %>% 
#   select(time) %>%
#   
#   mutate(clus = km$cluster) %>% 
#   right_join(tb, by = "time") %>% 
#   
#   group_by(clus, x, y) %>% 
#   summarise(pr = mean(pr)) %>% 
#   
#   ggplot(aes(x, y, fill = pr)) +
#   geom_raster() +
#   coord_equal() +
#   colorspace::scale_fill_continuous_divergingx(mid = 0.5,
#                                                rev = T,
#                                                #limits = c(-1,1), oob = scales::squish
#   ) +
#   facet_wrap(~clus)
# 
# km$cluster %>% table()
# 
# # save results
# tibble(
#   tb_for_clustering %>%
#     select(1:2) %>%
#     mutate(time =
#              str_glue("{time}-01-01") %>%
#              as_date() %>%
#              {. + months(10)}),
# 
#   clus = km$cluster
# 
# ) %>%
#   write_csv("gcm_clusters_4.csv")
# 
# 
# predict.kmeans <-
#   function(x, newdata) {
#     a <- apply(newdata, 1, function(r) which.min(colSums((t(x$centers) - r)^2)))
#     return(a)
#   }
# 
# predict.kmeans(km,
#                tb_for_clustering %>%
#                  filter(time == 1850,
#                         mem == "r1i1p1f2") %>%
#                  select(-time, -mem))

}



# APPROACH TWO: NORMALIZE ALL MEMBERS TOGETHER

s_all_mem_norm <- 
  s_all_mem_seas %>% 
  {do.call(c, c(., along = "time"))} %>% 
  st_apply(c(1,2), function(x) {
    
    if(all(is.na(x))){
      rep(NA, length(x))
    } else {
      
      # quantiles
      x %>%
        matrix(ncol = 12, byrow = T) %>%
        apply(2, function(xx) {
          ecdf(xx)(xx)
        }) %>%
        t() %>%
        as.vector()
      
      
      # # standarize mean
      # xx <- 
      #   x %>%
      #   matrix(ncol = 12, byrow = T)
      # 
      # x_mean <-
      #   xx %>%
      #   apply(2, mean, na.rm = T)
      # 
      # x_sd <-
      #   xx %>%
      #   apply(2, sd, na.rm = T)
      # 
      # (x-x_mean)/x_sd
      
    }
    
  },
  .fname = "time") %>% 
  aperm(c(2,3,1))    


tm_vector <- 
  
  tidyr::expand_grid(
    mem = 
      mems,
    time = s_all_mem[[1]] %>% st_get_dimension_values(3) %>% as_date()
  ) %>% 
  
  mutate(tm = str_glue("{mem}_{time}")) %>% 
  pull(tm)


tb <- 
  s_all_mem_norm %>% 
  st_set_dimensions(3, values = tm_vector) %>% 
  as_tibble() %>% 
  group_by(time) %>% 
  mutate(r = row_number()) %>% 
  ungroup() 
  

tb_for_clustering <- 
  tb %>% 
  drop_na() %>% 
  select(-x,-y) %>% 
  
  pivot_wider(names_from = r, 
              names_prefix = "g_",
              values_from = "pr")


set.seed(1)
km <- 
  tb_for_clustering %>%
  select(-time) %>% 
  kmeans(4, iter.max = 50)

km$cluster %>% table()

res <- 
  tibble(
    tb_for_clustering %>% 
      select(time) %>% 
      separate_wider_delim(time, "_", names = c("mem", "time")),
    
    clus = km$cluster,
    
    dist = distance_centroid(km,      # see function below
                             tb_for_clustering %>%
                               select(-time))) %>% 
  group_by(clus) %>% 
  slice_min(dist, n = 800) %>% # closest 800 per cluster
  ungroup()

res %>%
  write_csv("gcm_clusters_4_v2.csv")




# PLOTS

s_all_mem_norm %>% 
  st_set_dimensions(3, values = tm_vector) %>% 
  as_tibble() %>% 
  separate_wider_delim(time, "_", names = c("mem", "time")) %>% 
  left_join(res, by = c("time", "mem")) %>% 
  group_by(time, mem, clus) %>% 
  nest() %>% 
  filter(clus == 3) %>% 
  ungroup() %>% 
  slice_sample(n = 1) %>% 
  .$data %>% 
  .[[1]] %>% 
  ggplot(aes(x,y, fill = pr)) +
  geom_raster() +
  coord_equal() +
  colorspace::scale_fill_continuous_divergingx(mid = 0.5, rev = T, na.value = "transparent")


s_all_mem_norm %>% 
  st_set_dimensions(3, values = tm_vector) %>% 
  as_tibble() %>% 
  separate_wider_delim(time, "_", names = c("mem", "time")) %>% 
  left_join(res, by = c("time", "mem")) %>% 
  group_by(clus, x, y) %>% 
  summarise(pr = median(pr)) %>% 
  drop_na() %>% 
  ggplot(aes(x, y, fill = pr)) +
  geom_raster() +
  coord_equal() +
  colorspace::scale_fill_continuous_divergingx(#mid = 0.5,
                                               rev = T,
                                               na.value = "transparent") +
  facet_wrap(~clus)
  


# Predict to what cluster a new obs belongs

x <- km
r <- unname(unlist(tb_for_clustering[10, -1])) # a single row
which.min(colSums((t(x$centers) - r)^2))

predict_cluster <-
  function(x, newdata) {
    a <- apply(newdata, 1, function(r) which.min(colSums((t(x$centers) - r)^2)))
    return(a)
  }

predict_cluster(km,
                tb_for_clustering %>%
                  filter(time == "r1i1p1f2_1940-10-16") %>%
                  select(-time))



# Determine distance from centroid

distance_centroid <-
  function(x, newdata) {
    a <- apply(newdata, 1, function(r) min(colSums((t(x$centers) - r)^2)))
    return(a)
  }

distance_centroid(km,
                  tb_for_clustering %>%
                    filter(time == "r1i1p1f2_1940-10-16") %>%
                    select(-time))


inner_join(
  s_all_mem_norm %>% 
    st_set_dimensions(3, values = tm_vector) %>% 
    as_tibble() %>% 
    separate_wider_delim(time, "_", names = c("mem", "time")),
  
  res %>% 
    mutate(dist = distance_centroid(km,
                               tb_for_clustering %>%
                                 select(-time)),
            id = row_number()) %>% 
    group_by(clus) %>% 
    slice_min(dist, n = 800) %>%  # closest
    ungroup(),
  
  by = c("mem", "time")
  
) %>%
  group_by(clus, x, y) %>% 
  summarise(pr = mean(pr)) %>% 
  drop_na() %>% 
  ggplot(aes(x, y, fill = pr)) +
  geom_raster() +
  coord_equal() +
  colorspace::scale_fill_continuous_divergingx(#mid = 0.5,
                                               rev = T,
                                               na.value = "transparent") +
  facet_wrap(~clus)


# 
# 










