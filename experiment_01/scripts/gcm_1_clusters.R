

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
# jsonlite::fromJSON(ret)$dimensions
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
  st_as_stars(st_bbox(gcm_proxy_ext), dx = 1.25, values = NA)

land_r <- 
  land %>% 
  st_rasterize(new_grid, align = T) %>% 
  st_warp(new_grid)




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
      st_set_crs(4326)
    
    # change units
    pr_gcm <- 
      pr_gcm %>% 
      mutate(pr = 
               pr %>% 
               units::set_units(kg/m^2/d) %>% 
               units::set_units(NULL))
    
    # subset SON and aggregate
    pr_gcm <- 
      pr_gcm %>% 
      filter(month(time) %in% c(9,10,11)) %>% 
      aggregate(by = "1 year", sum) %>% 
      aperm(c(2,3,1))
    
    
    # regrid
    pr_gcm <- 
      pr_gcm %>%
      st_warp(new_grid,
              use_gdal = T, 
              method = "cubic") %>%
      setNames("pr") %>% 
      st_set_dimensions(which = c(1,2,3), names = c("lon", "lat", "time")) %>% 
      st_set_dimensions(3, values = seq(year(time_vector[time_i]), last(year(time_vector))))
    
    pr_gcm[is.na(land_r)] <- NA
    
    return(pr_gcm)
    
  })


# de-trend ?

# s_all_mem[[2]] %>%
#   pull() %>%
#   .[20,20,] %>%
#   plot(type = "l")


# s_all_mem_detrended <- 
#   s_all_mem %>% 
#   map(function(s) {
#     
#     r <- 
#       s %>% 
#       st_apply(c(1,2), function(x) {
#         
#         if(any(is.na(x))){
#           rep(NA, length(x))
#         } else {
#           
#           m <- lm(x ~ seq(x))
#           trend <- coef(m)[2] * seq(x) + coef(m)[1]
#           x - trend
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


# normalize
s_all_mem_norm <- 
  # s_all_mem_detrended %>%
  s_all_mem %>% 
  map(function(s) {
    
    r <- 
      s %>% 
      st_apply(c(1,2), function(x) {
        
        if(any(is.na(x))){
          rep(NA, length(x))
        } else {
          
          ecdf(x)(x)
          
        }
        
      },
      .fname = "time") %>% 
      aperm(c(2,3,1))
    
    st_dimensions(r) <- st_dimensions(s)
    return(r)
  })




# tb for clustering
tb <- 
  s_all_mem_norm %>% 
  imap_dfr(function(s, i) {
    
    s %>% 
      as_tibble() %>% 
      group_by(time) %>% 
      mutate(r = row_number()) %>% 
      ungroup() %>% 
      mutate(mem = i) %>% 
      na.omit()
    
  })


tb_for_clustering <- 
  tb %>%
  select(-lon, -lat) %>% 
  pivot_wider(names_from = r, 
              names_prefix = "g_",
              values_from = "pr")


set.seed(1)
km <- 
  tb_for_clustering %>%
  select(-time, -mem) %>% 
  kmeans(4, iter.max = 50)

# plot
tb_for_clustering %>% 
  select(time, mem) %>% 
  mutate(clus = km$cluster) %>% 
  right_join(tb, by = c("time", "mem")) %>% 
  
  group_by(clus, lon, lat) %>% 
  summarise(pr = mean(pr)) %>% 
  
  ggplot(aes(lon, lat, fill = pr)) +
  geom_raster() +
  coord_equal() +
  colorspace::scale_fill_continuous_divergingx(mid = 0.5,
                                               rev = T,
                                               #limits = c(-1,1), oob = scales::squish
  ) +
  facet_wrap(~clus)

km$cluster %>% table()



# save results
tibble(
  tb_for_clustering %>% 
    select(1:2) %>% 
    mutate(time = 
             str_glue("{time}-01-01") %>% 
             as_date() %>% 
             {. + months(10)}),
  
  clus = km$cluster
  
) %>% 
  write_csv("gcm_clusters_4.csv")




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
# 









# tb_anom <-
#   tb_pr_all_mem %>%
#   group_by(lon, lat) %>%
#   nest() %>%
#   mutate(pr_std = map(data, function(df) {
#     if (any(is.na(df$pr))) {
#       rep(NA, nrow(df))
#     } else {
#       ecdf(df$pr)(df$pr)
#     }
#   })) %>%
#   unnest(c(data, pr_std)) %>%
#   ungroup() %>%
#   
#   na.omit() %>%
#   group_by(time, mem) %>%
#   mutate(gcell = row_number()) %>%
#   ungroup()
# 
# 
# tb_for_clustering <-
#   tb_anom %>%
#   select(gcell, pr_std, time, mem) %>%
#   pivot_wider(names_from = gcell, names_prefix = "loc_", values_from = pr_std)
# 
# 
# 
# set.seed(111)
# km <-
#   kmeans(tb_for_clustering %>% select(-time, -mem),
#          4,
#          iter.max = 50)
# 
# bind_cols(tb_for_clustering %>%
#             select(time, mem),
#           cluster = km$cluster) %>%
#   right_join(tb_anom, by = c("time", "mem")) -> tb_c
# 
# tb_c %>%
#   group_by(cluster, lon, lat) %>%
#   summarize(pr_std = mean(pr_std)) %>%
#   ungroup() -> tb_d
# 
# tb_d %>%
#   ggplot(aes(lon, lat, fill = pr_std)) +
#   geom_raster() +
#   coord_equal() +
#   colorspace::scale_fill_continuous_divergingx(mid = 0.5, rev = T) +
#   facet_wrap(~cluster)
# 
# km$cluster %>% table()










