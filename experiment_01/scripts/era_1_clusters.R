
library(tidyverse)
library(stars)
library(furrr)

plan(multisession)


dir_data <- "/mnt/pers_disk/data"

land <- 
  "/mnt/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp" %>% 
  read_sf(quiet = T) %>% 
  mutate(a = 1) %>% 
  select(a)


# DOWNLOAD ----

# fs::dir_create(dir_data)
# 
# "gsutil -m cp gs://clim_data_reg_useast1/era/monthly/mean-daily-precip/* /mnt/pers_disk/data/" %>%
#   system()


# LOAD PRECIP DATA ----

s <- 
  dir_data %>% 
  fs::dir_ls(regexp = "precip") %>%
  str_subset(str_flatten(1950:2023, "|")) %>% 
  future_map(read_ncdf, 
             # drc extension
             ncsub = cbind(start = c(45, 335, 1),
                           count = c(90, 85, NA)),
             downsample = c(3,3,0)) %>% 
  suppressMessages() %>% 
  do.call(c, .)

s <- 
  s %>% 
  mutate(tp = tp %>% units::set_units("mm"))



# MASK LAND ----

# rasterize
land <- 
  land %>% 
  st_rasterize(s, align = T) %>% 
  st_warp(s %>% slice(time, 1))
  
# mask
s[is.na(land)] <- NA



# PREPARE DATA ----

# subset SON and aggregate
s_son <- 
  s %>% 
  filter(month(time) %in% c(9,10,11)) %>% 
  aggregate(by = "1 year", sum) %>% 
  aperm(c(2,3,1))



# calculate anomalies

# de-trend first
s_son_notrend <-
  s_son %>%
  st_apply(c(1,2), function(x) {

    if(any(is.na(x))){
      rep(NA, length(x))
    } else {

      m <- lm(x ~ seq(x))
      trend <- coef(m)[2] * seq(x) + coef(m)[1]
      x - trend

    }

  },
  .fname = "time") %>%
  aperm(c(2,3,1))

st_dimensions(s_son_notrend) <- st_dimensions(s_son)


# normalize
s_norm <- 
  s_son_notrend %>%
  st_apply(c(1,2), function(x) {
    
    if(any(is.na(x))){
      rep(NA, length(x))
    } else {
      # (x-mean(x))/sd(x)
      ecdf(x)(x)
      
    }
      
  }, .fname = "time") %>% 
  aperm(c(2,3,1))
  
st_dimensions(s_norm) <- st_dimensions(s_son)



# CLUSTER ----

tb <- 
  s_norm %>%
  as_tibble() %>% 
  group_by(time = year(time)) %>% 
  mutate(r = row_number()) %>% 
  ungroup() %>% 
  na.omit()
  
tb_for_clustering <- 
  tb %>% 
  select(-longitude, -latitude) %>% 
  pivot_wider(names_from = r, 
              names_prefix = "g_",
              values_from = "tp")

set.seed(1)
km <- 
  tb_for_clustering %>%
  select(-time) %>% 
  kmeans(4, iter.max = 50)

# plot
tb_for_clustering %>% 
  select(time) %>% 
  mutate(clus = km$cluster) %>% 
  right_join(tb, by = "time") %>% 
  
  group_by(clus, longitude, latitude) %>% 
  summarise(tp = mean(tp)) %>% 
  
  ggplot(aes(longitude, latitude, fill = tp)) +
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
  time = st_get_dimension_values(s_anom, "time") %>% 
    as_date() %>% 
    {. + months(10)},
  
  clus = km$cluster

  ) %>% 
  write_csv("era_clusters_4.csv")
  


