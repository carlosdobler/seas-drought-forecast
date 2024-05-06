
library(tidyverse)
library(stars)
library(furrr)

plan(multisession)


dir_data <- "/mnt/pers_disk/data"


# DOWNLOAD ----

# fs::dir_create(dir_data)
# 
# "gsutil -m cp gs://clim_data_reg_useast1/era/monthly/mean-daily-precip/* /mnt/pers_disk/data/" %>% 
#   system()


# LOAD DATA ----

s <- 
  dir_data %>% 
  fs::dir_ls(regexp = "precip") %>%
  str_subset(str_flatten(1940:2020, "|")) %>% 
  future_map(read_ncdf, 
             # drc extension
             ncsub = cbind(start = c(45, 335, 1),
                           count = c(90, 85, NA))) %>% 
  suppressMessages() %>% 
  do.call(c, .)


tb_time <- 
  tibble(date = st_get_dimension_values(s, 3),
         yr = str_sub(date, end = 4) %>% as.numeric(),
         mon = str_sub(date, 6,7) %>% as.numeric(),
         dy = str_sub(date, 9,10) %>% as.numeric(),
         seas = case_when(mon %in% c(12, 1, 2) ~ "1", # DJF
                          mon %in% c(3:5) ~ "2", # MAM
                          mon %in% c(6:8) ~ "3", # JJA
                          mon %in% c(9:11) ~ "4"), # SON
         yr_shift = if_else(mon %in% c(1, 2), yr - 1, yr))


seas <- 
  str_glue("{tb_time$yr_shift}_{tb_time$seas}") %>% 
  as.vector()

s_seas <- 
  s %>% 
  st_apply(c(1,2), function(x) {
    if(any(is.na(x))) {
      rep(NA, length(unique(seas)))
    } else {
      aggregate(x, by = list(seas), mean)$x
    }
  }, 
  .fname = "time") %>% 
  aperm(c(2,3,1))

time_seas <- 
  tb_time %>% 
  group_by(yr_shift, seas) %>% 
  mutate(s = str_glue("{first(yr_shift)}-{first(mon)}-{first(dy)}")) %>% 
  pull(s) %>% 
  as.vector() %>% 
  unique() %>% 
  PCICt::as.PCICt(cal = "gregorian")

s_seas <- 
  s_seas %>% 
  st_set_dimensions("time", values = time_seas)

s_winter <- 
  s_seas %>%
  filter(month(time) == 12)
  
s_winter %>% 
  st_apply(c(1,2), function(x) {
    
    # ecdf(x)(x)
    (x - mean(x))/sd(x)
    
  }, .fname = "time") %>% 
  aperm(c(2,3,1)) %>% 
  st_set_dimensions("time", values = st_get_dimension_values(s_winter, 3)) -> s_winter_anom
  
s_winter_anom %>% 
  as_tibble() %>% 
  group_by(yr = year(time)) %>% 
  select(-time) %>% 
  mutate(r = row_number()) %>% 
  ungroup() -> tb_a

tb_a %>% 
  select(-longitude, -latitude) %>% 
  pivot_wider(names_from = r, names_prefix = "loc_", values_from = tp) -> tb_b

kmeans(tb_b %>% select(-yr),
       6,
       iter.max = 50) -> km

bind_cols(yr = tb_b$yr,
          cluster = km$cluster) %>% 
  right_join(tb_a, by = "yr") -> tb_c

tb_c %>% 
  group_by(cluster, longitude, latitude) %>% 
  summarize(tp = mean(tp)) %>% 
  ungroup() -> tb_d

tb_d %>% 
  ggplot(aes(longitude, latitude, fill = tp)) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradient2() +
  facet_wrap(~cluster)

km$cluster %>% table()



  


