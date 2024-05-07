
library(tidyverse)
library(stars)
library(furrr)

plan(multisession)

dir_data <- "/mnt/pers_disk/data"



# DOWNLOAD ----

# fs::dir_create(dir_data)
# 
# "gsutil -m cp gs://clim_data_reg_useast1/era/monthly/mean-sst/* /mnt/pers_disk/data/" %>% 
#   system()



# LOAD DATA ----

# only TP
s <- 
  dir_data %>% 
  fs::dir_ls(regexp = "sst") %>% 
  future_map(read_ncdf,
             # TP extent
             ncsub = cbind(start = c(400, 250, 1),
                           count = c(750, 300, NA)),
             downsample = c(3,3,0)) %>% 
  suppressMessages() %>% 
  do.call(c, .)

s <- 
  s %>% 
  filter(year(time) >= 1950,
         year(time) <= 2020)

# PREPARE DATA ----

# 3-month rolling mean

s_3mon <- 
  s %>% 
  st_apply(c(1,2), function(x) {
    
    zoo::rollmean(x, k = 3, align = "right", na.pad = T)
    # slider::slide_dbl(x, mean, .before = 2, .complete = T)
    
  },
  .fname = "time",
  FUTURE = T) %>% 
  aperm(c(2,3,1))

st_dimensions(s_3mon) <- st_dimensions(s)


# no de-trending necessary (no significant trend identified)
# e.g.

# s_3mon %>% 
#   pull() %>% 
#   .[300,300,] %>% 
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


# ts at 130.5, 0 is incomplete
# move this up!
s_anom <- 
  s_anom %>% 
  st_apply(c(1,2), function(x) {
    
    n_na <- mean(is.na(x))
    
    if (n_na > 0.1) {
      rep(NA, length(x))
    } else {
      x
    }
  },
  .fname = "time") %>% 
  aperm(c(2,3,1))

st_dimensions(s_anom) <- st_dimensions(s)




# PREPARE DATA FOR PCA ----

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
  mutate(time = as.character(time))
  

# columns are timesteps (not grid cells)
tb_for_pca <- 
  tb_anom %>% 
  pivot_wider(names_from = "time", names_prefix = "t_", values_from = sst)




# PCA ----

pca <- 
  tb_for_pca %>% 
  select(-c(1,2)) %>% 
  prcomp(scale. = F,
         center = T)



# plot (spatial)
pca$x[,"PC2"] %>% 
  as_tibble() %>% 
  bind_cols(
    s_anom %>% 
      slice(time, 3) %>% 
      as_tibble() %>% 
      na.omit() %>% 
      select(1:2)
  ) %>% 
  ggplot(aes(longitude, latitude, fill = value)) +
  geom_raster() +
  scale_fill_gradient2()


# plot (temporal)
tibble(r = pca$rotation[, "PC2"],
       time = tail(st_get_dimension_values(s, 3), -2) %>% as_date()) %>% 
  
  ggplot(aes(time, r)) +
  geom_line() +
  # scale_y_reverse() +
  scale_x_date(date_breaks = "5 years", 
               date_labels = "%Y", 
               limits = c(as_date("1980-01-01"), 
                          NA)) +
  geom_hline(yintercept = 0, linetype = "2222")

  

# test prdiction
tb_for_pca %>% 
  slice(1) %>% 
  select(-c(1,2)) %>%
  mutate(t_1950_04 = t_1950_03) %>% 
  {predict(pca, .)} %>% .[,"PC2"]

pca$x[1,"PC1"]




# REVERSE !!!! ----
# now columns will be grid cells (not timesteps)

# PREPARE DATA FOR PCA ----

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
  mutate(time = as.character(time))


tb_for_pca <- 
  tb_anom %>%
  group_by(time) %>% 
  mutate(gcell = row_number()) %>% 
  ungroup() %>%
  select(time, sst, gcell) %>% 
  pivot_wider(names_from = "gcell", names_prefix = "g_", values_from = sst)

pca <- 
  tb_for_pca %>% 
  select(-1) %>%
  apply(1, scale, scale = F) %>% 
  aperm(c(2,1)) %>% 
  prcomp(scale. = T,
         center = F) # virtually no diff



# plot (spatial)
pca$rotation[,"PC2"] %>% 
  as_tibble() %>% 
  bind_cols(
    s_anom %>% 
      slice(time, 3) %>% 
      as_tibble() %>% 
      na.omit() %>% 
      select(1:2)
  ) %>% 
  ggplot(aes(longitude, latitude, fill = value)) +
  geom_raster() +
  scale_fill_gradient2()


# plot (temporal)
tibble(r = pca$x[, "PC1"],
       time = tail(st_get_dimension_values(s, 3), -2) %>% as_date()) %>% 
  
  ggplot(aes(time, r)) +
  geom_line() +
  # scale_y_reverse() +
  scale_x_date(date_breaks = "5 years", 
               date_labels = "%Y", 
               # limits = c(as_date("1980-01-01"), 
               #            NA)
               ) +
  geom_hline(yintercept = 0, linetype = "2222")



# test prediction
tb_for_pca %>% 
  slice(1) %>% 
  select(-1) %>% 
  apply(1, scale, scale = F) %>%
  aperm(c(2,1)) %>% 
  {predict(pca, .)} %>% .[,"PC1"]

pca$x[1,"PC1"]







# ***********


tb_anom_f <- 
  tb_anom %>% 
  ungroup() %>% 
  na.omit() %>% 
  group_by(time, mem) %>% 
  mutate(gcell = row_number()) %>% 
  ungroup()
  
  

s_07 <- 
  s %>% 
  filter(month(time) == 7) %>% 
  
  st_apply(c(1,2), \(x) x - mean(x), .fname = "yr") %>% 
  aperm(c(2,3,1)) %>% 
  st_set_dimensions(3, values = 1940:2023)


s_07 %>% 
  # units::drop_units() %>% 
  as_tibble() %>% 
  # mutate(yr = year(time)) %>% 
  # select(-time) %>% 
  pivot_wider(names_from = yr, values_from = sst, names_prefix = "x_") %>% 
  na.omit() -> tb

tb %>% 
  select(-c(1,2)) %>% 
  prcomp(scale. = F, 
         center = F
  ) -> pca_model

tb %>% 
  select(-c(1,2)) %>% 
  svd() -> svd_model

pca_model$sdev^2
screeplot(pca_model)

svd_model$d %>% plot()

pca_model$x[,"PC1"] %>% 
  as_tibble() %>% 
  bind_cols(tb %>% select(1:2)) %>% 
  ggplot(aes(longitude, latitude, fill = value)) +
  geom_raster() +
  scale_fill_gradient2()

svd_model$u[,1] %>% 
  as_tibble() %>% 
  bind_cols(tb %>% select(1:2)) %>% 
  ggplot(aes(longitude, latitude, fill = value)) +
  geom_raster() +
  scale_fill_gradient2()

tibble(r = pca_model$rotation[, "PC1"],
       time = 1940:2023) %>% 
  
  ggplot(aes(time, r)) +
  geom_line() +
  scale_y_reverse()

arrange(r) 

s_07 %>% 
  filter(yr == 2023) %>% 
  plot(breaks = "equal")


pca_model$rotation[, "PC1"] -> a
tb %>% select(-c(1,2)) %>% .[1,] %>% as.matrix() %>% as.vector() -> b
sum(a*b)
pca_model$x[1, "PC1"]


