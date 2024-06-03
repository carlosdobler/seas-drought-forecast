library(tidyverse)
library(stars)
library(furrr)

plan(multisession)

source("https://raw.github.com/carlosdobler/spatial-routines/master/cell_pos_and_import_subset.R")

# Load ERA data

pr_files <- 
  "gcloud storage ls gs://clim_data_reg_useast1/era/monthly/mean-daily-precip/*nc" %>% 
  system(intern = T)

f <- pr_files[1]

era_proxy <- 
  read_mdim(str_replace(f, "gs:/", "/vsigs"), 
            count = c(NA,NA,1)) %>% 
  adrop()

off_x <- fn_get_cell_pos(era_proxy, 1, 5)
count_x <- round((fn_get_cell_pos(era_proxy, 1, 37) - off_x)/5)
off_y <- fn_get_cell_pos(era_proxy, 2, 13)
count_y <- round((fn_get_cell_pos(era_proxy, 2, -20) - off_y)/5)

# era_test <- read_mdim(str_replace(f, "gs:/", "/vsigs"),
#                       offset = c(off_x, off_y, 0),
#                       count = c(count_x, count_y, NA),
#                       step = c(5,5,1))

pr_era <- 
  pr_files %>% 
  str_subset(seq(1940, 2014) %>% str_flatten("|")) %>%  
  future_map(function(f) {
    
    read_mdim(str_replace(f, "gs:/", "/vsigs"), 
              offset = c(off_x, off_y, 0),
              count = c(count_x, count_y, NA),
              step = c(5,5,1))
    
  }) 

pr_era <- 
  pr_era %>% 
  do.call(c, .)

pr_era <- 
  pr_era %>% 
  mutate(era = 
           tp %>% 
           units::set_units(mm) %>% 
           units::set_units(NULL)) %>% 
  select(era)

# r <- 
#   st_bbox(c(xmin = 0+1.2,
#             xmax = 360-1.2*2,
#             ymin = -90+1.2,
#             ymax = 90-1.2),
#           crs = 4326) %>% 
#   st_as_stars(dx = 1.2)
# 
# pr_era_bl <- 
#   pr_era %>% 
#   st_warp(r, use_gdal = T, method = "bilinear")




# Load CNRM-CM6-1 data

gcm <- "CNRM-CERFACS/CNRM-CM6-1"
mem <- "r1i1p1f2"
var <- "pr"


root_dsn <- str_glue('ZARR:"/vsicurl/https://storage.googleapis.com/cmip6/CMIP6/CMIP/{gcm}/historical/{mem}/Amon/{var}/gr/v20180917"/')

# Time vector
dsn <- str_glue('{root_dsn}:/time')

ori <- 
  gdal_utils("info", dsn, quiet = T) %>% 
  str_extract("time_origin.*\n") %>% 
  str_extract("....-..-..")

time_vector <- 
  read_stars(dsn) %>% 
  pull() %>% 
  as.vector() %>% 
  as_date(origin = ori)

# Data
dsn <- str_glue('{root_dsn}:zg.zarr/')

gcm_proxy <- 
  read_mdim(dsn,
            count = c(NA,NA,1)) %>% 
  adrop()

off_x <- fn_get_cell_pos(gcm_proxy, 1, 5)
count_x <- fn_get_cell_pos(gcm_proxy, 1, 37) - off_x
off_y <- fn_get_cell_pos(gcm_proxy, 2, -20)
count_y <- fn_get_cell_pos(gcm_proxy, 2, 13) - off_y

# gcm_test <- read_mdim(dsn,
#                       offset = c(off_x, off_y, 0),
#                       count = c(count_x, count_y, 1))
# 
# gcm_test %>% 
#   st_set_dimensions(2, values = seq(st_bbox(gcm_test)[2], 
#                                     st_bbox(gcm_test)[4], 
#                                     length.out = 24)) %>% 
#   st_set_crs(4326) %>% 
#   mapview::mapview()


y <- first(which(year(time_vector) == 1940))-1



pr_gcm <- 
  read_mdim(dsn, 
            offset = c(off_x, off_y, y), 
            count = c(count_x, count_y, NA))

pr_gcm <- 
  pr_gcm %>%
  st_set_dimensions(2, values = seq(st_bbox(pr_gcm)[2],
                                    st_bbox(pr_gcm)[4],
                                    length.out = dim(pr_gcm)[2])) %>%
  st_set_crs(4326)

pr_gcm <- 
  pr_gcm %>% 
  mutate(pr = 
           pr %>% 
           units::set_units(kg/m^2/d) %>% 
           units::set_units(NULL))

pr_gcm_bl <- 
  pr_gcm %>% 
  st_warp(pr_era, use_gdal = T, method = "bilinear") %>% 
  setNames("gcm")

st_dimensions(pr_gcm_bl) <- st_dimensions(pr_era)



s <- 
  c(pr_era, pr_gcm_bl, along = "source")





s %>% pull() %>% .[5,25,,] -> x

s %>% 
  st_apply(c(1,2), function(x) {
    
    if (any(is.na(x[,2]))) {
      
      c(dif_ = rep(NA, 12), t_ = rep(NA, 12)) 
      
    } else {
      
      tb <- 
        tibble(era = x[,1],
               gcm = x[,2],
               time = tail(time_vector, -y)) %>% 
        pivot_longer(1:2) %>% 
        group_by(mon = month(time)) %>%
        nest() %>% #.$data %>% .[[1]] %>% 
        mutate(d = map(data, function(df) {
          
          dif <- 
            df %>% 
            group_by(name) %>% 
            summarize(value = mean(value)) %>% 
            pull(value) %>% 
            diff()
          
          t_test <- 
            t.test(df$value ~ df$name)$p.value
          
          tibble(dif = dif, t = t_test)
          
        })) %>% 
        unnest(d) 
      
      c(dif_ = tb$dif, t_ = tb$t)
      
    }
    
  }, 
  .fname = "s") -> a

a %>% 
  split("s") %>% 
  as_tibble() %>% 
  pivot_longer(dif_1:t_12) %>% 
  separate_wider_delim(name, "_", names = c("v", "mon")) %>% 
  mutate(mon = as.integer(mon)) %>% 
  pivot_wider(names_from = v, values_from = value) %>% 
  na.omit() %>% 
  mutate(t = if_else(t <= 0.05 & abs(dif) > 1, "1", NA)) %>% 
  
  ggplot(aes(longitude, latitude)) +
  geom_raster(aes(fill = dif)) +
  geom_point(aes(shape = t), show.legend = F) +
  colorspace::scale_fill_binned_diverging(n.breaks = 7, limits = c(-6,6), oob = scales::squish) +
  scale_shape_manual(values = c("1" = 4)) +
  coord_equal() +
  facet_wrap(~mon, ncol = 4)
  









s %>% pull() %>% .[15,15,,] -> x

# seasonal cycle
s %>% 
  st_apply(c(1,2), function(x) {
    
    if (any(is.na(x[,2]))) {
      
      r_ <- NA
      p_ <- NA
      
    } else {
      
      era_seas <- stl(ts(x[,1],frequency = 12), s.window = "periodic")$time.series[1:12, "seasonal"]
      gcm_seas <- stl(ts(x[,2],frequency = 12), s.window = "periodic")$time.series[1:12, "seasonal"]
      
      r <- cor.test(era_seas, gcm_seas)
      
      r_ <- r$estimate
      p_ <- r$p.value
      
    }
    
    res <- c(r_, p_) %>% set_names(c("r", "p"))
    return(res)
    
  },
  .fname = "s") -> r

r %>% 
  split("s") %>% 
  as_tibble() %>%
  na.omit() %>% 
  mutate(p = if_else(p > 0.05, "1", NA)) %>% 
  
  ggplot(aes(longitude, latitude)) +
  geom_raster(aes(fill = r)) +
  geom_point(aes(shape = p), show.legend = F) +
  colorspace::scale_fill_continuous_diverging(n.breaks = 7,
                                          rev = T,
                                          limits = c(-1,1)
  ) +
  scale_shape_manual(values = c("1" = 4)) +
  coord_equal()



# wet and dry season bias

s %>% 
  st_apply(c(1,2), function(x) {
    
    if (any(is.na(x[,2]))) {
      
      bias_dry <- NA
      cor_dry <- NA
      bias_wet <- NA
      cor_wet <- NA
      
    } else { 
      
      seasonal_comp <- 
        ts(x[,1],frequency = 12) %>% 
        stl(s.window = "periodic") %>% 
        .$time.series %>% 
        .[,"seasonal"] %>% 
        as.vector()
      
      wet_seas <-
        slider::slide_dbl(seasonal_comp, 
                          sum, 
                          .before = 5,
                          .complete = T) %>% 
        .[13:24] %>% 
        {which(. == max(.))}
      
      wet_seas <- 
        rep(1:12, 2) %>% 
        .[(wet_seas+12-5):(wet_seas+12)]
      
      tb <-
        tibble(era = x[,1],
               gcm = x[,2],
               time = tail(time_vector, -y)) %>% 
        
        mutate(season = if_else(month(time) %in% wet_seas, "wet", "dry")) %>% 
        group_by(yr = year(time), season) %>% 
        summarize(era = sum(era),
                  gcm = sum(gcm)) %>% 
        suppressMessages() %>% 
        
        pivot_wider(names_from = season, values_from = c(era, gcm))
      
      bias_dry <- (mean(tb$gcm_dry) - mean(tb$era_dry))/mean(tb$era_dry)*100
      cor_dry <- cor(sort(tb$era_dry), sort(tb$gcm_dry))
      
      bias_wet <- (mean(tb$gcm_wet) - mean(tb$era_wet))/mean(tb$era_wet)*100
      cor_wet <- cor(sort(tb$era_wet), sort(tb$gcm_wet))
      
    }
    
    res <- 
      c(bias_dry, cor_dry, bias_wet, cor_wet) %>% 
      set_names(c("bias_dry", "cor_dry", "bias_wet", "cor_wet"))
    
    return(res)
    
  },
  .fname = "s") -> rr


rr %>% 
  as_tibble() %>%
  separate_wider_delim(s, "_", names = c("var", "seas")) %>% 
  pivot_wider(names_from = var, values_from = era) %>% 
  na.omit() %>% 
  mutate(cor = if_else(cor < 0.9, "1", NA)) %>%
  
  ggplot(aes(longitude, latitude)) +
  geom_raster(aes(fill = bias)) +
  geom_point(aes(shape = cor), show.legend = F) +
  colorspace::scale_fill_continuous_diverging(limits = c(-100,100), oob = scales::squish) +
  scale_shape_manual(values = c("1" = 4)) +
  coord_equal() +
  facet_wrap(~seas, ncol = 2)
