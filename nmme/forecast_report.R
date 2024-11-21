
# THIS SCRIPT STILL USES DATA FROM NCEI:
# NOT BIAS CORRECTED


# scrip parameters:
inic <- "202410" # initial forecast 
dir_data <- "/mnt/pers_disk/tmp2" # temporary local directory 
dir_tables <- "/mnt/bucket_mine/misc_data/temporary" # directory form where to pull tables




# SETUP -----------------------------------------------------------------------

library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
options(future.rng.onMisuse = "ignore")
plan(multicore)

sf_use_s2(F)

# delete temp dir if exists (clean slate)
if (exists(dir_data)) {
  fs::dir_delete(dir_data)
}


dir_root_fc <- "https://ftp.cpc.ncep.noaa.gov/NMME/realtime_anom"


# copy, import, delete countries
str_glue("gsutil -m cp gs://clim_data_reg_useast1/misc_data/admin_units/ne_110m_admin_0/* .") %>% 
  system()

countries <-
  "ne_110m_admin_0_countries.shp" %>%
  read_sf() %>% 
  select(ADMIN) %>% 
  mutate(id = row_number())

fs::dir_ls(regexp = "ne_110m_admin") %>% 
  fs::file_delete()


# target dates
target_dates <- 
  str_glue("{inic}01") %>% 
  as_date() %>% 
  seq(as_date(. + months(5)), by = "1 month")



# DOWNLOAD ERA CLIMATOLOGIES --------------------------------------------------

fs::dir_create(dir_data)

# file names for each month's gamma parameters
ff <- 
  "gsutil ls gs://clim_data_reg_useast1/era5/climatologies/" %>% 
  system(intern = T) %>% 
  str_subset("gamma-params") %>% 
  str_subset("1982")

# download files
ff %>% 
  future_walk(\(f){
    
    str_glue("gsutil cp {f} {dir_data}") %>% 
      system(ignore.stdout = T, ignore.stderr = T)
    
  })

# update files' path
ff <- 
  ff %>% 
  str_replace("gs://clim_data_reg_useast1/era5/climatologies", dir_data)



# INGEST ANOMALIES FORECASTS --------------------------------------------------

models <- c("CFSv2", 
            "CanESM5", 
            "GEM5.2_NEMO", 
            "GFDL_SPEAR", 
            "NASA_GEOS5v2", 
            "NCAR_CCSM4", 
            "NCAR_CESM1")


s_fcst <- 
  models %>% 
  map(\(mod){ # loop through models
    
    print(str_glue("IMPORTING {mod}"))
    
    # read data
    s <- 
      read_mdim(
        str_glue("/vsicurl/{dir_root_fc}/{mod}/{inic}0800/{mod}.prate.{inic}.anom.nc"),
        count = c(NA,NA,6,10) # 6 lead time steps, 10 members
      ) %>% 
      
      suppressWarnings() %>% 
      
      # fix dates
      st_set_dimensions("target",
                        values = target_dates) %>% 
      
      # convert units
      mutate(fcst = fcst %>% units::set_units(mm/d)) %>% 
      units::drop_units()
    
    return(s)
    
  })

# merge all models
s_fcst <- 
  do.call(c, c(s_fcst, along = "ensmem")) %>% 
  setNames("anom")



# IMPORT CLIMATOLOGIES --------------------------------------------------------

# months to import
mons <- 
  target_dates %>% 
  str_sub(6,7)

s_clim <- 
  mons %>%
  set_names() %>% 
  future_map(\(mon){ # loop through months
    
    # read data
    s <- 
      ff %>% 
      str_subset(mon) %>% 
      read_ncdf() %>% 
      suppressMessages()
    
    # resample to match grid
    s <- 
      s %>% 
      st_warp(s_fcst %>% slice(target, 1) %>% slice(ensmem, 1))
    
    return(s)
    
  })




# LAND MASK -------------------------------------------------------------------

# Rasterize the countries polygon and regrid to 0:360

countries_r <- 
  countries %>% 
  select(id) %>% 
  st_rasterize(
    st_as_stars(st_bbox(),
                xlim = c(-180, 180),
                ylim = c(st_bbox(s_fcst)[2], st_bbox(s_fcst)[4]),
                dx = 0.5,
                values = 0)
  ) %>%
  st_warp(st_as_stars(st_bbox(),
                      xlim = c(st_bbox(s_fcst)[1], st_bbox(s_fcst)[3]),
                      ylim = c(st_bbox(s_fcst)[2], st_bbox(s_fcst)[4]),
                      dx = 0.5,
                      values = 0))


# resample to match grid of forecast
countries_r <- 
  countries_r %>% 
  st_warp(s_fcst %>% slice(target, 1) %>% slice(ensmem, 1),
          method = "mode",
          use_gdal = T) %>% 
  setNames("id")

# remove oceans
countries_r[countries_r == 0] <- NA

# homogeneize dimensions
st_dimensions(countries_r) <- st_dimensions(s_fcst)[1:2]





# FIRST CALCULATIONS ----------------------------------------------------------

s_fcstexp <- 
  seq(mons) %>% 
  set_names(mons) %>% 
  map(\(mon_in){ # loop through months
    
    # select month's gamma params (climatology)
    s_clim_1mon <- s_clim[[mon_in]]
    s_clim_1mon[s_clim_1mon == -9999] <- NA
    
    # calculate mean precip
    s_mean_pr <- 
      s_clim_1mon %>% 
      mutate(mean = shape * 1/rate)
    
    # select month's forecasted anomaly (all members)
    s_fcst_1mon <- 
      s_fcst %>% 
      slice(target, mon_in)
    
    
    # calculate stats for each member
    s_fcstexp_1mon <- 
      map(seq(dim(s_fcst_1mon)[3]), \(mem){
        
        # subset 1 member from month's forecasted anomaly
        s_1mem <- 
          s_fcst_1mon %>% 
          slice(ensmem, mem)
        
        # calculations
        s <- 
          s_mean_pr %>% 
          c(s_1mem) %>% 
          
          # forecasted total precip amount
          mutate(fcst = mean + anom) %>% 
          
          # quintile in relation to climatology
          mutate(quint = pgamma(fcst, shape, rate) %>% 
                   cut(seq(0,1,0.2), labels = F, include.lowest = T)) %>% 
          
          # anomaly as percentage relative to climatology
          mutate(anom_pct = anom/mean) %>% 
          
          select(anom, anom_pct, fcst, quint)
        
        return(s)
        
      })
    
    # merge all members
    s_fcstexp_1mon <- 
      do.call(c, c(s_fcstexp_1mon, along = "ensmember"))
    
    return(s_fcstexp_1mon)
    
  })




# SECOND CALCULATIONS ---------------------------------------------------------
# (generates data for final plots)

data_f_plots <- 
  seq(mons) %>% 
  set_names(mons) %>% 
  map(\(mon_in){ # loop through months
    
    
    # ******************
    
    # 1. MEAN PRECIP
    
    # select month's gamma params (climatology)
    s_clim_1mon <- s_clim[[mon_in]]
    s_clim_1mon[s_clim_1mon == -9999] <- NA
    
    # calculate mean precip
    s_mean_pr <- 
      s_clim_1mon %>% 
      mutate(mean = shape * 1/rate,
             mean = days_in_month(target_dates[mon_in]) * mean) %>% # whole month
      select(mean)
    
    # mask land
    s_mean_pr[is.na(countries_r)] <- NA
    
    
    # ******************
    
    # 2. MEAN ANOMALY PERCENT
    
    # select month's forecast calculations (from previous section)
    s_fcst_1mon <- 
      s_fcstexp[[mon_in]]
    
    # mask land
    s_fcst_1mon[is.na(countries_r)] <- NA
    
    # calculate ensemble mean and agreement
    s_mean_anom <- 
      s_fcst_1mon %>% 
      select(anom_pct) %>%
      st_apply(c(1,2), \(x){
        
        if(all(is.na(x))){
          
          mean_anom <- NA
          agree <- NA
          
        } else {
          
          x <- x[!is.na(x)]
          mean_anom <- mean(x)*100
          agree <- sd(x)*100
          
        }
        
        c(mean_anom = mean_anom, agree = agree)
        
      },
      .fname = "q") %>% 
      split("q")
    
    
    # ******************
    
    # 3. QUINTILES
    
    # calculate mode of quintiles and agreement
    s_quint <- 
      s_fcst_1mon %>% 
      select(quint) %>% 
      st_apply(c(1,2), \(x){
        
        if(all(is.na(x))){
          
          mode <- NA
          agree <- NA
          
        } else {
          
          x <- x[!is.na(x)]
          mode <- seq(5)[which.max(tabulate(match(x, seq(5))))]
          agree <- unname(sum(x == mode)/length(x)*100)
          
        }
        
        c(mode = mode, agree = agree)
        
      },
      .fname = "q") %>% 
      split("q")
    
    
    # ******************
    
    # results
    r <- list(mean_pr = s_mean_pr,
              mean_anom = s_mean_anom,
              quint = s_quint)
    
    return(r)
    
    
  })





fs::dir_delete(dir_data)








# PART 2: FIGURES AND TABLES --------------------------------------------------

country_name <- "India" 
# run `sort(countries$ADMIN)` to see names 


# *****************

# country raster
country_r <- 
  countries_r[countries_r == countries[countries$ADMIN == country_name, ]$id] %>% 
  mutate(a = if_else(is.na(id), NA, 1)) %>% 
  select(a)

# country_polygon
country_pol <- 
  country_r %>% 
  st_as_sf(as_points = F, merge = T) %>% 
  summarize()

country_bbox <- 
  country_pol %>% 
  st_bbox()


if(country_bbox[3]-country_bbox[1] < 60 & country_bbox[4]-country_bbox[2] < 60) {
  
  country_centroid <- 
    country_pol %>% 
    st_centroid() %>% 
    st_coordinates() %>% 
    .[1,] %>% 
    unname()
  
  country_bbox <- 
    st_bbox(c(xmin = country_centroid[1]-30,
              ymin = country_centroid[2]-30,
              xmax = country_centroid[1]+30,
              ymax = country_centroid[2]+30),
            crs = 4326)
  
} else {
  
  x_length <- country_bbox[3]-country_bbox[1]
  y_length <- country_bbox[4]-country_bbox[2]
  longest_side <- which.max(c(x_length, y_length))
  dif_length <- c(x_length, y_length)[longest_side] - c(x_length, y_length)[c(1,2) != longest_side]

  if(longest_side == 1){
    country_bbox[2] <- country_bbox[2] - dif_length/2
    country_bbox[4] <- country_bbox[4] + dif_length/2
  } else {
    country_bbox[1] <- country_bbox[1] - dif_length/2
    country_bbox[3] <- country_bbox[3] + dif_length/2
  }
  
}


# fortified poligon of land (for background)
tb_land_pol <- 
  countries_r[country_bbox] %>% 
  mutate(a = if_else(is.na(id), NA, 1)) %>% 
  select(a) %>% 
  st_as_sf(as_points = F, merge = T) %>%
  st_coordinates() %>%
  as_tibble() %>% 
  mutate(g = str_c(L1, "_", L2))


figs <- 
  seq(mons) %>% 
  set_names(mons) %>% 
  map(\(mon_in){ # loop through months
    
    # crop country
    data_aoi <- 
      data_f_plots %>% 
      pluck(mon_in) %>% 
      map(\(s) {
        
        s[is.na(country_r)] <- NA
        st_crop(s, country_bbox)
        
      })
    
    
    # *************
    
    # FIG 1: mean precipitation
    
    # legend scale limit
    q90 <- 
      data_aoi[[1]] %>% 
      pull() %>% 
      quantile(0.9, na.rm = T) %>% 
      {round(./20)*20}
    
    
    p1 <-
      data_aoi[[1]] %>%
      as_tibble() %>% 
      ggplot() +
      
      geom_polygon(data = tb_land_pol, aes(X,Y,group = g), fill = "grey80") +
      
      geom_raster(aes(lon, lat, fill = mean)) +
      scale_fill_viridis_b(name = "mm",
                           direction = 1,# rev = F,
                           na.value = "transparent",
                           n.breaks = 7,
                           limits = c(0,unname(q90)),
                           oob = scales::squish,
                           guide = guide_colorsteps(barheight = 0.6,
                                                    barwidth = 12,
                                                    title.position = "top")) +
      coord_equal(expand = F) +
      theme(axis.title = element_blank(),
            legend.position = "bottom") +
      labs(subtitle = "Mean total precipitation")
    
    
    # *************
    
    # FIG 2: forecasted percent anomaly
    
    # legend scale limits
    qq <-
      data_aoi[[2]] %>%
      select(mean_anom) %>%
      pull() %>%
      quantile(c(0.1, 0.9), na.rm = T) %>%
      {round(./5)*5} %>%
      unname() %>% 
      abs() %>% 
      max()
    
    
    tb2 <- 
      data_aoi[[2]] %>% 
      as_tibble()
    
    p2 <-
      tb2 %>% 
      mutate(
        
        agree_label = 
          agree %>%
          cut(breaks = c(0, 30, 60, Inf),
              labels = c("< 30", "30-60", "> 60")) %>%
          factor(levels = c("< 30", "30-60", "> 60"))
        
      ) %>% 
      
      ggplot() +
      
      geom_polygon(data = tb_land_pol, aes(X,Y,group = g), fill = "grey80") +
      
      geom_point(aes(lon, lat, size = agree_label, color = mean_anom)) +
      
      scale_color_fermenter(palette = "RdBu",
                            name = "%",
                            direction = 1,
                            limits = c(-qq, qq),
                            oob = scales::squish,
                            guide = guide_colorsteps(barheight = 0.6,
                                                     barwidth = 10,
                                                     title.position = "top",
                                                     order = 1),
                            n.breaks = 7,
                            show.limits = F) +
      
      scale_size_manual(values = c("< 30" = 2.0, "30-60" = 1, "> 60" = 0.5),      
                        na.translate = F,
                        name = "std. dev. (%)",
                        guide = guide_legend(title.position = "top",
                                             theme = theme(legend.text.position = "bottom",
                                                           legend.key = element_blank(),
                                                           legend.key.height  = unit(3, "mm"),
                                                           legend.key.spacing.x = unit(2.5, "mm"))), 
                        drop = F) +
      
      coord_equal(expand = F) +
      theme(axis.title = element_blank(),
            legend.position = "bottom") +
      labs(subtitle = "Predicted mean percent anomaly")
    
    
    # results as table
    
    tb2_f <- 
      tb2 %>%
      filter(!is.na(mean_anom)) %>% 
      reframe(bin = cut(mean_anom, c(-Inf, seq(-80,80,20), Inf))) %>% 
      count(bin) %>% 
      mutate(percent_area = round(n/sum(n)*100)) %>%
      select(-n) %>% 
      mutate(bin = as.character(bin)) %>% 
      mutate(month = mons[mon_in])
    
    
    # *************
    
    # FIG 3: quintiles
    
    tb3 <- 
      data_aoi[[3]] %>% 
      as_tibble()
    
    p3 <- 
      tb3 %>% 
      mutate(agree_label = 
               agree %>% 
               cut(breaks = c(0, 33,66,100), 
                   labels = c("poor", "fair", "high")) %>% 
               factor(levels = c("poor", "fair", "high"))
      ) %>% 
      
      ggplot() +
      
      geom_polygon(data = tb_land_pol, aes(X,Y,group = g), fill = "grey80") +
      
      geom_point(aes(lon,lat,size = agree_label, color = factor(mode))) +
      
      scale_color_brewer(palette = "Spectral",
                         name = "quintile",
                         direction = 1,
                         na.translate = F,
                         guide = guide_legend(title.position = "top",
                                              order = 1,
                                              theme = theme(legend.text.position = "bottom",
                                                            legend.key = element_blank(),
                                                            legend.key.spacing = unit(0.1, "mm")),
                                              override.aes = list(size = 5, shape = 15))
      ) +
      
      scale_size_manual(values = c("poor" = 0.5, "fair" = 1, "high" = 2.0),                 
                        na.translate = F,
                        name = "agreement",
                        guide = guide_legend(title.position = "top",
                                             theme = theme(legend.text.position = "bottom",
                                                           legend.key = element_blank(),
                                                           legend.key.height  = unit(0.2, "mm"))), 
                        drop = F,
      ) +
      coord_equal(expand = F) +
      theme(axis.title = element_blank(),
            legend.position = "bottom",
            panel.background = element_rect(fill = 'grey90')) +
      labs(subtitle = "Predicted precipitation quintile")
    
    # results as table
    
    tb3_f <- 
      tb3 %>%
      filter(!is.na(mode)) %>% 
      count(mode) %>% 
      mutate(percent_area = round(n/sum(n)*100)) %>%
      select(-n) %>% 
      rename(quintile = mode) %>% 
      mutate(month = mons[mon_in])
    
    
    # *************
    
    # merge 3 plots
    pp <- 
      patchwork::wrap_plots(p1,p2,p3, nrow = 1) + 
      patchwork::plot_annotation(title = str_glue("{month.name[month(target_dates[mon_in])]} {year(target_dates[mon_in])}"))
    
    
    # *************
    
    # return all results (merged figs and tables)
    
    r <- list(pp, tb2_f, tb3_f)
    
    return(r)
    
    
  })


# RESULTS ---------------------------------------------------------------------

# maps

figs[[1]][[1]]
figs[[2]][[1]]
figs[[3]][[1]]
figs[[4]][[1]]
figs[[5]][[1]]
figs[[6]][[1]]


# tables

figs %>% 
  map(pluck,2) %>% 
  bind_rows() %>% 
  write_csv(str_glue("{dir_tables}/{country_name}_anom-precip-percentage_{str_sub(inic, end = 4)}-{str_sub(inic, start = 5, end = 6)}-01_plus_5months.csv"))

figs %>% 
  map(pluck,3) %>% 
  bind_rows() %>% 
  write_csv(str_glue("{dir_tables}/{country_name}_anom-precip-quintile_{str_sub(inic, end = 4)}-{str_sub(inic, start = 5, end = 6)}-01_plus_5months.csv"))




# THE END

