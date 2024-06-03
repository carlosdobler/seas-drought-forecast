
library(tidyverse)
library(tidymodels)




clusters <- 
  read_csv("gcm_clusters_4.csv")

eof <- 
  read_csv("gcm_tb_eof.csv")

var_names <- 
  eof %>% 
  names() %>% 
  .[-c(1:2)]

tb <- 
  
  clusters %>% 
  mutate(time = year(time),
         clus = factor(clus)) %>% 
  
  inner_join(eof %>% 
               filter(month(time) == 6) %>% # 3 month lead (9 = SEPT)
               mutate(time = year(time)) %>% 
               rename_with(~str_glue("{.x}_lead3"), .cols = -c(1:2)),
             by = c("time", "mem")) %>% 
  
  inner_join(eof %>% 
               filter(month(time) == 4) %>%
               mutate(time = year(time)) %>% 
               rename_with(~str_glue("{.x}_lead5"), .cols = -c(1:2)),
             by = c("time", "mem")) %>% 
  
  inner_join(eof %>% 
               filter(month(time) == 2) %>%
               mutate(time = year(time)) %>% 
               rename_with(~str_glue("{.x}_lead7"), .cols = -c(1:2)),
             by = c("time", "mem"))

    
   



set.seed(1)
folds <- 
  tb %>% 
  vfold_cv(v = 10)



spec <- 
  rand_forest(mode = "classification",
              trees = 5000) %>% 
  set_engine("randomForest")


set.seed(2)
f <- 
  fit_resamples(spec, clus ~ ., folds)

collect_metrics(f)



