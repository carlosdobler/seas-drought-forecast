
library(tidyverse)
library(tidymodels)




clusters <- 
  read_csv("era_clusters_4.csv")

eof <- 
  read_csv("era_eof_ts.csv")

var_names <- 
  eof %>% 
  names() %>% 
  .[-1]

tb <- 
  bind_cols(
    clusters %>% 
      mutate(time = year(time),
             clus = factor(clus)) %>% 
      filter(time < 2023),
    
    eof %>% 
      filter(month(time) == 6) %>% # 3 month lead (9 = SEPT)
      select(-time) %>%  
      rename_with(~str_glue("{.x}_lag3")),
    
    eof %>% 
      filter(month(time) == 5) %>% # 3 month lead (9 = SEPT)
      select(-time) %>%  
      rename_with(~str_glue("{.x}_lag4")),
    
    eof %>% 
      filter(month(time) == 4) %>% # 3 month lead (9 = SEPT)
      select(-time) %>% 
      rename_with(~str_glue("{.x}_lag5")),
    
    by = "time"
  ) %>% 
  select(-time)



set.seed(1)
folds <- 
  tb %>% 
  vfold_cv(v = 6, strata = clus)



spec <- 
  rand_forest(mode = "classification",
              trees = 5000) %>% 
  set_engine("randomForest")
  

set.seed(2)
f <- 
  fit_resamples(spec, clus ~ ., folds)

collect_metrics(f)



