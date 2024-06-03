
library(tidyverse)
library(tidymodels)
library(doFuture)

registerDoFuture()
options(future.fork.enable = T)
plan(multicore)


tbs_predictors <- read_rds("experiment_02/predictors_data_gcm.rds")
tbs_response <- read_rds("experiment_02/response_data_gcm.rds")

tb <-
  map2(tbs_predictors, tbs_response, function(p, r) {
    
    inner_join(
      r %>% 
        mutate(time = time - months(2)),
      
      p,
      
      by = c("time", "member")
    )
    
  }) %>% 
  
  bind_rows()
  
  

tb <-
  tb %>%
  mutate(yr = year(time))

tb <- 
  tb %>% 
  select(-loc20_ua850,
         -loc20_va850,
         -loc20_hus850,
         -loc20_zg850)


spec <- 
  rand_forest(trees = 3000) %>% 
  set_mode("regression") %>% 
  set_engine("randomForest")

rec <- 
  recipe(pr ~ ., data = tb) %>% 
  update_role(c(time, yr, 
                member), new_role = "id_var")

wf <- 
  workflow() %>% 
  add_model(spec) %>% 
  add_recipe(rec)

set.seed(1)
folds <- 
  group_vfold_cv(tb, group = yr, v = 10)

f_1 <- 
  fit_resamples(wf, folds, control = control_resamples(save_pred = T))


collect_metrics(f_1, summarize = F)

collect_predictions(f_1) %>% 
  arrange(.row) %>%
  bind_cols(tb %>% select(c(1,3))) %>%
  
  filter(member == "r6i1p1f2") %>%
  
  filter(year(time) > 1970) %>% 
  
  # filter(month(time) == 2) %>% 
  
  
  ggplot(aes(time)) +
  geom_line(aes(y = .pred), color = "red") +
  geom_line(aes(y = pr))


