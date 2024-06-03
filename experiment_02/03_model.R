
library(tidyverse)
library(tidymodels)
library(doFuture)

registerDoFuture()
options(future.fork.enable = T)
plan(multicore)


tb_predictors <- read_csv("experiment_02/predictors_data.csv")
tb_response <- read_csv("experiment_02/response_data.csv")

tb <- 
  inner_join(
    tb_response %>% 
      mutate(time = time - months(2)),
    
    tb_predictors,
    
    by = "time"
  )

tb <- 
  tb %>% 
  mutate(yr = year(time))


spec <- 
  rand_forest(trees = 2000) %>% 
  set_mode("regression") %>% 
  set_engine("randomForest")

rec <- 
  recipe(tp_anom ~ ., data = tb) %>% 
  update_role(c(tp_anom_std, tp, time, yr), new_role = "id_var")

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
  bind_cols(tb %>% select(c(1,2,4))) %>% 
  
  ggplot(aes(time)) +
  # geom_line(aes(y = .pred), color = "red") +
  geom_line(aes(y = tp))


