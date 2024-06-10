
library(tidyverse)
library(tidymodels)
library(doFuture)

registerDoFuture()
options(future.fork.enable = T)
plan(multicore)


tb <- 
  read_csv("experiment_05/data/gcm_tb_eof.csv") %>% 
  select(-pr_drc_eof2, -pr_drc_eof3)


tb_lags <- 
  tb %>% 
  group_by(mem) %>%
  group_split() %>%
  map(mutate, across(-c(1:2) & !contains("lag"), lag, n = 1, .names = "{.col}_lag1")) %>%
  map(mutate, across(-c(1:2) & !contains("lag"), lag, n = 2, .names = "{.col}_lag2")) %>%
  map(mutate, across(-c(1:2) & !contains("lag"), lag, n = 3, .names = "{.col}_lag3")) %>% 
  map_dfr(mutate, pr_drc_eof1_lead3 = lead(pr_drc_eof1, n = 3)) %>% 
  drop_na() %>% 
  mutate(mon = factor(month(time)))
  

train <- tb_lags %>% filter(mem != "r1i1p1f2")
test <- tb_lags %>% filter(mem == "r1i1p1f2")

# spec <- 
#   boost_tree() %>% 
#   set_mode("regression") %>% 
#   set_engine("xgboost")


spec <- 
  boost_tree(tree_depth = tune(),
             trees = 500,
             learn_rate = tune(),
             mtry = tune(),
             loss_reduction = tune(),
             min_n = tune(),
             sample_size = tune()) %>% 
  set_mode("regression") %>% 
  set_engine("xgboost")



r_grid <- 
  grid_latin_hypercube(
    tree_depth(),
    learn_rate(),
    finalize(mtry(), train),
    loss_reduction(),
    min_n(),
    sample_size = sample_prop(),
    size = 20) 

 
rec <- 
  recipe(pr_drc_eof1_lead3 ~ ., data = train) %>%
  step_dummy(mon) %>%
  update_role(c(mem, time), new_role = "id_var")

wf <- 
  workflow() %>% 
  add_model(spec) %>% 
  add_recipe(rec)

folds <- 
  group_vfold_cv(train, group = mem)


set.seed(1)
tune_res <- 
  tune_grid(wf,
            resamples = folds,
            grid = r_grid,
            control = control_grid(save_pred = T))




#  f <- 
#   fit(wf, data = train)
# 
# train %>%
#   select(mem, time, pr_drc_eof1_lead3) %>% 
#   bind_cols(predict(f, train)) %>%
#   
#   filter(mem == "r6i1p1f2") %>%  
#   
#   ggplot(aes(x = time)) +
#   geom_line(aes(y = pr_drc_eof1_lead3)) +
#   geom_line(aes(y = .pred), color = "red")
# 
# 
# test %>%
#   select(mem, time, pr_drc_eof1_lead3) %>% 
#   bind_cols(predict(f, test)) %>%
#   
#   # filter(mem == "r6i1p1f2") %>%  
#   
#   ggplot(aes(x = time)) +
#   geom_line(aes(y = pr_drc_eof1_lead3)) +
#   geom_line(aes(y = .pred), color = "red")
