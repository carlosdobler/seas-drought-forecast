
library(tidyverse)
library(tidymodels)
library(doFuture)

registerDoFuture()
options(future.fork.enable = T)
plan(multicore)



# response

r <- 
  read_csv("experiment_03/data/pr_1loc_gcm.csv") 

r <- 
  r %>% 
  group_by(mem) %>%
  group_split() %>%
  map(mutate, pr3 = zoo::rollmean(pr, 3, fill = NA, align = "right")) %>%
  bind_rows()
  

lead_t <- 0

rr <- 
  r %>% 
  filter(month(time) == 11) %>% 
  select(time, mem, pr_son = pr3) %>% 
  mutate(pr_son_perc = (ecdf(pr_son)(pr_son))*100, # try another norm?
         time = time - months(3+lead_t))
         
# discretize?
  
rr %>%
  mutate(cc = case_when(pr_son_perc < 25 ~ -1,
                        pr_son_perc < 75 ~ 0,
                        TRUE ~ 1)) %>%
  filter(mem == "r1i1p1f2") %>%

  ggplot(aes(time, pr_son)) +
  geom_line(alpha = 0.3) +
  geom_point(aes(color = factor(cc)))
  


# pred

p <- 
  read_csv("experiment_03/data/predictors_gcm.csv")


pp <- 
  r %>% 
  left_join(p, by = c("mem", "time")) %>% 
  
  # not needed for now
  select(-pr) %>%
  group_by(mem) %>%
  group_split() %>%
  
  map(mutate, across(-c(time, mem) & !contains("lag"), lag, n = 2, .names = "{.col}_lag2")) %>%
  map(mutate, across(-c(time, mem) & !contains("lag"), lag, n = 4, .names = "{.col}_lag4")) %>%
  map(mutate, across(-c(time, mem) & !contains("lag"), lag, n = 6, .names = "{.col}_lag6")) %>%

  # .[[1]] %>% filter(year(time) == 1941) %>% select(time, pr3, pr3_lag2, pr3_lag4)
  
  bind_rows()




tb_f <- 
  
  inner_join(
    rr %>% 
      select(mem, time, r = pr_son_perc),                                     # ***************
    
    pp,
    
    by = c("mem", "time")
  ) %>% 
  select(-contains("eof3")) %>% 
  drop_na()


tb_f <- 
  tb_f %>% 
  select(mem, time, r, pr3, tos_io_eof1)



train <- tb_f %>% filter(mem != "r11i1p1f2")
test <- tb_f %>% filter(mem == "r11i1p1f2")

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
  # grid_latin_hypercube(
  grid_regular(
    tree_depth(),
    learn_rate(),
    finalize(mtry(), train),
    loss_reduction(),
    min_n(),
    sample_size = sample_prop(),
    # size = 20
    levels = 3
    )


rec <- 
  recipe(r ~ ., data = train) %>%
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


show_best(tune_res, "rsq")
show_best(tune_res, "rmse")


final_xgb <- finalize_workflow(
  wf,
  select_best(tune_res, "rmse")
)



f <- fit(final_xgb, train)


# spec <- 
#   rand_forest(trees = 3000, min_n = 50) %>% 
#   set_mode("regression") %>% 
#   set_engine("ranger", num.threads = 10)
# 
# rec <- 
#   recipe(r ~ ., data = train) %>% 
#   update_role(c(mem, time), new_role = "id_var")
# 
# wf <- 
#   workflow() %>% 
#   add_model(spec) %>% 
#   add_recipe(rec)
# 
# f <- 
#   fit(wf, data = train)
# 
train %>%
  select(mem, time, r) %>%
  bind_cols(predict(f, train)) %>%

  filter(mem == "r6i1p1f2") %>%

  ggplot(aes(x = time)) +
  geom_line(aes(y = r)) +
  geom_line(aes(y = .pred), color = "red")


test %>%
  select(mem, time, r) %>%
  bind_cols(predict(f, test)) %>%

  ggplot(aes(x = time)) +
  geom_line(aes(y = r)) +
  geom_line(aes(y = .pred), color = "red")



