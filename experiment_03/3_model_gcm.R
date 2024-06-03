
library(tidyverse)
library(tidymodels)
library(doFuture)

registerDoFuture()
options(future.fork.enable = T)
plan(multicore)



# response

r <- 
  read_csv("experiment_03/response_gcm.csv")

r_all <- 
  r %>% 
  group_by(mem) %>% 
  group_split() %>% 
  map(mutate, pr3 = zoo::rollsum(pr, 3, fill = NA, align = "left")) %>% 
  bind_rows() %>% 
  drop_na() %>% 
  group_by(mon = month(time)) %>% 
  mutate(pr3_anom = pr3 - mean(pr3)) %>% 
  ungroup() %>% 
  select(-mon)



# pred

p <- 
  read_csv("experiment_03/predictors_gcm.csv")


p_all <-
  p %>% 
  left_join(r, by = c("mem", "time")) %>%
  group_by(mem) %>%
  group_split() %>%
  map(mutate, across(-c(1:2), lag, .names = "{.col}_lag1")) %>%
  map(mutate, across(-c(1:2) & !contains("lag"), lag, n = 2, .names = "{.col}_lag2")) %>%
  map(mutate, across(-c(1:2) & !contains("lag"), lag, n = 3, .names = "{.col}_lag3")) %>%

  # map(mutate, pr_lag1 = lag(pr, n = 1)) %>%  
  # map(mutate, pr_lag2 = lag(pr, n = 2)) %>% 
  # 
  bind_rows()




tb_f <- 
  
  inner_join(
    r_all %>% 
      mutate(time = time - months(3)) %>% 
      select(mem, time, r = pr3_anom),                                                  # ***************
    
    p_all,
    
    by = c("mem", "time")
  ) %>% 
  drop_na() %>% 
  mutate(mon = factor(month(time))) %>% 
  select(!contains("zg500"))



train <- tb_f %>% filter(mem != "r11i1p1f2")
test <- tb_f %>% filter(mem == "r11i1p1f2")



spec <- 
  rand_forest(trees = 3000, min_n = 50) %>% 
  set_mode("regression") %>% 
  set_engine("ranger", num.threads = 10)

rec <- 
  recipe(r ~ ., data = train) %>% 
  update_role(c(mem, time), new_role = "id_var")

wf <- 
  workflow() %>% 
  add_model(spec) %>% 
  add_recipe(rec)

f <- 
  fit(wf, data = train)

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



