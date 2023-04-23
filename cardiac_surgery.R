# packag imports ----
library(tidyverse) # for data wrangling and plotting
library(here) # convenience package for locating files
library(brms) # bayesian multilevel  modelling using stan in the backround
library(tidybayes) # for manipulating draws from the models in a tidy/ggplot framwork
library(rio) #for import of datasets
library(janitor) # for cleaning datasets

# data import and cleaning ----
cs_data <- import(here("data/cardiac_surgery.xlsx")) |> #original cardiac surgery dataset
  clean_names() |> #cleans variable names to lowercase and _ seperation
  mutate(year = str_extract(study, "\\d{4}"), #create year variable from study variable
         across(study,as.factor),
         across(contains(c("events","total","years")),as.numeric)) |>
  pivot_longer(cols = contains(c("events","total")),cols_vary = "fastest", names_to = c("group",".value"),names_pattern = "(.*)_(.*)") |> # convert to long format
  mutate(group12 = ifelse(group == "control",-0.5,0.5)) # create new treatment variable (for Smith's model)

cs_data_2 <- import(here("data/cardiac_surgery_2.xlsx")) |> # # corrected Likhvantsev to 30-day outcome
  clean_names() |>
  mutate(year = str_extract(study, "\\d{4}"),
         across(study,as.factor),
         across(contains(c("events","total","years")),as.numeric)) |>
  pivot_longer(cols = contains(c("events","total")),cols_vary = "fastest", names_to = c("group",".value"),names_pattern = "(.*)_(.*)")|>
  mutate(group12 = ifelse(group == "control",-0.5,0.5))

cs_data_3 <- import(here("data/cardiac_surgery_3.xlsx")) |> # corrected to Likhvantsev correct denominator
  clean_names() |>
  mutate(year = str_extract(study, "\\d{4}"),
         across(study,as.factor),
         across(contains(c("events","total","years")),as.numeric)) |>
  pivot_longer(cols = contains(c("events","total")),cols_vary = "fastest", names_to = c("group",".value"),names_pattern = "(.*)_(.*)")|>
  mutate(group12 = ifelse(group == "control",-0.5,0.5))

cs_data_4 <- import(here("data/cardiac_surgery_4.xlsx")) |> # added Jia-Li 2023 study
  clean_names() |>
  mutate(year = str_extract(study, "\\d{4}"),
         across(study,as.factor),
         across(contains(c("events","total","years")),as.numeric)) |>
  pivot_longer(cols = contains(c("events","total")),cols_vary = "fastest", names_to = c("group",".value"),names_pattern = "(.*)_(.*)")|>
  mutate(group12 = ifelse(group == "control",-0.5,0.5))

# check the datasets for Likhvantsev denominators
cs_data |>
  mutate(dataset = "original") |>
  bind_rows (cs_data_2 |>
               mutate(dataset = "30-day")
  ) |>
  bind_rows(cs_data_3 |>
              mutate (dataset = "Correct Denominator")) |>
  filter(grepl( "Likhvantsev",study))



# model fitting ----

model1 <- brm(events | trials(total) ~ 0 + Intercept + group + (1 + group || study), # multilevel (AKA random effects) logistic model, allowing intercept (base event rate) and difference to vary between studies (|| rather than | excludes covariance between the two, which can't be estimated here)
              family = binomial,
              prior = prior(normal(0,1.5),class = "b",coef="Intercept")+ #prior for base rate, different than original analysis (more sensible)
                prior(normal(0,2.82),class = "b") + # prior for difference propofol - comparator on log-odds scale, same as original analysis
                prior(normal(0,0.5), class = "sd"), # prior for standard deviation for random effects, same as original analysis
              data = cs_data, #original data
              control = list( adapt_delta = 0.95),
              file = here("fits/model1.rds"), #saves the model fit
              file_refit ="on_change") # only refits if model or data is changed from what is saved

summary(model1)

plot(model1)

hypothesis(model1, "grouppropofol>0")

model2 <- brm(events | trials(total) ~ 0 + Intercept + group + (1 + group || study),
              family = binomial,
              prior = prior(normal(0,1.5),class = "b",coef="Intercept")+
                prior(normal(0,2.82),class = "b") +
                prior(normal(0,0.5), class = "sd"),
              data = cs_data_2, #correct denominator data
              control = list( adapt_delta = 0.95),
              file = here("fits/model2.rds"),
              file_refit ="on_change")

summary(model2)

plot(model2)

hypothesis(model2, "grouppropofol>0")

model3 <- brm(events | trials(total) ~ 0 + Intercept + group + (1 + group || study),
              family = binomial,
              prior = prior(normal(0,1.5),class = "b",coef="Intercept")+
                prior(normal(0,2.82),class = "b") +
                prior(normal(0,0.5), class = "sd"),
              data = cs_data_3, #30-day outcome data
              control = list( adapt_delta = 0.95),
              file = here("fits/model3.rds"),
              file_refit ="on_change")

summary(model3)

plot(model3)

hypothesis(model3, "grouppropofol>0")

# compare the relative risks estimated from the three datasets ----
# samples are in the rvars forma:
comparison_data <- model1 |>
  epred_rvars(newdata = cs_data |> modelr::data_grid(group, total = 1),
              # adds draws from
              re_formula = NA) |> #means random effects are excluded, ie it takes the population mean not study deviations
  mutate(model = "Original") |>
  bind_rows(
    #adds the same from model 2
    model2 |>
      epred_rvars(cs_data_2 |> modelr::data_grid(group, total = 1), re_formula = NA) |>
      mutate(model = "30-day outcome")

  ) |>
  bind_rows( #now adding model 3
    model3 |>
      epred_rvars(cs_data_3 |> modelr::data_grid(group, total = 1),
                  re_formula = NA) |>
      mutate(model = "Correct numerator")
  ) |>
  group_by(model) |>
  compare_levels(.epred, #computes relative risk, really neat function for computing all kinds of contrasts, first chose which variable
                 by = group, #the by which grouping variable the comparision/contrast should be made
                 fun = `/`) #which function to use, / gives rr in this instance, - would give absolute differences

# create label data
comparison_labels <-
  comparison_data |>
  mutate(pp = mean(.epred > 1)) |> #compute posterior probability for harm from propofol by simply calculating what proportion of the draws have an rr above 1
  median_hdci(.epred) |> # takes .epred (the rr) and computes the median and 95% highest density interval (credible interval)
  mutate(across(where(is.numeric),  ~ round(.x, 2))) #round the numeric variable

# plot
  comparison_data|>
  ggplot()+
  aes(dist=.epred,y=model, # with draws in the rvars format, we specify which variable contains them as "dist", not "x"
      fill=after_stat(x>1)) + #after_stat to give the fill colour difference
  stat_halfeye()+ # nice geom for baysian samples, density + interval plot
  geom_text(aes(label = glue::glue("{.epred} [{.lower}, {.upper}]\n Post. prob. {pp}"),x=Inf),data=comparison_labels,hjust = "inward") +
  scale_x_continuous(name = "Relative risk") +
  scale_y_discrete(name = "") + # remove y-axis title
  scale_fill_manual( values = c("red","steelblue")) +
  theme_linedraw() +
  theme(legend.position = "none") + #remove legend
  labs(title = "Bayesian meta(re)analysis of propofol trials in cardiac surgery",subtitle = "Based on https://doi.org/10.1186/s13054-023-04431-8", caption = "@load_dependent")


  # added model with the Jia-Li 2023 study
  model4 <- brm(events | trials(total) ~ 0 + Intercept + group + (1 + group || study),
                family = binomial,
                prior = prior(normal(0,1.5),class = "b",coef="Intercept")+
                  prior(normal(0,2.82),class = "b") +
                  prior(normal(0,0.5), class = "sd"),
                data = cs_data_4, #30-day outcome data including Jia-Li 2023 study
                control = list( adapt_delta = 0.95),
                file = here("fits/model4.rds"),
                file_refit ="on_change")

  summary(model4)

  plot(model4)

  hypothesis(model4, "grouppropofol>0")
