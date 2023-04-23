library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(rio)
library(janitor)

df_full_30 <- import(here("data/propofol_metaanalysis_data3.csv")) |> #dataset from Ben Prytherch, Likhvantsev corrected to 30-day outcome
  clean_names() |>
  mutate(year = str_extract(study, "\\d{4}"),
         across(study,as.factor),
         across(contains(c("events","total","years")),as.numeric),
         across("setting",as.factor)) |>
  pivot_longer(cols = contains(c("events","total")),cols_vary = "fastest", names_to = c("group",".value"),names_pattern = "(.*)_(.*)") |>
  mutate(group12 = ifelse(group == "control",-0.5,0.5))

df_full_30_jl <- import(here("data/propofol_metaanalysis_data3.csv")) |> # added Jia-Li study
  clean_names() |>
  mutate(year = str_extract(study, "\\d{4}"),
         across(study,as.factor),
         across(contains(c("events","total","years")),as.numeric),
         across("setting",as.factor)) |>
  left_join(
import(here("data/cardiac_surgery_4.xlsx")) |>
  clean_names() |>
  mutate(year = str_extract(study, "\\d{4}"),
         across(study,as.factor),
         across(contains(c("events","total","years")),as.numeric)) |>
  filter(year==2023) |>
  mutate(setting="car_surg") |>
  rename_with(~ str_replace(., "^control", "comparator"), starts_with("control"))) |> #harmonise naming
  pivot_longer(cols = contains(c("events","total")),cols_vary = "fastest", names_to = c("group",".value"),names_pattern = "(.*)_(.*)")



# model fitting ----

model4_full <-  brm(events | trials(total) ~ 0 + Intercept + group + (1 + group || study), # overall model without interaction with setting
               family = binomial,
               prior = prior(normal(0,1.5),class = "b",coef="Intercept")+
                 prior(normal(0,2.82),class = "b") +
                 prior(normal(0,0.5), class = "sd"),
               data = df_full_30,
               control = list( adapt_delta = 0.95),
               file = here("fits/model4_full.rds"),
               file_refit ="on_change")

summary(model4_full)

plot(model4_full, ask = FALSE)

hypothesis(model4_full, "grouppropofol>0")

model4_full_interaction <-  brm(events | trials(total) ~ 0 + Intercept + group*setting + (1 + group || study), # overall model without interaction with setting
                    family = binomial,
                    prior = prior(normal(0,1.5),class = "b",coef="Intercept")+
                      prior(normal(0,2.82),class = "b") +
                      prior(normal(0,0.5), class = "sd"),
                    data = df_full_30,
                    control = list( adapt_delta = 0.95),
                    file = here("fits/model4_full_interaction.rds"),
                    file_refit ="on_change")

summary(model4_full_interaction)

plot(model4_full_interaction,ask = FALSE)

hypothesis(model4_full_interaction, "grouppropofol>0")

# model including Jie-Li study (and interactions with setting)
model4_full_jl <-  brm(events | trials(total) ~ 0 + Intercept + group*setting + (1 + group || study),
                    family = binomial,
                    prior = prior(normal(0,1.5),class = "b",coef="Intercept")+
                      prior(normal(0,2.82),class = "b") +
                      prior(normal(0,0.5), class = "sd"),
                    data = df_full_30_jl,
                    control = list( adapt_delta = 0.95),
                    file = here("fits/model4_full_jl.rds"),
                    file_refit ="on_change")

summary(model4_full_jl)

plot(model4_full_jl)


# plot relative risks ----

labels_data <- model4_full_jl |>
  epred_rvars(
    newdata = df_full_30_jl |> modelr::data_grid(group, setting, total = 1),
    re_formula = NA
  ) |>
  group_by(setting) |>
  compare_levels(.epred, by = group, fun = `/`) |>
  mutate(pp = mean(.epred > 1)) |>
  median_hdi(.epred) |>
  mutate(across(where(is.numeric),  ~ round(.x, 2)))

model4_full_jl |>
  epred_rvars(newdata = df_full_30_jl |> modelr::data_grid(group,setting,total=1,study="new_study"), allow_new_levels = TRUE ) |>
  group_by(setting) |>
  compare_levels(.epred,by=group, fun = `/`) |>
  ggplot()+
  aes(dist=.epred,y=setting,fill=after_stat(x>1)) +
  stat_halfeye() +
  geom_vline(xintercept = 1,linetype=2) +
  scale_fill_manual(values = c("steelblue","pink")) +
  scale_x_continuous(name = "") +
  scale_y_discrete(name = "setting") +
  theme_linedraw() +
  theme(legend.position = "none") +
  geom_text(aes(label = glue::glue("{.epred} [{.lower}, {.upper}]\n Post. prob. {pp}"),x=Inf),data=labels_3fjl,hjust = "inward")


