# make tables for models
library(brms)
library(tidyverse)
library(gt)
library(broom.mixed)

weight <- readRDS("output/haps_weight_mod.RDS")
surv <- readRDS("output/haps_fitness_mod.RDS")

# table for weight
tidy_weight <- tidy(weight) %>% 
                select(term:conf.high) %>% 
                .[c(1,4:7,2,3,10,9,8,11,12:14), ] %>% 
                mutate(term = c("Intercept",
                                "SEL05 (1 copy)",
                                "SEL05 (2 copies)",
                                "SEL07 (1 copy)",
                                "SEL07 (2 copies)",
                                "SEL18 (1 copy)",
                                "SEL18 (2 copies)",
                                "F<sub>ROH</sub>",
                                "Hindleg length",
                                "Sex",
                                "Twin",
                                "Birth Year",
                                "Mother ID",
                                "Residual"
                                )) %>% 
                mutate(add_info = c("", 
                                    rep("categorical", 6),
                                    "z-transformed (x-mean(x))/sd(x)",
                                    "z-transformed (x-mean(x))/sd(x)",
                                    "categorical (0=female, 1=male)",
                                    "categorical (0=singleton, 1=twin)",
                                    "n = 30",
                                    "n = 819",
                                    "")) %>% 
                mutate(effect = c("", rep("Population level/fixed effects", 10),
                          rep("Group level/random effects (standard deviation)", 2),
                          "")) %>% 
                select(7, 1:6) %>% 
                mutate_if(is.numeric, function(x) round(x, 3)) %>% 
                select(-std.error) %>% 
                setNames(c("effect", "Term", "Post.Mean", "CI (2.5%)", "CI (97.5%)", "Info"))

tidy_weight %>% gt(
        rowname_col = "term",
        groupname_col = "effect") %>% 
        tab_style(
                style = cell_text( weight = "bold"),
                locations = cells_column_labels(columns = everything())
        ) %>% 
        fmt_markdown(columns = everything()) %>% 
        gtsave("weight_model_table.png", path = "tables/")


# survival model

tidy_surv <- tidy(surv) %>% 
        select(term:conf.high) %>% 
        .[c(1,4:7,2,3,9, 11,8,10,12,13), ] %>% 
        mutate(term = c("Intercept",
                        "SEL05 (1 copy)",
                        "SEL05 (2 copies)",
                        "SEL07 (1 copy)",
                        "SEL07 (2 copies)",
                        "SEL18 (1 copy)",
                        "SEL18 (2 copies)",
                        "F<sub>ROH</sub>",
                        "Hindleg length",
                        "Sex",
                        "Twin",
                        "Birth Year",
                        "Mother ID"
        )) %>% 
        mutate(add_info = c("", 
                            rep("categorical", 6),
                            "z-transformed (x-mean(x))/sd(x)",
                            "z-transformed (x-mean(x))/sd(x)",
                            "categorical (0=female, 1=male)",
                            "categorical (0=singleton, 1=twin)",
                            "n = 30",
                            "n = 819")) %>% 
        mutate(effect = c("", rep("Population level/fixed effects", 10),
                          rep("Group level/random effects (standard deviation)", 2))) %>% 
        select(7, 1:6) %>% 
        mutate_if(is.numeric, function(x) paste0(round(x, 3)," (",round(exp(x), 3), ")")) %>% 
        select(-std.error) %>% 
        setNames(c("effect", "Term", "Post.Mean", "CI (2.5%)", "CI (97.5%)", "Info"))

tidy_surv %>% gt(
        rowname_col = "term",
        groupname_col = "effect") %>% 
        tab_style(
                style = cell_text( weight = "bold"),
                locations = cells_column_labels(columns = everything())
        ) %>% 
        fmt_markdown(columns = everything()) %>% 
        gtsave("surv_model_table.png", path = "tables/")

