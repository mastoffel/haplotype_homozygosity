library(lubridate)
library(tidyverse)
library(janitor)
library(here)
library(lme4)
library(broom.mixed)
library(sjPlot)
source("theme_simple.R")

# load fitness data
load("data/fitness_roh.RData") 
fitness <- fitness_data %>% 
        select(id, survival, sheep_year, age, birth_year, sex, mum_id,
               twin, offspring_born, offspring_survived, weight,
               hindleg, horn, father, froh_all, birthwt) %>% 
        group_by(id) %>% 
        mutate(lifespan = max(age)) %>% 
        mutate(id = as.numeric(id),
               focal_year =  as.numeric(as.character(sheep_year)),
               year_mating = focal_year - 1) %>% 
        #filter(age == 0) %>% 
        select(id, froh_all, year_mating, focal_year, sex, offspring_born, mum_id, age) %>% 
        mutate(mum_id = as.numeric(as.character(mum_id)))

# consorts
consorts <- read_delim("data/consorts.txt", " ")
# pregnancies
pregs <- read_delim("data/pregnancies.txt") %>% 
        rename(ewe_id = MumID) %>% 
        mutate(birth_date = mdy(paste(BirthMonth, BirthDay, BirthYear, "|"))) %>% 
        mutate(year_conc = BirthYear - 1) %>% 
        select(ewe_id, birth_date, year_conc)

# -1 means mounting observed, 0 means mounting not observed
cons <- consorts %>% 
        as_tibble() %>% 
        clean_names() %>% 
        separate(col = date, into = c("date", "out"), sep = " ") %>% 
        separate(col = time, into = c("out2", "time"), sep = " ") %>% 
        select(-out, -out2) %>% 
        mutate(date = ymd(date),
               time = hms(time)) %>% 
        mutate(month = month(date),
               year = year(date))
table(cons$mounts)

# rational failed matings:
# 1) females which were seen with only one male and didn't give birth
# the next year
# this is a little bit boring as it doesn't bring us closer to the mechanism
# pot_failed_matings <- cons %>% 
#         group_by(year, ewe_id) %>% 
#         # this finds ewes mating with only one male
#         filter(n_distinct(tup_id) == 1) %>% 
#         arrange(ewe_id) %>% 
#         # who didn't have offspring next year
#         left_join(fitness, by = c("ewe_id" ="id", "year" = "year_mating")) %>% 
#         filter(!is.na(tup_id)) %>% 
#         filter(offspring_born == 0) %>% 
#         slice(1) %>% 
#         #filter(!is.na(obs)) %>% 
#         #filter(mounts == -1) %>% 
#         select(date, tup_id, ewe_id, obs, age) %>%
#         rename(age_ewe = age) %>% 
#         mutate(failed = 1)


# first principles for this analysis
# chose all females with a "gap"
second_oestrus <- cons %>% 
        group_by(year, ewe_id) %>% 
        arrange(year, ewe_id, date) %>% 
        mutate(time_lag = date - lag(date)) %>% 
        #filter(any(time_lag >= 10)) %>% 
        mutate(time_lag = ifelse(is.na(time_lag), 0, time_lag)) %>% 
        ungroup() %>% 
        filter(time_lag >= 10) %>% 
        group_by(year, ewe_id) %>% 
        left_join(fitness %>% select(id, year_mating, age, offspring_born), by = c("ewe_id" ="id", "year" = "year_mating")) %>% 
        select(date, tup_id, ewe_id, obs, age, offspring_born) %>% 
        #filter(offspring_born > 0) %>% 
        rename(age_ewe = age) %>% 
        mutate(failed = 1)
    
# take consorts where a female only consorted once (though it might be over
# multiple days) with a single male
one_oestrus <- cons %>% 
        group_by(year, ewe_id) %>% 
        arrange(year, ewe_id, date) %>% 
        mutate(time_lag = date - lag(date)) %>% 
        #filter(any(time_lag >= 10)) %>% 
        mutate(time_lag = ifelse(is.na(time_lag), 0, time_lag)) %>% 
        # filter ewes which didn't have two oestrusses
        filter(all(time_lag < 10)) %>% 
        # just take one row of each ewe within a year
        slice(1) %>% 
        left_join(fitness, by = c("ewe_id" ="id", "year" = "year_mating")) %>% 
        # only take ewes which had offspring the next year
        filter(offspring_born > 0) %>% 
        select(date, tup_id, ewe_id, obs, age) %>% 
        rename(age_ewe = age) %>% 
        mutate(failed = 0)


# combine        
all_matings <- bind_rows(second_oestrus, one_oestrus) %>% 
                left_join(fitness %>% select(id, froh_all, year_mating), 
                          by = c("ewe_id"="id", "year"="year_mating")) %>% 
                rename(froh_ewe = froh_all) %>% 
                left_join(fitness %>% select(id, froh_all, year_mating, age), 
                          by = c("tup_id"="id", "year"="year_mating")) %>% 
                rename(froh_tup = froh_all, age_tup = age) %>% 
                filter(!is.na(year) & !is.na(froh_ewe) & !is.na(ewe_id) & !is.na(age_ewe))


# run models
# time saver function for modeling
nlopt <- function(par, fn, lower, upper, control) {
        .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                                  opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                              maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
        list(par = res$solution,
             fval = res$objective,
             conv = if (res$status > 0) 0 else res$status,
             message = res$message
        )
}


# ewes <- sort(table( carrier_matings_sub$ewe_id))
# keep_ewes <- names(ewes)[ewes < 6]
# carrier_matings_sub <- carrier_matings_sub %>% filter(ewe_id %in% keep_ewes)
all_matings <- all_matings %>%
                ungroup() %>% 
                mutate(age_ewe_std = as.numeric(scale(age_ewe)),
                       age_ewe_std2 = age_ewe_std^2)

fit <- glmer(failed ~ froh_ewe + froh_tup + age_ewe_std + age_ewe_std2 + (1|year) + (1|ewe_id) + (1|tup_id), #+ (1|obs2),
             family = binomial, data = all_matings,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

library(performance)
check_model(fit)
binned_residuals(fit)
tidy(fit, conf.int = TRUE)
p1 <- plot_model(fit, type = "pred", terms = "froh_ewe [all]", show.data = TRUE) +
   # geom_point(alpha = 0.4) +
    theme_simple(grid_lines = FALSE) +
    xlab(expression(F[ROH])) +
    ylab("Probability of\nfailed conception") +
    theme(axis.line.x = element_line(size = 0.5),
          axis.ticks.y = element_line(size = 0.5)) +
    ggtitle("")
p1
ggsave("figs/failed_conception.jpg", width=3.5, height=3)

# inbreeding depression in breeding success?
abs <- fitness %>% 
        filter(sex == "F") %>% 
        mutate(age_std = as.numeric(scale(age)),
           age_std2 = age_std^2) %>% 
        mutate(offspring_born = ifelse(offspring_born > 0, 1, 0)) %>% 
        # remove individuals which came into second oestrus from this
        # analysis
        filter(!(id %in% second_oestrus))

fit <- glmer(offspring_born ~ froh_all + age_std + age_std2 + (1|year_mating) + (1|id), #+ (1|obs2),
             family = binomial, data = abs,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(fit)
tidy(fit, conf.int = TRUE)

p2 <- plot_model(fit, type = "pred", terms = "froh_all [all]", show.data = TRUE) +
    theme_simple(grid_lines = FALSE) +
    xlab(expression(F[ROH])) +
    theme(axis.line.x = element_line(size = 0.5),
          axis.ticks.y = element_line(size = 0.5)) +
    ylab("Predicted\nannual breeding success") +
    ggtitle("")
p2
ggsave("figs/abs.jpg", p2, width=3.5, height=3)

library(patchwork)
p <- p1 + p2 + plot_annotation(tag_levels = 'A')

ggsave("figs/froh_breeding.jpg", width = 6, height = 2.8)

# another approach: duration from first mating to offspring

cons %>% 
    left_join(fitness, by )



