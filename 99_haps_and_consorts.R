library(dplyr)
library(dbplyr)
library(lubridate)
library(tidyverse)
library(janitor)
library(here)
library(lme4)
library(broom.mixed)
library(sjPlot)
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


# first approach: take females which mated a second time after two week
# assuming that the pregnancy failed. Take only
pot_failed_matings <- cons %>% 
        group_by(year, ewe_id) %>% 
        #filter(year > 2000) %>% 
        #filter(mounts == -1) %>% 
        # this finds ewes mating with several males
        filter(n() > 1) %>% 
        arrange(desc(ewe_id)) %>% 
        # ewe being consorted by at least two males
        #filter(n_distinct(tup_id) > 1) %>% 
        arrange(year, ewe_id, date) %>% 
        mutate(time_lag = date - lag(date)) %>% 
        # filter ewes where any two consorts are at least 10 days apart
        filter(any(time_lag >= 10)) %>% 
        # first mating gets NA, so rather give it 0
        mutate(time_lag = ifelse(is.na(time_lag), 0, time_lag)) %>% 
        # filter matings before the gap (potentially failed pregnancies)
        slice(1:((match(TRUE, time_lag == max(time_lag)))-1)) %>% 
        #slice(1:((match(TRUE, time_lag >= 7)-1))) %>%
        # 
        # filter all which were mated/seen mated by only only one tup in previous estrous
        # remove unknown tups
        filter(n_distinct(tup_id) == 1) %>% 
        filter(!is.na(tup_id)) %>% 
        left_join(fitness %>% select(id, year_mating, age), by = c("ewe_id" ="id", "year" = "year_mating")) %>% 
        #filter(year > 1990) %>% 
        # mounting observed
        #filter(mounts == -1) %>% 
        select(date, tup_id, ewe_id, obs, age) %>% 
        rename(age_ewe = age) %>% 
        mutate(failed = 1)

# take consorts where a female only consorted once (though it might be over
# multiple days) with a single male
pot_succ_matings <- cons %>% 
        mutate(month = month(date),
               year = year(date)) %>% 
        # this finds ewes mating with several males
        group_by(year, ewe_id) %>% 
        filter(length(unique(tup_id))==1) %>% 
        # mounting observed
        #filter(mounts == -1) %>% 
        #filter(year > 2000) %>% 
        arrange(year,ewe_id) %>%
        select(date, tup_id, ewe_id, obs) %>% 
        filter(!is.na(tup_id)) %>% 
        # if ewe mated with same tup on multiple days in a year, take only
        # one
        slice(1) %>% 
        left_join(fitness, by = c("ewe_id" ="id", "year" = "year_mating")) %>% 
        filter(offspring_born == 1) %>% 
        #mutate(preg_time = birth_date - date) %>% 
        # filter where pregnancy time is odd or no birth
        #filter(!is.na(preg_time)) %>% 
        select(date, tup_id, ewe_id, obs, age) %>% 
        rename(age_ewe = age) %>% 
        mutate(failed = 0)

# combine        
all_matings <- bind_rows(pot_failed_matings, pot_succ_matings) %>% 
                # combine rare observers into other group
                group_by(obs) %>% 
                mutate(obs2 = ifelse(n() < 5, "other", obs)) %>% 
                ungroup()

# load haplotypes
haps <- read_delim(here("output", "sheep_top_haps_500.txt"), " ") %>% 
                select(region, id, gt)

# test:
carrier_matings <- all_matings %>% 
        left_join(haps, by = c("tup_id" = "id")) %>% 
        left_join(haps, by = c("ewe_id" = "id", "region" = "region")) %>% 
        filter(!is.na(gt.x) & !is.na(gt.y)) %>% 
        mutate(mating_type = case_when(
                gt.x > 0 & gt.y > 0 ~ "c_c",
                gt.x == 0 & gt.y > 0 ~ "c_nonc",
                gt.x > 0 & gt.y == 0 ~ "c_nonc",
                gt.x == 0 & gt.y == 0 ~ "nonc_nonc"),
               mating_type = factor(mating_type, 
                                    levels = c("nonc_nonc", "c_nonc", "c_c"))) %>% 
        mutate(mating_type2 = case_when(
                gt.x > 0 & gt.y > 0 ~ "c_c",
                TRUE ~ "aa"
        )) %>% 
        left_join(fitness, by = c("ewe_id" = "id", "year" = "year_mating")) %>% 
        rename(froh_ewe = froh_all) %>% 
        left_join(fitness[c("id", "froh_all", "year_mating")], by = c("tup_id" = "id", "year" = "year_mating")) %>% 
        rename(froh_tup = froh_all)


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

options(pillar.sigfig = 4)

ins_effects <- function(hap, carrier_matings) {
        carrier_matings_sub <- carrier_matings %>% 
                                        filter(region == hap) %>% 
                                        mutate(froh_ewe = as.numeric(scale(froh_ewe)),
                                               froh_tup = as.numeric(scale(froh_tup))) %>% 
                                        group_by(year, ewe_id) %>%
                                        slice(1) %>%
                                        ungroup()
        
        # ewes <- sort(table( carrier_matings_sub$ewe_id))
        # keep_ewes <- names(ewes)[ewes < 6]
        # carrier_matings_sub <- carrier_matings_sub %>% filter(ewe_id %in% keep_ewes)
        
        fit <- glmer(failed ~ mating_type2 + froh_ewe + age_ewe + (1|year) + (1|ewe_id), #+ (1|obs2),
                     family = binomial, data = carrier_matings_sub,
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        fit
}

all_regions <- map(unique(carrier_matings$region), ins_effects, carrier_matings)
names(all_regions) <- unique(carrier_matings$region)
#saveRDS(all_regions, file = "output/ins_effects_500.rds")

out <- map(all_regions, tidy, conf.int = TRUE)

plot_model(all_regions[[4]], type = "pred", terms = c("mating_type"))

all_regions <- map_dfr(unique(carrier_matings$region), 
                       ins_effects, carrier_matings, .id = "region")        

print(all_regions, n = 40)

        

library(brms)
fit <- brm(failed ~ mating_type + scale(froh_all) + (1|year) + (1|ewe_id), 
           data = carrier_matings_sub, family = bernoulli)

# second approach: which matings did not lead to offspring the next year?
matings_offsp <- cons %>% 
        mutate(month = month(date),
               year = as.numeric(year(date))) %>% 
        # this finds ewes mating with several males
        group_by(year, ewe_id) %>% 
        #filter(year > 2000) %>% 
       # filter(n() > 1) %>% 
        # mounting observed
        #filter(mounts == -1) %>% 
        arrange(desc(ewe_id)) %>% 
        # ewe being consorted by only one male
        filter(n_distinct(tup_id) == 1) %>% 
       # filter(mounts == -1) %>%
        arrange(year, ewe_id, date) %>% 
        filter(!is.na(tup_id)) %>% 
        # filter all which were mated/seen mated only once in the previous estrous
        slice(1) %>% 
       # rename(id = ewe_id, year_mating  = year) %>% 
        #ungroup() %>% 
        inner_join(fitness, by = c("ewe_id" = "id", "year" = "year_mating")) %>% 
        rename(year_mating = year) %>% 
        #filter(year_mating > 1992) %>% 
       # filter(offspring_born == 0) %>% 
        select(date, tup_id, ewe_id, offspring_born, froh_all) %>% 
        left_join(haps, by = c("tup_id" = "id")) %>% 
        left_join(haps, by = c("ewe_id" = "id", "region" = "region")) %>% 
        filter(!is.na(gt.x) & !is.na(gt.y)) %>% 
        mutate(mating_type = case_when(
                gt.x > 0 & gt.y > 0 ~ "c_c",
                gt.x == 0 & gt.y > 0 ~ "c_nonc",
                gt.x > 0 & gt.y == 0 ~ "c_nonc",
                gt.x == 0 & gt.y == 0 ~ "nonc_nonc"),
               mating_type = factor(mating_type, 
                                    levels = c("nonc_nonc", "c_nonc", "c_c"))) %>% 
        mutate(mating_type2 = case_when(
                gt.x > 0 & gt.y > 0 ~ "c_c",
                TRUE ~ "aa"
        )) %>% 
        mutate(rs = ifelse(offspring_born > 0, 1, 0)) %>% 
        ungroup()


offspring_mod <- function(hap, matings_offsp) {
        matings_offsp_sub <- matings_offsp %>% 
                                filter(region == hap) %>% 
                                mutate(froh_std = as.numeric(scale(froh_all))) %>% 
                                mutate(year_mating = as.numeric(year_mating))
        fit <- glmer(rs ~ mating_type + froh_std + (1|year_mating) + (1|ewe_id) + (1|tup_id),
                     family = binomial, data = matings_offsp_sub)
        fit
}

all_regions_rs <- map(unique(carrier_matings$region), 
                      offspring_mod, matings_offsp)
map(all_regions_rs, tidy, conf.int = TRUE)

all_regions <- map_dfr(unique(carrier_matings$region), 
                       offspring_mod, matings_offsp, .id = "region")        


print(all_regions, n = 24)
