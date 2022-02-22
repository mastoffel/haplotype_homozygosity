# model haplotype effects on postnatal fitness
library(tidyverse)
library(lme4)
library(performance)
library(sjPlot)
library(brms)
library(here)

# 
haps_fit <- read_delim(here("output", "haps_and_fitness.txt"))

mod_df <- haps_fit %>% 
        filter(#region == location,
                #sex == s,
                age == 0) %>% 
        mutate(gt = as.factor(gt),
               froh_std = as.numeric(scale(froh_all)),
               weight_std = as.numeric(scale(weight)),
               mum_age_std = as.numeric(scale(mum_age)),
               mum_age_std2 = mum_age_std^2) %>% 
        select(survival, region, gt, sex, froh_std, twin, weight_std, mum_age_std,
               birth_year, mum_id) %>% 
        pivot_wider(names_from = region, values_from = gt) %>% 
        # mutate(gt9 = ifelse(chr9_6571 == 2, 1, 0),
        #        gt18 = ifelse(chr18_267 == 2, 1, 0),
        #        gt5 = ifelse(chr5_6193 == 2, 1, 0),
        #        gt7 = ifelse(chr7_12119 == 2, 1, 0))
        mutate(gt9 = as.factor(chr9_6571),
               gt18 = as.factor(chr18_267),
               gt5 = as.factor(chr5_6193),
               gt7 = as.factor(chr7_12119))

# lme4
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

fit_glmer <- glmer(survival ~ gt9 + gt18 + gt5 + gt7 + sex + froh_std + twin + weight_std + mum_age_std + (1|birth_year) + (1|mum_id), #gt + sex + weight_std +  twin + froh_std + (1|birth_year) + (1|mum_id)
             data = dat, family = binomial(link = "logit"),
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(fit_glmer, conf.int = TRUE)


# brms
fit <- brm(survival ~ gt9 + gt18 + gt5 + gt7 + sex + froh_std + twin + weight_std + mum_age_std + (1|birth_year) + (1|mum_id),
           data = dat, family = bernoulli(), 
           iter = 10000,
           set_prior("normal(0,5)", class = "b"))
#saveRDS(fit, "output/haps_fitness_mod.RDS")
fit <- readRDS("output/haps_fitness_mod.RDS")
prior_summary(fit)
summary(fit)
plot(fit)
loo(fit)
pp_check(fit, nsamples = 100)



# check both haplotypes on chr 9 || not going anywhere here I think
dat <- haps_all %>% 
        filter(region == locations[1],
               #sex == s,
               age == 0) %>% 
        mutate(gt_hap1 = case_when(
                gt == 1 & gt_hap2 != 1 ~ 1,
                #gt == 2 & gt_hap2 == 2 ~ 0,
                gt == 2 & gt_hap2 == 1 ~ 1,
                gt == 2 & gt_hap2 == 0 ~ 2,
                TRUE ~ 0
        )) %>%
        mutate(gt_hap1 = case_when(
                gt_hap2 == 2 ~ NaN,
                gt_hap2 == 1 ~ NaN,
                TRUE ~ gt_hap1
        )) %>% 
        mutate(gt_hap2 = case_when(
                gt_hap1 == 2 ~ NaN,
                gt_hap1 == 1 ~ NaN,
                TRUE ~ gt_hap2
        )) %>% 
        select(region, id, gt, gt_hap1, gt_hap2, everything()) %>% 
        mutate(gt = as.factor(gt),
               gt_hap1 = as.factor(gt_hap1),
               gt_hap2 = as.factor(gt_hap2),
               froh_std = as.numeric(scale(froh_all)),
               weight_std = as.numeric(scale(weight)),
               mum_age_std = as.numeric(scale(mum_age)),
               mum_age_std2 = mum_age_std^2)
#  gt_alt = ifelse(gt == 0 | gt ==1, 0, 1))
fit <- glmer(survival ~ gt_hap2 + sex + froh_std + twin + weight_std + mum_age_std + (1|birth_year) + (1|mum_id), #gt + sex + weight_std +  twin + froh_std + (1|birth_year) + (1|mum_id)
             data = dat, family = binomial(link = "logit"),
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

tidy(fit, conf.int=TRUE)





# haplotype effects on annual survival?
survival_mod <- function(location, haps_all) {
        
        dat <- haps_all %>% 
                ungroup() %>% 
                filter(region == location) %>% 
                mutate(gt = as.factor(gt),
                       #sex = relevel(sex, "M"),
                       age_std = as.numeric(scale(age)),
                       age_std2 = age_std^2,
                       froh_std = scale(froh_all),
                       weight_std = as.numeric(scale(weight)),
                       lamb = ifelse(age == 0, 1, 0),
                       lamb = as.factor(lamb),
                       life_stage = case_when(
                               age == 0 ~ "lamb",
                               age > 0 & age <= 2 ~ "early_life",
                               age > 2 & age <= 4 ~ "mid_life",
                               age > 4 ~ "late_life",
                       )) 
        
        fit <- glmer(survival ~ gt + twin + froh_std + life_stage * sex + weight_std + (1|sheep_year) + (1|birth_year) + (1|id),
                     data = dat, family = binomial(link = "logit"), 
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        fit
        # fit %>% 
        #   tidy(conf.int = TRUE) %>% 
        #   mutate_if(is.numeric, round, 5)
}

locations <- unique(haps_all$region)
survival <- map(locations, survival_mod, haps_all)
names(survival) <- locations
map(survival, tidy, conf.int=TRUE)
survival

saveRDS(survival, file = "output/lifetimesurv_mods_hap500.rds")

#survival_f <- map(locations, survival_mod, haps_all, "F")
#survival_m <- map(locations, survival_mod, haps_all, "M")
library(performance)
library(sjPlot)
map(survival, tidy, conf.int=TRUE)
plot_model(survival[[3]], type = "pred", terms = "gt")
binned_residuals(survival[[2]])


library(brms)
fit_brm <- brm(survival ~ gt + sex + twin + froh_std + life_stage * sex + (1|sheep_year) + (1|id), # (1|birth_year) +
               data = dat, family = bernoulli(link = "logit"))




# reproduction
lrs_mod <- function(location, haps_all, s) {
        
        dat <- haps_all %>% 
                filter(region == location,
                       sex == s) %>% 
                mutate(gt = as.factor(gt),
                       olre = 1:nrow(.),
                       froh_scaled = scale(froh_all),
                       age_scaled = scale(age),
                       age_scaled_2 = age_scaled ^ 2)
        fit <- glmer(offspring_born ~ gt + froh_scaled + twin + age_scaled + age_scaled_2 + (1|sheep_year) + (1|birth_year) + (1|id) + (1|olre),
                     data = dat, family = poisson,
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        fit %>% 
                tidy(conf.int = TRUE) %>% 
                mutate_if(is.numeric, round, 5)
}

locations <- unique(haps_all$region)
lrs_f <- map(locations, lrs_mod, haps_all, "F")
lrs_m <- map(locations, lrs_mod, haps_all, "M")

# check homozygous sex ratio
haps_all %>% 
        filter(gt == 2, age == 0) %>% 
        group_by(region, sex) %>% 
        tally()


# traits: first august weight
weight_mod <- function(location, haps_all, s) {
        
        dat <- haps_all %>% 
                filter(region == location,
                       sex == s,
                       age == 0) %>% 
                mutate(gt = as.factor(gt),
                       olre = 1:nrow(.))
        fit <- lmer(weight ~ gt + froh_all + twin + (1|sheep_year) + (1|birth_year) + (1|mum_id),
                    data = dat)
        fit %>% 
                tidy(conf.int = TRUE) %>% 
                mutate_if(is.numeric, round, 5)
}

locations <- unique(haps_all$region)
weight_f <- map(locations, weight_mod, haps_all, "F")
weight_m <- map(locations, weight_mod, haps_all, "M")


# traits: first hindleg
hindleg_mod <- function(location, haps_all, s) {
        
        dat <- haps_all %>% 
                filter(region == location,
                       sex == s,
                       age == 0) %>% 
                mutate(gt = as.factor(gt),
                       olre = 1:nrow(.))
        fit <- lmer(hindleg ~ gt + froh_all + twin + (1|sheep_year) + (1|birth_year) + (1|mum_id),
                    data = dat)
        fit %>% 
                tidy(conf.int = TRUE) %>% 
                mutate_if(is.numeric, round, 5)
}

locations <- unique(haps_all$region)
hindleg_f <- map(locations, hindleg_mod, haps_all, "F")
hindleg_m <- map(locations, hindleg_mod, haps_all, "M")







# lifespan
haps_all %>% 
        group_by(region, id) %>% 
        summarise(lifespan = max(lifespan),
                  gt = first(gt),
                  sex = first(sex))


haps_all %>% 
        group_by(region, id) %>% 
        summarise(lifespan = max(lifespan),
                  gt = first(gt),
                  sex = first(sex)) %>%
        #filter(sex == "M") %>% 
        ggplot(aes(gt, lifespan, color = sex)) +
        facet_wrap(~region) +
        geom_boxplot(outlier.color=NA) +
        geom_point(pch = 21, size = 0.5, alpha = 0.3,
                   position = position_jitterdodge()) +
        scale_y_log10() +
        scale_fill_viridis_d() +
        theme_simple()

pal <- wes_palette("Darjeeling2", 2)

# first year survival
haps_all %>% 
        filter(age == 0) %>% 
        group_by(region, sex, gt) %>% 
        summarise(survival = mean(survival, na.rm = TRUE)) %>% 
        ggplot(aes(gt, survival, fill = sex)) +
        facet_wrap(~region) +
        geom_bar(position="dodge", stat="identity") +
        theme_minimal() +
        scale_fill_manual(values = pal) 
# scale_y_log10()

haps_all %>% 
        group_by(region, id) %>% 
        summarise(LRS = sum(offspring_survived, na.rm = TRUE),
                  sex = first(sex),
                  gt = first(gt)) %>%  
        drop_na() %>% 
        ggplot(aes(gt, LRS, fill = sex)) +
        facet_wrap(~region) +
        geom_boxplot(outlier.color=NA) +
        # geom_point(pch = 21, size = 0.5, alpha = 0.3, 
        #           position = position_jitterdodge()) +
        #scale_y_log10() +
        
        scale_fill_viridis_d()

mod_dat <- haps_all %>% 
        
        
        
        options(scipen=999)
library(broom.mixed)
library(lme4)
mod_dat <- haps_all %>% 
        filter(region == "chr22_788") %>% 
        mutate(gt = as.factor(gt)) %>% 
        filter(sex == "M")

mod <- glmer(survival ~ gt + twin + froh_all + (1|sheep_year) + (1|id),
             data = mod_dat, family = binomial)

tidy(mod, conf.int = TRUE)
