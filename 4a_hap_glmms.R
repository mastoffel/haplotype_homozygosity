# model haplotype effects on postnatal fitness
library(tidyverse)
library(lme4)
library(performance)
library(sjPlot)
library(brms)
#library(modelbased)
library(here)
source("theme_simple.R")
library(broom.mixed)

haps_fit <- read_delim(here("output", "haps400_and_fitness.txt"))

mod_df <- haps_fit %>% 
        filter(#region == location,
                #sex == s,
                age == 0) %>% 
        select(survival, region, gt, sex, froh_all, twin, weight, hindleg, 
               birth_year, mum_id, mum_age) %>% 
        # drop all rows with missing values
        #drop_na() %>% 
        mutate(gt = as.factor(gt),
               froh_std = as.numeric(scale(froh_all)),
               weight_std = as.numeric(scale(weight)),
               hindleg_std = as.numeric(scale(hindleg)),
               mum_age_std = as.numeric(scale(mum_age))) %>% 
        select(survival, region, gt, sex, froh_std, twin, weight_std,  weight, hindleg_std, 
               birth_year, mum_id, mum_age_std, hindleg) %>% 
        pivot_wider(names_from = region, values_from = gt) %>% 
        mutate(gt18 = as.factor(chr18_267),
                gt5 = as.factor(chr5_6293),
                gt7 = as.factor(chr7_12196))

mod_df$hindleg

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

fit_glmer <- glmer(survival ~ gt5 + gt18 + gt7 + hindleg_std + froh_std + sex + twin + (1|birth_year) + (1|mum_id), 
                   data = mod_df, family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
out <- tidy(fit_glmer, conf.int = TRUE)
out
binned_residuals(fit_glmer)

fit_glmer <- lmer(weight ~ gt5 + gt18 + gt7 + froh_std + sex + twin + (1|birth_year) + (1|mum_id), 
                   data = mod_df)
out <- tidy(fit_glmer, conf.int = TRUE)
summary(fit_glmer)

fit_lmer <- lmer(weight_std ~ gt5 + gt18 + gt7 + sex + froh_std + hindleg_std + twin + mum_age_std +  (1|birth_year) + (1|mum_id), #gt + sex + weight_std +  twin + froh_std + (1|birth_year) + (1|mum_id)
                   data = mod_df)
tidy(fit_lmer, conf.int = TRUE)

fit_glmer <- glmer(survival ~ gt18 + gt5 + gt7 + sex + froh_std + twin + (1|birth_year) + (1|mum_id), #gt + sex + weight_std +  twin + froh_std + (1|birth_year) + (1|mum_id)
             data = mod_df, family = binomial(link = "probit"))
tidy(fit_glmer, conf.int = TRUE)
binned_residuals(fit_glmer)

estimate_contrasts(fit_glmer, transform = "response", contrast = "gt9")

# brms #  + gt18 + gt5 + gt7 +
fit <- brm(survival ~  gt18 + gt5 + gt7 + sex + froh_std + twin + hindleg_std + (1|birth_year) + (1|mum_id),
           data = mod_df, family = bernoulli(), 
           iter = 10000, thin = 1,
           set_prior("normal(0,5)", class = "b"))
saveRDS(fit, "output/haps_fitness_mod.RDS")
fit <- readRDS("output/haps_fitness_mod.RDS")
prior_summary(fit)
summary(fit)
plot(fit)
loo(fit)
pp_check(fit, nsamples = 100)
plot_model(fit)
binned_residuals(fit, term = "sex")

# get marginal differences on probability scale 
# posterior estimates for all fixed effects
post <- as.matrix(posterior_samples(fit)[c(1:8, 10)])
dim(post)
preds <- rbind(diag(1, nrow = 6, ncol = 6), rep(0, 6)) %>% 
                as.matrix() %>% 
                as_tibble() %>% 
                setNames(colnames(post)[2:7]) %>% 
        add_column(intercept = 1, .before = 1) %>% 
        # average marginal effect
        mutate(sex = mean(ifelse(mod_df$sex == "F", 0, 1))) %>% 
        # standardised measure have mean of 0, so not necessary to include here
        mutate(twin = mean(mod_df$twin))
                #froh = mean(mod_df$froh_std),
               #weight = mean(mod_df$weight_std))

# get predictions
out <- map(1:3, function(i) as_tibble( (t(post) * preds[i, ])))

# vals = vector 
# post = matrix (post. probs)
get_surv <- function(vals, post){
        effs <- sweep(post, MARGIN = 2, vals, "*")
        logodds <- rowSums(effs)
        # inv logit
        probs <- plogis(logodds)
}

diff_probs <- apply(preds, 1, get_surv, post) %>% 
        as_tibble() %>% 
        setNames(c("gt181", "gt182", "gt51", "gt52", 
                   "gt71", "gt72", "nogt")) %>% 
        map_df(function(x) (x - .$nogt)*100) %>% 
        select(-nogt) %>% 
        pivot_longer(everything(), names_to = "predictor", values_to = "estimate") %>% 
        mutate(hap = str_sub(predictor, end = -2))

write_delim(diff_probs, here("output", "survival_marginal_effects.txt"))


# pleiotropy? Fit model with weigh
# brms #  + gt18 + gt5 + gt7 +
fit <- brm(weight ~  gt18 + gt5 + gt7 + sex + hindleg_std + froh_std + twin + (1|birth_year) + (1|mum_id),
           data = mod_df, family = "gaussian",
           iter = 10000, thin = 1,
           set_prior("normal(0,5)", class = "b"))
saveRDS(fit, "output/haps_weight_mod.RDS")
fit <- readRDS("output/haps_weight_mod.RDS")
tidy(fit)
prior_summary(fit)
summary(fit)
plot(fit)
loo(fit)
pp_check(fit, nsamples = 100)
plot_model(fit)


library(tidybayes)
diff_probs %>% 
        ggplot(aes(x = estimate, y = predictor)) +
        stat_halfeye() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = c("gray80", "skyblue")) +
        facet_wrap(~hap ,scales = "free_y", nrow = 1) +
        scale_x_continuous(breaks = seq(-30, 30, 10)) +
        theme(axis.text.y = element_blank()) +
        theme_simple()

hist(out$gt92 - out$nogt)
           
hist(invlogit(rowSums(post[, c(1, 2, 12)])) - invlogit(rowSums(post[, c(1, 12)])))

emm2 <- lsmeans(fit, pairwise ~ gt9)
summary(emm2, type = "response")

emm3 <- lsmeans(fit, ~gt9)
pairs(regrid(emm3))

estimate_contrasts()
# design matrix
mod_mat <- model.matrix(~ gt9 + gt18 + gt5 + gt7 + sex + froh_std + twin + weight_std, data = mod_df)

dat <- standata(fit)


library(tidybayes)
library(stringr)
post %>% 
        select(b_gt91:b_gt72) %>% 
        pivot_longer(everything(), names_to = "predictor", values_to = "estimate") %>% 
        mutate(hap = str_sub(predictor, end = -2)) %>% 
        ggplot(aes(x = plogis(estimate), y = predictor, fill = stat(x < 1))) +
        stat_halfeye() +
        geom_vline(xintercept = 1, linetype = "dashed") +
        scale_fill_manual(values = c("gray80", "skyblue")) +
        facet_grid(hap~. ,scales = "free_y") +
        scale_x_continuous(limits = c(0, 2)) +
        theme_simple()

invlogit <- function(x) exp(x)/(1+exp(x))
summary(fit)
dim(post)
post %>% 
        select(b_gt91:b_gt72) %>% 
        pivot_longer(everything(), names_to = "predictor", values_to = "estimate") %>% 
        mutate(hap = str_sub(predictor, end = -2)) %>% 
        ggplot(aes(x = invlogit(estimate), y = predictor)) +
        stat_halfeye() +
        geom_vline(xintercept = 1, linetype = "dashed") +
        scale_fill_manual(values = c("gray80", "skyblue")) +
        facet_grid(hap~. ,scales = "free_y") +
        scale_x_continuous(limits = c(0, 2)) +
        theme_simple()

marginal_effects(fit, effects = "gt9")


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
