library(tidyverse)
library(brms)
library(here)
library(nadiv)
library(asreml)
library(broom.mixed)
library(MCMCglmm)
load(here("data", "sheep_ped.RData"))

# add rel mat
Amat <- as.matrix(nadiv::makeA(sheep_ped))
ainv <- ainverse(sheep_ped)
#Ainv <- inverseA(sheep_ped)$Ainv

# fitness
haps_fit <- read_delim(here("output", "haps400_and_fitness.txt"))

mod_df <- haps_fit %>% 
        filter(#region == location,
                #sex == s,
                age == 0) %>% 
        select(id, survival, region, gt, sex, froh_all, twin, weight, weight, hindleg, 
               birth_year, mum_id) %>% 
        drop_na() %>% 
        mutate(gt = as.factor(gt),
               froh_std = as.numeric(scale(froh_all)),
               weight_std = as.numeric(scale(weight)),
               hindleg_std = as.numeric(scale(hindleg))) %>% 
        select(id, survival, region, gt, sex, froh_std, twin, weight_std, weight, hindleg_std, 
               birth_year, mum_id) %>% 
        pivot_wider(names_from = region, values_from = gt) %>% 
        mutate(gt18 = as.factor(chr18_267),
               gt5 = as.factor(chr5_6293),
               gt7 = as.factor(chr7_12196)) %>% 
        mutate(id = as.factor(id),
               sex = as.factor(sex),
               twin = as.factor(twin)) %>% 
        rename(animal = id) %>% 
        mutate(animal_pe = animal) %>% 
        as.data.frame()


fit <- brm(
        weight ~ 1 + gt18 + gt5 + gt7 + sex + hindleg_std + froh_std + twin + (1 | gr(animal, cov = Amat)) + (1 | birth_year) + (1 | mum_id),
        data = mod_df,
        family = gaussian(), data2 = list(Amat = Amat),
        chains = 2, cores = 2, iter = 1000
)

summary(fit)
tidy(fit)


fit <- asreml(
        fixed =  weight ~ 1 + gt18 + gt5 + gt7 + sex + hindleg_std + froh_std + twin,
        random = ~  vm(animal, ainv) + birth_year + mum_id, # 
        residual = ~idv(units),
        data = mod_df,
        na.action = na.method(x = "omit", y = "omit")
)
summary(fit)
summary(fit, coef = TRUE)$coef.fixed
wald.asreml(fit, ssType = "conditional", denDF = "numeric")
summary(fit)$varcomp
plot(fit)


### MCMCglmm ####
prior1.1 <- list(
        G = list(G1 = list(V = 1, nu = 0.002),
                 G2 = list(V = 1, nu = 0.002),
                 G3 = list(V = 1, nu = 0.002)),
        R = list(V = 1, nu = 0.002)
)

fit <- MCMCglmm(weight ~ gt18 + gt5 + gt7 + sex + hindleg_std + froh_std + twin,
                random = ~animal + birth_year + mum_id, ginv = list(animal = Ainv),
                data = mod_df, prior = prior1.1,
                nitt = 65000, thin = 50, burnin = 15000
)
