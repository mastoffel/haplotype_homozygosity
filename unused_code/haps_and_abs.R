library(tidyverse)
library(here)
library(glue)

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

# load haplotypes
haps <- read_delim(here("output", "sheep_top_haps_500.txt"), " ") %>% 
        select(region, id, gt)

abs <- fitness %>% 
        left_join(haps) %>% 
        filter(sex == "F",
               region == "chr9_6571") %>% 
        mutate(offspring_born = ifelse(offspring_born > 0, 1, 0),
               gt = ifelse(gt > 0, 1, 0)) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age_std2 = age_std^2)


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

options(pillar.sigfig = 2)

fit <- glmer(offspring_born ~ gt + froh_all + age_std + age_std2 + (1|year_mating) + (1|id), #+ (1|obs2),
             family = binomial, data = abs,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(fit)
