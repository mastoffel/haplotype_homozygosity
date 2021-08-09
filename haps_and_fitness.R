library(dplyr)
library(readr)
library(here)
library(glue)
library(tidyr)
library(purrr)
library(tidyverse)
library(stringr)
library(data.table)
library(broom.mixed)
library(wesanderson)
source("theme_simple.R")
library(lme4)
# fitness
load("data/fitness_roh.RData") 
fitness <- fitness_data %>% 
                select(id, survival, sheep_year, age, birth_year, sex, mum_id,
                       twin, offspring_born, offspring_survived, weight,
                       hindleg, horn, father, froh_all, birthwt) %>% 
                group_by(id) %>% 
                mutate(lifespan = max(age))

# snp map
snp_map <- fread(here("data", "plink", "sheep.bim")) %>% 
  setNames(c("chr", "snp", "cM", "pos", "a", "b")) %>% 
  select(chr, snp, pos) %>% 
  group_by(chr) %>% 
  mutate(number = 1,
         snp_num = cumsum(number)) %>% 
  select(-number)

# target haplotypes:
# read results from haplotype homozygosity scan
all_files <- list.files(here("output", "hap_results_imputed", "hap_len_500"), full.names = TRUE)
results <- map(all_files, read_delim, delim = "\t") %>%  #read_delim, delim = "\t"
        bind_rows()

# note: sometimes there are two significant haplotypes starting with the same 
# snp -> need to be merged eventually
top_haps <- results %>% 
        # get haplotypes which are significant because there are fewer homozyous
        # haplotypes as expected
        filter(obs < exp) %>% 
        filter(p_val < (0.05/39184)) %>% # 417373 # 39184
        arrange(chr, snp_start) %>% 
        # workflow to get only one haplotype per genome-region
        group_by(chr) %>% 
        # work with lag, the default gets around the a problem
        # when the haplotype starts at the very beginning of a chromosome
        # (snp numbers below 10)
        mutate(lag1 = lag(snp_start),
               diff = snp_start - lag1,
               diff = ifelse(is.na(diff), snp_start, diff)) %>% 
        mutate(region_start = ifelse(diff < 100, NA, snp_start)) %>% 
        fill(region_start) %>% 
        # mutate(lag1 = lag(snp_start, default = -1000),
        #        diff = snp_start - lag1) %>% 
        # # count new region when further than 10 snps away
        # mutate(diff = ifelse(diff < 10, NA, diff-1000)) %>% 
        # fill(diff) %>% 
        # give regions name: chr + start_snp
        mutate(region = paste0("chr", chr, "_", region_start)) %>% 
        group_by(region, snp_start) %>% 
        # if multiple haplotypes with same starting snp are significant, 
        # add them to "haps" column
        mutate(haps = paste0(hap, collapse = "_"),
               obs = paste0(obs, collapse = ","),
               exp = paste0(exp, collapse = ","),
               num_carriers = paste0(num_carriers, collapse = ",")) %>% 
        #mutate(haps = hap) %>% 
        # now haplotypes with lowest pvalues per regions
        group_by(region) %>% 
        arrange(p_val) %>% 
        filter(row_number()==1) 
      
  
# (1) get haplotypes for all individuals 
#hap_length <- 100
#chr <- top_haps$chr[1]
#snp_start <- top_haps$snp_start[1]
#chr <- 17
#top_hap <- "0010010010"
#snp_start <- 1798

# take top haplotype and get genotypes for hom/het/alt_hom 
hap_to_geno <- function(i, top_haps) {
        
        top_hap <- top_haps[i, ]
        
        chr <- top_hap$chr
        
        # get haplotypes, remove everything but genotypes to make it a matrix
        haps_raw <- fread(here("output", "phased_matrix", glue("sheep_phased_chr_{chr}.hap.gz"))) %>% 
                select(-V1, -V2, -V3, -V4, -V5) %>% 
                as.matrix()
        
        # get individual names, should be 5952
        # here is a bit clunky here so stay with fixed path
        ind_names <- system(glue("zgrep '^#CHROM*' output/phased/sheep_phased_chr_{chr}.vcf.gz"),
                            intern = TRUE) %>% 
                str_split("\t") %>%
                unlist() %>% 
                .[-(1:9)] %>% 
                str_split("_") %>% 
                map_chr(2)
        
        # double each name and add _a _b for haplotype 1 / haplotype 2
        ind_names_hap <- rep(ind_names, each = 2) %>% 
                paste0(., c("_a", "_b"))
        # add to matrix
        colnames(haps_raw) <- ind_names_hap
        
        snp_start <- top_hap$snp_start
        hap_length <- str_count(top_hap$hap)
        
        # combine genotypes to haplotypes
        haps <- apply(haps_raw[snp_start:(snp_start+hap_length-1), ],
                      paste, collapse = "", MARGIN = 2) %>% 
                enframe(name = "id_hap", value = "hap") %>% 
                mutate(id = str_remove(id_hap, "_[a-z]")) %>%
                select(-id_hap) %>% 
                mutate(hap_pos = rep(c("hap_a", "hap_b"), nrow(.)/2), .before = 1) %>% 
                mutate(id = as.numeric(id)) %>% 
                pivot_wider(names_from = hap_pos, values_from = hap)
        
        hap_seq <- unlist(str_split(top_hap$haps, "_"))
        # genotype for haplotype
        # in case of multiple significant haplotypes, this assumes that all
        # significant backgrounds contain the target mutation, and are 
        # therefore clustered as one haplotype
        hap_gts <- haps %>% 
                rowwise() %>% 
                mutate(gt = case_when(
                        (hap_a %in% hap_seq) & (hap_b %in% hap_seq) ~ 2,
                        (hap_a %in% hap_seq) & !(hap_b %in% hap_seq) ~ 1,
                        !(hap_a %in% hap_seq) & (hap_b %in% hap_seq) ~ 1,
                        !(hap_a %in% hap_seq) & !(hap_b %in% hap_seq) ~ 0,
                        TRUE ~ NA_real_
                )) %>% 
                ungroup() %>% 
                mutate(snp_start = snp_start, chr = chr)
}

# get haplotype genotypes for top haplotypes
haps <- map(1:nrow(top_haps), hap_to_geno, top_haps)
#saveRDS(haps, file="output/haps_600.rds")

# data.frame with all haplotypes and regions
haps_all <- haps %>% 
        setNames(top_haps$region) %>% 
        bind_rows(.id = "region") %>% 
        mutate(id = as.character(id)) %>% 
        left_join(fitness) #%>% 
        #mutate(gt = as.factor(gt))

# top haplotype genotypes per individual
haps_ind <- haps %>% 
            setNames(top_haps$region) %>% 
            bind_rows(.id = "region") %>% 
            write_delim(here("output", "sheep_top_haps_500.txt"), " ")

# haplotype frequency over time
#write_delim(haps_all, here("output", "haps_all_500.txt"), " ")

focal_hap <- top_haps$hap[1]
haps_all %>% 
  filter(as.numeric(as.character(birth_year)) > 1990) %>% 
  group_by(birth_year, region) %>% 
  summarise(freq = sum(gt) / (n()*2)) %>% 
  ungroup() %>% 
  ggplot(aes(birth_year, freq, group = 1)) +
    geom_line() +
    facet_wrap(~region, scales = "free_y")

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
                   lamb = ifelse(age == 0, 1, 0),
                   lamb = as.factor(lamb),
                   life_stage = case_when(
                     age == 0 ~ "lamb",
                     age > 0 & age <= 2 ~ "early_life",
                     age > 2 & age <= 4 ~ "mid_life",
                     age > 4 ~ "late_life",
                   )) 
  
  fit <- glmer(survival ~ gt + sex + twin + froh_std + lamb + age_std * sex + age_std2 * sex + (1|sheep_year) + (1|birth_year) + (1|id),
               data = dat, family = binomial(link = "probit"), 
               control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
  fit
  # fit %>% 
  #   tidy(conf.int = TRUE) %>% 
  #   mutate_if(is.numeric, round, 5)
}

locations <- unique(haps_all$region)
survival <- map(locations, survival_mod, haps_all)
names(survival) <- locations
survival

saveRDS(survival, file = "output/surv_mods_hap500.rds")

#survival_f <- map(locations, survival_mod, haps_all, "F")
#survival_m <- map(locations, survival_mod, haps_all, "M")
library(performance)
map(survival, tidy, conf.int=TRUE)
binned_residuals(survival[[2]])


library(brms)
fit_brm <- brm(survival ~ gt + sex + twin + froh_std + life_stage * sex + (1|sheep_year) + (1|id), # (1|birth_year) +
               data = dat, family = bernoulli(link = "logit"))


# haplotype effects on first year survival?
# haplotype effects on annual survival?
survival_mod2 <- function(location, haps_all, s) {
  
  dat <- haps_all %>% 
    filter(region == location,
           sex == s,
           age == 0) %>% 
    mutate(gt = as.factor(gt),
           froh_scaled = scale(froh_all))
  fit <- glmer(survival ~ gt + twin + froh_scaled + (1|birth_year) + (1|mum_id),
               data = dat, family = binomial(link = "probit"),
               control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
  fit %>% 
    tidy(conf.int = TRUE) %>% 
    mutate_if(is.numeric, round, 5)
}

locations <- unique(haps_all$region)
survival <- map(locations, survival_mod2, haps_all, "F")
survival2_f <- map(locations, survival_mod2, haps_all, "F")
survival2_m <- map(locations, survival_mod2, haps_all, "M")

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
