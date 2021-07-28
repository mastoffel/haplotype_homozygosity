library(RJDBC)
library(dplyr)
library(dbplyr)
library(lubridate)
library(tidyverse)
library(janitor)

# load fitness data
load("data/fitness_roh.RData") 
fitness <- fitness_data %>% 
        select(id, survival, sheep_year, age, birth_year, sex, mum_id,
               twin, offspring_born, offspring_survived, weight,
               hindleg, horn, father, froh_all, birthwt) %>% 
        group_by(id) %>% 
        mutate(lifespan = max(age)) %>% 
        mutate(id = as.numeric(id),
               year_mating = as.numeric(as.character(sheep_year)) - 1)

dbname <- "../sheep/data/db/StKilda_Data.accdb"
driver <- "net.ucanaccess.jdbc.UcanloadDriver"
driverpath <- "../sheep/data/db/UCanAccess/loader/ucanload.jar"
options <- paste0("jdbc:ucanaccess://", dbname, ";memory=false")

con <- DBI::dbConnect(JDBC(driver, driverpath), options)
# src <- src_dbi(con)

tbls <- dbGetTables(con)
flds <- dbGetFields(con, "tblPregnancies")
consorts <- dbGetQuery(con, "Select * from Consorts")
dbDisconnect(con)

# -1 means mounting observed, 0 means mounting not observed
cons <- consorts %>% 
        as_tibble() %>% 
        clean_names() %>% 
        separate(col = date, into = c("date", "out"), sep = " ") %>% 
        separate(col = time, into = c("out2", "time"), sep = " ") %>% 
        select(-out, -out2) %>% 
        mutate(date = ymd(date),
               time = hms(time))
table(cons$mounts)

# first approach: take females which mated a second time after two week
# assuming that the pregnancy failed. Take only
pot_failed_matings <- cons %>% 
        mutate(month = month(date),
               year = year(date)) %>% 
        # this finds ewes mating with several males
        group_by(year, ewe_id) %>% 
        #filter(mounts == -1) %>% 
        #filter(year > 2000) %>% 
        filter(n() > 1) %>% 
        arrange(desc(ewe_id)) %>% 
        # ewe being consorted by at least two males
        filter(n_distinct(tup_id) > 1) %>% 
        arrange(year, ewe_id, date) %>% 
        mutate(time_lag = date - lag(date)) %>% 
        # filter ewes where any two consorts are at least 10 days apart
        filter(any(time_lag > 7)) %>% 
        # first mating gets NA, so rather give it 0
        mutate(time_lag = ifelse(is.na(time_lag), 0, time_lag)) %>% 
        # filter matings before the gap (potentially failed pregnancies)
        slice(1:(match(TRUE, time_lag == max(time_lag)))-1) %>% 
        # remove unknown tups
        filter(!is.na(tup_id)) %>% 
        # filter all which were mated/seen mated only once in the previous estrous
        filter(n() == 1) %>% 
        select(date, tup_id, ewe_id) %>% 
        mutate(failed = 1)

pot_succ_matings <- cons %>% 
        mutate(month = month(date),
               year = year(date)) %>% 
        # this finds ewes mating with several males
        group_by(year, ewe_id) %>% 
        #filter(mounts == -1) %>% 
        #filter(year > 2000) %>% 
        filter(length(unique(tup_id))==1) %>% 
        arrange(ewe_id) %>%
        select(date, tup_id, ewe_id) %>% 
        mutate(failed = 0) 

# combine        
all_matings <- bind_rows(pot_failed_matings, pot_succ_matings)

# load haplotypes
haps <- read_delim(here("output", "sheep_top_haps.txt"), " ") %>% 
                select(region, id, gt)

gts <- 


# join failed matings \ 2 or 1 means individual carries focal haplotype
gts <- pot_failed_matings %>% 
        left_join(haps, by = c("tup_id" = "id")) %>% 
        left_join(haps, by = c("ewe_id" = "id", "region" = "region")) %>% 
        filter(!is.na(gt.x) & !is.na(gt.y))

exp_homs <- gts %>% 
        group_by(region) %>% 
        mutate(offspr_hom_prob = case_when(
                gt.x == 1 & gt.y == 1 ~ 0.25,
                gt.x == 2 & gt.y == 1 ~ 0.5,
                gt.x == 1 & gt.y == 2 ~ 0.5,
                TRUE ~ 0
        )) %>% 
        summarise(exp_hom_offspr = sum(offspr_hom_prob)/n())
        

# second approach: which matings did not lead to offspring the next year?

# first approach: take females which mated a second time after two week
# assuming that the pregnancy failed
failed_preg <- cons %>% 
        mutate(month = month(date),
               year = year(date)) %>% 
        # this finds ewes mating with several males
        group_by(year, ewe_id) %>% 
        #filter(mounts == -1) %>% 
        #filter(year > 2000) %>% 
       # filter(n() > 1) %>% 
        arrange(desc(ewe_id)) %>% 
        # ewe being consorted by at least two males
        filter(n_distinct(tup_id) == 1) %>% 
        arrange(year, ewe_id, date) %>% 
        filter(!is.na(tup_id)) %>% 
        # filter all which were mated/seen mated only once in the previous estrous
        slice(1) %>% 
        left_join(fitness, by = c("ewe_id" = "id", "year" = "year_mating")) %>% 
        rename(year_mating = year) %>% 
        filter(offspring_born == 0) %>% 
        select(date, tup_id, ewe_id, offspring_born, survival) %>% 
        # chose only individuals which survived to the next year
        filter(survival == 1) %>% 
        left_join(haps, by = c("tup_id" = "id")) %>% 
        left_join(haps, by = c("ewe_id" = "id", "region" = "region")) %>% 
        filter(!is.na(gt.x) & !is.na(gt.y))

failed_preg %>% 
        group_by(region) %>% 
        summarise(freq = sum(gt.x + gt.y)/(n()*4))

