# format pedigree for alpha peel

library(tidyverse)

load("data/sheep_ped.RData")
ped_alphapeel <- sheep_ped %>% 
        select(ID, FATHER, MOTHER) %>% 
        mutate(
                across(everything(), ~replace_na(.x, 0))
        ) 
write_delim(ped_alphapeel, "data/ped.txt", " ", col_names = FALSE)
