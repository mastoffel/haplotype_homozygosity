library(data.table)
library(tidyverse)
library(here)
library(glue)
library(slider)
library(windowscanr)
sheep_genos <- here("data", "plink", "sheep")
roh <- here("output", "roh")

# run ROH from plink (with data filtered for low call SNPs)
system(glue("plink --bfile {sheep_genos} --sheep --out {roh}/sheep ",
            "--homozyg --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 1000 ",
            "--homozyg-gap 300 --homozyg-density 200 --homozyg-window-missing 2 ",
            "--homozyg-het 2 --homozyg-window-het 2 --keep data/id.txt"))

homs <- fread(glue("{roh}/sheep.hom"))
homs %>% 
        filter(IID == 7695) %>% 
        select(CHR, POS1, POS2)







out <- winScan(homs, 
        groups = "CHR",
        position = "BP",
        values = "UNAFF",
        win_size = 1000000,
        win_step = 500000,
        funs = c("mean"))

frollmean(homs, .)

slide_period(homs$UNAFF, )

?windowscanr
