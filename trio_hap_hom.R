library(data.table)
library(tidyverse)
library(furrr)
library(glue)
# carrier mating tests
load("data/sheep_ped.RData")
ped <- as_tibble(sheep_ped)
genos <- fread(glue("data/genos/chr23.txt")) %>% 
                rename(id = V1)

# filter ped for individuals with genotypes
ped <- ped %>% 
        filter(ID %in% genos$id) %>% 
        drop_na() %>% 
        filter(MOTHER %in% genos$id,
               FATHER %in% genos$id) %>% 
        rename(id = ID, mum = MOTHER, dad = FATHER) %>% 
        mutate_if(is.character, as.numeric) %>% 
        as.data.table()


# at a given SNP, test for transmission distortion for minor allele homozygosity
# this is done by selecting heterozygous parents
# probability for offspring to have hh is thus 25%

hom_test <- function(snp, ped, genos_chr) {
        
        cols <- c("id", snp)
        geno_sub <- genos_chr[, ..cols]
        
        ped_geno <- ped |> 
               # as.data.table() |>
                merge(geno_sub, by = "id") |>
                setnames(old = snp, new = "snp_id") |>
                merge(geno_sub, by.x = "mum", by.y = "id") |>
                setnames(old = snp, new = "snp_mum") |>
                merge(geno_sub, by.x = "dad", by.y = "id") |> 
                setnames(old = snp, new = "snp_dad")
        
        gt <- ped_geno[snp_mum == 1 & snp_dad == 1 & snp_id != 9, .N, by =.(snp_id)]
        
        # if there are no heterozygote matings
        if (is.data.frame(gt) && nrow(gt) == 0) return(NA_real_)
        
        # get count of homozygous for minor allele
        hh <- gt[snp_id == 2, N]
        if (length(hh) == 0) hh <- 0
        # get count of other genotypes
        non_hh <- sum(gt[, N]) - hh
        #all_g <- sum(gt[, N])
        #hH <- gt[snp_id == 1, N]
        #HH <- gt[snp_id == 0, N]
        # chisqu
        cst <- chisq.test(x = c(hh, non_hh), p = c(0.25, 0.75))$p.value
        
        # fisher exact
        #fet <- rbind(c(hh, non_hh), c(0.25 * (hh + non_hh), 0.75 * (hh + non_hh)))
        #fet <- fisher.test(x, alternative = "less")$p.value 
        
}

hap_hom_gwas <- function(chr_num) {
        
        gc()
        print(glue("running transmission distortion tests on chromosome {chr_num}"))
        # get all genos
        genos <- fread(glue("data/genos/chr{chr_num}.txt")) %>% rename(id = V1)
        
        snp_names <- names(genos[2:ncol(genos)])
        #snp_names <- snp_names[seq(1, length(snp_names), by = 10)]
        hom_test_safely <- safely(hom_test, otherwise = NA_real_)
        
        plan(multiprocess, workers = 4)
        out <- future_map(snp_names, hom_test_safely, ped, genos, .progress = TRUE) # ,
        
        all_p <- out %>% 
                transpose() %>% 
                simplify_all() %>% 
                .$result
        
        all_p

}


out <- map(1:26, hap_hom_gwas)
saveRDS(out, file = "output/trio_transmission_distortion2.RDS")
out

#out2 <- readRDS("output/trio_transmission_distortion.RDS")

test <- rbind(c(29, 99), c(0.25*128, 0.75*128))
fisher.test(x = test, alternative = "less")


plot(-log10(unlist(out)))


max(-log10(out[[5]]), na.rm=TRUE)

which.max(-log10(out[[5]]))

snp <- "V44883"

chr_num <- 20

plot(-log10(all_p))

which.max(-log(all_p))
which(-log(all_p) == max(-log(all_p), na.rm = TRUE))

all <-  out %>% 
        transpose() %>% 
        simplify_all()
 
which.max(-log(all$result))

out <- hom_test("V8662", ped, genos)

# population test
# chop off ids
genos <- as.matrix(genos[, 2:ncol(genos)])
genos[genos == 9] <- NA

# find complete lack of hh animals
no_hh <- colSums(genos == 2, na.rm = TRUE) == 0
num_geno <- colSums(!is.na(genos))
freq_Hh <- colSums(genos == 1, na.rm = TRUE) / num_geno

Phh <- (1 - (freq_Hh^2 / 4)) ^ num_geno
Phh_where_no_hh <- Phh[no_hh]
Phh_where_no_hh
hist(-log(Phh_where_no_hh))