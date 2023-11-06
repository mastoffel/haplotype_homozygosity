# Main analysis to find embryonic lethals
# Idea: For a focal haplotype, find triplets where mum and dad are carriers
# and check whether there are fewer homozygous offspring than expected 
# -> Figure 1

library(data.table)
library(tidyverse)
#library(collapse)
library(glue)
library(here)
library(furrr)
library(data.table)

args_in <- commandArgs(trailingOnly=TRUE)
hap_length <- as.numeric(args_in[1])

# carrier mating tests
load(here("data", "sheep_ped.RData"))
ped <- as_tibble(sheep_ped) %>% 
        setNames(c("id", "mum", "dad"))

# get individuals with genotypes
inds <- fread(here("data", "plink", "sheep.fam")) %>% 
                select(V2) %>% 
                rename(id = V2) %>% 
                .$id

# filter ped for individuals with genotypes (id, mum and dad must have genotypes)
ped <- ped %>% 
        filter(id %in% inds) %>% 
        drop_na() %>% 
        filter(mum %in% inds,
               dad %in% inds) %>% 
        mutate_if(is.character, as.numeric) %>% 
        as.data.table()

# chromosome
# chr <- 26

run_hap_hom_by_chr <- function(chr) {
        
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
        
        
        # calculate homozygote deficiency
        calc_hom_def <- function(hap, hap_mat) {
                
                # who carries the haplotype
                hap_num <- (hap_mat == hap) * 1
                # calculate probability of transmitting haplotype for id, mum, dad
                hap_prob <- (hap_num[, c(1, 3, 5)] + hap_num[, c(2, 4, 6)])/2
                
                # filter parents where both mum and dad are carriers
                hap_prob_carriers <- hap_prob[(hap_prob[, "mum_hap_a"] > 0) & (hap_prob[, "dad_hap_a"] > 0),,drop=FALSE]
                num_carrier_matings <- sum((hap_prob[, "mum_hap_a"] > 0) & (hap_prob[, "dad_hap_a"] > 0))
                
                # when these are 0, return NA
                if (num_carrier_matings <= 1) return(tibble(hap = hap, p_val = NA_real_,
                                                            obs = NA_real_, exp = NA_real_,
                                                            num_carriers = 0))
                
                # number offspring w homozygous haplotype 
                O_k_hom <- sum(hap_prob_carriers[, "id_hap_a"] == 1)
                # expected offspring w homozygous haplotype
                E_k_hom <- sum(hap_prob_carriers[, "mum_hap_a"] * hap_prob_carriers[, "dad_hap_a"])
                
                # number of offspring wo homozygous focal haplotype (i.e. het or hom for alternative haplotype)
                O_k_nonhom <- num_carrier_matings - O_k_hom
                E_k_nonhom <- num_carrier_matings - E_k_hom
                
                # calculate chisq statistic and p-value
                chisq <- ((O_k_hom - E_k_hom)^2)/E_k_hom + ((O_k_nonhom - E_k_nonhom)^2)/E_k_nonhom
                chi_p <- pchisq(chisq, df = 1, lower.tail = FALSE)
                
                cst_p <- tibble(hap = hap, p_val = chi_p, obs = O_k_hom, 
                                exp = E_k_hom, num_carriers = num_carrier_matings)
                cst_p
        }
        
        # running window haplotypes
        test_hap_hom <- function(start_snp, haps_raw, hap_length, calc_hom_def) {
                
                # reshape 
                haps <- apply(haps_raw[start_snp:(start_snp+hap_length-1), ], 
                               paste, collapse = "", MARGIN = 2) %>% 
                        enframe(name = "id_hap", value = "hap") %>% 
                        mutate(id = str_remove(id_hap, "_[a-z]")) %>%
                        select(-id_hap) %>% 
                        mutate(hap_pos = rep(c("hap_a", "hap_b"), nrow(.)/2), .before = 1) %>% 
                        mutate(id = as.numeric(id)) %>% 
                        pivot_wider(names_from = hap_pos, values_from = hap)
                
                # list haplotypes
                haps_tab <- table(as.matrix(haps[, c("hap_a", "hap_b")]))
                
                # haplotypes with high frequencies hf > 1%
                haps_hf <- haps_tab[haps_tab/sum(haps_tab) > 0.005]
                
                # get haplotypes matched for parents and offspring
                haps_all <- ped %>% 
                        merge(haps, by = "id") %>% 
                        setnames(old = c("hap_a", "hap_b"), new = paste0("id_", c("hap_a", "hap_b"))) |>
                        merge(haps, by.x = "mum", by.y = "id") |>
                        setnames(old = c("hap_a", "hap_b"), new = paste0("mum_", c("hap_a", "hap_b"))) |>
                        merge(haps, by.x = "dad", by.y = "id") |>
                        setnames(old = c("hap_a", "hap_b"), new = paste0("dad_", c("hap_a", "hap_b"))) 
                
                # make a matrix 
                hap_mat <- as.matrix(haps_all[, 4:9])
                
                # get transmission probabilities 
               # hap_one <- names(haps_hf)[1]
                
                hap_tests <- map_dfr(names(haps_hf), calc_hom_def, hap_mat) %>% 
                        mutate(chr = chr, snp_start = start_snp) 
                
                # sometimes a genotyping error causes a haplotype to be 
                # 2 haplotypes, both of which are significant. If this happens
                # let's combine these haplotypes 
                
                # check whether more than one haplotype is fairly significant
                hap_sub <- hap_tests %>% 
                        filter(p_val < (0.05/1000))
                # if so, check whether these are similar
                if (nrow(hap_sub) == 2){
                        hap1 <- hap_sub$hap[1]
                        hap2 <- hap_sub$hap[2]
                        diff <- str_split(hap1, "")[[1]] == str_split(hap2, "")[[1]]
                        # if they differ only by 1 
                        if (sum(!diff) < 2) {
                                hap_mat[hap_mat==hap2] <- hap1
                                haps_hf2 <- names(haps_hf)[names(haps_hf) != hap2]
                                hap_tests <- map_dfr(haps_hf2, calc_hom_def, hap_mat) %>% 
                                        mutate(chr = chr, snp_start = start_snp) 
                        }
                }
                
                hap_tests
                
        }
        
        # number of snps
        num_snps <- nrow(haps_raw)
        # loop safely
        test_hap_hom_safe <- safely(test_hap_hom)
        # run in parallel
        plan(multiprocess, workers = 4)
        
        out <- future_map(1:num_snps, test_hap_hom_safe, # num_snps
                          haps_raw = haps_raw, hap_length = hap_length,
                          calc_hom_def = calc_hom_def,
                          .progress = TRUE,
                          .options = furrr_options(seed = 123))
        
        results <- out %>% 
                transpose() %>% 
                simplify_all() %>% 
                .$result %>% 
                bind_rows()
        
        # write to folder by haplotype length
        file_name <- glue("hap_res_chr_{chr}.txt")
        dir_name <- here("output",glue("hap_len_{hap_length}"))
        dir.create(dir_name)
        write_delim(results, glue("{dir_name}/{file_name}"), "\t")
        
}

walk(1:26, run_hap_hom_by_chr)
