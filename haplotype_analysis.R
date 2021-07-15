library(data.table)
library(tidyverse)
#library(collapse)
library(glue)
library(here)
library(furrr)


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
                
                # number parents where mum and dad have are carriers
                num_carrier_matings <- sum((hap_prob[, "mum_hap_a"] > 0) & (hap_prob[, "dad_hap_a"] > 0))
                # when these are 0, return NA
                if (num_carrier_matings == 0) return(tibble(hap = hap, p_val = NA_real_))
                
                # number offspring w homozygous haplotype 
                O_k_hom <- sum(hap_prob[, "id_hap_a"] == 1)
                E_k_hom <- sum(hap_prob[, "mum_hap_a"] * hap_prob[, "dad_hap_a"])
                
                # number of offspring wo homozygous focal haplotype (i.e. het or hom for alternative haplotype)
                O_k_nonhom <- num_carrier_matings - O_k_hom
                E_k_nonhom <- num_carrier_matings - E_k_hom
                
                # calculate chisq statistic and p-value
                chisq <- ((O_k_hom - E_k_hom)^2)/E_k_hom + ((O_k_nonhom - E_k_nonhom)^2)/E_k_nonhom
                chi_p <- pchisq(chisq, df = 1, lower.tail = FALSE)
                
                # 
                #O_k_nonhom <-  num_carrier_matings - O_k_hom
                
                # Fritz et al. expected offspring carrying homozygous haplotype / proportion for chi squ
                #E_k_hom <- sum(hap_prob[, "mum_hap_a"] * hap_prob[, "dad_hap_a"])/num_carrier_matings
                #E_k_nonhom <- 1 - E_k_hom 
                #chi <- chisq.test(x = c(O_k_hom, O_k_nonhom), p = c(E_k_hom, E_k_nonhom))
                
                cst_p <- tibble(hap = hap, p_val = chi_p)
                cst_p
        }
        
        # running window haplotypes
        test_hap_hom <- function(start_snp, haps_raw, hap_length, calc_hom_def) {
                
                # reshape 
                haps <- apply(haps_raw[start_snp:(start_snp+hap_length), ], 
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
                haps_hf <- haps_tab[haps_tab/sum(haps_tab) > 0.01]
                names(haps_hf)
                
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
                hap_one <- names(haps_hf)[1]
                
                out <- map_dfr(names(haps_hf), calc_hom_def, hap_mat) %>% 
                        mutate(chr = chr, snp_start = start_snp) 
                
                out
        }
        
        # number of snps
        num_snps <- nrow(haps_raw)
        # loop safely
        test_hap_hom_safe <- safely(test_hap_hom)
        # run in parallel
        plan(multiprocess, workers = 4)
        
        out <- future_map(1:num_snps, test_hap_hom_safe, # num_snps
                          haps_raw = haps_raw, hap_length = 20,
                          calc_hom_def = calc_hom_def,
                          .progress = TRUE,
                          .options = furrr_options(seed = 123))
        
        results <- out %>% 
                transpose() %>% 
                simplify_all() %>% 
                .$result %>% 
                bind_rows()

        file_name <- glue("hap_res_chr_{chr}.txt")
        write_delim(results, here("output", file_name), "\t")
        
}

walk(26:1, run_hap_hom_by_chr)
