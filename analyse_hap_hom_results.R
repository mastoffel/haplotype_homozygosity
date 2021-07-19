library(tidyverse)
library(here)

# read results from haplotype homozygosity scan
all_files <- list.files(here("output", "hap_len_100"), full.names = TRUE)
results <- map(all_files, read_delim, delim = "\t") %>% 
                bind_rows()

# any potential lethals?
results %>% 
        filter(obs == 0 & exp > 9)

res2 <- results %>% 
        filter(obs < exp) %>% 
        filter(p_val < 0.05) #%>% 
        #filter(chr == 2)

ggplot(res2, aes(snp_start, -log10(p_val))) + 
        geom_point(size = 1, alpha = 1) +
        facet_wrap(~chr, scales = "free_x") +
        geom_hline(yintercept = -log10(0.05/38000)) #+
       # ylim(c(0, 8))

-log10(0.05/38000)

results %>% 
        filter(p_val < (0.05/40000)) %>% 
        filter(chr == 4) %>% 
        arrange(p_val)

snp_map <- read_delim(here("data", "plink", "sheep.bim"), delim = "\t",
                      col_names = FALSE) %>% 
                        setNames(c("chr", "snp", "cM", "bp", "a1", "a2")) %>% 
                        filter(chr == 17)



