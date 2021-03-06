library(tidyverse)
source("theme_simple.R")
library(ggdist)
library(ggrepel)
library(ggridges)
library(tidybayes)

# 10 simulation runs 
files <- list.files("slim_sim/sims/muts_shape05_scale01/", full.names = TRUE)
full <- map_dfr(files, read_delim, .id = "run")
n_sample <- 200
muts <- full %>% 
        #sample_frac(0.1) %>% 
        group_by(run, mut_id, pos, s, originG) %>% 
        tally() %>% 
        mutate(freq = n/(n_sample*2)) %>% 
        ungroup() #%>% 
        #mutate(classes = cut(s, breaks = seq(0, -1, -0.05))) %>% 
       # mutate(classes = cut(s, breaks = 30)) %>% 
        #group_by(classes) %>% 
        #summarise(n = mean(n), freq = mean(freq))

# get proportion moderate (10 < NeS < 100) muts
# and strong (NeS > 100)
muts %>% 
  mutate(s_class = case_when(
      200*s > -1 ~ "neutral",
      200*s > -10 & 200*s < -1 ~ "small",
      (200*s < -10) & (200*s > -100) ~ "moderate",
      200*s <= -100 ~ "strong")) %>% 
  group_by(run, s_class) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

full %>% 
  group_by(run) %>% 
  mutate(s_class = case_when(
    200*s > -1 ~ "neutral",
    200*s > -10 & 200*s < -1 ~ "small",
    (200*s < -10) & (200*s > -100) ~ "moderate",
    200*s <= -100 ~ "strong")) %>%  
  group_by(s_class) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

  
  
ggplot(muts, aes(s, freq)) +
  geom_hex(bins = 30) + 
  scale_fill_viridis_c(trans = "log", direction = 1, option = "D",
                       breaks = c(10, 100, 1000, 5000),
                       labels = c(10, 100, 1000, 5000)/10) +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE)# + 

#ggsave(filename = "figs/s_vs_freq.jpg",p, width = 5.5, height = 4)

# expected
get_p <- function(freq, s, n = 5952) {
  E_k_hom <- freq^2 * n
  O_k_hom <- (1-s) * E_k_hom
  O_k_nonhom <- n-O_k_hom
  E_k_nonhom <- n-E_k_hom
  chisq <- ((O_k_hom - E_k_hom)^2)/E_k_hom + ((O_k_nonhom - E_k_nonhom)^2)/E_k_nonhom
  chi_p <- pchisq(chisq, df = 1, lower.tail = FALSE)
  out <- chi_p * 39149
  out <- ifelse(out>1, 1, out)
  out
}

pars <- expand_grid(freq = seq(0, 1, 0.001), 
                    s = seq(0, -1, -0.001))

pars$p <- map2(pars$freq,pars$s, get_p) %>% 
  unlist()

pars_thres <- pars %>% 
  group_by(s) %>% 
  filter(p <= 0.05) %>% 
  filter(freq == min(freq)) %>% 
  ungroup()
#pars_thres[8, 2] <- -0.3

p1 <- ggplot(muts, aes(s, freq)) +
  geom_hex(bins = 30) + 
  scale_fill_viridis_c("# mutations\nper simulation",
                       trans = "log", direction = 1, option = "D",
                       breaks = c(1, 10, 100, 1000, 5000),
                       labels = c(1, 10, 100, 1000, 5000)/10) +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, -0.2, -0.4, -0.6, -0.8, -1.0)) +
  geom_line(data = pars_thres, color = "#D8DEE9",
            linetype = 1, size = 0.8) +
  theme(axis.line.y = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="#D8DEE9"),
        axis.line.x = element_blank(),
        axis.ticks = element_blank()) +
  annotate("text", x = -0.8, y = 0.18, label = "Detection\nthreshold",
           family = "Lato", color = "#4C566A")+
  xlab("Selection coefficient") +
  ylab("Population frequency") 
p1

# plot 2
n <- 200
gamma <- data.frame(s = c(-rgamma(n*0.95, shape = 0.5, scale = 0.1), rep( -1,n*0.05)))
mean(rgamma(1000000, shape = 0.5, scale = 0.1))

df <- full %>% 
  #filter(run == 1) %>% 
  bind_rows(gamma, .id = "type") 

ggplot(full) +
  geom_histogram(
                 aes(x = s, y = -stat(count) / sum(count)), fill="#D8DEE9", 
                 bins = 20, color = "black", size = 0.1)


# h = 0.5e-13s
p2 <- ggplot(df %>% filter(type == 2), aes(x=s, color = type)) +
  geom_histogram( aes(x = s, y =  stat(count) / sum(count)), fill="#4C566A", 
                  bins = 20, size = 0.1, color = "black") +
  geom_function(fun = function(s) 0.5*exp(13 * s), mapping = aes(x=s),
                color = "#81A1C1", size = 0.8) +
  annotate("text", x = -0.15, y = 0.2, label = "h", color = "#81A1C1") +
  geom_histogram(data = df %>% filter(type == 1), 
                 aes(x = s, y = -stat(count) / sum(count)), fill="#D8DEE9", 
                 bins = 20, color = "black", size = 0.1) +
  theme_simple(grid_lines = FALSE) +
  scale_y_continuous(name = "Proportion of mutations", breaks = seq(-0.6, 0.6, 0.2), 
                     labels = paste0(c(seq(60, 10, -20), seq(0, 60, 20)), "%")) +
                     #labels = scales::percent) +
  xlab("Selection coefficient") +
  theme(panel.grid.major.y = element_line( size=.1, color="#D8DEE9")) +
  annotate("text", label = "New mutations", x = -0.5, y = 0.3, color = "#4C566A") +
  annotate("text", label = "Simulated Soay sheep", x = -0.5, y = -0.3,
           color = "#4C566A") +
  scale_x_continuous(breaks = c(0, -0.2, -0.4, -0.6, -0.8, -1.0))
# coord_flip() 
#scale_fill_manual(name="Bar",values=c("#4C566A", "#D8DEE9")) +
p2

#ggsave("figs/dfe_vs_sim_hist.jpg",p, width = 4.5, height = 4)
#pd <- position_dodge(0.1)

library(patchwork)
p_final <- p2 + p1 +
  plot_layout(widths = c(1, 1.1)) +
  plot_annotation(tag_levels = "A")

ggsave("figs/sims.jpg", p_final, width = 9, height = 3.5)






df <- mut_classes %>% 
  mutate(s_class = str_extract(s_class, "(?<=^.{1})([^,])+"),
         s = as.numeric(s_class)) %>% 
  left_join(pars_thres)

p <- df %>% 
  ggplot(aes(x = mf, y = s_class, size = n)) + 
  geom_linerange(aes(xmin = mf, xmax = hf2), size = 0.5, color = "#88C0D0") +
  geom_linerange(aes(xmin = mf, xmax = hf), size = 1.3, color = "#81A1C1") +
  geom_point(fill= "#5E81AC", size = 3, shape = 21) +
  theme_simple(grid_lines = FALSE) +
  #geom_label(label = mut_classes$n,
  #           nudge_x = 0.05, nudge_y = 0.4, 
  #            size = 3,
  #            color = "#5E81AC") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_discrete(labels = c("-0.01", "-0.05", seq(-0.1, -0.6, -0.1), "-1")) +
  xlab("Mutation frequency\n(mean, 90th, 99th percentile)") +
  theme(legend.position = "none") +
  ylab("Selection coefficient s") #+
#geom_point(aes(x = freq, y = s_class), 
#           size = 0.8, color = "#4C566A") +
#geom_line(aes(x = freq, y = s_class, group = 1), size = 0.3, color = "#4C566A")
p
ggsave(filename = "figs/s_freq_nothresh.jpg",p, width = 4.5, height = 4)


df_sum <- df %>% 
  group_by(s_class) %>% 
  summarise(threshold = mean(freq))

ggplot(df, aes(s_class, mut_freq)) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_point(aes(s_class, threshold), data = df_sum, color = "blue", size = 2) +
  theme_simple()





s = 0.3
freq = 0.2
E_k_hom <- 500
O_k_hom <- (1-s) * 500
n <- 5952
O_k_nonhom <- n-O_k_hom
E_k_nonhom <- n-E_k_hom

# 5952
chisq <- ((O_k_hom - E_k_hom)^2)/E_k_hom + ((O_k_nonhom - E_k_nonhom)^2)/E_k_nonhom
chi_p <- pchisq(chisq, df = 1, lower.tail = FALSE)
chi_p 


p <- mut_classes %>% 
  ggplot(aes(x = mf, y = s_class)) +  # size = n
  
  # geom_point(aes(x = hf, y = s_class), 
  #           size = 2, shape = "|", stroke = 1) +
  geom_linerange(aes(xmin = mf, xmax = hf2), size = 0.5, color = "#88C0D0") +
  geom_linerange(aes(xmin = mf, xmax = hf), size = 1.3, color = "#81A1C1") +
  geom_point(fill= "#5E81AC", size = 3, shape = 21) +
  theme_simple(grid_lines = FALSE) +
  #geom_label(label = mut_classes$n,
  #           nudge_x = 0.05, nudge_y = 0.4, 
  #            size = 3,
  #            color = "#5E81AC") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_discrete(labels = c("-0.01", "-0.05", seq(-0.1, -0.6, -0.1), "-1")) +
  xlab("Mutation frequency\n(mean, 90th, 99th percentile)") +
  theme(legend.position = "none") +
  ylab("Selection coefficient s")
#xlim(c(0, 0.5))
p
ggsave(filename = "figs/s_freq.jpg",p, width = 4.5, height = 4)


# option 2 plotting

ggplot(df %>% filter(type == 2), aes(x=s, color = type)) +
  geom_histogram(data = df %>% filter(type == 1), 
                 aes(x = s, y = stat(count) / sum(count)), fill="#4C566A", 
                 bins = 50, color = "black", size = 0.1, alpha = 0.5) +
  geom_histogram( aes(x = s, y =  stat(count) / sum(count)), fill="white",  # #D8DEE9
                   bins = 50, size = 0.1, color = "black", alpha = 0.5) +
 
  theme_simple(grid_lines = FALSE) +
  scale_y_continuous(name = "% deleterious mutations", breaks = seq(-0.6, 0.6, 0.2), 
                     labels = c(seq(60, 10, -20), seq(0, 60, 20))) +
  xlab("selection coefficient") +
  theme(panel.grid.major.y = element_line( size=.1, color="#D8DEE9")) +
  annotate("text", label = "New mutations", x = -0.5, y = 0.3) +
  annotate("text", label = "Simulated Soay sheep", x = -0.5, y = -0.3) +
  scale_x_continuous(breaks = c(0, -0.2, -0.4, -0.6, -0.8, -1.0))



   

p <- ggplot(mut_classes, aes(x = freq, y=s_class)) +
        stat_interval(
                     position = position_nudge(y = -0.2),
                     orientation = "horizontal",
                     .width = c(.5,.8,.99),
                     stroke = 2,
                     size = c(6)) +
        theme_simple(grid_lines = FALSE) +
        xlim(c(0, 0.3)) +
        scale_colour_viridis_d("Percentile", direction = -1) +
        stat_dots(size = 2,
                  scale = 0.65,
                  color = "#2E3440",
                  fill = "#4C566A",
                  shape=23,
                  stroke = 0.001,
                  position = "dodge",
                  alpha = 0.1)
p
        
ggsave("figs/sims.jpg",p, width = 6, height = 6)

stat_dots(position = "dodge", orientation = "horizontal", 
                  layout = "swarm",
                  scale = 0.6,
                  size = 3,
                  color = "#6f4e37",
                  fill = "#6f4e37")
        
#stat_slab()
        stat_dots(
                quantiles = NA,
                orientation = "horizontal",
                #normalize = "none",
                #scale = .87,
                scale = 1,
                dotsize = 0.01,
                stackratio = 0.1,
               # layout = "swarm",
                color = "#6f4e37",
                fill = "#6f4e37"
        )


        

