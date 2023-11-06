# Back of the envelope popgen calculations
# Question: How do we expect haplotype frequencies to change over time,
# assuming single locus models?

# Calculate expected frequencies using single locus pop gen
# assumes random mating, large population size etc.
# also assumes deleterious allele to be fully recessive
# Falconer and Mackay (Introduction to Quantitative Genetics), p28


# expected decline for SEL05 and SEL07 -----------------------------------------

# let's assume a soay sheep generation time of 3 years
# and time frame 1990 - 2018, then we look at 18/3 = 6 generations
# 

update_freqs <- function(s, q, num_gen = 6, qs = c(), max_gen = 6) {
        #delta_q <- - (q^2 * (1-q)) / (1 - s*(1-q^2))
        delta_q <- -(s*q^2*(1-q) / (1-s*q^2)) 
        new_q <- q + delta_q
        if (num_gen == max_gen) {
                qs <- c(q, new_q) # if first iteration, add initial q value
        } else {
                qs <- c(qs, new_q)  # Store the current value of q  
        }
        
        # Termination condition: specified number of iterations reached
        if (num_gen == 1) return(qs)
        
        return(update_freqs(s, new_q, num_gen - 1, qs))
}


# sel07 (emp: 19 to 7%)
s_sel07 <- 0.47 
q_init_sel07 <- 0.19
fs07 <- update_freqs(s_sel07, q_init_sel07)
# -> freq after six generations: 12.8%
fs07

# sel18 (emp: 32 to 18%)
s_sel18 <- 0.31
q_init_sel18 <- 0.32 # 0.32
fs18 <- update_freqs(s_sel18, q_init_sel18)
# -> freq after six generations: 21.7%
fs18


s_sel18 <- 0.31
q_init_sel18 <- 0.20
fs18 <- update_freqs(s_sel07, q_init_sel07)


# equilibrium frequency for SEL18 ----------------------------------------------
# Formula for Overdominance from Falconer and Mackay (Table 2.2 (5))
# Fitnesses of s1 (A1A1) and s2 (A2A2) expressed relative to A1A2

# SEL05: s = -0.27, emp. freq from 20 - 23%, het higher survival prob: 6.6%

# Define the equation for the change in allele frequency
delta_q <- function(q) {
        p <- 1 - q
        s1 <- -0.066
        s2 <- -0.33
        return(abs(p*q*(s1*p - s2*q) / (1 - s1*p^2 - s2*q^2)))
}

# find the minimum of the function,
# where delta_q is 0 (i.e. equilibrium)
solution <- optim(0.5, delta_q, method = "Brent", lower = 0, upper = 1)

# Equilibrium frequency of the A2 allele
cat("Equilibrium frequency of A2 (q):", solution$par, "\n")

# Equilibrium frequency of the A1 allele
cat("Equilibrium frequency of A1 (p):", 1 - solution$par, "\n")


