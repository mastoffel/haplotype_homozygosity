gd_slopes <- function (genedrop_object_summary, n_founder_cohorts = NULL, 
          remove_founders = T, method = "lm", obs_line_col = "red") 
{
        sim.slopes <- NULL
        genedrop_object_summary$simulated_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Simulation))
        genedrop_object_summary$simulated_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Cohort))
        genedrop_object_summary$observed_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Simulation))
        genedrop_object_summary$observed_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Cohort))
        if (remove_founders) {
                if (length(n_founder_cohorts) == 1) {
                        genedrop_object_summary$simulated_frequencies <- subset(genedrop_object_summary$simulated_frequencies, 
                                                                                !Cohort %in% unique(sort(genedrop_object_summary$simulated_frequencies$Cohort))[1:n_founder_cohorts])
                        genedrop_object_summary$observed_frequencies <- subset(genedrop_object_summary$observed_frequencies, 
                                                                               !Cohort %in% unique(sort(genedrop_object_summary$observed_frequencies$Cohort))[1:n_founder_cohorts])
                }
        }
        if ("Allele" %in% names(genedrop_object_summary$simulated_frequencies)) {
                Allele = unique(genedrop_object_summary$simulated_frequencies$Allele)
        }
        else {
                Allele = "p"
        }
        for (i in 1:max(genedrop_object_summary$simulated_frequencies$Simulation)) {
                if (i %in% seq(1, max(genedrop_object_summary$simulated_frequencies$Simulation), 
                               100)) 
                        print(paste("Calculating slope", i, "in", max(genedrop_object_summary$simulated_frequencies$Simulation)))
                for (j in Allele) {
                        x <- subset(genedrop_object_summary$simulated_frequencies, 
                                    Simulation == i & Allele == j)
                        x1 <- lm(p ~ Cohort, data = x)$coefficients[[2]]
                        sim.slopes <- rbind(sim.slopes, data.frame(Iteration = i, 
                                                                   Allele = j, Slope = x1))
                        rm(x1)
                }
        }
        true.slopes <- NULL
        for (i in 0) {
                for (j in Allele) {
                        x <- subset(genedrop_object_summary$observed_frequencies, 
                                    Simulation == i & Allele == j)
                        x1 <- lm(p ~ Cohort, data = x)$coefficients[[2]]
                        true.slopes <- rbind(true.slopes, data.frame(Iteration = i, 
                                                                     Allele = j, Slope = x1))
                        rm(x1)
                }
        }
        print(ggplot(sim.slopes, aes(Slope)) + geom_histogram() + 
                      facet_wrap(~Allele) + geom_vline(data = true.slopes, 
                                                       aes(xintercept = Slope), col = obs_line_col) + ggtitle(paste0("Distribution of Regression Slopes: Nsim = ", 
                                                                                                                     max(genedrop_object_summary$simulated_frequencies$Simulation))))
        true.slopes$Slopes.Lower <- NA
        true.slopes$Slopes.Higher <- NA
        for (i in 1:nrow(true.slopes)) {
                true.slopes$Slopes.Lower[i] <- length(which(sim.slopes$Allele == 
                                                                    true.slopes$Allele[i] & sim.slopes$Slope < true.slopes$Slope[i]))
                true.slopes$Slopes.Higher[i] <- length(which(sim.slopes$Allele == 
                                                                     true.slopes$Allele[i] & sim.slopes$Slope > true.slopes$Slope[i]))
        }
        out <- list(sim.slopes, true.slopes)
        out
}