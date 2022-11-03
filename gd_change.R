gd_change <- function (genedrop_object_summary, n_founder_cohorts = NULL, 
          remove_founders = T, obs_line_col = "red") 
{
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
                genedrop_object_summary$simulated_frequencies$Allele <- "p"
                genedrop_object_summary$observed_frequencies$Allele <- "p"
        }
        Allele <- as.character(Allele)
        cumu.func <- function(x) {
                x <- diff(x)
                x <- ifelse(x < 0, x * -1, x)
                sum(x, na.rm = F)
        }
        sim.changes <- tapply(genedrop_object_summary$simulated_frequencies$p, 
                              list(genedrop_object_summary$simulated_frequencies$Simulation, 
                                   genedrop_object_summary$simulated_frequencies$Allele), 
                              cumu.func)
        sim.changes <- melt(sim.changes)
        true.changes <- tapply(genedrop_object_summary$observed_frequencies$p, 
                               list(genedrop_object_summary$observed_frequencies$Simulation, 
                                    genedrop_object_summary$observed_frequencies$Allele), 
                               cumu.func)
        true.changes <- melt(true.changes)
        print(ggplot(sim.changes, aes(value)) + geom_histogram() + 
                      facet_wrap(~Var2) + geom_vline(data = true.changes, aes(xintercept = value), 
                                                     col = obs_line_col) + ggtitle(paste0("Distribution of Cumulative Change: Nsim = ", 
                                                                                          max(genedrop_object_summary$simulated_frequencies$Simulation))))
        true.changes$Cumulative.Change.Lower <- NA
        true.changes$Cumulative.Change.Higher <- NA
        for (i in 1:nrow(true.changes)) {
                true.changes$Cumulative.Change.Lower[i] <- length(which(sim.changes$Var2 == 
                                                                                true.changes$Var2[i] & sim.changes$value < true.changes$value[i]))
                true.changes$Cumulative.Change.Higher[i] <- length(which(sim.changes$Var2 == 
                                                                                 true.changes$Var2[i] & sim.changes$value > true.changes$value[i]))
        }
        out <- list(sim.changes, true.changes)
        out
}