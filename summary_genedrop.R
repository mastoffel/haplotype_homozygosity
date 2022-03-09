summary_genedrop <- function (genedrop_object, genotype_delim = "") 
{
        genedrop_object <- subset(genedrop_object, !is.na(True.Geno))
        if (is.numeric(genedrop_object$Simulated.Geno)) {
                x1 <- reshape2::melt(tapply(genedrop_object$Simulated.Geno, list(genedrop_object$cohort, 
                                                                       genedrop_object$Simulation), sum))
                x2 <- reshape2::melt(tapply(genedrop_object$Simulated.Geno, list(genedrop_object$cohort, 
                                                                       genedrop_object$Simulation), length))
                names(x2)[3] <- "Count"
                suppressMessages(x1 <- plyr::join(x1, x2))
                head(x1)
                names(x1) <- c("Cohort", "Simulation", "Sum", "Count")
                x1$p <- 1 - (x1$Sum/(2 * x1$Count))
                x1$Sum <- NULL
                x3 <- subset(genedrop_object, Simulation == 1) %>% subset(select = c(True.Geno, 
                                                                                     cohort)) %>% na.omit()
                x3 <- data.frame(Sum = tapply(x3$True.Geno, x3$cohort, 
                                              sum), Count = tapply(x3$True.Geno, x3$cohort, length))
                x3$p <- 1 - (x3$Sum/(2 * x3$Count))
                x3$Simulation <- 0
                x3$Cohort <- as.numeric(row.names(x3))
                x3$Sum <- NULL
                head(x1)
                head(x3)
                x3 <- x3[, c("Cohort", "Simulation", "Count", "p")]
        }
        else {
                head(genedrop_object)
                x1 <- melt(genedrop_object[, c("cohort", "Simulation", 
                                               "Mum.Allele", "Dad.Allele")], id.vars = c("cohort", 
                                                                                         "Simulation"))
                x2a <- data.frame(table(x1$cohort, x1$Simulation, x1$value))
                x2b <- data.frame(table(x1$cohort, x1$Simulation))
                names(x2b)[3] <- "Count"
                suppressMessages(x1 <- join(x2a, x2b))
                x1$Freq <- x1$Freq/x1$Count
                head(x1)
                names(x1) <- c("Cohort", "Simulation", "Allele", "p", 
                               "Count")
                x1 <- x1[, c("Cohort", "Simulation", "Count", "p", "Allele")]
                rm(x2a, x2b)
                x3 <- subset(genedrop_object, Simulation == 1) %>% subset(select = c(True.Geno, 
                                                                                     cohort)) %>% na.omit
                x3$Allele1 <- sapply(x3$True.Geno, function(foo) strsplit(foo, 
                                                                          split = genotype_delim, fixed = T)[[1]][1])
                x3$Allele2 <- sapply(x3$True.Geno, function(foo) strsplit(foo, 
                                                                          split = genotype_delim, fixed = T)[[1]][2])
                x4 <- melt(x3[, c("cohort", "Allele1", "Allele2")], id.vars = c("cohort"))
                x2a <- data.frame(table(x4$cohort, x4$value))
                x2b <- data.frame(table(x4$cohort))
                names(x2b)[2] <- "Count"
                suppressMessages(x3 <- join(x2a, x2b))
                x3$Freq <- x3$Freq/x3$Count
                x3$Simulation <- 0
                head(x3)
                names(x3) <- c("Cohort", "Allele", "p", "Count", "Simulation")
                x3 <- x3[, c("Cohort", "Simulation", "Count", "p", "Allele")]
                rm(x2a, x2b, x4)
        }
        return(list(observed_frequencies = x3, simulated_frequencies = x1))
}
