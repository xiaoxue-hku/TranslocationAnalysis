args <- commandArgs(TRUE)

# data import

mut <- read.table("args[1]", sep='\t', header=TRUE)

# factor linear regression for normalized translocation number on replication timing and HR for regions with 1Mb bin size across the whole genome

mut_model_replication <- lm(translocation.Mb.sample ~ replication_timing, data = mut)

# summarize results

mut_model_replication_sum <- summary(mut_model_replication)

exp(cbind(OR=coef(mut_model_replication), confint(mut_model_replication, level=0.9)))

summary(mut_model_replication)$coefficients[,"Pr(>|t|)"]