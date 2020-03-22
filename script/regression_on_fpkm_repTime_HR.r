args <- commandArgs(TRUE)

# data import

mut <- read.table("args[1]", sep='\t', header=TRUE)

# factor linear regression for normalized translocation number on gene_expression, replication timing and HR for all / non-silent genes

mut_model <- lm(translocation.Mb.sample ~ gene_expression + replication_timing + factor(HR), data = mut)

mut_model_transcription <- lm(translocation.Mb.sample ~ gene_expression, data = mut)

mut_model_replication <- lm(translocation.Mb.sample ~ replication_timing, data = mut)

mut_model_hr <- lm(translocation.Mb.sample ~ factor(HR), data = mut)

# summarize results

mut_model_sum <- summary(mut_model)

mut_model_transcription_sum <- summary(mut_model_transcription)

mut_model_replication_sum <- summary(mut_model_replication)

mut_model_hr_sum <- summary(mut_model_hr)

exp(cbind(OR=coef(mut_model_transcription), confint(mut_model_transcription, level=0.9)))

exp(cbind(OR=coef(mut_model_replication), confint(mut_model_replication, level=0.9)))

exp(cbind(OR=coef(mut_model_hr), confint(mut_model_hr, level=0.9)))

exp(cbind(OR=coef(mut_model), confint(mut_model, level=0.9)))

summary(mut_model_transcription)$coefficients[,"Pr(>|t|)"]

summary(mut_model_replication)$coefficients[,"Pr(>|t|)"]

summary(mut_model_hr)$coefficients[,"Pr(>|t|)"]

summary(mut_model)$coefficients[,"Pr(>|t|)"]
