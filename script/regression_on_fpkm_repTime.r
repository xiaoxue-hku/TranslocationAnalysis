args <- commandArgs(TRUE)

# data import

mut <- read.table(args[1], sep='\t', header=TRUE)

# linear regression for normalized translocation number on gene_expression and replication timing for all / silent / non-silent genes

mut_model <- lm(translocation.Mb.sample ~ gene_expression + replication_timing, data = mut)

mut_model_transcription <- lm(translocation.Mb.sample ~ gene_expression, data = mut)

mut_model_replication <- lm(translocation.Mb.sample ~ replication_timing, data = mut)

# summarize results

mut_model_sum <- summary(mut_model)

mut_model_transcription_sum <- summary(mut_model_transcription)

mut_model_replication_sum <- summary(mut_model_replication)

exp(cbind(OR=coef(mut_model_transcription), confint(mut_model_transcription, level=0.9)))

exp(cbind(OR=coef(mut_model_replication), confint(mut_model_replication, level=0.9)))

exp(cbind(OR=coef(mut_model), confint(mut_model, level=0.9)))

summary(mut_model_transcription)$coefficients[,"Pr(>|t|)"]

summary(mut_model_replication)$coefficients[,"Pr(>|t|)"]

summary(mut_model)$coefficients[,"Pr(>|t|)"]