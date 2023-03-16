args <- commandArgs(trailingOnly = TRUE)

covar_fname <- args[1]
varexp_cutoff <- 0.9

covars <- read.delim(covar_fname, check.names = F, header = T, row.names = 1)
covars <- covars[,colSums(covars) > 0]
# From covariates get PCs explaining up to 90% of variance
pca <- prcomp(covars, center = TRUE, scale = TRUE)
exp_var <- summary(pca)$importance[3,]
print(summary(pca)$importance)
main_pcs <- pca$x[,1:min(which(exp_var > varexp_cutoff))]
write.table(main_pcs, file = paste0( covar_fname, ".PCs_", varexp_cutoff, ".txt"), sep = "\t", quote = F, col.names = NA)
cat("Number of PCs explaining up to 90% of variance: ", ncol(main_pcs), "\n")

