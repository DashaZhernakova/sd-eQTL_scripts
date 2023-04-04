args <- commandArgs(trailingOnly = TRUE)

library(granulator)
library(dplyr)
library(tidyverse)

bulk_expr_tpm <- as.matrix(read.delim(args[1], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
lm22 <- as.matrix(read.delim(args[2], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
abis <- as.matrix(read.delim(args[3], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
out_basepath <- args[4]

sigList = list(as.matrix(lm22), as.matrix(abis))
abis <- na.omit(abis)
#plot_similarity(sigMatrix=sigList)

decon_abis <- deconvolute(m = bulk_expr_tpm, sigMatrix = abis, methods = "nnls")
decon_lm22 <- deconvolute(m = bulk_expr_tpm, sigMatrix = lm22, methods = "nnls")

decon_dtangle_abis <- deconvolute(m = bulk_expr_tpm, sigMatrix = abis, methods = "dtangle")
decon_dtangle_lm22 <- deconvolute(m = bulk_expr_tpm, sigMatrix = lm22, methods = "dtangle")

decon_abis_res <- decon_abis$proportions$nnls_sig1
decon_abis_res <- decon_abis_res[,colSums(decon_abis_res) > 0]
decon_lm22_res <- decon_lm22$proportions$nnls_sig1
decon_lm22_res <- decon_lm22_res[,colSums(decon_lm22_res) > 0]

decon_dtangle_abis_res <- decon_dtangle_abis$proportions$dtangle_sig1
decon_dtangle_abis_res <- decon_dtangle_abis_res[,colSums(decon_dtangle_abis_res) > 0]

decon_dtangle_lm22_res <- decon_dtangle_lm22$proportions$dtangle_sig1
decon_dtangle_lm22_res <- decon_dtangle_lm22_res[,colSums(decon_dtangle_lm22_res) > 0]

write.table(decon_abis_res, file = paste0(out_basepath, ".TPM.nnls.ABIS.txt"), sep = "\t", quote = F, col.names = NA)
write.table(decon_lm22_res, file = paste0(out_basepath, ".TPM.nnls.LM22.txt"), sep = "\t", quote = F, col.names = NA)

write.table(decon_dtangle_abis_res, file = paste0(out_basepath, ".TPM.dtangle.ABIS.txt"), sep = "\t", quote = F, col.names = NA)
write.table(decon_dtangle_lm22_res, file = paste0(out_basepath, ".TPM.dtangle.LM22.txt"), sep = "\t", quote = F, col.names = NA)

#
# EPIC
#
library(EPIC)
out <- EPIC(bulk = bulk_expr_tpm)
write.table(out$cellFraction, file = paste0(out_basepath, ".TPM.EPIC.txt"), sep = "\t", quote = F, col.names = NA)

