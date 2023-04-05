args <- commandArgs(trailingOnly = TRUE)

library(granulator)
library(dplyr)
library(tidyverse)

bulk_expr_tpm <- as.matrix(read.delim(args[1], row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
out_basepath <- args[2]

bulk_expr_tpm <- as.matrix(read.delim("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v3/data/LLS_660Q.gene_read_counts_BIOS_and_LLD_passQC.tsv.TPM.txt.gz", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
out_basepath <-"/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v3/cell_counts/LLS_660Q"
lm22 <- as.matrix(read.delim("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/data/LM22.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
abis <- as.matrix(read.delim("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/data/ABIS_sigmatrixRNAseq.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))


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

