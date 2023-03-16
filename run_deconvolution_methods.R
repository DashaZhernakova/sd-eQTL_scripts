library(granulator)
library(dplyr)
library(tidyverse)

setwd("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/")
bulk_expr_tpm <- as.matrix(read.delim("data/gene.counts-LLDBIOSSamples.LLD_subset.txt.gz.TPM.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
cell_counts1 <- na.omit(as.matrix(read.delim("data/LLD_cell_counts.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t")))
cell_counts2 <- read.delim("data/LLD_scrnaseq_cell_counts_perc.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t")
lm22 <- as.matrix(read.delim("data/LM22.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))
abis <- as.matrix(read.delim("data/ABIS_sigmatrixRNAseq.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t"))

sigList = list(as.matrix(lm22), as.matrix(abis))
abis <- na.omit(abis)
#plot_similarity(sigMatrix=sigList)

decon_abis <- deconvolute(m = bulk_expr_tpm, sigMatrix = abis, methods = "nnls")
decon_lm22 <- deconvolute(m = bulk_expr_tpm, sigMatrix = lm22, methods = "nnls")

decon_abis <- decon_abis[,colSums(decon_abis) > 0]
decon_lm22 <- decon_lm22[,colSums(decon_lm22) > 0]
write.table(decon_abis$proportions$nnls_sig1, file = "cell_counts/LLD.TPM.nnls.ABIS.txt", sep = "\t", quote = F, col.names = NA)
write.table(decon_lm22$proportions$nnls_sig1, file = "cell_counts/LLD.TPM.nnls.LM22.txt", sep = "\t", quote = F, col.names = NA)

#
# EPIC
#
library(EPIC)
out <- EPIC(bulk = bulk_expr_tpm)
write.table(out$cellFraction, file = "cell_counts/LLD.TPM.EPIC.txt", sep = "\t", quote = F, col.names = NA)


nnls_abis <- read.delim("cell_counts/LLD.TPM.nnls.ABIS.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t")
nnls_lm22 <- read.delim("cell_counts/LLD.TPM.nnls.LM22.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t")
epic <- out$cellFraction

nnls_abis_fmt <- rename_cell_types_abis(nnls_abis)
nnls_lm22_fmt <- rename_cell_types_lm22(nnls_lm22)
epic_fmt <- rename_cell_types_epic(epic)

res1 <- benchmark_deconvolution(nnls_abis_fmt, cell_counts2, "nnls_abis_vs_sc")
res2 <- benchmark_deconvolution(nnls_abis_fmt, cell_counts1, "nnls_abis_vs_facs")
res3 <- benchmark_deconvolution(nnls_lm22_fmt, cell_counts2, "nnls_lm22_vs_sc")
res4 <- benchmark_deconvolution(nnls_lm22_fmt, cell_counts1, "nnls_lm22_vs_facs")
res5 <- benchmark_deconvolution(epic_fmt, cell_counts2, "epic_vs_sc")
#res6 <- benchmark_deconvolution(epic_fmt, cell_counts1, "epic_vs_facs")

res_all <- list(res1,res2,res3, res4, res5) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="cell_type"), .)

res_all_long <- res_all %>% pivot_longer(!cell_type, names_to = "method", values_to = "pcc")
res_all_long$truth <- "scRNA-seq"
res_all_long[endsWith(res_all_long$method, "facs"), "truth"] <- "FACS"
res_all_long$method_short <- gsub("_vs.*", "",res_all_long$method) 
pdf("plots/pcc_nnls_epic.pdf", height = 5, width = 5)
ggplot(res_all_long, aes(x = method_short, y = pcc, color = truth))  + 
  geom_point() + 
  facet_wrap(~cell_type) + 
  geom_segment( aes(x=method_short, xend=method_short, y=0, yend=pcc)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
benchmark_deconvolution <- function(decon, truth, plot_path){
  #decon <- decon_lm22$proportions$nnls_sig1
  #truth <- cell_counts2
  #plot_path <- "nnls_abis_vs_sc.pdf"
  ids <- intersect(row.names(decon), row.names(truth))
  cell_types<- intersect(colnames(decon), colnames(truth))

  pdf(paste0("plots/", plot_path, "pdf"), width = 9, height = 6, useDingbats = F)
  par(mfrow = c(2, 3))
  res <- data.frame(matrix(nrow = length(cell_types), ncol = 1))
  
  row.names(res) <- cell_types
  colnames(res) <- plot_path
  for (c in cell_types){
    correl <- cor(decon[ids, c], truth[ids,c])
    res[c,1] <- correl
    plot(decon[ids, c], truth[ids,c], pch = 16, main = paste0(c,  "\nr = ",formatC(correl, digits = 3)) , xlab = "deconvoluted", ylab = "true")
  }
  dev.off()
  return(rownames_to_column(res, var="cell_type"))
}



rename_cell_types_epic <- function(epic_res){
  colnames(epic_res) <- gsub("Bcells", "B", colnames(epic_res))
  colnames(epic_res) <- gsub("CD4_Tcells", "CD4T", colnames(epic_res))
  colnames(epic_res) <- gsub("CD8_Tcells", "T cells CD8", colnames(epic_res))
  colnames(epic_res) <- gsub("NKcells", "NK", colnames(epic_res))
  return(epic_res)
}

rename_cell_types_lm22 <- function(lm22_res){
  lm22_res$B <- lm22_res[,"B cells memory"] + lm22_res[,"B cells naive"]
  lm22_res$DC <- lm22_res[,"Dendritic cells activated"] + lm22_res[,"Dendritic cells resting"]
  lm22_res$NK <- lm22_res[,"NK cells activated"] + lm22_res[,"NK cells resting"]
  lm22_res$CD4T <- lm22_res[,"T cells CD4 memory activated"] + lm22_res[,"T cells CD4 memory resting"] + lm22_res[,"T cells CD4 naive"]
  return(lm22_res)
  }
  
  
rename_cell_types_abis <- function(abis_res){
  abis_res$B <- abis_res[,"B Memory"] + abis_res[,"B Naive"]
  abis_res$DC <- abis_res[,"mDCs"] + abis_res[,"pDCs"]
  abis_res$Monocytes <- abis_res[,"Monocytes C"] + abis_res[,"Monocytes NC+I"]
  abis_res$CD4T <- abis_res[,"T CD4 Memory"] + abis_res[,"T CD4 Naive"]
  abis_res[,"T cells CD8"] <- abis_res[,"T CD8 Memory"] + abis_res[,"T CD8 Naive"]
  colnames(abis_res) <- gsub("Neutrophils LD", "Neutrophils", colnames(abis_res))
  colnames(abis_res) <- gsub("Basophils LD", "Basophils", colnames(abis_res))
  return(abis_res)
}


#cell_counts_full <- read.delim("data/LLD_scrnaseq_cell_counts_full.txt", row.names = 1, header = T, as.is = T, check.names = F, sep = "\t")
#cell_counts_perc <- t(apply(cell_counts_full,1, function(x) x/sum(x)))
#write.table(cell_counts_perc, file = "data/LLD_scrnaseq_cell_counts_perc.txt", sep = "\t", quote = F, col.names = NA)
