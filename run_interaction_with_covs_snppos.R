args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(patchwork)

# setwd("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/")
# 
expr_fname <- "data/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt.table.gz.TMM.txt"
geno_fname <- "data/500FG.eqtl_genotypes.dosages.txt"
gte_fname <- "data/500FG_samples_3cols.txt"
covar_fname <- "cell_counts/tmp/LLD.TPM.EPIC+sex.txt"
eqtls_fname <- "data/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.top_hits.txt"

# 
# expr_fname <- "data/LLD_gene_read_counts.TMM.ProbesCentered.SamplesZTransformed.txt.gz"
# geno_fname <- "data/LL.eqtl_genotypes.filtered.top.dosages.txt.gz"
# gte_fname <- "data/LLD_samples_3cols.txt"
# covar_fname <- "cell_counts/tmp/LLD.TPM.nnls.LM22+sex.txt"
# out_fname <- "eQTLs/LL_inter_res.nnls.LM22.txt"
#   
expr_fname <- args[1]
geno_fname <- args[2]
gte_fname <- args[3]
covar_fname <- args[4]
eqtls_fname <- args[5]
out_fname <- args[6]

# Read in the data
expr_noadj <- as.data.frame(t(read.delim(gzfile(expr_fname), check.names = F, header = T, row.names = 1)))
geno <- as.data.frame(t(read.delim(gzfile(geno_fname), check.names = F, header = T, row.names = 1)))
gte <- read.delim(gte_fname, check.names = F, header = F)
sex <- read.delim(covar_fname, check.names = F, header = T, row.names = 1)
sex <- na.omit(sex)
# rename expression sample ids to genotype sample ids
geno_m <- geno[match(gte[,2], row.names(geno), nomatch = 0),]
gte_m <- gte[match(row.names(geno_m), gte[,2], nomatch = 0),]
row.names(geno_m) <- gte_m[,3]
cat("Coverted genotype sample ids to expression sample ids.\nN samples before: ", nrow(geno), "\nN samples after: ", nrow(geno_m), "\n")

# get ids for samples shared between expr_noadjession,  genotypes and sex
ids <- intersect(intersect(row.names(expr_noadj), row.names(geno_m)), row.names(sex))
cat ("Number of overlapping samples between genotype and expr_noadjession data: ", length(ids), "\n")
expr_noadj2 <- expr_noadj[ids,]
geno2 <- geno_m[ids, ]
sex2 <- sex[ids,]

# Read eQTLs to test:
eqtls_full <- read.delim(eqtls_fname, check.names = F, header = T)
eqtls_full$chrpos <- paste0(eqtls_full$SNPChr, ":", eqtls_full$SNPPos)
eqtls <- eqtls_full[eqtls_full$chrpos %in% colnames(geno),]
cat("Read ", nrow(eqtls), " eQTLs")

# Run analysis
res <- data.frame(matrix(nrow = nrow(eqtls), ncol = 9))
colnames(res) <- c("gene", "snp", "b_geno", "b_sex", "b_geno:sex", "p_geno", "p_sex", "p_geno:sex", "Wilcox_p_sex-af")
cnt <- 1
for (e in 1:nrow(eqtls)){
  gene <- eqtls[e, "Gene"]
  snp <- eqtls[e, "chrpos"]
  
  if (gene %in% colnames(expr_noadj2) & snp %in% colnames(geno2)){
    m <- as.data.frame(cbind(expr_noadj2[, gene], geno2[, snp], sex2[ids,]))
    colnames(m) <- c("gene", "dosage", colnames(sex))
    colnames(m) <- gsub("gender_F1M2", "sex", colnames(m))
    row.names(m) <- ids
    m$sex <- as.factor(m$sex)
    lm_fit <- lm(gene ~ . + sex*(.) + dosage*(.) , data = m)
    #summary(lm_fit)
    af_sex <- wilcox.test(dosage ~ sex, data = m)$p.value
    coef <- summary(lm_fit)$coefficients
    if ("dosage:sex2" %in% row.names(coef)){
      res[cnt,] <- c(gene, snp, coef["dosage",1], coef["sex2",1], coef["dosage:sex2",1], coef["dosage",4], coef["sex2",4], coef["dosage:sex2",4], af_sex)
    }
    cnt <- cnt + 1
  }
  if (cnt %% 10000 == 0){
    cat("processed ", cnt, " eQTLs\n")
  }
}

res[,3:ncol(res)] = apply(res[,3:ncol(res)], 2, function(x) as.numeric(as.character(x)))

res_srt <- res[order(res[,"p_geno:sex"]),]
res_srt[,"p_geno:sex_BH"] <- p.adjust(res_srt[,"p_geno:sex"], method = "BH")
write.table(res_srt, file = out_fname, sep = "\t", quote = F, col.names = NA)

# 
# plot_interaction <- function(snp, gene, inter_p = NA){
#   #snp = "rs885070"
#   #gene = "ENSG00000107317"
#   #snp = "rs12151742"
#   #gene = "ENSG00000115523"
#   #gene="ENSG00000100116"
#   #snp = "rs6000897"
#   gene_name <- eqtls[eqtls$SNP == snp & eqtls$Gene == gene,"GeneSymbol"]
#   eqtl_p <- eqtls[eqtls$SNP == snp & eqtls$Gene == gene,"Pvalue"]
#   
#   m <- as.data.frame(cbind(expr_noadj2[, gene], geno2[, snp], sex[ids,]))
#   colnames(m) <- c("gene", "dosage", colnames(sex))
#   colnames(m) <- gsub("gender_F1M2", "sex", colnames(m))
#   row.names(m) <- ids
#   
#   m$sex <- as.factor(m$sex)
#   
#   
#   m <- na.omit(m)
#   
#   counts <- m %>% 
#     group_by(interaction(sex, genotype)) %>% 
#     tally()
#   
#   counts$yloc <- max(m$gene) - 0.005
#   counts$genotype <- c(0,0,1,1,2,2)
#   counts$sex <- c(1,2,1,2,1,2)
#   
#   lm_covars <- lm(gene ~ . + sex*(.) - dosage - dosage:sex, data = m)
#   m$residuals <- resid(lm_covars)
#   m$genotype <- round(m$dosage)
#   
#   p <- ggplot(data = m , aes(y = residuals, x = genotype))  + 
#     geom_boxplot(aes(group = paste(sex, genotype), color = sex)) +
#     geom_smooth(aes(x = as.integer(genotype), y = gene, color = sex), method = 'lm', formula = y ~ x) +
#     ylab(paste0(gene_name, " adjusted expression")) + xlab(snp) +
#     scale_color_manual(values = c("#F65C00", "#006DB3"), labels=c("female", "male")) +
#     theme_bw() +
#     theme(plot.title = element_text(size=12)) +
#     ggtitle(paste0("eQTL P = ", format(eqtl_p, digits = 3), "\nGxS interaction P = ", format(inter_p, digits = 3))) +
#     geom_text(data = counts, aes(x = genotype,y =  yloc, label = n, color = as.factor(sex)), position = position_dodge(width = .75)) 
#   return (p)
# }

# 
# library(patchwork)
# n4plot <- 6
# plots <- list()
# cnt <- 1
# for (i in 1:n4plot){
#   plots[[cnt]] <- plot_interaction(res_srt[i,"snp"], res_srt[i,"gene"], res_srt[i,"p_geno:sex"])
#   cnt <- cnt + 1
# }
# pdf(paste0(out_fname, ".top_interactions.pdf"), width = 10, height = 15, useDingbats = F)
# (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) / (plots[[5]] + plots[[6]])
# dev.off()
# 
# 

