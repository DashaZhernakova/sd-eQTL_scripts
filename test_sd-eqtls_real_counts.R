args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(patchwork)

setwd("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v3/")
 
cohort="LLS_660Q"
expr_fname <- paste0("data/", cohort, ".gene_read_counts_BIOS_and_LLD_passQC.TMM.ProbesCentered.SamplesZTransformed.txt.gz")
geno_fname <- paste0("data/", cohort, ".eqtl_genotypes.filtered.dosages.txt.gz")
gte_fname <- paste0("data/", cohort, "_samples.txt")
covar_fname <- "data/BIOS_cell_percentages+sex.txt"


#out_fname <- paste0(eqtls_fname, ".tested_LL_cellcounts.txt")
 
# expr_fname <- args[1]
# geno_fname <- args[2]
# gte_fname <- args[3]
# covar_fname <- args[4]
# eqtls_fname <- args[5]
# out_fname <- args[6]

# Read in the data
expr_noadj <- as.data.frame(t(read.delim(gzfile(expr_fname), check.names = F, header = T, row.names = 1)))
geno <- as.data.frame(t(read.delim(gzfile(geno_fname), check.names = F, header = T, row.names = 1)))
gte <- read.delim(gte_fname, check.names = F, header = F)
sex <- read.delim(covar_fname, check.names = F, header = T, row.names = 1)

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

filenames <- list.files(path="/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v3/eQTLs/", pattern=glob2rx(paste0(cohort, "*txt")), full.names=TRUE)

for (eqtls_fname in filenames){
  
# Read eQTLs to test:
eqtls_full <- read.delim(eqtls_fname, check.names = F, header = T)
#cat("Read ", nrow(eqtls), " eQTLs")
#eqtls <- eqtls_full
eqtls <- eqtls_full[eqtls_full[,"p_geno:sex"] < 0.001,]
eqtls <- na.omit(eqtls)
res <- data.frame(gene = character(), snp = character(), decon_inter_p = double(), decon_inter_b = double(), counts_inter_p = double(), counts_inter_b = double(), stringsAsFactors = F)
cnt = 1
for (e in 1:nrow(eqtls)){

  gene <- eqtls[e, "gene"]
  snp <- eqtls[e, "snp"]
  #cat(gene, snp, "\n")
  m <- as.data.frame(cbind(expr_noadj2[, gene], geno2[, snp], sex2[ids,]))
  colnames(m) <- c("gene", "dosage", colnames(sex))
  colnames(m) <- gsub("gender_F1M2", "sex", colnames(m))
  row.names(m) <- ids
  m <- Filter(function(x)!all(is.na(x)), m)
  m$sex <- as.factor(m$sex)
  lm_fit <- lm(gene ~ . + sex*(.) + dosage*(.) , data = m)
  summary(lm_fit)
  coef <- summary(lm_fit)$coefficients
  
  res[cnt, "gene"] <- gene
  res[cnt, "gene"] <- snp
  res[cnt, "decon_inter_p"] <- eqtls[e, "p_geno:sex"]
  res[cnt, "decon_inter_b"]  <- eqtls[e, "b_geno:sex"]
  res[cnt, "counts_inter_p"] <- coef["dosage:sex2",4]
  res[cnt, "counts_inter_b"] <- coef["dosage:sex2",1]
  cnt <- cnt + 1
}
cat(eqtls_fname, nrow(eqtls), nrow(res[res$counts_inter_p < 0.05,]), nrow(res[res$counts_inter_p < 0.05,])/nrow(eqtls), formatC(cor(res$decon_inter_b, res$counts_inter_b), digits = 3), "\n")

}

plot(res$decon_inter_p, res$counts_inter_p)



