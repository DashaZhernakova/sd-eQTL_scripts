args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(patchwork)

setwd("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v2/")


expr_fname <- "data/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt.table.gz.TMM.txt"
geno_fname <- "data/500FG.eqtl_genotypes.dosages.txt"
gte_fname <- "data/500FG_samples_3cols.txt"
covar_fname <- "data/500FG_sex+cellcounts.txt"
eqtls_fname <- "eQTLs/LL_inter_res.nnls.LM22.txt"
#out_fname <- paste0(eqtls_fname, ".tested_LL_cellcounts.txt")

# expr_fname <- args[1]
# geno_fname <- args[2]
# gtxx?? xe_fname <- args[3]
# covar_fname <- args[4]
# eqtls_fname <- args[5]
# out_fname <- args[6]

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

filenames <- list.files(path="/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v2/eQTLs/", pattern="*txt", full.names=TRUE, recursive=FALSE)

for (eqtls_fname in filenames){
  
  # Read eQTLs to test:
  eqtls_full <- read.delim(eqtls_fname, check.names = F, header = T)
  #cat("Read ", nrow(eqtls), " eQTLs")
  #eqtls <- eqtls_full
  eqtls <- na.omit(eqtls_full[eqtls_full[,"p_geno:sex"] < 0.001,])

  res <- data.frame(gene = character(), snp = character(), decon_inter_p = double(), decon_inter_b = double(), counts_inter_p = double(), counts_inter_b = double(), stringsAsFactors = F)
  cnt = 1
  for (e in 1:nrow(eqtls)){
    gene <- eqtls[e, "gene"]
    snp <- eqtls[e, "snp"]
    #cat(gene, snp, "\n")
    max_p = 0
    for (ct in colnames(sex2)[2:ncol(sex2)]){
      #cat (e, gene, snp, ct, "\n")
      m <- as.data.frame(cbind(expr_noadj2[, gene], geno2[, snp], sex2[ids,c("gender_F1M2", ct)]))
      colnames(m) <- c("gene", "dosage", "gender_F1M2", ct)
      colnames(m) <- gsub("gender_F1M2", "sex", colnames(m))
      row.names(m) <- ids
      m$sex <- as.factor(m$sex)
      
      lm_fit <- lm(gene ~ . + sex*(.) + dosage*(.) + sex*dosage, data = m)
      summary(lm_fit)
      coef <- summary(lm_fit)$coefficients
      max_p <- max(max_p, coef["dosage:sex2",4])
      #cat(ct, coef["dosage:sex2",4], "\n")
    }
    res[cnt, "gene"] <- gene
    res[cnt, "gene"] <- snp
    res[cnt, "decon_inter_p"] <- eqtls[e, "p_geno:sex"]
    res[cnt, "decon_inter_b"]  <- eqtls[e, "b_geno:sex"]
    res[cnt, "counts_inter_p"] <- max_p
    
    cnt <- cnt + 1
  }
  cat(eqtls_fname, nrow(eqtls), nrow(res[res$counts_inter_p < 0.05,]), nrow(res[res$counts_inter_p < 0.05,])/nrow(eqtls),  "\n")
  
}

plot(res$decon_inter_p, res$counts_inter_p)



