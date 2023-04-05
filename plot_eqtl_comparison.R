res <- data.frame()
filenames <- list.files(path="/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/eQTLs/", pattern="*txt", full.names=TRUE, recursive=FALSE)

for (fname in filenames){
  f <- read.delim(fname, check.names = F, header = T, row.names = 1)
  f$id <- paste0(f$snp, ":", f$gene)
  f$method <- gsub(".txt", "", gsub("LL_inter_res.", "", tail(unlist(strsplit(fname, "/")), n = 1)))
  f <-f[,c("id", "method", "b_geno:sex", "p_geno:sex")] 
  res <- rbind(res, f)
}
res$sign <-ifelse(res$`p_geno:sex` < 0.05, 1, 0)
colnames(res)[3] <- "beta"
colnames(res)[4] <- "P"

res_wide = res %>% spread(method, P)

res_wide = res %>%   pivot_wider(names_from = method, values_from = beta, id_cols = -c(P,sign))
pdf("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/eQTLs/beta_cmp.pdf", width = 7, height =7)
column_to_rownames(res_wide, var = "id") %>% select(-ends_with("PCs")) %>% ggpairs(diag = NULL) + theme_bw() + xlim(-7,7) + ylim(-7,7)
dev.off()

pdf("/Users/Dasha/work/UMCG/data/gender_differences/eQTLgen/v1/eQTLs/cc_vs_pcs.pdf", width = 7, height =7)
p1 <- ggplot(res_wide, aes(x = nnls.LM22, y = nnls.LM22.PCs)) + geom_point()
p2 <- ggplot(res_wide, aes(x = nnls.ABIS, y = nnls.ABIS.PCs)) + geom_point()
p3 <- ggplot(res_wide, aes(x = EPIC, y = EPIC.PCs)) + geom_point()
p4 <- ggplot(res_wide, aes(x = deconcell, y = deconcell.PCs)) + geom_point()
(p1 + p2) / (p3 + p4)
dev.off()


res_wide[res_wide$nnls.ABIS < -4,]

unlist(strsplit(res[res$P < 0.0001, "id"], ":"))



d <- read.delim("benchmark.txt", sep = "\t", as.is = T, check.names = F)
pdf("BIOS_benchmark.pdf")
ggplot(d, aes(x = method)) + geom_boxplot(aes(y = fraction)) + facet_grid(rows = vars(cohort)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(d, aes(x = method, y = fraction)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(d, aes(x = method)) + geom_boxplot(aes(y = num_replicated)) + facet_grid(rows = vars(cohort)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(d, aes(x = method, y = num_replicated)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

