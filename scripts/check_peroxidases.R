res <- read.csv("results/deseq2/differential_expression_results.csv", row.names = 1)
res$gene <- sub("T[0-9]+$", "", sub("^transcript:", "", rownames(res)))
ref <- read.csv("annotation/peroxidase_genes.csv", stringsAsFactors = FALSE)
m <- merge(ref, res[, c("gene", "baseMean", "log2FoldChange", "padj")], by = "gene", all.x = TRUE)
m <- m[order(m$group, -m$log2FoldChange), ]
print(m, digits = 3, row.names = FALSE)
cat("\n== Resumo por grupo ==\n")
for (g in unique(m$group)) {
  s <- m[m$group == g, ]
  cat(sprintf("%s: %d genes na lista, %d no resultado, %d mais expressos em 96h (padj<0.05, log2FC>=1), %d mais expressos em 40h (padj<0.05, log2FC<=-1), mediana log2FC = %.2f\n",
    g, nrow(s), sum(!is.na(s$log2FoldChange)),
    sum(s$padj < 0.05 & s$log2FoldChange >= 1, na.rm = TRUE),
    sum(s$padj < 0.05 & s$log2FoldChange <= -1, na.rm = TRUE),
    median(s$log2FoldChange, na.rm = TRUE)))
}
write.csv(m, "results/deseq2/peroxidase_genes_results.csv", row.names = FALSE)
