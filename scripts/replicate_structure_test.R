library(tximport); library(DESeq2); library(ggplot2)

dirs <- list.dirs("results/quant_test", recursive = FALSE, full.names = FALSE)
files <- file.path("results/quant_test", dirs, "quant.sf"); names(files) <- dirs
txi <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)

parts <- do.call(rbind, strsplit(dirs, "_"))
coldata <- data.frame(row.names = dirs, time = parts[,1], rep = parts[,2], bc = parts[,3])

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
vsd <- vst(dds, blind = TRUE)

for (t in c("40h", "96h")) {
  v <- vsd[, vsd$time == t]
  pca <- plotPCA(v, intgroup = c("rep", "bc"), returnData = TRUE, ntop = 1000)
  p <- ggplot(pca, aes(PC1, PC2, color = bc, shape = rep)) +
    geom_point(size = 4) + theme_minimal() + ggtitle(paste("PCA dentro de", t))
  ggsave(paste0("results/deseq2/pca_test_", t, ".png"), p, width = 6, height = 5)
}

