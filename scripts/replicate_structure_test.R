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
  pv <- round(100 * attr(pca, "percentVar"))
  p <- ggplot(pca, aes(PC1, PC2, color = bc, shape = rep)) +
    geom_point(size = 4, alpha = 0.85) +
    scale_color_manual(values = c(bc1 = "#0072B2", bc2 = "#D55E00", bc3 = "#009E73")) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance"),
         color = "GEO barcode label", shape = "GEO rep label",
         title = paste0("Run-level PCA at ", sub("h", " h", t), " (5 M reads per run)")) +
    theme_minimal(base_size = 12)
  ggsave(paste0("results/deseq2/pca_test_", t, ".png"), p, width = 6.5, height = 5, dpi = 150)
}
