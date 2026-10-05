library(tximport)
library(DESeq2)
library(ggplot2)
library(ggrepel)

samples <- c("40h_bc1","40h_bc2","40h_bc3","96h_bc1","96h_bc2","96h_bc3")
condition <- factor(c("h40","h40","h40","h96","h96","h96"), levels = c("h40","h96"))

files <- file.path("results/quant", samples, "quant.sf")
names(files) <- samples
txi <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
coldata <- data.frame(row.names = samples, condition = condition)

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)
saveRDS(dds, "results/deseq2/dds.rds")

# log2FC positivo = mais expresso em 96 h
res <- results(dds, contrast = c("condition", "h96", "h40"))
res <- res[order(res$padj), ]
df <- as.data.frame(res)
write.csv(df, "results/deseq2/differential_expression_results.csv")

# PCA
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "condition", ntop = 500, returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), size = 3.5, box.padding = 0.6,
                  min.segment.length = 0, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.25)) +
  scale_color_manual(values = c(h40 = "#0072B2", h96 = "#D55E00"),
                     labels = c(h40 = "40 h", h96 = "96 h")) +
  labs(x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance"),
       color = "Time point", title = "PCA of variance-stabilized counts (top 500 genes)") +
  theme_minimal(base_size = 12)
ggsave("results/deseq2/pca_plot.png", p_pca, width = 7, height = 5.5, dpi = 200)

# Volcano (padj abaixo de 1e-300 limitado em 300 no eixo)
d <- df[!is.na(df$padj), ]
d$nlp <- -log10(pmax(d$padj, 1e-300))
d$class <- "Not significant"
d$class[d$padj < 0.05 & d$log2FoldChange >= 1] <- "Higher at 96 h"
d$class[d$padj < 0.05 & d$log2FoldChange <= -1] <- "Higher at 40 h"
p_vol <- ggplot(d, aes(log2FoldChange, nlp, color = class)) +
  geom_point(size = 0.9, alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  scale_color_manual(values = c("Higher at 96 h" = "#D55E00", "Higher at 40 h" = "#0072B2",
                                "Not significant" = "grey75")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = "log2 fold change (96 h vs 40 h)", y = "-log10 adjusted p-value (capped at 300)",
       color = NULL, title = "Differential expression, 96 h vs 40 h") +
  theme_minimal(base_size = 12)
ggsave("results/deseq2/volcano_plot.png", p_vol, width = 7.5, height = 5.5, dpi = 200)

# Numeros usados no README
sig <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) >= 1
cat("Genes testados:", sum(!is.na(df$padj)), "\n")
cat("Up em 96h:", sum(sig & df$log2FoldChange > 0), "\n")
cat("Down em 96h:", sum(sig & df$log2FoldChange < 0), "\n")
cat("Mediana da dispersao por gene:", signif(median(mcols(dds)$dispGeneEst, na.rm = TRUE), 3), "\n")
res_strict <- results(dds, contrast = c("condition", "h96", "h40"), lfcThreshold = 1, altHypothesis = "greaterAbs")
cat("Teste contra |log2FC| > 1, padj < 0.05:", sum(res_strict$padj < 0.05, na.rm = TRUE), "\n")
cat("Up em 96h, log2FC >= 2, padj < 0.05:", sum(df$log2FoldChange >= 2 & df$padj < 0.05, na.rm = TRUE), "\n")
cat("   o mesmo com baseMean >= 100:", sum(df$log2FoldChange >= 2 & df$padj < 0.05 & df$baseMean >= 100, na.rm = TRUE), "\n")
