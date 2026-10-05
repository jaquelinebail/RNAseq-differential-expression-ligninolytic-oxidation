library(ggplot2)
m <- read.csv("results/deseq2/peroxidase_genes_results.csv")
m$group <- factor(m$group, levels = c("classII_lignin_mn", "heme_other", "other_peroxidases"),
                  labels = c("Class II peroxidases (lignin/Mn)", "Heme peroxidases (unspecified)", "Other peroxidases"))
m <- m[order(m$group, m$log2FoldChange), ]
m$gene <- factor(m$gene, levels = m$gene)
m$class <- "Not significant or < 2-fold"
m$class[!is.na(m$padj) & m$padj < 0.05 & m$log2FoldChange >= 1] <- "Higher at 96 h"
m$class[!is.na(m$padj) & m$padj < 0.05 & m$log2FoldChange <= -1] <- "Higher at 40 h"
p <- ggplot(m, aes(log2FoldChange, gene, color = class, size = log10(baseMean))) +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point() +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("Higher at 96 h" = "#D55E00", "Higher at 40 h" = "#0072B2",
                                "Not significant or < 2-fold" = "grey65")) +
  labs(x = "log2 fold change (96 h vs 40 h)", y = NULL, color = NULL,
       size = "log10 mean\nnormalized counts", title = "Peroxidase genes, 96 h vs 40 h") +
  theme_minimal(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))
ggsave("results/deseq2/peroxidase_validation.png", p, width = 7.5, height = 7.5, dpi = 200)
