library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ── Load data ──────────────────────────────────────────────────────────────────

metadata <- read.csv("F:/CRC_Tumour_WGS/SraRunTable.csv")
counts_path <- "F:/CRC_Tumour_WGS/counts"

# WGS (metagenomic DNA) samples only
meta <- metadata[metadata$Assay.Type == "WGS", ]
meta$host_phenotype <- factor(meta$host_phenotype, levels = c("healthy", "cancer"))

# ── Load and aggregate counts to phage (contig) level ─────────────────────────

read_counts <- function(f) {
  df <- read.table(f, header = TRUE, skip = 1, sep = "\t")
  sample_id <- sub("\\.cds\\.counts\\.txt$", "", basename(f))
  agg <- aggregate(df[[ncol(df)]], by = list(phage = df$Chr), FUN = sum)
  colnames(agg)[2] <- sample_id
  agg
}

files <- list.files(counts_path, pattern = "\\.cds\\.counts\\.txt$", full.names = TRUE)
count_list <- lapply(files, read_counts)
counts <- Reduce(function(a, b) merge(a, b, by = "phage"), count_list)
rownames(counts) <- counts$phage
counts$phage <- NULL

# Keep WGS samples only and align order with metadata
counts <- counts[, colnames(counts) %in% meta$Run]
meta <- meta[match(colnames(counts), meta$Run), ]

# ── DESeq2 ────────────────────────────────────────────────────────────────────

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ host_phenotype
)
dds <- dds[rowSums(counts(dds)) > 0, ]  # drop zero-count phages
dds <- DESeq(dds)

res <- results(dds, contrast = c("host_phenotype", "cancer", "healthy"))
res_df <- as.data.frame(res)
res_df$phage <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# Save results table
write.csv(res_df, "F:/CRC_Tumour_WGS/differential_abundance.csv", row.names = FALSE)

# ── Volcano plot ───────────────────────────────────────────────────────────────

res_df$sig <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 2, alpha = 0.7) +
  geom_label_repel(
    data = subset(res_df, sig),
    aes(label = phage), size = 3, max.overlaps = 20
  ) +
  scale_color_manual(values = c("grey70", "firebrick"),
                     labels = c("NS", "FDR < 0.05 & |LFC| > 1")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  theme_bw() +
  labs(
    title = "Phage differential abundance: cancer vs healthy",
    x = "log2 Fold Change (cancer / healthy)",
    y = "-log10(adjusted p-value)",
    color = NULL
  )

ggsave("F:/CRC_Tumour_WGS/volcano_plot.pdf", width = 8, height = 6)

# ── Normalised count boxplots for significant phages ──────────────────────────

sig_phages <- res_df$phage[res_df$sig & !is.na(res_df$sig)]

if (length(sig_phages) > 0) {
  norm_counts <- counts(dds, normalized = TRUE)
  plot_data <- as.data.frame(norm_counts[sig_phages, , drop = FALSE])
  plot_data$phage <- rownames(plot_data)
  plot_long <- tidyr::pivot_longer(plot_data, -phage, names_to = "Run", values_to = "count")
  plot_long <- merge(plot_long, meta[, c("Run", "host_phenotype")], by = "Run")

  ggplot(plot_long, aes(x = host_phenotype, y = count + 1, fill = host_phenotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    facet_wrap(~ phage, scales = "free_y") +
    scale_y_log10() +
    scale_fill_manual(values = c("healthy" = "steelblue", "cancer" = "firebrick")) +
    theme_bw() +
    labs(
      title = "Normalised counts — significant phages",
      x = NULL, y = "Normalised count (log10 + 1)", fill = NULL
    )

  ggsave("F:/CRC_Tumour_WGS/significant_phages_boxplots.pdf", width = 10, height = 8)
}
