library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

pdf(NULL)  # suppress Rplots.pdf

# ── Load data ──────────────────────────────────────────────────────────────────

metadata <- read.csv("F:/CRC_Tumour_WGS/SraRunTable.csv")
counts_path <- "F:/CRC_Tumour_WGS/counts"

# WGS (metagenomic DNA) samples only
meta <- metadata[metadata$Assay.Type == "WGS", ]
meta$host_phenotype <- factor(meta$host_phenotype, levels = c("healthy", "cancer"))

# ── Build contig -> phage name mapping (Bakta renames in input order) ──────────

fasta_names <- sub("^>", "", grep("^>", readLines("F:/CRC_Tumour_WGS/phages/combined_revised_phages.fasta"), value = TRUE))
fna_names   <- sub(" .*", "", sub("^>", "", grep("^>", readLines("F:/CRC_Tumour_WGS/phages/combined_revised_phages.fna"), value = TRUE)))
name_map <- setNames(fasta_names, fna_names)  # contig_N -> phage name

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
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("host_phenotype", "cancer", "healthy"))
res_df <- as.data.frame(res)
res_df$contig <- rownames(res_df)
res_df$phage  <- name_map[res_df$contig]
res_df <- res_df[order(res_df$pvalue, na.last = TRUE), ]

# Save results table
write.csv(res_df, "F:/CRC_Tumour_WGS/differential_abundance.csv", row.names = FALSE)

# ── Volcano plot (nominal p < 0.05 labelled) ───────────────────────────────────

res_df$nominal_sig <- !is.na(res_df$pvalue) & res_df$pvalue < 0.05

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = nominal_sig), size = 2, alpha = 0.7) +
  geom_label_repel(
    data = subset(res_df, nominal_sig),
    aes(label = phage), size = 3, max.overlaps = 20
  ) +
  scale_color_manual(values = c("grey70", "firebrick"),
                     labels = c("NS", "Nominal p < 0.05")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  theme_bw() +
  labs(
    title = "Phage differential abundance: cancer vs healthy",
    subtitle = "Labelled: nominal p < 0.05 (n=10 per group; FDR correction too stringent)",
    x = "log2 Fold Change (cancer / healthy)",
    y = "-log10(p-value)",
    color = NULL
  )

ggsave("F:/CRC_Tumour_WGS/volcano_plot.pdf", width = 8, height = 6)

# ── Assign species groups ──────────────────────────────────────────────────────

assign_group <- function(name) {
  if (grepl("cryptic", name))          return("Cryptic phages")
  if (grepl("_FU",     name))          return("FU phages")
  if (grepl("_ODE",    name))          return("ODE phages")
  if (grepl("_TAA",    name))          return("TAA phages")
  sub("Bacteroides_phage_", "", name)  # singletons: use short name
}

res_df$group <- sapply(res_df$phage, assign_group)

# ── Normalised counts long format ──────────────────────────────────────────────

norm_counts <- counts(dds, normalized = TRUE)
all_contigs <- res_df$contig[!is.na(res_df$pvalue)]

plot_data <- as.data.frame(norm_counts[all_contigs, , drop = FALSE])
plot_data$contig <- rownames(plot_data)
plot_long <- tidyr::pivot_longer(plot_data, -contig, names_to = "Run", values_to = "count")
plot_long <- merge(plot_long, meta[, c("Run", "host_phenotype")], by = "Run")
plot_long$phage <- name_map[plot_long$contig]
plot_long$group <- sapply(plot_long$phage, assign_group)

# ── Main PDF: one boxplot per species group (counts summed within group) ───────

group_long <- plot_long |>
  dplyr::group_by(Run, host_phenotype, group) |>
  dplyr::summarise(count = sum(count), .groups = "drop")

ggplot(group_long, aes(x = host_phenotype, y = count + 1, fill = host_phenotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ group, scales = "free_y", ncol = 4) +
  scale_y_log10() +
  scale_fill_manual(values = c("healthy" = "steelblue", "cancer" = "firebrick")) +
  theme_bw() +
  theme(strip.text = element_text(size = 9)) +
  labs(
    title = "Normalised counts by phage group",
    x = NULL, y = "Normalised count (log10 + 1)", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/phage_groups_boxplots.pdf", width = 16, height = 12)

# ── Separate PDF: cryptic phages individually ──────────────────────────────────

cryptic_long <- plot_long[plot_long$group == "Cryptic phages", ]
cryptic_long$short_name <- sub("Bacteroides_cryptic_phage_", "", cryptic_long$phage)

ggplot(cryptic_long, aes(x = host_phenotype, y = count + 1, fill = host_phenotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ short_name, scales = "free_y", ncol = 4) +
  scale_y_log10() +
  scale_fill_manual(values = c("healthy" = "steelblue", "cancer" = "firebrick")) +
  theme_bw() +
  theme(strip.text = element_text(size = 10)) +
  labs(
    title = "Normalised counts — Bacteroides cryptic phages",
    x = NULL, y = "Normalised count (log10 + 1)", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/cryptic_phages_boxplots.pdf", width = 14, height = 10)

# ── Separate PDF: FU phages individually ──────────────────────────────────────

fu_long <- plot_long[plot_long$group == "FU phages", ]
fu_long$short_name <- sub("Bacteroides_phage_", "", fu_long$phage)

ggplot(fu_long, aes(x = host_phenotype, y = count + 1, fill = host_phenotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ short_name, scales = "free_y", ncol = 4) +
  scale_y_log10() +
  scale_fill_manual(values = c("healthy" = "steelblue", "cancer" = "firebrick")) +
  theme_bw() +
  theme(strip.text = element_text(size = 10)) +
  labs(
    title = "Normalised counts — Bacteroides FU phages",
    x = NULL, y = "Normalised count (log10 + 1)", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/FU_phages_boxplots.pdf", width = 14, height = 10)

# ── Separate PDF: ODE phages individually ─────────────────────────────────────

ode_long <- plot_long[plot_long$group == "ODE phages", ]
ode_long$short_name <- sub("Bacteroides_phage_", "", ode_long$phage)

ggplot(ode_long, aes(x = host_phenotype, y = count + 1, fill = host_phenotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ short_name, scales = "free_y", ncol = 2) +
  scale_y_log10() +
  scale_fill_manual(values = c("healthy" = "steelblue", "cancer" = "firebrick")) +
  theme_bw() +
  theme(strip.text = element_text(size = 10)) +
  labs(
    title = "Normalised counts — Bacteroides ODE phages",
    x = NULL, y = "Normalised count (log10 + 1)", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/ODE_phages_boxplots.pdf", width = 10, height = 8)
