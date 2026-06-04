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
  dplyr::filter(group != "Cryptic phages") |>
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

# ── CDS-level read counts along genome ────────────────────────────────────────

read_counts_cds <- function(f) {
  df <- read.table(f, header = TRUE, skip = 1, sep = "\t")
  data.frame(
    gene_id = df$Geneid,
    phage   = df$Chr,
    start   = df$Start,
    end     = df$End,
    count   = df[[ncol(df)]],
    Run     = sub("\\.cds\\.counts\\.txt$", "", basename(f))
  )
}

cds_long <- do.call(rbind, lapply(files, read_counts_cds))
cds_long <- cds_long[cds_long$Run %in% meta$Run, ]

sf <- sizeFactors(dds)
cds_long$norm_count <- cds_long$count / sf[cds_long$Run]
cds_long$phage_name <- name_map[cds_long$phage]
cds_long$mid        <- (cds_long$start + cds_long$end) / 2
cds_long <- merge(cds_long, meta[, c("Run", "host_phenotype")], by = "Run")
cds_long <- cds_long[order(cds_long$phage_name, cds_long$Run, cds_long$mid), ]

# ── ODE phages: reads per CDS along genome ────────────────────────────────────

ode_cds <- cds_long[grepl("_ODE", cds_long$phage_name), ]
ode_cds$short_name <- sub("Bacteroides_phage_", "", ode_cds$phage_name)
n_ode <- length(unique(ode_cds$short_name))

ggplot(ode_cds, aes(x = mid, y = norm_count + 1, colour = host_phenotype, group = Run)) +
  geom_line(alpha = 0.4, linewidth = 0.4) +
  geom_point(alpha = 0.7, size = 1.2) +
  facet_wrap(~ short_name, scales = "free", ncol = 1) +
  scale_y_log10() +
  scale_colour_manual(values = c("healthy" = "steelblue", "cancer" = "firebrick")) +
  theme_bw() +
  labs(
    title = "Read counts per CDS along genome — ODE phages",
    x = "Genomic position (bp)", y = "Normalised count (log10 + 1)", colour = NULL
  )

ggsave("F:/CRC_Tumour_WGS/ODE_phages_genome_coverage.pdf", width = 10, height = 4 * n_ode)

# ── FU phages: reads per CDS along genome ─────────────────────────────────────

fu_cds <- cds_long[grepl("_FU", cds_long$phage_name), ]
fu_cds$short_name <- sub("Bacteroides_phage_", "", fu_cds$phage_name)
n_fu <- length(unique(fu_cds$short_name))

ggplot(fu_cds, aes(x = mid, y = norm_count + 1, colour = host_phenotype, group = Run)) +
  geom_line(alpha = 0.4, linewidth = 0.4) +
  geom_point(alpha = 0.7, size = 1.2) +
  facet_wrap(~ short_name, scales = "free", ncol = 1) +
  scale_y_log10() +
  scale_colour_manual(values = c("healthy" = "steelblue", "cancer" = "firebrick")) +
  theme_bw() +
  labs(
    title = "Read counts per CDS along genome — FU phages",
    x = "Genomic position (bp)", y = "Normalised count (log10 + 1)", colour = NULL
  )

ggsave("F:/CRC_Tumour_WGS/FU_phages_genome_coverage.pdf", width = 10, height = 4 * n_fu)

# ── Phage prevalence barplots (cumulative: ≥N genes with ≥1 read) ──────────────

gene_presence <- cds_long |>
  dplyr::group_by(phage_name, Run, host_phenotype) |>
  dplyr::summarise(n_genes_detected = sum(count > 0), .groups = "drop")

n_cancer  <- sum(meta$host_phenotype == "cancer")
n_healthy <- sum(meta$host_phenotype == "healthy")

thresholds <- 1:12

# All phage × phenotype × threshold combinations, defaulting to 0
all_combos <- expand.grid(
  phage_name     = unique(gene_presence$phage_name),
  host_phenotype = c("cancer", "healthy"),
  threshold      = thresholds,
  stringsAsFactors = FALSE
)

prevalence_cum <- do.call(rbind, lapply(thresholds, function(thr) {
  gene_presence |>
    dplyr::filter(n_genes_detected >= thr) |>
    dplyr::group_by(phage_name, host_phenotype) |>
    dplyr::summarise(n_positive = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(threshold = thr)
}))

prevalence_cum <- merge(all_combos, prevalence_cum,
                        by = c("phage_name", "host_phenotype", "threshold"),
                        all.x = TRUE)
prevalence_cum$n_positive[is.na(prevalence_cum$n_positive)] <- 0L
prevalence_cum$threshold <- factor(prevalence_cum$threshold, levels = as.character(thresholds))
prevalence_cum$group     <- sapply(prevalence_cum$phage_name, assign_group)
prevalence_cum$short_name <- sub("Bacteroides_cryptic_phage_", "", prevalence_cum$phage_name)
prevalence_cum$short_name <- sub("Bacteroides_phage_",         "", prevalence_cum$short_name)

# Order panels by group then cancer count at threshold 8 (descending)
order_at_8 <- prevalence_cum[prevalence_cum$threshold == "8" &
                               prevalence_cum$host_phenotype == "cancer", ]
order_at_8 <- order_at_8[order(order_at_8$group, -order_at_8$n_positive), ]
prevalence_cum$short_name <- factor(prevalence_cum$short_name,
                                    levels = unique(order_at_8$short_name))

n_phages <- length(unique(prevalence_cum$short_name))

ggplot(prevalence_cum, aes(x = threshold, y = n_positive, fill = host_phenotype)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ short_name, ncol = 4, scales = "free_y") +
  scale_fill_manual(
    values = c("healthy" = "steelblue", "cancer" = "firebrick"),
    labels = c(
      "healthy" = paste0("Healthy (n=", n_healthy, ")"),
      "cancer"  = paste0("Cancer (n=",  n_cancer,  ")")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  theme_bw() +
  theme(strip.text = element_text(size = 9)) +
  labs(
    title    = "Phage prevalence by sample type",
    subtitle = "Each bar = samples with ≥N genes having at least 1 mapped read",
    x = "Minimum genes detected (N)", y = "Number of positive samples", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/phage_prevalence_barplots.pdf",
       width = 14, height = ceiling(n_phages / 4) * 3)

# ── ODE phages: distribution of genes-detected per sample ─────────────────────

ode_per_sample <- cds_long[grepl("_ODE", cds_long$phage_name), ] |>
  dplyr::group_by(phage_name, Run, host_phenotype) |>
  dplyr::summarise(n_genes_detected = sum(count > 0), .groups = "drop")

ode_pooled_per_sample <- cds_long[grepl("_ODE", cds_long$phage_name), ] |>
  dplyr::group_by(Run, host_phenotype) |>
  dplyr::summarise(n_genes_detected = sum(count > 0), .groups = "drop") |>
  dplyr::mutate(phage_name = "All ODE phages (pooled)")

ode_all_combos <- expand.grid(
  phage_name     = c(unique(ode_per_sample$phage_name), "All ODE phages (pooled)"),
  host_phenotype = c("cancer", "healthy"),
  threshold      = thresholds,
  stringsAsFactors = FALSE
)

ode_dist <- do.call(rbind, lapply(thresholds, function(thr) {
  rbind(
    ode_per_sample |>
      dplyr::filter(n_genes_detected >= thr) |>
      dplyr::group_by(phage_name, host_phenotype) |>
      dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(threshold = thr),
    ode_pooled_per_sample |>
      dplyr::filter(n_genes_detected >= thr) |>
      dplyr::group_by(phage_name, host_phenotype) |>
      dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(threshold = thr)
  )
}))

ode_dist <- merge(ode_all_combos, ode_dist,
                  by = c("phage_name", "host_phenotype", "threshold"),
                  all.x = TRUE)
ode_dist$n_samples[is.na(ode_dist$n_samples)] <- 0L
ode_dist$threshold <- factor(ode_dist$threshold, levels = as.character(thresholds))

ode_dist$short_name <- sub("Bacteroides_phage_", "", ode_dist$phage_name)
ode_dist$short_name <- factor(ode_dist$short_name,
                               levels = c("All ODE phages (pooled)",
                                          setdiff(unique(ode_dist$short_name),
                                                  "All ODE phages (pooled)")))

n_ode_phages <- length(unique(ode_dist$short_name))

ggplot(ode_dist, aes(x = threshold, y = n_samples, fill = host_phenotype)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ short_name, ncol = 1) +
  scale_fill_manual(
    values = c("healthy" = "steelblue", "cancer" = "firebrick"),
    labels = c(
      "healthy" = paste0("Healthy (n=", n_healthy, ")"),
      "cancer"  = paste0("Cancer (n=",  n_cancer,  ")")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  theme_bw() +
  theme(strip.text = element_text(size = 10)) +
  labs(
    title    = "ODE phages — gene detection distribution per sample",
    subtitle = "Each bar = samples with ≥N genes having at least 1 mapped read",
    x = "Minimum genes detected (N)", y = "Number of samples", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/ODE_phages_gene_detection_distribution.pdf",
       width = 8, height = 4 * n_ode_phages)

# ── ODE phages: line plot ──────────────────────────────────────────────────────

ode_lines <- ode_dist[as.character(ode_dist$short_name) != "All ODE phages (pooled)", ]
ode_lines <- droplevels(ode_lines)
ode_lines$threshold_num <- as.integer(as.character(ode_lines$threshold))

ggplot(ode_lines, aes(x = threshold_num, y = n_samples,
                       color = host_phenotype, shape = short_name,
                       group = interaction(short_name, host_phenotype))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  scale_color_manual(values = c("cancer"  = "firebrick", "healthy" = "steelblue"),
                     labels = c("cancer"  = paste0("CRC (n=",  n_cancer,  ")"),
                                "healthy" = paste0("CTR (n=", n_healthy, ")"))) +
  scale_shape_manual(values = setNames(seq_along(levels(ode_lines$short_name)),
                                       levels(ode_lines$short_name))) +
  theme_bw() +
  labs(
    title    = "ODE phages — gene detection distribution",
    subtitle = "Each line = samples with ≥N genes having at least 1 mapped read",
    x = "Minimum genes detected (N)", y = "Number of samples",
    color = NULL, shape = "Phage"
  )

ggsave("F:/CRC_Tumour_WGS/ODE_phages_gene_detection_lines.pdf", width = 8, height = 5)

# ── FU phages: distribution of genes-detected per sample ──────────────────────

fu_per_sample <- cds_long[grepl("_FU", cds_long$phage_name), ] |>
  dplyr::group_by(phage_name, Run, host_phenotype) |>
  dplyr::summarise(n_genes_detected = sum(count > 0), .groups = "drop")

fu_pooled_per_sample <- cds_long[grepl("_FU", cds_long$phage_name), ] |>
  dplyr::group_by(Run, host_phenotype) |>
  dplyr::summarise(n_genes_detected = sum(count > 0), .groups = "drop") |>
  dplyr::mutate(phage_name = "All FU phages (pooled)")

fu_all_combos <- expand.grid(
  phage_name     = c(unique(fu_per_sample$phage_name), "All FU phages (pooled)"),
  host_phenotype = c("cancer", "healthy"),
  threshold      = thresholds,
  stringsAsFactors = FALSE
)

fu_dist <- do.call(rbind, lapply(thresholds, function(thr) {
  rbind(
    fu_per_sample |>
      dplyr::filter(n_genes_detected >= thr) |>
      dplyr::group_by(phage_name, host_phenotype) |>
      dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(threshold = thr),
    fu_pooled_per_sample |>
      dplyr::filter(n_genes_detected >= thr) |>
      dplyr::group_by(phage_name, host_phenotype) |>
      dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(threshold = thr)
  )
}))

fu_dist <- merge(fu_all_combos, fu_dist,
                 by = c("phage_name", "host_phenotype", "threshold"),
                 all.x = TRUE)
fu_dist$n_samples[is.na(fu_dist$n_samples)] <- 0L
fu_dist$threshold <- factor(fu_dist$threshold, levels = as.character(thresholds))

fu_dist$short_name <- sub("Bacteroides_phage_", "", fu_dist$phage_name)
fu_dist$short_name <- factor(fu_dist$short_name,
                              levels = c("All FU phages (pooled)",
                                         setdiff(unique(fu_dist$short_name),
                                                 "All FU phages (pooled)")))

n_fu_phages <- length(unique(fu_dist$short_name))

ggplot(fu_dist, aes(x = threshold, y = n_samples, fill = host_phenotype)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ short_name, ncol = 1) +
  scale_fill_manual(
    values = c("healthy" = "steelblue", "cancer" = "firebrick"),
    labels = c(
      "healthy" = paste0("Healthy (n=", n_healthy, ")"),
      "cancer"  = paste0("Cancer (n=",  n_cancer,  ")")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  theme_bw() +
  theme(strip.text = element_text(size = 10)) +
  labs(
    title    = "FU phages — gene detection distribution per sample",
    subtitle = "Each bar = samples with ≥N genes having at least 1 mapped read",
    x = "Minimum genes detected (N)", y = "Number of samples", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/FU_phages_gene_detection_distribution.pdf",
       width = 8, height = 4 * n_fu_phages)

# ── FU phages: line plot ───────────────────────────────────────────────────────

fu_lines <- fu_dist[as.character(fu_dist$short_name) != "All FU phages (pooled)", ]
fu_lines <- droplevels(fu_lines)
fu_lines$threshold_num <- as.integer(as.character(fu_lines$threshold))

ggplot(fu_lines, aes(x = threshold_num, y = n_samples,
                      color = host_phenotype, shape = short_name,
                      group = interaction(short_name, host_phenotype))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  scale_color_manual(values = c("cancer"  = "firebrick", "healthy" = "steelblue"),
                     labels = c("cancer"  = paste0("CRC (n=",  n_cancer,  ")"),
                                "healthy" = paste0("CTR (n=", n_healthy, ")"))) +
  scale_shape_manual(values = setNames(seq_along(levels(fu_lines$short_name)),
                                       levels(fu_lines$short_name))) +
  theme_bw() +
  labs(
    title    = "FU phages — gene detection distribution",
    subtitle = "Each line = samples with ≥N genes having at least 1 mapped read",
    x = "Minimum genes detected (N)", y = "Number of samples",
    color = NULL, shape = "Phage"
  )

ggsave("F:/CRC_Tumour_WGS/FU_phages_gene_detection_lines.pdf", width = 8, height = 5)

# ── Prevalence at ≥50% gene coverage per phage ────────────────────────────────

total_genes <- cds_long |>
  dplyr::group_by(phage_name) |>
  dplyr::summarise(total_genes = dplyr::n_distinct(gene_id), .groups = "drop") |>
  dplyr::mutate(threshold_50pct = ceiling(total_genes / 2))

prev_50pct <- merge(gene_presence, total_genes, by = "phage_name") |>
  dplyr::mutate(phage_positive = n_genes_detected >= threshold_50pct) |>
  dplyr::group_by(phage_name, host_phenotype, total_genes, threshold_50pct) |>
  dplyr::summarise(n_positive = sum(phage_positive), .groups = "drop")

prev_50pct$group      <- sapply(prev_50pct$phage_name, assign_group)
prev_50pct$short_name <- sub("Bacteroides_cryptic_phage_", "", prev_50pct$phage_name)
prev_50pct$short_name <- sub("Bacteroides_phage_",         "", prev_50pct$short_name)
prev_50pct$label      <- paste0(prev_50pct$short_name,
                                 " (", prev_50pct$threshold_50pct, "/",
                                 prev_50pct$total_genes, ")")

cancer_order_50 <- prev_50pct[prev_50pct$host_phenotype == "cancer", ]
cancer_order_50 <- cancer_order_50[order(cancer_order_50$group, -cancer_order_50$n_positive), ]
prev_50pct$label <- factor(prev_50pct$label, levels = rev(unique(cancer_order_50$label)))

ggplot(prev_50pct, aes(y = label, x = n_positive, fill = host_phenotype)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c("healthy" = "steelblue", "cancer" = "firebrick"),
    labels = c(
      "healthy" = paste0("Healthy (n=", n_healthy, ")"),
      "cancer"  = paste0("Cancer (n=",  n_cancer,  ")")
    )
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  theme_bw() +
  theme(strip.text = element_text(size = 10)) +
  labs(
    title    = "Phage prevalence at ≥50% gene coverage",
    subtitle = "Positive = ≥ ceiling(total_genes / 2) genes with ≥1 mapped read; labels show threshold/total",
    y = NULL, x = "Number of positive samples", fill = NULL
  )

ggsave("F:/CRC_Tumour_WGS/phage_prevalence_50pct_threshold.pdf", width = 11, height = 14)

# ── All phages: gene detection line plot ───────────────────────────────────────
# FU/ODE/TAA kept as group panels; cryptic phages and singletons get individual panels

all_lines <- prevalence_cum
all_lines$threshold_num <- as.integer(as.character(all_lines$threshold))

group_panels <- c("FU phages", "ODE phages", "TAA phages")
all_lines$facet_var <- ifelse(
  all_lines$group %in% group_panels,
  all_lines$group,
  as.character(all_lines$short_name)
)

# Panel order: group panels first, then individual phages in existing factor order
indiv_ordered <- levels(all_lines$short_name)[
  levels(all_lines$short_name) %in%
    unique(all_lines$facet_var[!all_lines$facet_var %in% group_panels])
]
all_lines$facet_var <- factor(all_lines$facet_var,
                               levels = c(group_panels, indiv_ordered))

# Shapes for phages within group panels (only needed there)
group_phages <- unique(as.character(
  all_lines$short_name[all_lines$facet_var %in% group_panels]
))
group_shape_vals <- setNames(seq_along(group_phages), group_phages)

# End-of-line labels for group panels (one label per phage, cancer line only)
line_labels <- all_lines[
  all_lines$threshold_num == 12 &
  all_lines$host_phenotype == "cancer" &
  all_lines$facet_var %in% group_panels, ]

n_panels <- nlevels(all_lines$facet_var)

ggplot(all_lines, aes(x = threshold_num, y = n_positive,
                       color = host_phenotype,
                       group = interaction(short_name, host_phenotype))) +
  geom_line(linewidth = 0.7) +
  geom_point(
    data = all_lines[all_lines$facet_var %in% group_panels, ],
    aes(shape = short_name), size = 2.5
  ) +
  geom_point(
    data = all_lines[!all_lines$facet_var %in% group_panels, ],
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data    = line_labels,
    aes(label = short_name),
    color   = "black", size = 2.5, hjust = 0,
    direction = "y", nudge_x = 0.3,
    segment.size = 0.3, show.legend = FALSE
  ) +
  facet_wrap(~ facet_var, ncol = 4) +
  scale_x_continuous(breaks = 1:12, expand = expansion(mult = c(0.05, 0.2))) +
  scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
  scale_color_manual(values = c("cancer"  = "firebrick", "healthy" = "steelblue"),
                     labels = c("cancer"  = paste0("CRC (n=",  n_cancer,  ")"),
                                "healthy" = paste0("CTR (n=", n_healthy, ")"))) +
  scale_shape_manual(values = group_shape_vals, name = "Phage") +
  theme_bw() +
  theme(strip.text      = element_text(size = 8),
        legend.position = "right") +
  labs(
    title    = "All phages — gene detection distribution",
    subtitle = "Each line = samples with ≥N genes having at least 1 mapped read",
    x = "Minimum genes detected (N)", y = "Number of positive samples",
    color = NULL
  )

ggsave("F:/CRC_Tumour_WGS/all_phages_gene_detection_lines.pdf",
       width = 16, height = ceiling(n_panels / 4) * 4)

# ── All phages: gene detection line plot extended to 30 genes ─────────────────

thresholds_30 <- 1:30

all_combos_30 <- expand.grid(
  phage_name     = unique(gene_presence$phage_name),
  host_phenotype = c("cancer", "healthy"),
  threshold      = thresholds_30,
  stringsAsFactors = FALSE
)

prevalence_cum_30 <- do.call(rbind, lapply(thresholds_30, function(thr) {
  gene_presence |>
    dplyr::filter(n_genes_detected >= thr) |>
    dplyr::group_by(phage_name, host_phenotype) |>
    dplyr::summarise(n_positive = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(threshold = thr)
}))

prevalence_cum_30 <- merge(all_combos_30, prevalence_cum_30,
                            by = c("phage_name", "host_phenotype", "threshold"),
                            all.x = TRUE)
prevalence_cum_30$n_positive[is.na(prevalence_cum_30$n_positive)] <- 0L
prevalence_cum_30$threshold  <- factor(prevalence_cum_30$threshold,
                                        levels = as.character(thresholds_30))
prevalence_cum_30$group      <- sapply(prevalence_cum_30$phage_name, assign_group)
prevalence_cum_30$short_name <- sub("Bacteroides_cryptic_phage_", "",
                                     prevalence_cum_30$phage_name)
prevalence_cum_30$short_name <- sub("Bacteroides_phage_", "",
                                     prevalence_cum_30$short_name)

# Reuse same panel ordering as the 1-12 plot (by group then cancer count at thr 8)
prevalence_cum_30$short_name <- factor(prevalence_cum_30$short_name,
                                        levels = levels(prevalence_cum$short_name))

all_lines_30 <- prevalence_cum_30[prevalence_cum_30$group != "Cryptic phages", ]
all_lines_30$threshold_num <- as.integer(as.character(all_lines_30$threshold))
all_lines_30$facet_var <- ifelse(
  all_lines_30$group %in% group_panels,
  all_lines_30$group,
  as.character(all_lines_30$short_name)
)
all_lines_30$facet_var <- factor(all_lines_30$facet_var,
                                  levels = levels(all_lines$facet_var))

line_labels_30 <- all_lines_30[
  all_lines_30$threshold_num == 30 &
  all_lines_30$host_phenotype == "cancer" &
  all_lines_30$facet_var %in% group_panels, ]

ggplot(all_lines_30, aes(x = threshold_num, y = n_positive,
                          color = host_phenotype,
                          group = interaction(short_name, host_phenotype))) +
  geom_line(linewidth = 0.7) +
  geom_point(
    data = all_lines_30[all_lines_30$facet_var %in% group_panels, ],
    aes(shape = short_name), size = 2.5
  ) +
  geom_point(
    data = all_lines_30[!all_lines_30$facet_var %in% group_panels, ],
    size = 2
  ) +
  ggrepel::geom_text_repel(
    data      = line_labels_30,
    aes(label = short_name),
    color     = "black", size = 2.5, hjust = 0,
    direction = "y", nudge_x = 0.5,
    segment.size = 0.3, show.legend = FALSE
  ) +
  facet_wrap(~ facet_var, ncol = 4) +
  scale_x_continuous(breaks = seq(0, 30, 5),
                     expand = expansion(mult = c(0.05, 0.2))) +
  scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
  scale_color_manual(values = c("cancer"  = "firebrick", "healthy" = "steelblue"),
                     labels = c("cancer"  = paste0("CRC (n=",  n_cancer,  ")"),
                                "healthy" = paste0("CTR (n=", n_healthy, ")"))) +
  scale_shape_manual(values = group_shape_vals, name = "Phage") +
  theme_bw() +
  theme(strip.text      = element_text(size = 8),
        legend.position = "right") +
  labs(
    title    = "All phages — gene detection distribution (up to 30 genes)",
    subtitle = "Each line = samples with ≥N genes having at least 1 mapped read",
    x = "Minimum genes detected (N)", y = "Number of positive samples",
    color = NULL
  )

ggsave("F:/CRC_Tumour_WGS/all_phages_gene_detection_lines_30.pdf",
       width = 16, height = ceiling(n_panels / 4) * 4)

# ── Prevalence table (≥8 genes) exported to Excel ─────────────────────────────

prev_table <- gene_presence |>
  dplyr::mutate(positive = n_genes_detected >= 8) |>
  dplyr::group_by(phage_name, host_phenotype) |>
  dplyr::summarise(n_positive = sum(positive), .groups = "drop") |>
  tidyr::pivot_wider(names_from = host_phenotype, values_from = n_positive,
                     values_fill = 0) |>
  dplyr::rename(n_CTR = healthy, n_CRC = cancer)

prev_table$group      <- sapply(prev_table$phage_name, assign_group)
prev_table$short_name <- sub("Bacteroides_cryptic_phage_", "", prev_table$phage_name)
prev_table$short_name <- sub("Bacteroides_phage_",         "", prev_table$short_name)

prev_table <- prev_table[order(prev_table$group, -prev_table$n_CRC), ]
prev_table <- prev_table[, c("short_name", "phage_name", "group", "n_CRC", "n_CTR")]

writexl::write_xlsx(prev_table,
                    "F:/CRC_Tumour_WGS/phage_prevalence_8genes.xlsx")
