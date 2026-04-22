##############################################################
# title: "Beta diversity in R - rumen solids by age and calf diet"
# version: alpha-matched styling for ordination plots
##############################################################

# Set this to the folder containing:
# - prewean_calf_metadata_simple.txt
# - 8w/core-metrics-results/
# - 1yr/core-metrics-results/
# - 2yr/core-metrics-results/
setwd("C:/ANSC516_Project")

library(tidyverse)
library(qiime2R)
library(vegan)

if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/beta_styled")) dir.create("output/beta_styled", recursive = TRUE)

# ---------------------------
# Metadata
# ---------------------------
meta <- read.delim(
  "metadata.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(meta)[colnames(meta) == "sample-id"] <- "SampleID"

meta <- meta %>%
  mutate(
    calf_diet = trimws(calf_diet),
    calf_diet = gsub("[ -]+", "_", calf_diet),
    calf_diet = tolower(calf_diet),
    Host_Age = trimws(Host_age)
  )

# Standardize diet labels
meta$calf_diet[meta$calf_diet %in% c("traditional_calf_starter", "trad_calf_starter", "traditional_starter")] <- "calf_starter"
meta$calf_diet[meta$calf_diet %in% c("cornsilage", "corn_silage")] <- "corn_silage"
meta$calf_diet[meta$calf_diet %in% c("mixed", "mixed_starter")] <- "mix"

meta <- meta %>%
  mutate(
    calf_diet = factor(
      calf_diet,
      levels = c("calf_starter", "mix", "corn_silage"),
      labels = c("Calf starter", "Mixed starter", "Corn silage")
    ),
    Host_Age = factor(Host_Age, levels = c("8w", "1yr", "2yr")),
    host_subject_id = as.factor(host_subject_id),
    host_sex = as.factor(host_sex)
  )

diet_colors <- c(
  "Calf starter" = "#E67E22",
  "Mixed starter" = "#7B75C8",
  "Corn silage" = "#2CB17E"
)

diet_shapes <- c(
  "Calf starter" = 21,   # circle
  "Mixed starter" = 24,  # triangle
  "Corn silage" = 22     # square
)

age_inputs <- tibble(
  Host_Age = c("8w", "1yr", "2yr"),
  core_dir = c(
    "8w-core-metrics-results",
    "1yr-core-metrics-results",
    "2yr-core-metrics-results"
  )
)

beta_files <- tribble(
  ~metric,                ~metric_label,           ~pcoa_file,                                 ~dist_file,
  "jaccard",              "Jaccard",               "jaccard_pcoa_results.qza",                 "jaccard_distance_matrix.qza",
  "bray_curtis",          "Bray-Curtis",           "bray_curtis_pcoa_results.qza",             "bray_curtis_distance_matrix.qza",
  "unweighted_unifrac",   "Unweighted UniFrac",    "unweighted_unifrac_pcoa_results.qza",      "unweighted_unifrac_distance_matrix.qza",
  "weighted_unifrac",     "Weighted UniFrac",      "weighted_unifrac_pcoa_results.qza",        "weighted_unifrac_distance_matrix.qza"
)

pairwise.adonis2 <- function(x, data, nperm = 999, ...) {
  lhs <- x[[2]]
  x1 <- x
  x1[[2]] <- NULL
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)

  group_var <- rhs.frame[, 1]
  group_levels <- unique(as.character(group_var))
  if (length(group_levels) < 2) {
    return(list(message = "Fewer than 2 groups available for pairwise PERMANOVA"))
  }

  pairs <- combn(group_levels, 2, simplify = FALSE)
  out_list <- vector("list", length(pairs))
  names(out_list) <- sapply(pairs, function(z) paste(z, collapse = "_vs_"))

  for (i in seq_along(pairs)) {
    keep <- group_var %in% pairs[[i]]

    if (inherits(eval(lhs), "dist")) {
      full_mat <- as.matrix(eval(lhs))
      reduced <- as.dist(full_mat[keep, keep, drop = FALSE])
    } else {
      reduced <- eval(lhs)[keep, , drop = FALSE]
    }

    sub_data <- data[keep, , drop = FALSE]
    sub_formula <- as.formula("reduced ~ calf_diet")
    out_list[[i]] <- adonis2(sub_formula, data = sub_data, permutations = nperm, ...)
  }

  out_list
}

make_styled_ordination <- function(pcoa_obj, plot_df, metric_label, age_label, permanova_p) {
  axis1 <- round(100 * pcoa_obj$data$ProportionExplained[1], 2)
  axis2 <- round(100 * pcoa_obj$data$ProportionExplained[2], 2)
aes(x = PC1, y = PC2, fill = calf_diet, shape = calf_diet)
ggplot(plot_df, aes(x = PC1, y = PC2, fill = calf_diet, shape = calf_diet)) +
  geom_point(size = 3, alpha = 0.9, color = "black", stroke = 0.8) +
  scale_fill_manual(values = diet_colors, drop = FALSE) +
  scale_shape_manual(values = diet_shapes, drop = FALSE) +
  labs(
    title = paste0(age_label, " ", metric_label),
    subtitle = paste0("PERMANOVA p = ", signif(permanova_p, 3)),
    x = paste0("PC1 (", axis1, "%)"),
    y = paste0("PC2 (", axis2, "%)"),
    fill = "Preweaning diet",
    shape = "Preweaning diet"
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
  ) +
    theme(
      panel.background = element_rect(fill = "#EBEBEB", color = NA),
      plot.background = element_rect(fill = "#EBEBEB", color = NA),
      panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.7),
      panel.grid.minor = element_line(color = "#DDDDDD", linewidth = 0.5),
      legend.position = "bottom",
      axis.line = element_line(color = "#555555", linewidth = 0.6),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(hjust = 0, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0, size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 12)
    )
}

all_stats <- list()

for (i in seq_len(nrow(age_inputs))) {
  age_label <- age_inputs$Host_Age[i]
  core_dir  <- age_inputs$core_dir[i]

  meta_age <- meta %>%
    filter(Host_Age == age_label)

  if (nrow(meta_age) == 0) {
    warning(paste("No metadata rows for age group:", age_label))
    next
  }

  for (j in seq_len(nrow(beta_files))) {
    metric_name  <- beta_files$metric[j]
    metric_label <- beta_files$metric_label[j]
    pcoa_path <- file.path(core_dir, beta_files$pcoa_file[j])
    dist_path <- file.path(core_dir, beta_files$dist_file[j])

    if (!file.exists(pcoa_path)) {
      warning(paste("Missing PCoA file:", pcoa_path))
      next
    }
    if (!file.exists(dist_path)) {
      warning(paste("Missing distance matrix file:", dist_path))
      next
    }

    pcoa_obj <- read_qza(pcoa_path)
    dist_obj <- read_qza(dist_path)

    plot_df <- pcoa_obj$data$Vectors %>%
      select(any_of(c("SampleID", "PC1", "PC2", "PC3"))) %>%
      inner_join(meta_age, by = "SampleID") %>%
      filter(!is.na(PC1), !is.na(PC2), !is.na(calf_diet))

    if (nrow(plot_df) == 0) {
      warning(paste("No overlapping samples for", age_label, metric_name))
      next
    }

    dist_mat <- as.matrix(dist_obj$data)
    common_ids <- intersect(rownames(dist_mat), meta_age$SampleID)

    if (length(common_ids) < 3) {
      warning(paste("Not enough samples for PERMANOVA:", age_label, metric_name))
      next
    }

    dist_mat_sub <- dist_mat[common_ids, common_ids, drop = FALSE]
    meta_sub <- meta_age %>%
      filter(SampleID %in% common_ids) %>%
      arrange(match(SampleID, common_ids)) %>%
      filter(!is.na(calf_diet))

    dist_mat_sub <- dist_mat_sub[meta_sub$SampleID, meta_sub$SampleID, drop = FALSE]

    if (!all(meta_sub$SampleID == rownames(dist_mat_sub))) {
      stop(paste("Sample order mismatch in", age_label, metric_name))
    }

    if (length(unique(meta_sub$calf_diet)) < 2) {
      warning(paste("Fewer than 2 diet groups for PERMANOVA:", age_label, metric_name))
      next
    }

    permanova_out <- adonis2(as.dist(dist_mat_sub) ~ calf_diet, data = meta_sub, permutations = 999)
    permanova_p <- permanova_out$`Pr(>F)`[1]

    stats_df <- tibble(
      Host_Age = age_label,
      metric = metric_name,
      permanova_F = permanova_out$F[1],
      p_value = permanova_p,
      r2 = permanova_out$R2[1]
    )
    all_stats[[paste(age_label, metric_name, sep = "_")]] <- stats_df

    ord_plot <- make_styled_ordination(pcoa_obj, plot_df, metric_label, age_label, permanova_p)

    write.csv(
      plot_df,
      file = file.path("output/beta_styled", paste0(age_label, "_", metric_name, "_ordination_data.csv")),
      row.names = FALSE
    )

    write.csv(
      as.data.frame(permanova_out),
      file = file.path("output/beta_styled", paste0(age_label, "_", metric_name, "_PERMANOVA.csv"))
    )

    pairwise_out <- pairwise.adonis2(as.dist(dist_mat_sub) ~ calf_diet, data = meta_sub, nperm = 999)
    pairwise_file <- file.path("output/beta_styled", paste0(age_label, "_", metric_name, "_pairwise_PERMANOVA.txt"))
    sink(pairwise_file)
    print(pairwise_out)
    sink()

    ggsave(
      filename = file.path("output/beta_styled", paste0(age_label, "_", metric_name, "_ordination_styled.png")),
      plot = ord_plot, width = 10, height = 8, dpi = 300
    )

    ggsave(
      filename = file.path("output/beta_styled", paste0(age_label, "_", metric_name, "_ordination_styled.pdf")),
      plot = ord_plot, width = 10, height = 8
    )
  }
}

if (length(all_stats) > 0) {
  combined_stats <- bind_rows(all_stats)
  write.csv(
    combined_stats,
    file = file.path("output/beta_styled", "all_beta_permanova_stats.csv"),
    row.names = FALSE
  )
}

message("Styled beta diversity plots completed.")
