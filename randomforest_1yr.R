# ============================================================
# Random Forest analysis for microbiome data by age group
# Publication-quality figures for manuscript / class submission
# ============================================================

# ---------------------------
# Install packages if needed
# ---------------------------
# install.packages(c("randomForest", "tidyverse", "reshape2"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")
# install.packages("remotes")
# remotes::install_github("jbisanz/qiime2R")

# ---------------------------
# Load libraries
# ---------------------------
library(qiime2R)
library(phyloseq)
library(randomForest)
library(tidyverse)
library(reshape2)

# ---------------------------
# User settings
# ---------------------------
project_dir <- "C:/ANSC516_Project"

# Change these for each script/run:
# age_group <- "8w"
# feature_table_file <- "8w-core-metrics-results/rarefied_table.qza"

# Example alternatives:
age_group <- "1yr"
feature_table_file <- "1yr-core-metrics-results/rarefied_table.qza"

# age_group <- "2yr"
# feature_table_file <- "2yr-core-metrics-results/rarefied_table.qza"

taxonomy_file <- "taxonomy.qza"
tree_file <- "rooted-tree_simple.qza"
metadata_file <- "metadata.txt"

# Output folder
out_dir <- file.path(project_dir, paste0("output_", age_group), "random_forest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setwd(project_dir)

# ---------------------------
# Import data
# ---------------------------
physeq <- qza_to_phyloseq(
  features = feature_table_file,
  tree = tree_file,
  taxonomy = taxonomy_file,
  metadata = metadata_file
)

# ---------------------------
# Clean metadata
# ---------------------------
meta <- data.frame(sample_data(physeq), check.names = FALSE)

meta$calf_diet <- factor(
  meta$calf_diet,
  levels = c("calf_starter", "corn_silage", "mix")
)

sample_data(physeq) <- sample_data(meta)

# ---------------------------
# Transform to relative abundance
# ---------------------------
physeq.rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# ---------------------------
# Filter low-abundance ASVs
# Keep ASVs with mean relative abundance > 0.001
# ---------------------------
physeq.filt <- filter_taxa(
  physeq.rel,
  function(x) mean(x) > 0.001,
  prune = TRUE
)

# Remove taxa with zero total abundance after filtering
physeq.filt <- prune_taxa(taxa_sums(physeq.filt) > 0, physeq.filt)

# ---------------------------
# Prepare model matrix
# ---------------------------
otu <- as.data.frame(t(otu_table(physeq.filt)))
meta <- data.frame(sample_data(physeq.filt), check.names = FALSE)

meta$calf_diet <- droplevels(meta$calf_diet)

# Ensure row order matches
otu <- otu[rownames(meta), , drop = FALSE]

# Safety checks
stopifnot(all(rownames(otu) == rownames(meta)))
stopifnot(nrow(otu) == nrow(meta))

# ---------------------------
# Run Random Forest
# ---------------------------
set.seed(123)

rf_model <- randomForest(
  x = otu,
  y = meta$calf_diet,
  ntree = 500,
  importance = TRUE
)

# Print model in console
print(rf_model)

# Save model summary
capture.output(
  print(rf_model),
  file = file.path(out_dir, paste0("rf_model_summary_", age_group, ".txt"))
)

# ---------------------------
# Predictions and accuracy
# ---------------------------
rf_pred <- predict(rf_model, otu)

overall_accuracy <- mean(rf_pred == meta$calf_diet)

write.table(
  data.frame(
    Age = age_group,
    Accuracy = overall_accuracy
  ),
  file = file.path(out_dir, paste0("model_accuracy_", age_group, ".txt")),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# ---------------------------
# Confusion matrix
# ---------------------------
conf_mat <- table(
  Actual = meta$calf_diet,
  Predicted = rf_pred
)

write.csv(
  as.data.frame.matrix(conf_mat),
  file = file.path(out_dir, paste0("confusion_matrix_", age_group, ".csv"))
)

# ---------------------------
# Variable importance
# ---------------------------
imp <- as.data.frame(importance(rf_model))
imp$ASV <- rownames(imp)

# Sort by MeanDecreaseGini
imp_sorted <- imp %>%
  arrange(desc(MeanDecreaseGini))

# ---------------------------
# Add taxonomy labels
# ---------------------------
tax_df <- as.data.frame(tax_table(physeq.filt), check.names = FALSE)

# Make sure ASVs match
tax_df$ASV <- rownames(tax_df)

# Join taxonomy to importance table
imp_sorted <- imp_sorted %>%
  left_join(tax_df, by = "ASV")

# Create readable label priority:
# Genus -> Family -> Phylum -> ASV ID
imp_sorted <- imp_sorted %>%
  mutate(
    TaxonLabel = case_when(
      !is.na(Genus) & Genus != "" ~ paste0("g__", Genus),
      !is.na(Family) & Family != "" ~ paste0("f__", Family),
      !is.na(Phylum) & Phylum != "" ~ paste0("p__", Phylum),
      TRUE ~ ASV
    )
  )

# Optional: make labels unique if repeated
imp_sorted$TaxonLabel <- make.unique(imp_sorted$TaxonLabel)

# Save full importance table
write.csv(
  imp_sorted,
  file = file.path(out_dir, paste0("ASV_importance_", age_group, ".csv")),
  row.names = FALSE
)

# ---------------------------
# Top 15 importance plot
# Publication-quality barplot
# ---------------------------
top15 <- imp_sorted %>%
  slice(1:15)

p_importance <- ggplot(
  top15,
  aes(x = reorder(TaxonLabel, MeanDecreaseGini), y = MeanDecreaseGini)
) +
  geom_col(fill = "#2C7FB8", width = 0.8) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = paste0(age_group, " Random Forest: Top 15 Predictive Taxa"),
    x = "Taxa",
    y = "Variable Importance (Mean Decrease Gini)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

ggsave(
  filename = file.path(out_dir, paste0("top15_taxa_importance_", age_group, ".png")),
  plot = p_importance,
  width = 8,
  height = 6,
  dpi = 300
)

# ---------------------------
# Confusion matrix heatmap
# Publication-quality figure
# ---------------------------
conf_df <- as.data.frame(conf_mat)

p_conf <- ggplot(conf_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Freq), size = 5, color = "black") +
  scale_fill_gradient(low = "#C6DBEF", high = "#08519C") +
  theme_classic(base_size = 14) +
  labs(
    title = paste0(age_group, " Random Forest Classification"),
    x = "Predicted diet",
    y = "Actual diet",
    fill = "Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

ggsave(
  filename = file.path(out_dir, paste0("confusion_matrix_heatmap_", age_group, ".png")),
  plot = p_conf,
  width = 6,
  height = 5,
  dpi = 300
)

# ---------------------------
# Optional diagnostic OOB plot
# Useful for checking model, but not ideal for final paper
# ---------------------------
png(
  filename = file.path(out_dir, paste0("OOB_error_", age_group, ".png")),
  width = 800,
  height = 600
)
plot(rf_model, main = paste0(age_group, " Random Forest OOB Error"))
dev.off()

# ---------------------------
# Save top 15 table only
# ---------------------------
write.csv(
  top15,
  file = file.path(out_dir, paste0("top15_taxa_importance_", age_group, ".csv")),
  row.names = FALSE
)

# ---------------------------
# Console summary
# ---------------------------
cat("\n=============================\n")
cat("Random Forest analysis complete\n")
cat("Age group:", age_group, "\n")
cat("Overall accuracy:", round(overall_accuracy, 3), "\n")
cat("Output folder:", out_dir, "\n")
cat("=============================\n\n")