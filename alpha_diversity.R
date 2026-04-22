##############################################################
# Alpha diversity in R - qiime2R adapted for preweaning diet
# rumen solid samples from dairy cattle at 8 weeks, 1 year, 2 years
##############################################################

# -----------------------------
# WHAT TO EDIT
# -----------------------------
# 1) Set your working directory to the folder containing:
#    - prewean_calf_metadata_simple.txt
#    - one subfolder per age group with core-metrics results
# 2) Update the three paths in age_inputs below if your folders differ.
#
# Example expected structure:
# project_folder/
#   prewean_calf_metadata_simple.txt
#   8w-core-metrics-results/
#   1yr-core-metrics-results/
#   2yr-core-metrics-results/
getwd()
setwd("C:/ANSC516_Project")

# -----------------------------
# PACKAGES
# -----------------------------
libs <- c("tidyverse", "qiime2R", "ggpubr")
for (pkg in libs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/alpha2")) dir.create("output/alpha2", recursive = TRUE)

# -----------------------------
# METADATA
# -----------------------------
metadata <- read.delim(
  "metadata.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  rename(
    SampleID = `sample-id`,
    calf_diet = calf_diet,
    Host_Age = Host_age,
    host_sex = host_sex,
    host_subject_id = host_subject_id
  ) %>%
  mutate(
    Host_Age = factor(Host_Age, levels = c("8w", "1yr", "2yr")),
    calf_diet = factor(calf_diet,
                       levels = c("calf_starter", "mix", "corn_silage"),
                       labels = c("Calf starter", "Mixed starter", "Corn silage")),
    host_subject_id = as.factor(host_subject_id)
  )

rownames(metadata) <- metadata$SampleID
str(metadata)
table(metadata$Host_Age, metadata$calf_diet)

# -----------------------------
# INPUT PATHS FOR EACH AGE GROUP
# -----------------------------
age_inputs <- tibble(
  Host_Age = c("8w", "1yr", "2yr"),
  core_dir = c(
    "8w-core-metrics-results",
    "1yr-core-metrics-results",
    "2yr-core-metrics-results"
  )
)

# -----------------------------
# APPEARANCE
# -----------------------------
diet_colors <- c(
  "Calf starter" = "#d95f02",
  "Mixed starter" = "#7570b3",
  "Corn silage" = "#1b9e77"
)

metric_y_labels <- c(
  observed_features = "Observed features",
  shannon = "Shannon diversity",
  pielou_e = "Pielou's evenness",
  faith_pd = "Faith's PD"
)

# -----------------------------
# FUNCTIONS
# -----------------------------
read_alpha_metric <- function(core_dir, file_name, metric_name) {
  obj <- read_qza(file.path(core_dir, file_name))
  df <- obj$data %>% rownames_to_column("SampleID")
  colnames(df)[2] <- metric_name
  df
}

make_alpha_plot <- function(df, metric_name, age_name) {
  y_lab <- metric_y_labels[[metric_name]]

  kw <- kruskal.test(df[[metric_name]] ~ df$calf_diet)
  p_lab <- paste0("Kruskal-Wallis p = ", signif(kw$p.value, 3))

  p <- ggplot(df, aes(x = calf_diet, y = .data[[metric_name]], fill = calf_diet)) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85) +
    scale_fill_manual(values = diet_colors, name = "Preweaning diet") +
    labs(
      title = paste(age_name, y_lab),
      subtitle = p_lab,
      x = "Preweaning diet treatment",
      y = y_lab
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "bottom"
    )

  p
}

# -----------------------------
# MAIN LOOP: MAKE 4 ALPHA PLOTS PER AGE GROUP
# -----------------------------
for (i in seq_len(nrow(age_inputs))) {
  age_name <- age_inputs$Host_Age[i]
  core_dir <- age_inputs$core_dir[i]

  age_meta <- metadata %>% filter(Host_Age == age_name)

  evenness <- read_alpha_metric(core_dir, "evenness_vector.qza", "pielou_e")
  observed_features <- read_alpha_metric(core_dir, "observed_features_vector.qza", "observed_features")
  shannon <- read_alpha_metric(core_dir, "shannon_vector.qza", "shannon")
  faith_pd <- read_alpha_metric(core_dir, "faith_pd_vector.qza", "faith_pd")

  alpha_df <- age_meta %>%
    left_join(evenness, by = "SampleID") %>%
    left_join(observed_features, by = "SampleID") %>%
    left_join(shannon, by = "SampleID") %>%
    left_join(faith_pd, by = "SampleID")

  write.csv(alpha_df,
            file = file.path("output/alpha", paste0(age_name, "_alpha_diversity_values.csv")),
            row.names = FALSE)

  metrics <- c("observed_features", "shannon", "pielou_e", "faith_pd")

  for (metric_name in metrics) {
    plot_obj <- make_alpha_plot(alpha_df, metric_name, age_name)
    print(plot_obj)
    ggsave(
      filename = file.path("output/alpha2", paste0(age_name, "_", metric_name, "_alpha_plot.png")),
      plot = plot_obj,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
}
warnings()
# -----------------------------
# OPTIONAL: COMBINED SUMMARY TABLE
# -----------------------------
alpha_summary <- metadata %>%
  count(Host_Age, calf_diet, name = "n_samples")
write.csv(alpha_summary, "output/alpha2/sample_counts_by_age_and_diet.csv", row.names = FALSE)

message("Alpha-diversity plots saved in output/alpha2/")

