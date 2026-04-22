# Co-occurrence analysis for one age group
# Edited from demo script

library(Hmisc)
library(plyr)
library(reshape2)
library(qiime2R)
library(ggplot2)

getwd()
setwd("C:/ANSC516_Project")

dir.create("output", showWarnings = FALSE)
dir.create("output/cooccurrence", showWarnings = FALSE)

# ---------------------------
# Input files: change age group here
# ---------------------------
age_group <- "2yr"
table_file <- "2yr-table.qza"      # change to 1yr-table.qza or 2yr-table.qza later
metadata_file <- "metadata.txt"
taxonomy_file <- "taxonomy.qza"

ASVs <- read_qza(table_file)
ASV_table <- as.data.frame(ASVs$data)

# Keep ASV IDs in a key
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)]
ASV_table <- ASV_table[, -(ncol(ASV_table)-1):-ncol(ASV_table)]

dataset <- as.data.frame(t(ASV_table))

# ---------------------------
# Metadata
# ---------------------------
metadata <- read_tsv(metadata_file)
str(metadata)

# Keep only this age group in case metadata has all ages
metadata <- subset(metadata, Host_age == age_group)

# Merge metadata with ASV table
dataset <- merge(metadata, dataset, by.x = "sample-id", by.y = 0)

# Use diet as treatment
my_column <- "calf_diet"
treatments <- as.vector(unique(dataset[[my_column]]))

datasetn <- dataset
datasetn[datasetn == 0] <- NA

summary(metadata[[my_column]])

# At least this many non-zero values per treatment
# For 8w, sample sizes are smaller, so keep this lenient
n1 <- 2
n2 <- 2
n3 <- 2

num_metadata_columns <- ncol(metadata)
q_cutoff <- 0.05

final_results <- data.frame()

for (i in 1:length(treatments)) {
  print(paste("reading", treatments[i]))
  
  temp <- subset(dataset, get(my_column) == treatments[i])
  tempn <- subset(datasetn, get(my_column) == treatments[i])
  
  results <- rcorr(as.matrix(temp[, -c(1:num_metadata_columns)]), type = "spearman")
  resultsn <- rcorr(as.matrix(tempn[, -c(1:num_metadata_columns)]), type = "spearman")
  
  rhos <- results$r
  ps <- results$P
  ns <- resultsn$n
  
  ps_melt <- na.omit(melt(ps))
  ps_melt$qval <- p.adjust(ps_melt$value, method = "BH")
  names(ps_melt)[3] <- "pval"
  ps_sub <- subset(ps_melt, qval < q_cutoff)
  
  rhos_melt <- na.omit(melt(rhos))
  names(rhos_melt)[3] <- "rho"
  
  ns_melt <- melt(ns)
  names(ns_melt)[3] <- "n"
  
  merged <- merge(ps_sub, rhos_melt, by = c("Var1", "Var2"))
  
  if (treatments[i] == treatments[1]) {
    merged <- merge(merged, subset(ns_melt, n > n1), by = c("Var1", "Var2"))
  } else if (treatments[i] == treatments[2]) {
    merged <- merge(merged, subset(ns_melt, n > n2), by = c("Var1", "Var2"))
  } else if (treatments[i] == treatments[3]) {
    merged <- merge(merged, subset(ns_melt, n > n3), by = c("Var1", "Var2"))
  }
  
  if (nrow(merged) > 0) {
    merged$trt <- treatments[i]
    final_results <- rbind(final_results, merged)
  } else {
    print(paste("no correlations for", treatments[i]))
  }
  
  print(paste("finished", treatments[i]))
}

# Keep strong correlations only
strong_results <- subset(final_results, abs(rho) >= 0.9)

# ---------------------------
# Optional scatterplot for one ASV pair
# ---------------------------
# diet_ASVs <- subset(dataset, get(my_column) == treatments[1])
# ggplot(diet_ASVs, aes(x = ASV15, y = ASV5)) + geom_point()

# ---------------------------
# Taxonomy
# ---------------------------
taxonomy <- read_qza(taxonomy_file)
tax.clean <- parse_taxonomy(taxonomy$data)

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i,2] == "") {
    kingdom <- paste("unclassified_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == "") {
    phylum <- paste("unclassified_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == "") {
    class <- paste("unclassified_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == "") {
    order <- paste("unclassified_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == "") {
    family <- paste("unclassified_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == "") {
    tax.clean$Species[i] <- paste("unclassified_", tax.clean$Genus[i], sep = "_")
  }
}

strong_results_taxa <- merge(strong_results, ASVkey, by.x = "Var1", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, ASVkey, by.x = "Var2", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.x", by.y = 0)
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.y", by.y = 0)

write.csv(
  strong_results_taxa,
  paste0("output/cooccurrence/", age_group, "_strong_results_taxa.csv"),
  row.names = FALSE
)

write.csv(
  subset(strong_results_taxa, trt == "calf_starter"),
  paste0("output/cooccurrence/", age_group, "_strong_results_taxa_calf_starter.csv"),
  row.names = FALSE
)

write.csv(
  subset(strong_results_taxa, trt == "corn_silage"),
  paste0("output/cooccurrence/", age_group, "_strong_results_taxa_corn_silage.csv"),
  row.names = FALSE
)

write.csv(
  subset(strong_results_taxa, trt == "mix"),
  paste0("output/cooccurrence/", age_group, "_strong_results_taxa_mix.csv"),
  row.names = FALSE
)

library(dplyr)
getwd()

df <- read.csv("C:/ANSC516_Project/output/cooccurrence/2yr_strong_results_taxa.csv")

filtered <- df %>%
  filter(abs(correlation) >= 0.7 & p_value <= 0.05)

write.csv(filtered, "filtered_network.csv", row.names = FALSE)
