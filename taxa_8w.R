library(qiime2R)
library(phyloseq)
library(tidyverse)
library(DESeq2)

setwd("C:/ANSC516_Project")

# ---------------------------
# Choose age group here
# ---------------------------
age_group <- "8w"   # change to "1yr" or "2yr" in the other scripts
feature_table_file <- "8w-core-metrics-results/rarefied_table.qza"

# Output folders
out_dir <- paste0("output_", age_group)
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "taxa"), showWarnings = FALSE)

# ---------------------------
# Read phyloseq object
# ---------------------------
physeq <- qza_to_phyloseq(
  features = feature_table_file,
  tree = "rooted-tree_simple.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.txt"
)

asv_table <- data.frame(otu_table(physeq), check.names = FALSE)
metadata <- data.frame(sample_data(physeq), check.names = FALSE)
taxonomy <- data.frame(tax_table(physeq), check.names = FALSE)

# Keep only this age group
metadata <- metadata %>%
  filter(Host_age == age_group)

# Make sure sample order matches
asv_table <- asv_table[, rownames(metadata), drop = FALSE]

# Fix factor order
metadata$Host_age <- factor(metadata$Host_age, levels = c("8w", "1yr", "2yr"))
metadata$calf_diet <- factor(metadata$calf_diet, levels = c("calf_starter", "corn_silage", "mix"))

# ---------------------------
# Clean taxonomy
# ---------------------------
tax.clean <- taxonomy
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i,2] == "") {
    kingdom <- paste("uncl_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == "") {
    phylum <- paste("uncl_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == "") {
    class <- paste("uncl_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == "") {
    order <- paste("uncl_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == "") {
    family <- paste("uncl_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == "") {
    tax.clean$Species[i] <- paste("uncl_", tax.clean$Genus[i], sep = "_")
  }
}

# ---------------------------
# Rebuild phyloseq object
# ---------------------------
OTU.physeq <- otu_table(as.matrix(asv_table), taxa_are_rows = TRUE)
tax.physeq <- tax_table(as.matrix(tax.clean))
meta.physeq <- sample_data(metadata)

physeq_age <- phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_age <- prune_taxa(taxa_sums(physeq_age) > 0, physeq_age)

# ---------------------------
# Colors
# ---------------------------
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

#################################################################
## Taxa barplot
#################################################################

my_level <- c("Phylum", "Family", "Genus")
my_column <- "calf_diet"
my_column_ordered <- c("calf_starter", "corn_silage", "mix")
abund_filter <- 0.02

for (ml in my_level) {
  
  taxa.summary <- physeq_age %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    group_by(.data[[my_column]], .data[[ml]]) %>%
    summarise(Abundance.average = mean(Abundance), .groups = "drop")
  
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>%
    group_by(.data[[ml]]) %>%
    summarise(overall.max = max(Abundance.average), .groups = "drop")
  
  physeq.taxa.mean <- taxa.summary %>%
    group_by(.data[[ml]]) %>%
    summarise(overall.mean = mean(Abundance.average), .groups = "drop")
  
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  physeq_meta <- merge(physeq_meta, physeq.taxa.mean)
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max > abund_filter)
  
  physeq_meta_filtered$diet_ordered <- factor(
    physeq_meta_filtered[[my_column]],
    levels = my_column_ordered
  )
  
  physeq_meta_filtered[[ml]] <- factor(physeq_meta_filtered[[ml]])
  y <- tapply(
    physeq_meta_filtered$overall.mean,
    physeq_meta_filtered[[ml]],
    function(z) max(z)
  )
  y <- sort(y, decreasing = TRUE)
  
  physeq_meta_filtered[[ml]] <- factor(
    as.character(physeq_meta_filtered[[ml]]),
    levels = names(y)
  )
  
  p <- ggplot(
    physeq_meta_filtered,
    aes(x = diet_ordered, y = Abundance.average, fill = .data[[ml]])
  ) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my_colors) +
    ylim(c(0, 1)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 0.5, keyheight = 0.5, ncol = 1)) +
    theme(legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab("Diet") +
    ggtitle(paste0(age_group, " ", ml, " (>", abund_filter * 100, "% in at least 1 sample)"))
  
  ggsave(
    file.path(out_dir, "taxa", paste0(age_group, "_", ml, "_BarPlot_calf_diet.png")),
    plot = p,
    height = 5,
    width = 4
  )
}

#################################################################
### Differential Abundance with DESeq2
#################################################################

physeq_otu_table <- data.frame(otu_table(physeq_age), check.names = FALSE)
OTU.clean2 <- physeq_otu_table + 1

OTU.physeq2 <- otu_table(as.matrix(OTU.clean2), taxa_are_rows = TRUE)
physeq_deseq <- phyloseq(OTU.physeq2, tax.physeq, meta.physeq)

diagdds <- phyloseq_to_deseq2(physeq_deseq, ~ calf_diet)
diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")

alpha <- 0.05

run_deseq2 <- function(x, y) {
  
  my_contrast <- c("calf_diet", x, y)
  res <- results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)
  
  sigtab <- res[which(res$padj < alpha & !is.na(res$padj)), ]
  sigtab <- cbind(
    as(sigtab, "data.frame"),
    as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix")
  )
  
  if (nrow(sigtab) == 0) {
    message(paste("No significant taxa for", x, "vs", y, "in", age_group))
    return(NULL)
  }
  
  theme_set(theme_bw())
  
  x_order <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(z) max(z))
  x_order <- sort(x_order, decreasing = TRUE)
  sigtab$Genus <- factor(as.character(sigtab$Genus), levels = names(x_order))
  
  DESeq_fig <- ggplot(sigtab, aes(x = Genus, y = log2FoldChange, color = Phylum)) +
    geom_point(size = 3) +
    ylab(paste0("(", x, "/", y, ")\nlog2FoldChange")) +
    scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    ggtitle(paste0(age_group, ": ", x, " vs ", y))
  
  ggsave(
    file.path(out_dir, "taxa", paste0("DESeq2_", age_group, "_", x, "_vs_", y, ".png")),
    DESeq_fig,
    height = 5,
    width = 10
  )
  
  write.csv(
    sigtab,
    file.path(out_dir, "taxa", paste0("DESeq2_", age_group, "_", x, "_vs_", y, ".csv")),
    row.names = TRUE
  )
}

run_deseq2("calf_starter", "corn_silage")
run_deseq2("calf_starter", "mix")
run_deseq2("corn_silage", "mix")
