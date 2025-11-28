library("ggplot2")
library("dplyr")
library("DESeq2")
library("tximport")
library("tidyverse")
library("wesanderson")

# ============================================================================
# PARAMETERS - EDIT THIS SECTION TO CONFIGURE YOUR ANALYSIS
# ============================================================================

# Base directory
base_dir <- "/home/diegoj/bianca"

# Output directory
output_dir <- "/home/diegoj/bianca/rnaseq/full_run1_no_collapse/star_salmon/deseq2_qc"

# Design formula - using interaction term to allow all comparisons
design_formula <- ~ genotype + soil + genotype:soil

# Sample filtering - use "infected" from status column
sample_filter_status <- "infected"

# Gene filtering thresholds
min_count <- 1      # Minimum read count threshold
min_samples <- 7    # Minimum number of samples with >= min_count

# DEG thresholds
lfc_threshold <- 1      # Log2 fold change threshold
padj_threshold <- 0.05  # Adjusted p-value threshold

# PCA parameters
pca_top_genes <- 500  # Number of most variable genes for PCA

# Define contrasts to run
contrasts_to_run <- list(
  # Within genotype comparisons (soil effect)
#  list(
#    variable = "soil",
#    numerator = "clay",
#    denominator = "sandy",
#    genotype_subset = "5503",
#    name = "5503_clay_vs_sandy",
#    description = "5503: clay vs sandy"
#  ),
#  list(
#    variable = "soil",
#    numerator = "clay",
#    denominator = "sandy",
#    genotype_subset = "6007",
#    name = "6007_clay_vs_sandy",
#    description = "6007: clay vs sandy"
#  ),
#  
  # Within soil comparisons (genotype effect)
#  list(
#    variable = "genotype",
#    numerator = "5503",
#    denominator = "6007",
#    soil_subset = "clay",
#    name = "5503_vs_6007_clay",
#    description = "Clay: 5503 (susceptible) vs 6007 (resistant)"
#  ),
#  list(
#    variable = "genotype",
#    numerator = "5503",
#    denominator = "6007",
#    soil_subset = "sandy",
#    name = "5503_vs_6007_sandy",
#    description = "Sandy: 5503 (susceptible) vs 6007 (resistant)"
#  )
  
  # Overall main effects
  list(
    variable = "soil",
    numerator = "clay",
    denominator = "sandy",
    name = "overall_clay_vs_sandy",
    description = "Overall: clay vs sandy (across genotypes)"
  ),
  list(
    variable = "genotype",
    numerator = "5503",
    denominator = "6007",
    name = "overall_5503_vs_6007",
    description = "Overall: 5503 vs 6007 (across soil types)"
  )
)

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=================================================================\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS PIPELINE\n")
cat("=================================================================\n\n")

# ----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# ----------------------------------------------------------------------------

cat("[1/6] Loading metadata and count data...\n")

# Read metadata - filter for infected samples
metadata <- read.table(
  file.path(base_dir, "raw_reads/Info_pipeline.csv"),
  header = TRUE, 
  sep = ","
) %>%
  filter(status == sample_filter_status)

# Display sample distribution
cat("\nSample distribution:\n")
sample_table <- table(metadata$genotype, metadata$soil)
print(sample_table)
cat("\nTotal infected samples:", nrow(metadata), "\n")

# Load Salmon quantification files
sample_files <- file.path(
  base_dir, 
  "rnaseq/full_run1_no_collapse/star_salmon",
  metadata$sample, 
  "quant.sf"
)
names(sample_files) <- metadata$sample

# Load transcript-to-gene mapping
tx2gene <- read.table(
  file.path(base_dir, "rnaseq/full_run1_no_collapse/star_salmon/tx2gene.tsv"),
  header = TRUE
)

# Import count data
count_data <- tximport(
  files = sample_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE,
  ignoreAfterBar = TRUE
)

# Prepare metadata with proper factor levels
coldata <- metadata %>%
  mutate(
    genotype = factor(genotype, levels = c("6007", "5503")),  # 6007 as reference (resistant)
    soil = factor(soil, levels = c("sandy", "clay")),          # sandy as reference
    group = factor(paste(genotype, soil, sep = "_"))
  )

rownames(coldata) <- coldata$sample

# Verify alignment
stopifnot(all(colnames(count_data$counts) == rownames(coldata)))

cat("✓ Data loaded successfully\n")
cat(sprintf("  Samples: %d\n", ncol(count_data$counts)))
cat(sprintf("  Genes: %d\n", nrow(count_data$counts)))

# ----------------------------------------------------------------------------
# 2. CREATE DESEQ2 OBJECT AND FILTER GENES
# ----------------------------------------------------------------------------

cat("\n[2/6] Creating DESeq2 object and filtering genes...\n")

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(
  txi = count_data,
  colData = coldata,
  design = design_formula
)

cat(sprintf("  Initial genes: %d\n", nrow(dds)))

# Filter low-count genes
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds_filtered <- dds[keep, ]

cat(sprintf("  Genes after filtering: %d\n", nrow(dds_filtered)))
cat(sprintf("  Genes removed: %d (%.1f%%)\n", 
            sum(!keep), 100 * sum(!keep) / nrow(dds)))

# ----------------------------------------------------------------------------
# 3. RUN DESEQ2
# ----------------------------------------------------------------------------

cat("\n[3/6] Running DESeq2 analysis...\n")
dds_filtered <- DESeq(dds_filtered, parallel = TRUE)
cat("✓ DESeq2 analysis complete\n")

# ----------------------------------------------------------------------------
# 4. GENERATE QC PLOTS
# ----------------------------------------------------------------------------

cat("\n[4/6] Generating quality control plots...\n")

# VST transformation
vst <- varianceStabilizingTransformation(dds_filtered)

# PCA plot
colors <- wes_palette("Darjeeling1", 4, type = "discrete")
pca_data <- plotPCA(
  vst, 
  intgroup = c("genotype", "soil"),
  returnData = TRUE, 
  ntop = pca_top_genes
)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_data <- pca_data %>%
  dplyr::rename(Genotype = genotype, Soil = soil)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, shape = Genotype, color = Soil)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  scale_colour_manual(values = colors) +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 22),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

pca_file <- file.path(output_dir, "PCA_infected_samples.png")
ggsave(p, filename = pca_file, units = "cm", width = 19.5, height = 15, dpi = 320)
cat(sprintf("✓ PCA plot saved: %s\n", pca_file))

# ----------------------------------------------------------------------------
# 5. RUN DIFFERENTIAL EXPRESSION CONTRASTS
# ----------------------------------------------------------------------------

cat("\n[5/6] Running differential expression contrasts...\n")
cat(sprintf("  Number of contrasts: %d\n\n", length(contrasts_to_run)))

# Function to run a single contrast
run_contrast <- function(contrast_info, dds, output_dir, lfc_thresh, padj_thresh) {
  
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("Processing: %s\n", contrast_info$name))
  cat(sprintf("Description: %s\n", contrast_info$description))
  
  # Determine subset and appropriate design
  dds_subset <- dds
  needs_refit <- FALSE
  new_design <- NULL
  
  # Case 1: Subset by genotype (compare soils within one genotype)
  if (!is.null(contrast_info$genotype_subset)) {
    dds_subset <- dds[, dds$genotype == contrast_info$genotype_subset]
    cat(sprintf("Subsetting to genotype: %s\n", contrast_info$genotype_subset))
    new_design <- ~ soil
    needs_refit <- TRUE
  } 
  # Case 2: Subset by soil (compare genotypes within one soil)
  else if (!is.null(contrast_info$soil_subset)) {
    dds_subset <- dds[, dds$soil == contrast_info$soil_subset]
    cat(sprintf("Subsetting to soil: %s\n", contrast_info$soil_subset))
    new_design <- ~ genotype
    needs_refit <- TRUE
  }
  # Case 3: Overall comparison (use simpler additive model)
  else {
    cat("Using all samples\n")
    new_design <- ~ genotype + soil  # Additive model for main effects
    needs_refit <- TRUE
  }
  
  cat(sprintf("Samples in analysis: %d\n", ncol(dds_subset)))
  cat(sprintf("Design formula: %s\n", deparse(new_design)))
  
  # Refit model with appropriate design
  if (needs_refit) {
    cat("Refitting DESeq2 model...\n")
    design(dds_subset) <- new_design
    dds_subset <- DESeq(dds_subset, parallel = TRUE, quiet = TRUE)
  }
  
  # Extract results
  dea_contrast <- results(
    dds_subset,
    contrast = c(contrast_info$variable, 
                 contrast_info$numerator, 
                 contrast_info$denominator),
    alpha = padj_thresh,
    parallel = TRUE
  )
  
  # Calculate base means for each condition
  numerator_samples <- colData(dds_subset)[[contrast_info$variable]] == contrast_info$numerator
  denominator_samples <- colData(dds_subset)[[contrast_info$variable]] == contrast_info$denominator
  
  baseMean_numerator <- rowMeans(counts(dds_subset, normalized = TRUE)[, numerator_samples])
  baseMean_denominator <- rowMeans(counts(dds_subset, normalized = TRUE)[, denominator_samples])
  
  # Create results data frame
  res_df <- data.frame(
    gene_id = rownames(dea_contrast),
    baseMean_numerator = baseMean_numerator,
    baseMean_denominator = baseMean_denominator,
    as.data.frame(dea_contrast)
  ) %>%
    filter(complete.cases(.))
  
  # Rename baseMean columns
  colnames(res_df)[2:3] <- c(
    paste0("baseMean_", contrast_info$numerator),
    paste0("baseMean_", contrast_info$denominator)
  )
  
  # Filter for significant DEGs
  deg_all <- res_df %>%
    filter(abs(log2FoldChange) > lfc_thresh & padj < padj_thresh) %>%
    arrange(desc(abs(log2FoldChange)))
  
  deg_up <- deg_all %>% filter(log2FoldChange > lfc_thresh)
  deg_down <- deg_all %>% filter(log2FoldChange < -lfc_thresh)
  
  # Print summary
  cat(sprintf("\nResults:\n"))
  cat(sprintf("  Total genes tested: %d\n", nrow(res_df)))
  cat(sprintf("  Upregulated in %s: %d\n", contrast_info$numerator, nrow(deg_up)))
  cat(sprintf("  Downregulated in %s: %d\n", contrast_info$numerator, nrow(deg_down)))
  cat(sprintf("  Total DEGs: %d\n", nrow(deg_all)))
  
  # Save results
  contrast_dir <- file.path(output_dir, contrast_info$name)
  dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(res_df, 
            file.path(contrast_dir, paste0("all_genes_", contrast_info$name, ".csv")), 
            row.names = FALSE)
  write.csv(deg_all, 
            file.path(contrast_dir, paste0("DEGs_", contrast_info$name, ".csv")), 
            row.names = FALSE)
  write.csv(deg_up, 
            file.path(contrast_dir, paste0("DEGs_upregulated_", contrast_info$name, ".csv")), 
            row.names = FALSE)
  write.csv(deg_down, 
            file.path(contrast_dir, paste0("DEGs_downregulated_", contrast_info$name, ".csv")), 
            row.names = FALSE)
  
  cat(sprintf("✓ Results saved to: %s/\n", contrast_dir))
  
  return(list(
    contrast = contrast_info,
    results = res_df,
    degs = deg_all,
    up = deg_up,
    down = deg_down
  ))
}

# Run all contrasts
contrast_results <- lapply(contrasts_to_run, function(contrast) {
  tryCatch({
    run_contrast(contrast, dds_filtered, output_dir, lfc_threshold, padj_threshold)
  }, error = function(e) {
    cat(sprintf("✗ ERROR processing %s: %s\n", contrast$name, e$message))
    return(NULL)
  })
})

# ----------------------------------------------------------------------------
# 6. SAVE DESEQ2 OBJECTS AND SUMMARY
# ----------------------------------------------------------------------------

cat("\n[6/6] Saving DESeq2 objects and creating summary...\n")

save(dds_filtered, vst, coldata, contrast_results,
     file = file.path(output_dir, "deseq2_objects.RData"))

cat(sprintf("✓ DESeq2 objects saved\n"))

# Create summary table
summary_df <- data.frame(
  Contrast = character(),
  Description = character(),
  Upregulated = integer(),
  Downregulated = integer(),
  Total_DEGs = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(contrast_results)) {
  if (!is.null(contrast_results[[i]])) {
    summary_df <- rbind(summary_df, data.frame(
      Contrast = contrast_results[[i]]$contrast$name,
      Description = contrast_results[[i]]$contrast$description,
      Upregulated = nrow(contrast_results[[i]]$up),
      Downregulated = nrow(contrast_results[[i]]$down),
      Total_DEGs = nrow(contrast_results[[i]]$degs)
    ))
  }
}

write.csv(summary_df, 
          file.path(output_dir, "DEG_summary.csv"), 
          row.names = FALSE)

# ----------------------------------------------------------------------------
# FINAL SUMMARY
# ----------------------------------------------------------------------------

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=================================================================\n\n")
cat(sprintf("Output directory: %s\n\n", output_dir))

cat("Summary of DEG counts (|LFC|>1, padj<0.05):\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
print(summary_df, row.names = FALSE)

cat("\n=================================================================\n")
