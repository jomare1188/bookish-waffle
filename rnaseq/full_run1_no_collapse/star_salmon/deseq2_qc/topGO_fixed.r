# =============================================================================
# GO over-representation analysis (topGO) -- corrected version
#
# Replaces topGO.r. Fixes, relative to that script:
#   1. BH correction is now computed over ALL tested GO terms, not over the
#      subset that already passed p < 0.05. In the old script the FDR filter
#      was vacuous (100% of nominally significant terms "passed" FDR).
#   2. p-values are taken from score(result) as NUMERIC, instead of from the
#      character column returned by GenTable(). The old `filter(table,
#      Classic < 0.05)` was a lexicographic STRING comparison, which silently
#      discarded every p-value printed in scientific notation -- i.e. the 31
#      most significant terms of the up-regulated set.
#   3. nodeSize = 5 drops GO terms with almost no annotated genes.
#   4. numChar = 1000 so GO term names are not truncated to 40 characters.
#   5. up/down handled by a loop, so output names can no longer disagree with
#      the gene set that was actually tested.
#   6. weight01 (topology-aware) is reported alongside classic. Per the topGO
#      vignette, weight01 p-values are NOT FDR-corrected: they are already
#      adjusted for the dependency between GO terms.
# =============================================================================

suppressPackageStartupMessages({
  library(topGO)
  library(ggplot2)
})

# ---- configuration ----------------------------------------------------------
go_table   <- "/dados04/jorge/bianca/references/PEDRO_genome/annotation/gene_go_table.txt"
contrast   <- "overall_clay_vs_sandy"
base_dir   <- file.path("/dados04/jorge/bianca/rnaseq/full_run1_no_collapse",
                        "star_salmon/deseq2_qc", contrast)
out_dir    <- file.path(base_dir, "topgo_corrected")

ontology   <- "BP"
node_size  <- 5      # ignore GO terms with < 5 annotated genes
alpha      <- 0.05   # FDR threshold
ntop       <- 20     # terms shown in the figure

dir.create(out_dir, showWarnings = FALSE)

# ---- gene -> GO map ---------------------------------------------------------
GO <- read.table(go_table, header = TRUE, stringsAsFactors = FALSE,
                 colClasses = c("character", "character"))
colnames(GO) <- c("Gene", "GO_term")
GO <- unique(GO)                                    # drop duplicated annotations

gene2GO   <- split(GO$GO_term, GO$Gene)
geneNames <- names(gene2GO)

# ---- plotting ---------------------------------------------------------------
# Sequential ramp for fold enrichment; point area encodes the number of DE genes
# in the term. White ring separates overlapping points.
ink        <- "#0b0b0b"
ink_soft   <- "#52514e"

# Scientific colour map 'batlow' (Crameri 2018, doi:10.5281/zenodo.1243862),
# sampled with scico::scico() and truncated so the palest step still clears
# ~2.8:1 contrast on white and the darkest stops short of black (~14:1).
# Perceptually uniform, colour-vision-deficiency safe, and monotonic in
# lightness so it survives greyscale printing. Embedded as hex rather than
# calling scico() so the script has no dependency beyond topGO + ggplot2.
batlow_ramp <- c("#C6903A", "#94862B", "#687A3D", "#406E53",
                 "#205E61", "#134961", "#0A2C5C")
ramp_mid    <- "#406E53"   # size-legend keys (fill is unmapped in that legend)

plot_go <- function(df, title, xlab) {
  df <- df[order(df$p_plot), ]
  df <- head(df, ntop)
  df$Term <- factor(df$Term, levels = rev(df$Term))

  ggplot(df, aes(x = Term, y = -log10(p_plot))) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed",
               colour = ink_soft, linewidth = 0.4) +
    geom_point(aes(size = Significant, fill = FoldEnrichment),
               shape = 21, colour = "white", stroke = 0.8) +
    scale_size(range = c(3, 11), name = "DE genes") +
    scale_fill_gradientn(colours = batlow_ramp, name = "Fold\nenrichment") +
    # shape 21 has no default fill, so the size-legend keys would render as an
    # invisible white ring on a white key -- give them an explicit fill
    guides(size = guide_legend(override.aes = list(fill = ramp_mid,
                                                   colour = "white"))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    labs(title = title, x = NULL, y = xlab) +
    coord_flip() +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "grey80"),
      # align titles to the whole plot, not the panel: long GO term names push
      # the panel to the right and a panel-centred title overflows the canvas
      plot.title.position = "plot",
      plot.title         = element_text(colour = ink, face = "bold", hjust = 0),
      axis.text          = element_text(colour = ink),
      axis.title         = element_text(colour = ink_soft),
      legend.title       = element_text(colour = ink_soft, size = 12),
      legend.text        = element_text(colour = ink_soft)
    )
}

save_plot <- function(p, stem) {
  for (dev in c("pdf", "svg", "png")) {
    ggsave(file.path(out_dir, paste0(stem, ".", dev)), plot = p, device = dev,
           width = 34, height = 22, units = "cm", dpi = 300)
  }
}

# ---- run one direction ------------------------------------------------------
run_direction <- function(direction) {

  deg_file <- file.path(base_dir, sprintf("DEGs_%sregulated_%s.csv",
                                          direction, contrast))
  interesting <- read.csv(deg_file, stringsAsFactors = FALSE)$gene_id

  message("\n=== ", direction, "-regulated ===")
  message("DE genes: ", length(interesting),
          " | with GO annotation: ", sum(interesting %in% geneNames))

  geneList <- factor(as.integer(geneNames %in% interesting), levels = c(0, 1))
  names(geneList) <- geneNames

  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = gene2GO, nodeSize = node_size)

  n_all <- numGenes(GOdata)
  n_sig <- length(sigGenes(GOdata))
  message("universe (feasible genes): ", n_all, " | DE genes tested: ", n_sig)

  if (n_sig < 2) {
    message("!! fewer than 2 annotated DE genes -- nothing to test, skipping.")
    return(invisible(NULL))
  }

  res_classic  <- runTest(GOdata, algorithm = "classic",  statistic = "fisher")
  res_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  # NUMERIC p-values for every tested term -- never parse GenTable's strings
  p_classic  <- score(res_classic)
  p_weight01 <- score(res_weight01)[names(p_classic)]

  # BH over ALL tested terms (this is the fix)
  fdr_classic <- p.adjust(p_classic, method = "BH")

  message("GO terms tested: ", length(p_classic))
  message("raw p < 0.05: ", sum(p_classic < 0.05),
          " | BH < ", alpha, ": ", sum(fdr_classic < alpha),
          " | weight01 p < 0.05: ", sum(p_weight01 < 0.05))

  # term names / counts, full-length
  tab <- GenTable(GOdata, Classic = res_classic, topNodes = length(p_classic),
                  orderBy = "Classic", numChar = 1000)
  tab <- tab[match(names(p_classic), tab$GO.ID), ]

  results <- data.frame(
    GO.ID           = names(p_classic),
    Term            = tab$Term,
    Annotated       = tab$Annotated,
    Significant     = tab$Significant,
    Expected        = round(tab$Annotated * n_sig / n_all, 2),
    FoldEnrichment  = round((tab$Significant / n_sig) /
                            (tab$Annotated   / n_all), 2),
    p_classic       = p_classic,
    fdr_classic_BH  = fdr_classic,
    p_weight01      = p_weight01,
    stringsAsFactors = FALSE
  )
  results <- results[order(results$p_classic), ]

  # NB: keep write.csv's default quote = TRUE. Some GO term names contain commas
  # ("regulation of vacuole fusion, non-autophagic"), which shift every column
  # after Term if the text is written unquoted.
  write.csv(results, file.path(out_dir, sprintf("GO_%s_all_tested_terms.csv", direction)),
            row.names = FALSE)

  sig <- results[results$fdr_classic_BH < alpha, ]
  write.csv(sig, file.path(out_dir, sprintf("GO_%s_FDR%g.csv", direction, alpha)),
            row.names = FALSE)

  sig_w <- results[results$p_weight01 < alpha, ]
  sig_w <- sig_w[order(sig_w$p_weight01), ]
  write.csv(sig_w, file.path(out_dir, sprintf("GO_%s_weight01_p%g.csv", direction, alpha)),
            row.names = FALSE)

  # ---- figures --------------------------------------------------------------
  lbl <- if (direction == "up") "Up-regulated" else "Down-regulated"

  if (nrow(sig) > 0) {
    d <- sig; d$p_plot <- d$fdr_classic_BH
    p <- plot_go(d,
      title    = sprintf("GO biological processes - %s (clay vs sandy)", lbl),
      xlab     = expression(-log[10]~"(FDR)"))
    save_plot(p, sprintf("GO_%s_classic_BH", direction))
  } else {
    message("!! no term passes FDR < ", alpha,
            " -- plotting top terms by RAW p instead, clearly labelled.")
    d <- head(results, ntop); d$p_plot <- d$p_classic
    p <- plot_go(d,
      title    = sprintf("GO biological processes - %s (clay vs sandy)", lbl),
      xlab     = expression(-log[10]~"(uncorrected p)"))
    save_plot(p, sprintf("GO_%s_classic_rawp_NOT_significant", direction))
  }

  if (nrow(sig_w) > 0) {
    d <- sig_w; d$p_plot <- d$p_weight01
    p <- plot_go(d,
      title    = sprintf("GO biological processes - %s (clay vs sandy)", lbl),
      xlab     = expression(-log[10]~"(weight01 p)"))
    save_plot(p, sprintf("GO_%s_weight01", direction))
  }

  invisible(results)
}

for (d in c("up", "down")) run_direction(d)

message("\nDone. Output written to: ", out_dir)
