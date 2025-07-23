#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2); library(dplyr); library(cowplot)

# Set input/output paths
input_root <- "02.eRNA_quantification/03.count_normalized_with_RPKM"
output_dir <- "03.eRNA_visualization"; dir.create(output_dir, FALSE, TRUE)

# Define input filenames for intergenic and intragenic enhancers
files <- list(
  inter = c("ES_non_active_E_intergenic_RPKM.txt", "ES_active_E_intergenic_RPKM.txt"),
  intra = c("Strand_ES_non_active_E_intragenic_AS_RPKM.txt", "Strand_ES_active_E_intragenic_AS_RPKM.txt")
)

# Function to load and format RPKM data
load_rpkm <- function(p, l) {
  stopifnot(file.exists(p))
  df <- read.table(p, header = FALSE)[, 1:5]
  colnames(df) <- c("region", "count", "length", "RPK", "RPKM")
  df %>% mutate(log2RPKM = log2(RPKM + 1), Group = factor(paste0(l, " E"), levels = c("non-active E", "active E")))
}

# Function to generate plot and statistics
make_plot <- function(df, label) {
  pval <- wilcox.test(log2RPKM ~ Group, df)$p.value
  signif <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns")))
  pval_str <- if (pval < 2.2e-16) "< 2.2e-16" else sprintf("%.2e", pval)
  y_max <- max(df$log2RPKM, na.rm = TRUE); y_pos <- y_max * 1.05
  p <- ggplot(df, aes(Group, log2RPKM, fill = Group)) +
    geom_violin(trim = FALSE, width = 1) +
    geom_boxplot(width = 0.05, outlier.shape = NA, color = "black", fill = "white") +
    annotate("text", x = 1.5, y = y_pos, label = signif, size = 5) +
    scale_fill_manual(values = c("non-active E" = "#F8F6BC", "active E" = "#ECA5A4")) +
    labs(title = paste(label, "Enhancer RPKM"), y = "log2(RPKM + 1)", x = NULL) +
    theme_minimal(base_size = 14) + theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "none"
    )
  stat <- data.frame(Type = label, Group1 = "non-active", Group2 = "active", P_value = pval_str, Signif = signif)
  list(plot = p, stat = stat)
}

# Process both intergenic and intragenic enhancer types
results <- lapply(names(files), function(type) {
  df <- bind_rows(
    load_rpkm(file.path(input_root, type, files[[type]][1]), "non-active"),
    load_rpkm(file.path(input_root, type, files[[type]][2]), "active")
  )
  make_plot(df, paste0(toupper(substring(type, 1, 1)), substring(type, 2)))
}); names(results) <- names(files)

# Combine and save violin plots
ggsave(file.path(output_dir, "eRNA_RPKM_Inter_Intra_violin.pdf"),
       plot_grid(results$inter$plot, results$intra$plot, ncol = 1, align = "v"),
       width = 4, height = 8)

# Save p-value summary table
write.table(bind_rows(results$inter$stat, results$intra$stat),
            file.path(output_dir, "eRNA_RPKM_Inter_Intra_violin_pvalue.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("All processes completed")
