#####---------------------
#Tryouts to understand the model outcomes 
## Author JD

#---------------------




# ==============================================================================
# SETUP & LIBRARIES
# ==============================================================================

library(corrplot)      # Correlation visualization
library(ggplot2)       # Plotting
library(dplyr)         # Data manipulation
library(tidyr)         # Data reshaping
library(patchwork)     # Combine plots
library(vegan)         # Ordination (NMDS, envfit)
library(ggrepel)       # Text labels without overlap

# Create output directory
output_dir <- "Lichens/data_audit/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ==============================================================================
# LOAD DATA
# ==============================================================================

load("Lichens/model_outputs/02c_models_complete_standard.RData")


#═══════════════════════════════════════════════════════════
#ACTION 1: Correlation Matrix
#═══════════════════════════════════════════════════════════

# Extract predictor data
predictor_data <- data_scaled[, predictors_scaled]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")

# Identify high correlations (|r| > 0.7)
high_cor <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)

if (nrow(high_cor) > 0) {
  cat("⚠️  HIGH CORRELATIONS DETECTED (|r| > 0.7):\n\n")
  
  high_cor_pairs <- data.frame(
    Var1 = rownames(cor_matrix)[high_cor[, 1]],
    Var2 = colnames(cor_matrix)[high_cor[, 2]],
    Correlation = cor_matrix[high_cor],
    stringsAsFactors = FALSE
  )
  
  # Remove duplicates (A-B same as B-A)
  high_cor_pairs <- high_cor_pairs[!duplicated(t(apply(high_cor_pairs[,1:2], 1, sort))), ]
  high_cor_pairs <- high_cor_pairs[order(abs(high_cor_pairs$Correlation), decreasing = TRUE), ]
  
  print(high_cor_pairs)
  cat("\n")
  
  write.csv(high_cor_pairs, paste0(output_dir, "high_correlations.csv"), row.names = FALSE)
} else {
  cat("✓ No problematic correlations detected (all |r| < 0.7)\n\n")
}

cat("═══════════════════════════════════════════════════════════\n")
cat("ACTION 1: Correlation Matrix\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Extract predictor data
predictor_data <- data_scaled[, predictors_scaled]

# Calculate correlation matrix
cor_matrix <- cor(predictor_data, use = "complete.obs")

# Identify high correlations (|r| > 0.7)
high_cor <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)

if (nrow(high_cor) > 0) {
  cat("⚠️  HIGH CORRELATIONS DETECTED (|r| > 0.7):\n\n")
  
  high_cor_pairs <- data.frame(
    Var1 = rownames(cor_matrix)[high_cor[, 1]],
    Var2 = colnames(cor_matrix)[high_cor[, 2]],
    Correlation = cor_matrix[high_cor],
    stringsAsFactors = FALSE
  )
  
  # Remove duplicates (A-B same as B-A)
  high_cor_pairs <- high_cor_pairs[!duplicated(t(apply(high_cor_pairs[,1:2], 1, sort))), ]
  high_cor_pairs <- high_cor_pairs[order(abs(high_cor_pairs$Correlation), decreasing = TRUE), ]
  
  print(high_cor_pairs)
  cat("\n")
  
  write.csv(high_cor_pairs, paste0(output_dir, "high_correlations.csv"), row.names = FALSE)
} else {
  cat("✓ No problematic correlations detected (all |r| < 0.7)\n\n")
}

# Create correlation plot
png(paste0(output_dir, "correlation_matrix.png"), width = 10, height = 10, units = "in", res = 300)

corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.8,
         addCoef.col = "black",
         number.cex = 0.6,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         title = "Predictor Correlation Matrix\n(Hierarchical clustering order)",
         mar = c(0,0,2,0))

# Add reference line for |r| = 0.7
legend("bottomright", 
       legend = c("|r| > 0.7 = Problematic", "|r| < 0.7 = Acceptable"),
       cex = 0.8, bty = "n")

dev.off()

cat("✓ Correlation matrix saved: correlation_matrix.png\n")
cat("✓ High correlation pairs saved: high_correlations.csv\n\n")

# ==============================================================================
# ACTION 2: PREDICTOR DISTRIBUTIONS BY LICHEN PRESENCE
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("ACTION 2: Predictor Distributions by Lichen Presence\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Function to create comparison plots
plot_predictor_distribution <- function(lichen_var, predictor_var, data_df) {
  
  # Prepare data
  plot_data <- data.frame(
    Presence = factor(data_df[[lichen_var]], 
                      levels = c(0, 1),
                      labels = c("Absent", "Present")),
    Value = data_df[[predictor_var]]
  )
  
  # Remove NA
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  # Calculate means
  means <- plot_data %>%
    group_by(Presence) %>%
    summarise(Mean = mean(Value), .groups = "drop")
  
  # T-test
  if (length(unique(plot_data$Presence)) == 2) {
    t_result <- t.test(Value ~ Presence, data = plot_data)
    p_val <- t_result$p.value
    p_label <- ifelse(p_val < 0.001, "p < 0.001",
                      ifelse(p_val < 0.01, "p < 0.01",
                             ifelse(p_val < 0.05, "p < 0.05",
                                    sprintf("p = %.3f", p_val))))
  } else {
    p_label <- "N/A"
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Presence, y = Value, fill = Presence)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 1) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
                 fill = "red", color = "black") +
    scale_fill_manual(values = c("Absent" = "gray70", "Present" = "steelblue")) +
    labs(title = gsub("_scaled", "", predictor_var),
         subtitle = p_label,
         x = NULL,
         y = "Scaled value") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9))
  
  return(p)
}

# Create plots for core_ogf_presence (main model)
cat("Creating distribution plots for core_ogf_presence...\n")

plot_list <- list()

for (pred in predictors_scaled) {
  p <- plot_predictor_distribution("core_ogf_presence", pred, data_scaled)
  plot_list[[pred]] <- p
}

# Combine into multi-panel figure
combined_plot <- wrap_plots(plot_list, ncol = 3)

ggsave(paste0(output_dir, "distributions_core_ogf.png"),
       combined_plot, width = 14, height = 12, dpi = 300)

cat("✓ Distribution plots saved: distributions_core_ogf.png\n\n")

# Create summary table of differences
cat("Summary of predictor differences (Present vs Absent):\n\n")

summary_stats <- data.frame()

for (pred in predictors_scaled) {
  
  present_vals <- data_scaled[data_scaled$core_ogf_presence == 1, pred]
  absent_vals <- data_scaled[data_scaled$core_ogf_presence == 0, pred]
  
  if (length(present_vals) > 0 && length(absent_vals) > 0) {
    t_result <- t.test(present_vals, absent_vals)
    
    summary_stats <- rbind(summary_stats, data.frame(
      Predictor = gsub("_scaled", "", pred),
      Mean_Absent = round(mean(absent_vals, na.rm = TRUE), 3),
      Mean_Present = round(mean(present_vals, na.rm = TRUE), 3),
      Difference = round(mean(present_vals, na.rm = TRUE) - mean(absent_vals, na.rm = TRUE), 3),
      T_statistic = round(t_result$statistic, 3),
      P_value = round(t_result$p.value, 4),
      Significant = t_result$p.value < 0.05,
      stringsAsFactors = FALSE
    ))
  }
}

summary_stats <- summary_stats[order(summary_stats$P_value), ]
print(summary_stats)

write.csv(summary_stats, paste0(output_dir, "predictor_differences_core_ogf.csv"), 
          row.names = FALSE)

cat("\n✓ Summary statistics saved: predictor_differences_core_ogf.csv\n\n")


# ==============================================================================
# ACTION 2b: RED FLAG VARIABLE - TREE HEIGHT MEDIAN
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("RED FLAG ANALYSIS: tree_height_median\n")
cat("════════════════════════════════════��══════════════════════\n\n")

# Detailed examination of tree_height_median
cat("Distribution of tree_height_median:\n")
cat(sprintf("  Min:    %.2f\n", min(data_scaled$tree_height_median_scaled, na.rm = TRUE)))
cat(sprintf("  Q1:     %.2f\n", quantile(data_scaled$tree_height_median_scaled, 0.25, na.rm = TRUE)))
cat(sprintf("  Median: %.2f\n", median(data_scaled$tree_height_median_scaled, na.rm = TRUE)))
cat(sprintf("  Q3:     %.2f\n", quantile(data_scaled$tree_height_median_scaled, 0.75, na.rm = TRUE)))
cat(sprintf("  Max:    %.2f\n\n", max(data_scaled$tree_height_median_scaled, na.rm = TRUE)))

# Correlation with other variables
tree_height_cors <- cor_matrix["tree_height_median_scaled", ]
tree_height_cors <- sort(abs(tree_height_cors), decreasing = TRUE)

cat("Correlations with tree_height_median:\n")
print(head(tree_height_cors[tree_height_cors < 1], 5))
cat("\n")

# If tree_height_SD or similar exists, compare
if ("tree_height_sd_scaled" %in% colnames(data_scaled)) {
  cat("��️  tree_height_sd also available - consider using instead of median!\n\n")
}


# ==============================================================================
# ACTION 3: SPECIES-LEVEL ORDINATION (if data available)
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("ACTION 3: Species-Level Ordination\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Check if species-level data exists
# You need to modify this based on your actual data structure
# Assuming you have columns with individual species presence/absence or abundance

# Identify species columns (this is a guess - adjust to your data!)
# Common patterns: species names, or columns ending in _sp, or in separate file

# Example 1: If species data is in the same dataframe
species_cols <- grep("^[A-Z][a-z]+_[a-z]+", names(data), value = TRUE)  # Genus_species pattern

# Example 2: If you have a separate species matrix file
# species_data <- read.csv("species_matrix.csv")

if (length(species_cols) > 0) {
  
  cat(sprintf("Found %d potential species columns\n", length(species_cols)))
  cat("First few:", paste(head(species_cols, 5), collapse = ", "), "\n\n")
  
  # Create species matrix
  species_matrix <- data[, species_cols]
  
  # Remove species with < 3 occurrences (too rare for ordination)
  species_freq <- colSums(species_matrix > 0)
  species_matrix <- species_matrix[, species_freq >= 3]
  
  cat(sprintf("Using %d species (≥3 occurrences)\n\n", ncol(species_matrix)))
  
  # NMDS ordination
  cat("Running NMDS ordination...\n")
  
  set.seed(123)  # For reproducibility
  nmds <- metaMDS(species_matrix, distance = "bray", k = 2, trymax = 100)
  
  cat(sprintf("  Stress: %.3f ", nmds$stress))
  if (nmds$stress < 0.1) {
    cat("(excellent)\n")
  } else if (nmds$stress < 0.2) {
    cat("(good)\n")
  } else {
    cat("(acceptable)\n")
  }
  cat("\n")
  
  # Environmental fit
  cat("Testing environmental variables...\n")
  env_fit <- envfit(nmds, predictor_data, permutations = 999)
  
  # Extract significant variables
  env_scores <- as.data.frame(scores(env_fit, display = "vectors"))
  env_pvals <- env_fit$vectors$pvals
  env_r2 <- env_fit$vectors$r
  
  env_summary <- data.frame(
    Variable = rownames(env_scores),
    NMDS1 = env_scores$NMDS1,
    NMDS2 = env_scores$NMDS2,
    R2 = round(env_r2, 3),
    P_value = round(env_pvals, 4),
    Significant = env_pvals < 0.05,
    stringsAsFactors = FALSE
  )
  
  env_summary <- env_summary[order(env_summary$P_value), ]
  
  cat("\nEnvironmental variable importance:\n")
  print(env_summary)
  cat("\n")
  
  write.csv(env_summary, paste0(output_dir, "ordination_envfit.csv"), 
            row.names = FALSE)
  
  # Create ordination plot
  png(paste0(output_dir, "NMDS_ordination.png"), 
      width = 10, height = 8, units = "in", res = 300)
  
  # Base plot
  plot(nmds, type = "n", main = "NMDS Ordination: Species Composition")
  
  # Add points colored by core_ogf_presence
  points(nmds, display = "sites", 
         pch = ifelse(data$core_ogf_presence == 1, 19, 1),
         col = ifelse(data$core_ogf_presence == 1, "blue", "gray"),
         cex = 1.2)
  
  # Add significant environmental vectors
  plot(env_fit, p.max = 0.05, col = "red", cex = 0.8)
  
  # Add legend
  legend("topright", 
         legend = c("Core OGF present", "Core OGF absent", "Env. variables (p<0.05)"),
         pch = c(19, 1, NA),
         col = c("blue", "gray", "red"),
         lty = c(NA, NA, 1),
         bty = "n")
  
  # Add stress value
  text(min(scores(nmds)[,1]), max(scores(nmds)[,2]), 
       paste("Stress =", round(nmds$stress, 3)), 
       pos = 4, cex = 0.9)
  
  dev.off()
  
  cat("✓ Ordination plot saved: NMDS_ordination.png\n")
  cat("✓ Environmental fit saved: ordination_envfit.csv\n\n")
  
} else {
  cat("⚠️  No species-level data detected in current workspace\n")
  cat("   To run ordination, you need:\n")
  cat("   • Species abundance/presence matrix (plots × species)\n")
  cat("   • Minimum 3 occurrences per species\n")
  cat("   • Load with: species_data <- read.csv('your_species_file.csv')\n\n")
  
  cat("   Example structure:\n")
  cat("      Plot_ID  Lecidea_fuscoatra  Mycoblastus_sanguinarius  ...\n")
  cat("      Plot_01              1                    0           ...\n")
  cat("      Plot_02              0                    1           ...\n\n")
}


# ==============================================================================
# SUMMARY & RECOMMENDATIONS
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  AUDIT SUMMARY & NEXT STEPS                                 ║\n")
cat("╚════════════════════════════════════════════���═════════════════╝\n\n")

cat("Files created in", output_dir, ":\n")
created_files <- list.files(output_dir)
for (f in created_files) {
  cat(sprintf("  ✓ %s\n", f))
}
cat("\n")

cat("REVIEW CHECKLIST:\n")
cat("─────────────────────────────────────────────────────────────\n")
cat("□ Check correlation_matrix.png for red cells (|r| > 0.7)\n")
cat("□ Review high_correlations.csv - which variables are redundant?\n")
cat("□ Examine distributions_core_ogf.png:\n")
cat("    • Do Present vs Absent differ as expected?\n")
cat("    • Any surprising patterns?\n")
cat("□ Review predictor_differences_core_ogf.csv:\n")
cat("    • Which predictors discriminate well? (low p-value)\n")
cat("    • Which show no difference? (high p-value)\n")
cat("□ If ordination ran: Check NMDS_ordination.png\n")
cat("    • Do current lichen groups separate in ordination space?\n")
cat("    • Which environmental variables are most important?\n\n")

cat("DECISION POINTS:\n")
cat("─────────────────────────────────────────────────────────────\n")
cat("1. Correlated predictors (|r| > 0.7):\n")
cat("   → Keep which one? Drop or combine others?\n\n")
cat("2. tree_height_median:\n")
cat("   → Does it discriminate core_ogf presence well?\n")
cat("   → Is it correlated with elevation/composition?\n")
cat("   → Consider replacing with height_SD or removing?\n\n")
cat("3. Non-discriminating predictors (high p-values):\n")
cat("   → Are these ecologically unimportant?\n")
cat("   → Or measurement issues?\n\n")
cat("4. Species composition:\n")
cat("   → Do grouped lichens cluster together in ordination?\n")
cat("   → Or are groups artificial?\n\n")

cat("═══════════════════════════════════════════════════════════\n")
cat("NEXT: Schedule meeting to discuss findings\n")
cat("═══���═══════════════════════════════════════════════════════\n\n")


# Save audit workspace
save.image(paste0(output_dir, "data_audit_workspace.RData"))
cat("✓ Workspace saved for later review\n\n")

################################################################################
# END OF DATA AUDIT SCRIPT
################################################################################
