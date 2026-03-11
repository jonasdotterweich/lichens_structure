################################################################################
# ŠUMAVA NATIONAL PARK LICHEN-STRUCTURE MODELING
# Script 04: Model Interpretation and Effect Plots
# Author: Jonas Dotterweich
# Date: 2026-02-14
# Purpose: Interpret coefficients and visualize lichen responses to predictors
################################################################################

# ==============================================================================
# CONTEXT: WHERE WE ARE
# ==============================================================================

# ALL 7 LICHEN MODELS ARE FINALIZED AND VALIDATED:
#   - 6 Binomial GLMs (presence/absence)
#   - 1 Negative Binomial GLM (calicioids richness)
#   - All diagnostics passed (no spatial autocorrelation, no overdispersion)
#   - !!!! However, watch out for mycoblastus model with high max coefficient (potential outlier influence) - interpret with caution
#
# NOW WE WILL:
#   1. Extract and summarize significant predictors
#   2. Create effect plots (marginal effects)
#   3. Compare patterns across lichen groups
#   4. Generate publication-ready figures
#   5. Identify key drivers of old-growth lichen occurrence

# ==============================================================================
# LOAD WORKSPACE AND LIBRARIES
# ==============================================================================

# Load finalized models from Script 02c
load("Lichens/model_outputs/02c_models_complete_standard.RData")

# Load additional libraries for visualization
library(tidyverse)
library(ggeffects)      # Marginal effects
library(sjPlot)         # Model tables
library(patchwork)      # Combine plots
library(scales)         # Formatting
library(broom)          # Tidy model outputs
library(broom.mixed)    # For glmmTMB models

cat("✅ Libraries loaded successfully\n\n")


# ==============================================================================
# SECTION 1: EXTRACT AND SUMMARIZE ALL MODEL COEFFICIENTS
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 1: COEFFICIENT SUMMARY ACROSS ALL MODELS              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Function to extract coefficients from models
extract_coefs <- function(model, model_name, model_type = "glm") {
  
  if (model_type == "glmmTMB") {
    # For negative binomial model from glmmTMB
    coef_table <- summary(model)$coefficients$cond
  } else {
    # For standard GLM
    coef_table <- summary(model)$coefficients
  }
  
  # Create tidy data frame
  coef_df <- data.frame(
    Model = model_name,
    Predictor = rownames(coef_table),
    Estimate = coef_table[, "Estimate"],
    SE = coef_table[, "Std. Error"],
    Z_value = coef_table[, grep("z value|t value", colnames(coef_table))],
    P_value = coef_table[, grep("Pr\\(", colnames(coef_table))],
    stringsAsFactors = FALSE
  )
  
  # Remove intercept
  coef_df <- coef_df[coef_df$Predictor != "(Intercept)", ]
  
  # Add significance stars
  coef_df$Sig <- cut(coef_df$P_value, 
                     breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                     labels = c("***", "**", "*", ".", "ns"))
  
  # Add direction
  coef_df$Direction <- ifelse(coef_df$Estimate > 0, "Positive", "Negative")
  
  # Add significance flag
  coef_df$Significant <- coef_df$P_value < 0.05
  
  return(coef_df)
}

# Extract coefficients from all models
all_coefs <- data.frame()

cat("Extracting coefficients from all models...\n\n")

# Binary models
for (lichen in lichen_binary) {
  coefs <- extract_coefs(models_glm[[lichen]], lichen, "glm")
  all_coefs <- rbind(all_coefs, coefs)
}

# Negative binomial model (calicioids)
coefs_nb <- extract_coefs(models_glm$calicioids_richness, "calicioids_richness", "glmmTMB")
all_coefs <- rbind(all_coefs, coefs_nb)

cat("✓ Coefficients extracted from 7 models\n\n")

# Clean up predictor names (remove "_scaled")
all_coefs$Predictor_clean <- gsub("_scaled", "", all_coefs$Predictor)

# Save full coefficient table
write.csv(all_coefs, "Lichens/interpret/04_all_model_coefficients.csv", row.names = FALSE)
cat("✓ Full coefficient table saved: 04_all_model_coefficients.csv\n\n")


# ==============================================================================
# SECTION 2: SUMMARY TABLE OF SIGNIFICANT PREDICTORS
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("SIGNIFICANT PREDICTORS BY LICHEN GROUP (p < 0.05)\n")
cat("═══════════════════════════════════════════════════��═══════\n\n")

for (lichen in all_lichen) {
  cat(sprintf("%-30s:\n", lichen))
  
  lichen_coefs <- all_coefs[all_coefs$Model == lichen & all_coefs$Significant, ]
  
  if (nrow(lichen_coefs) > 0) {
    # Sort by p-value
    lichen_coefs <- lichen_coefs[order(lichen_coefs$P_value), ]
    
    for (i in 1:nrow(lichen_coefs)) {
      direction_symbol <- ifelse(lichen_coefs$Direction[i] == "Positive", "↑", "↓")
      cat(sprintf("  %s %-25s %6.3f %3s (p = %.4f)\n",
                  direction_symbol,
                  lichen_coefs$Predictor_clean[i],
                  lichen_coefs$Estimate[i],
                  lichen_coefs$Sig[i],
                  lichen_coefs$P_value[i]))
    }
  } else {
    cat("  (no significant predictors)\n")
  }
  cat("\n")
}


# ==============================================================================
# SECTION 3: HEATMAP OF COEFFICIENTS ACROSS MODELS
# ==============================================================================

cat("\n╔══════════════════════════════════════════════��═══════════════╗\n")
cat("║  STEP 2: VISUALIZING PREDICTOR EFFECTS ACROSS MODELS       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Create wide-format matrix for heatmap
coef_matrix <- all_coefs %>%
  dplyr::select(Model, Predictor_clean, Estimate, P_value) %>%
  mutate(Model = factor(Model, levels = all_lichen)) %>%
  pivot_wider(names_from = Model, values_from = Estimate)

# Create significance matrix
sig_matrix <- all_coefs %>%
  dplyr::select(Model, Predictor_clean, P_value) %>%
  mutate(Significant = ifelse(P_value < 0.05, "*", "")) %>%
  dplyr::select(-P_value) %>%
  pivot_wider(names_from = Model, values_from = Significant)

# Merge for plotting
plot_data <- all_coefs %>%
  mutate(
    Model = factor(Model, levels = all_lichen),
    Predictor_clean = factor(Predictor_clean, levels = unique(Predictor_clean)),
    Sig_label = ifelse(P_value < 0.001, "***",
                       ifelse(P_value < 0.01, "**",
                              ifelse(P_value < 0.05, "*", "")))
  )

# Create coefficient heatmap
p_heatmap <- ggplot(plot_data, aes(x = Model, y = Predictor_clean, fill = Estimate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Sig_label), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                       midpoint = 0, 
                       limits = c(-max(abs(plot_data$Estimate)), 
                                  max(abs(plot_data$Estimate))),
                       name = "Coefficient\n(scaled)") +
  scale_x_discrete(position = "top") +
  labs(title = "Predictor Effects Across Lichen Groups",
       subtitle = "Significance: *** p<0.001, ** p<0.01, * p<0.05",
       x = NULL, y = "Predictor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 0.18, vjust = 0.28, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right",
        panel.grid = element_blank())

print(p_heatmap)
ggsave("Lichens/interpret/04_coefficient_heatmap.png", p_heatmap, width = 12, height = 8, dpi = 300)
cat("\n✓ Coefficient heatmap saved: 04_coefficient_heatmap.png\n\n")


# ==============================================================================
# SECTION 4: PREDICTOR IMPORTANCE SUMMARY
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("PREDICTOR IMPORTANCE: HOW MANY MODELS AFFECTED?\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Count significant associations per predictor
predictor_importance <- all_coefs %>%
  group_by(Predictor_clean) %>%
  summarize(
    N_models = n(),
    N_significant = sum(Significant),
    Pct_significant = round(N_significant / N_models * 100, 1),
    Avg_effect = mean(abs(Estimate)),
    Direction_consistency = sum(Estimate > 0) / N_models
  ) %>%
  arrange(desc(N_significant), desc(Avg_effect))

cat("Predictor ranking by number of significant associations:\n\n")
print(predictor_importance)
cat("\n")

# Save predictor importance table
write.csv(predictor_importance, "Lichens/interpret/04_predictor_importance.csv", row.names = FALSE)
cat("✓ Predictor importance table saved: 04_predictor_importance.csv\n\n")

# Identify key predictors
key_predictors <- predictor_importance$Predictor_clean[predictor_importance$N_significant >= 2]
cat(sprintf("Key predictors (significant in ≥2 models): %s\n\n", 
            paste(key_predictors, collapse = ", ")))


# ==============================================================================
# SECTION 5: MARGINAL EFFECT PLOTS FOR KEY PREDICTORS
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 3: MARGINAL EFFECTS PLOTS (KEY PREDICTORS)           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Focus on core_ogf_presence (main model) for detailed plots
cat("Creating detailed effect plots for core_ogf_presence...\n\n")

# Select top predictors for core_ogf
core_ogf_coefs <- all_coefs[all_coefs$Model == "core_ogf_presence", ]
core_ogf_sig <- core_ogf_coefs[core_ogf_coefs$Significant, ]
core_ogf_sig <- core_ogf_sig[order(-abs(core_ogf_sig$Estimate)), ]

if (nrow(core_ogf_sig) > 0) {
  
  cat(sprintf("Significant predictors for core_ogf_presence: %d\n", nrow(core_ogf_sig)))
  cat("  ", paste(core_ogf_sig$Predictor_clean, collapse = ", "), "\n\n")
  
  # Create effect plots for each significant predictor
  effect_plots <- list()
  
  for (i in 1:nrow(core_ogf_sig)) {
    pred <- core_ogf_sig$Predictor[i]
    pred_clean <- core_ogf_sig$Predictor_clean[i]
    
    cat(sprintf("  Creating effect plot for: %s\n", pred_clean))
    
    # Get marginal effects
    effects <- ggpredict(models_glm$core_ogf_presence, terms = pred)
    
    # Create plot
    p <- plot(effects) +
      labs(title = paste("Core OGF Presence:", pred_clean),
           x = paste(pred_clean, "(scaled)"),
           y = "Predicted Probability") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 11))
    
    effect_plots[[pred_clean]] <- p
  }
  
  # Combine plots
  if (length(effect_plots) > 0) {
    combined_plot <- wrap_plots(effect_plots, ncol = 2)
    ggsave("Lichens/interpret/04_core_ogf_effects.png", combined_plot, width = 10, height = 6 * ceiling(length(effect_plots)/2), dpi = 300)
    cat("\n✓ Combined effect plots saved: 04_core_ogf_effects.png\n\n")
  }
  
} else {
  cat("⚠️  No significant predictors for core_ogf_presence\n\n")
}


# ==============================================================================
# SECTION 6: DECAY STAGE COMPARISON (KEY ECOLOGICAL QUESTION)
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 4: DECAY STAGE EFFECTS (decay2 vs decay3)            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Extract decay2 and decay3 effects
decay_effects <- all_coefs %>%
  filter(Predictor_clean %in% c("decay2", "decay3")) %>%
  dplyr::select(Model, Predictor_clean, Estimate, SE, P_value, Sig, Significant)

cat("Decay stage effects across all lichen groups:\n\n")
print(decay_effects)
cat("\n")

# Create comparison plot
p_decay <- ggplot(decay_effects, 
                  aes(x = Model, y = Estimate, fill = Predictor_clean)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE),
                position = position_dodge(width = 0.8), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = Sig, y = Estimate + sign(Estimate) * 2 * SE),
            position = position_dodge(width = 0.8), size = 5) +
  scale_fill_manual(values = c("decay2" = "#8da0cb", "decay3" = "#fc8d62"),
                    labels = c("Decay stage 2 (early)", "Decay stage 3 (mid)"),
                    name = "Deadwood decay stage") +
  labs(title = "Effect of Deadwood Decay Stages on Lichen Occurrence",
       subtitle = "Coefficients from binomial/negative binomial GLMs (± 95% CI)",
       x = "Lichen Group",
       y = "Coefficient (scaled)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

print(p_decay)
ggsave("Lichens/interpret/04_decay_stage_comparison.png", p_decay, width = 10, height = 6, dpi = 300)
cat("\n✓ Decay stage comparison plot saved: 04_decay_stage_comparison.png\n\n")

# Summary statistics
cat("═══════════════════════════════════════════════════════════\n")
cat("DECAY STAGE SUMMARY:\n")
cat("═══════════════════════════════════════════════════════════\n\n")

decay2_sig <- sum(decay_effects$Predictor_clean == "decay2" & decay_effects$Significant)
decay3_sig <- sum(decay_effects$Predictor_clean == "decay3" & decay_effects$Significant)

cat(sprintf("Decay stage 2: Significant in %d/%d models\n", decay2_sig, 7))
cat(sprintf("Decay stage 3: Significant in %d/%d models\n", decay3_sig, 7))
cat("\n")


# ==============================================================================
# SECTION 7: ELEVATION GRADIENT EFFECTS
# ==============================================================================

cat("\n╔═════════════════════════════════════════════════════════════��╗\n")
cat("║  STEP 5: ELEVATION GRADIENT EFFECTS                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Extract elevation effects
elev_effects <- all_coefs %>%
  filter(Predictor_clean == "elevation") %>%
  dplyr::select(Model, Estimate, SE, P_value, Sig, Significant) %>%
  mutate(Model = factor(Model, levels = all_lichen))

cat("Elevation effects across all lichen groups:\n\n")
print(elev_effects)
cat("\n")

# Create elevation effect plot
p_elevation <- ggplot(elev_effects, 
                      aes(x = Model, y = Estimate, fill = Significant)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE),
                width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = Sig, y = Estimate + sign(Estimate) * 2 * SE),
            size = 5) +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70"),
                    name = "Significant (p<0.05)") +
  labs(title = "Elevation Effects on Lichen Occurrence and Richness",
       subtitle = "Coefficients from GLMs (± 95% CI)",
       x = "Lichen Group",
       y = "Coefficient (scaled)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

print(p_elevation)
ggsave("Lichens/interpret/04_elevation_effects.png", p_elevation, width = 10, height = 6, dpi = 300)
cat("\n✓ Elevation effect plot saved: 04_elevation_effects.png\n\n")


# ==============================================================================
# SECTION 8: CREATE MULTI-PANEL EFFECT PLOT FOR PUBLICATION
# ==============================================================================

cat("\n╔═══════════════���══════════════════════════════════════════════╗\n")
cat("║  STEP 6: PUBLICATION-READY MULTI-PANEL FIGURE              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Creating multi-panel figure with key predictors...\n\n")

# Select top 4 most important predictors overall
top_predictors <- predictor_importance$Predictor_clean[1:min(4, nrow(predictor_importance))]
cat("Top predictors for visualization:", paste(top_predictors, collapse = ", "), "\n\n")

# Create effect plots for core_ogf_presence with top predictors
pub_plots <- list()

for (pred_clean in top_predictors) {
  pred_scaled <- paste0(pred_clean, "_scaled")
  
  # Check if predictor exists in model
  if (pred_scaled %in% names(coef(models_glm$core_ogf_presence))) {
    
    cat(sprintf("  Panel: %s\n", pred_clean))
    
    # Get marginal effects
    effects <- ggpredict(models_glm$core_ogf_presence, terms = pred_scaled)
    
    # Get coefficient info
    coef_info <- all_coefs[all_coefs$Model == "core_ogf_presence" & 
                             all_coefs$Predictor_clean == pred_clean, ]
    
    # Create plot with consistent style
    p <- ggplot(effects, aes(x = x, y = predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                  fill = "steelblue", alpha = 0.3) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_rug(data = data_scaled, aes(x = get(pred_scaled), y = NULL),
               sides = "b", alpha = 0.3, inherit.aes = FALSE) +
      labs(title = pred_clean,
           subtitle = sprintf("β = %.2f %s", coef_info$Estimate, coef_info$Sig),
           x = paste(pred_clean, "(standardized)"),
           y = "P(Core OGF presence)") +
      ylim(0, 1) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9),
            axis.title = element_text(size = 9))
    
    pub_plots[[pred_clean]] <- p
  }
}

# Combine into publication figure
if (length(pub_plots) > 0) {
  pub_figure <- wrap_plots(pub_plots, ncol = 2) +
    plot_annotation(
      title = "Structural Drivers of Core Old-Growth Forest Lichen Presence",
      subtitle = "Marginal effects from binomial GLM (n=120 plots, Šumava NP)",
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 10))
    )
  
  print(pub_figure)
  ggsave("Lichens/interpret/04_figure_total_effects.png", pub_figure, width = 10, height = 8, dpi = 300)
  cat("\n✓ Publication figure saved: 04_publication_figure_effects.png\n\n")
}


# ==============================================================================
# SECTION 9: MODEL SUMMARY TABLE FOR PUBLICATION
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 7: PUBLICATION-READY SUMMARY TABLES                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Create formatted coefficient table
coef_table_pub <- all_coefs %>%
  mutate(
    Estimate_formatted = sprintf("%.3f", Estimate),
    SE_formatted = sprintf("%.3f", SE),
    Coef_CI = paste0(Estimate_formatted, " (", SE_formatted, ")"),
    P_formatted = ifelse(P_value < 0.001, "<0.001",
                         sprintf("%.3f", P_value))
  ) %>%
  select(Model, Predictor_clean, Coef_CI, P_formatted, Sig) %>%
  pivot_wider(names_from = Model, 
              values_from = c(Coef_CI, P_formatted),
              names_glue = "{Model}_{.value}")

write.csv(coef_table_pub, "04_coefficient_table_publication.csv", row.names = FALSE)
cat("✓ Publication coefficient table saved: 04_coefficient_table_publication.csv\n\n")

# Create HTML table using sjPlot
cat("Creating HTML model summary tables...\n\n")

# Binary models
tab_model(
  models_glm$core_ogf_presence,
  models_glm$mycoblastus_presence,
  models_glm$ochrolechia_presence,
  dv.labels = c("Core OGF", "Mycoblastus", "Ochrolechia"),
  show.ci = TRUE,
  show.se = TRUE,
  p.style = "stars",
  file = "04_model_table_binary_1.html"
)

tab_model(
  models_glm$xylographa_presence,
  models_glm$parmelia_agg_presence,
  models_glm$elite_rare_presence,
  dv.labels = c("Xylographa", "Parmelia agg.", "Elite rare"),
  show.ci = TRUE,
  show.se = TRUE,
  p.style = "stars",
  file = "04_model_table_binary_2.html"
)

cat("✓ HTML model tables saved:\n")
cat("    - 04_model_table_binary_1.html\n")
cat("    - 04_model_table_binary_2.html\n\n")


# ==============================================================================
# SECTION 10: KEY FINDINGS SUMMARY
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  INTERPRETATION COMPLETE: KEY FINDINGS                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("═══════════════════════════════════════════════════════════\n")
cat("KEY ECOLOGICAL FINDINGS:\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Count significant predictors per model
sig_summary <- all_coefs %>%
  group_by(Model) %>%
  summarize(
    N_predictors = n(),
    N_significant = sum(Significant),
    Significant_predictors = paste(Predictor_clean[Significant], collapse = ", ")
  )

for (i in 1:nrow(sig_summary)) {
  cat(sprintf("%s:\n", sig_summary$Model[i]))
  cat(sprintf("  - %d/%d predictors significant\n", 
              sig_summary$N_significant[i], sig_summary$N_predictors[i]))
  if (sig_summary$N_significant[i] > 0) {
    cat(sprintf("  - Key drivers: %s\n", sig_summary$Significant_predictors[i]))
  }
  cat("\n")
}

cat("═══════════════════════════════════════════════════════════\n")
cat("MOST IMPORTANT PREDICTORS (across all models):\n")
cat("═══════════════════════════════════════════════════════════\n\n")

print(predictor_importance[1:5, c("Predictor_clean", "N_significant", "Pct_significant", "Avg_effect")])
cat("\n")

cat("═══════════════════════════════════════════════════════════\n")
cat("FILES GENERATED:\n")
cat("══════���════════════════════════════════════════════════════\n\n")

output_files <- c(
  "04_all_model_coefficients.csv",
  "04_predictor_importance.csv",
  "04_coefficient_heatmap.png",
  "04_core_ogf_effects.png",
  "04_decay_stage_comparison.png",
  "04_elevation_effects.png",
  "04_publication_figure_effects.png",
  "04_coefficient_table_publication.csv",
  "04_model_table_binary_1.html",
  "04_model_table_binary_2.html"
)

for (file in output_files) {
  if (file.exists(file)) {
    cat(sprintf("  ✓ %s\n", file))
  }
}

cat("\n")

cat("═══════════════════════════════════════════════════════════\n")
cat("NEXT STEPS:\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("  → Review effect plots and identify key patterns\n")
cat("  → (Optional) Model selection to simplify predictor sets\n")
cat("  → (Optional) Threshold detection with segmented regression\n")
cat("  → Manuscript preparation\n")
cat("═══════════════════════════════════════════════════════════\n\n")


# Save workspace with all results
save.image("Lichens/to_model/04_interpretation_complete.RData")
cat("✅ Complete workspace saved: 04_interpretation_complete.RData\n\n")


################################################################################
# END OF SCRIPT 04
################################################################################