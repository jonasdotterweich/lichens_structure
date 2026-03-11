################################################################################
# ŠUMAVA NATIONAL PARK LICHEN-STRUCTURE MODELING
# Script 04b: Model Interpretation and Effect Plots
# Author: Jonas Dotterweich
# Date: 2026-02-23
# Purpose: Extract coefficients, visualize effects, summarize results
################################################################################

# ==============================================================================
# SETUP
# ==============================================================================

# Resolve namespace conflicts
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggeffects)
library(patchwork)
library(glmmTMB)

# Load results
load("Lichens/model_outputs/02c_models_complete_standard.RData")

load("Lichens/model_outputs/reduced/02c_models_complete_reduced.RData")

# Create output directory

interp_dir <-  "Lichens/interpret/"


# ==============================================================================
# SECTION 1: EXTRACT COEFFICIENTS
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("EXTRACTING COEFFICIENTS\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Extraction function
extract_coefs <- function(model, model_name) {
  if (inherits(model, "glmmTMB")) {
    coef_table <- summary(model)$coefficients$cond
  } else {
    coef_table <- summary(model)$coefficients
  }
  
  df <- data.frame(
    Model = model_name,
    Predictor = rownames(coef_table),
    Estimate = coef_table[, 1],
    Std_Error = coef_table[, 2],
    Z_value = coef_table[, 3],
    P_value = coef_table[, 4],
    stringsAsFactors = FALSE
  )
  
  df$Predictor_clean <- gsub("_scaled", "", df$Predictor)
  df$Significant <- df$P_value < 0.05
  df$Sig <- cut(df$P_value, 
                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                labels = c("***", "**", "*", "ns"))
  
  return(df)
}

# Extract all
all_coefs <- data.frame()
for (lichen in names(models_glm)) {
  all_coefs <- rbind(all_coefs, extract_coefs(models_glm[[lichen]], lichen))
}

# Remove intercepts
all_coefs <- all_coefs %>% dplyr::filter(Predictor != "(Intercept)")

cat(sprintf("Extracted: %d coefficients from %d models\n\n", 
            nrow(all_coefs), length(names(models_glm))))


# ==============================================================================
# SECTION 2: MODEL SUMMARY TABLE
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("MODEL SUMMARY\n")
cat("═══════════════════════════════════════════════════════════\n\n")

model_summary <- all_coefs %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    N_predictors = dplyr::n(),              # Count from extracted coefficients
    N_significant = sum(Significant),
    .groups = "drop"
  )

# Add model statistics
for (i in 1:nrow(model_summary)) {
  lichen <- model_summary$Model[i]
  model <- models_glm[[lichen]]
  
  model_summary$Family[i] <- ifelse(inherits(model, "glmmTMB"), 
                                    "Negative Binomial", "Binomial")
  model_summary$N[i] <- nrow(data)
  model_summary$AIC[i] <- round(AIC(model), 1)
  
  if (!inherits(model, "glmmTMB")) {
    model_summary$Pseudo_R2[i] <- round(1 - (model$deviance / model$null.deviance), 3)
  } else {
    model_summary$Pseudo_R2[i] <- NA
  }
}

# Reorder columns
model_summary <- model_summary[, c("Model", "Family", "N", "N_predictors", 
                                   "N_significant", "AIC", "Pseudo_R2")]

print(model_summary)
write.csv(model_summary, paste0(interp_dir, "04b_model_summary.csv"), row.names = FALSE)
cat("\n")


# ==============================================================================
# SECTION 3: VARIABLE IMPORTANCE
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("VARIABLE IMPORTANCE\n")
cat("═══════════════════════════════════════════════════════════\n\n")

var_importance <- all_coefs %>%
  dplyr::group_by(Predictor_clean) %>%
  dplyr::summarise(
    N_models = dplyr::n(),
    N_significant = sum(Significant),
    Pct_significant = round((N_significant / N_models) * 100, 1),
    Mean_abs_coef = round(mean(abs(Estimate)), 3),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(N_significant), desc(Mean_abs_coef))

print(var_importance)
write.csv(var_importance, paste0(interp_dir, "04b_variable_importance.csv"), row.names = FALSE)
cat("\n")


# ==============================================================================
# SECTION 4: COEFFICIENT HEATMAP
# ==============================================================================

cat("Creating coefficient heatmap...\n")

plot_data <- all_coefs %>%
  dplyr::mutate(
    Model = factor(Model, levels = names(models_glm)),
    Predictor_clean = factor(Predictor_clean, levels = rev(unique(Predictor_clean)))
  )

p_heatmap <- ggplot(plot_data, aes(x = Model, y = Predictor_clean, fill = Estimate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Significant, "*", "")), size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Coefficient") +
  labs(title = "Predictor Effects Across Models",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave(paste0(interp_dir, "04b_coefficient_heatmap.png"), 
       p_heatmap, width = 12, height = 8, dpi = 300)
cat("✓ Heatmap saved\n\n")


# ==============================================================================
# SECTION 5: VARIABLE IMPORTANCE PLOT
# ==============================================================================

cat("Creating variable importance plot...\n")

p_importance <- ggplot(var_importance, 
                       aes(x = reorder(Predictor_clean, N_significant),
                           y = N_significant)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = N_significant), hjust = -0.3) +
  coord_flip() +
  labs(title = "Variable Importance (# significant models)",
       x = NULL, y = "Number of models") +
  theme_minimal() +
  ylim(0, length(names(models_glm)) + 0.5)

ggsave(paste0(interp_dir, "04b_variable_importance.png"),
       p_importance, width = 8, height = 6, dpi = 300)
cat("✓ Importance plot saved\n\n")


# ==============================================================================
# SECTION 6: ODDS RATIOS (BINOMIAL MODELS)
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("ODDS RATIOS (BINOMIAL MODELS)\n")
cat("════════════════════════════════��══════════════════════════\n\n")

odds_ratios <- all_coefs %>%
  dplyr::filter(Model != "calicioids_richness") %>%
  dplyr::mutate(
    OR = exp(Estimate),
    OR_lower = exp(Estimate - 1.96 * Std_Error),
    OR_upper = exp(Estimate + 1.96 * Std_Error),
    Pct_change = round((OR - 1) * 100, 1)
  ) %>%
  dplyr::select(Model, Predictor_clean, Estimate, OR, OR_lower, OR_upper, 
                Pct_change, P_value, Significant)

# Show significant only
sig_or <- odds_ratios %>% 
  dplyr::filter(Significant) %>%
  dplyr::arrange(desc(abs(Pct_change)))

cat("Top 10 significant effects:\n")
print(head(sig_or %>% dplyr::select(Model, Predictor_clean, OR, Pct_change), 10))

write.csv(odds_ratios, paste0(interp_dir, "04b_odds_ratios.csv"), row.names = FALSE)
cat("\n✓ Odds ratios saved\n\n")


# ==============================================================================
# SECTION 7: MARGINAL EFFECTS (TOP 4 PREDICTORS)
# ==============================================================================

cat("Creating marginal effects plots...\n")

top_4 <- var_importance$Predictor_clean[1:4]

for (pred in top_4) {
  pred_scaled <- paste0(pred, "_scaled")
  
  if (pred_scaled %in% names(coef(models_glm$core_ogf_presence))) {
    
    effects <- ggpredict(models_glm$core_ogf_presence, terms = pred_scaled)
    
    p <- plot(effects) +
      labs(title = paste("Core OGF:", pred),
           x = pred, y = "Predicted Probability") +
      theme_minimal()
    
    ggsave(paste0(interp_dir, "marginal_", pred, ".png"),
           p, width = 6, height = 4, dpi = 300)
  }
}

cat("✓ Marginal effects saved\n\n")


# ==============================================================================
# SECTION 8: DECAY STAGE COMPARISON
# ==============================================================================

cat("Creating decay stage comparison...\n")

decay_effects <- all_coefs %>%
  dplyr::filter(Predictor_clean %in% c("decay2", "decay3"))             # here add decay stages if you have more than 3 stages

p_decay <- ggplot(decay_effects, 
                  aes(x = Model, y = Estimate, fill = Predictor_clean)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Estimate - 1.96*Std_Error, 
                    ymax = Estimate + 1.96*Std_Error),
                position = position_dodge(0.8), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("decay2" = "#8da0cb", "decay3" = "#fc8d62"),
                    name = "Decay Stage") +
  labs(title = "Decay Stage Effects",
       x = NULL, y = "Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(interp_dir, "04b_decay_comparison.png"),
       p_decay, width = 10, height = 6, dpi = 300)
cat("✓ Decay comparison saved\n\n")


# ==============================================================================
# SECTION 9: ELEVATION EFFECTS
# ==============================================================================

cat("Creating elevation effects plot...\n")

elev_effects <- all_coefs %>%
  dplyr::filter(Predictor_clean == "elevation")

p_elev <- ggplot(elev_effects, aes(x = Model, y = Estimate, fill = Significant)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Estimate - 1.96*Std_Error, 
                    ymax = Estimate + 1.96*Std_Error), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70")) +
  labs(title = "Elevation Effects",
       x = NULL, y = "Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(interp_dir, "04b_elevation_effects.png"),
       p_elev, width = 10, height = 6, dpi = 300)
cat("✓ Elevation effects saved\n\n")


# ==============================================================================
# SECTION 10: PUBLICATION FIGURE (MULTI-PANEL)
# ==============================================================================

cat("Creating publication figure...\n")

plot_list <- list()

for (pred in top_4) {
  pred_scaled <- paste0(pred, "_scaled")
  
  if (pred_scaled %in% names(coef(models_glm$core_ogf_presence))) {
    
    effects <- ggpredict(models_glm$core_ogf_presence, terms = pred_scaled)
    
    coef_info <- all_coefs %>% 
      dplyr::filter(Model == "core_ogf_presence", Predictor_clean == pred)
    
    p <- ggplot(effects, aes(x = x, y = predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                  fill = "steelblue", alpha = 0.3) +
      geom_line(color = "steelblue", linewidth = 1) +
      labs(title = pred,
           subtitle = sprintf("β=%.2f %s", coef_info$Estimate, coef_info$Sig),
           x = paste(pred, "(scaled)"),
           y = "P(presence)") +
      theme_minimal(base_size = 10) +
      theme(plot.subtitle = element_text(size = 8))
    
    plot_list[[pred]] <- p
  }
}

pub_fig <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "Core Old-Growth Forest Lichen Drivers")

ggsave(paste0(interp_dir, "04b_publication_figure.png"),
       pub_fig, width = 10, height = 8, dpi = 300)
cat("✓ Publication figure saved\n\n")


# ==============================================================================
# SECTION 11: SUMMARY TABLES
# ==============================================================================

cat("Creating summary tables...\n")

# Significant predictors per model
sig_summary <- all_coefs %>%
  dplyr::filter(Significant) %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    N_significant = dplyr::n(),
    Predictors = paste(Predictor_clean, collapse = ", "),
    .groups = "drop"
  )

write.csv(sig_summary, paste0(interp_dir, "significant_predictors.csv"), 
          row.names = FALSE)

# Full coefficient table
write.csv(all_coefs, paste0(interp_dir, "04b_all_coefficients.csv"), 
          row.names = FALSE)

cat("✓ Summary tables saved\n\n")


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n╔══════════════��═══════════════════════════════════════════════╗\n")
cat("║  INTERPRETATION COMPLETE                                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Files created:\n")
cat(sprintf("  • %s files in %s\n", 
            length(list.files(interp_dir)), interp_dir))

cat("\nTop 3 predictors:\n")
print(var_importance[1:3, c("Predictor_clean", "N_significant", "Pct_significant")])

cat("\nModels with most significant predictors:\n")
print(sig_summary %>% dplyr::arrange(desc(N_significant)) %>% head(3))

cat("\n")

################################################################################
# END OF SCRIPT 04
################################################################################