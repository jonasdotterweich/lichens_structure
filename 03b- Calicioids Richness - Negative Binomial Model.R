################################################################################
# ŠUMAVA NATIONAL PARK LICHEN-STRUCTURE MODELING
# Script 03b: Calicioids Richness - Negative Binomial Model
# Author: Jonas Dotterweich
# Date: 2026-02-14
# Purpose: Address overdispersion using negative binomial model
################################################################################

# ==============================================================================
# SETUP
# ==============================================================================

load("Lichens/to_model/02_GLM_models_and_diagnostics.RData")

library(glmmTMB)
library(DHARMa)
library(ggplot2)


# ==============================================================================
# SECTION 1: PROBLEM RECAP
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  CALICIOIDS RICHNESS: DIAGNOSTIC ISSUES                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

poisson_diag <- diagnostics_summary[diagnostics_summary$Model == "calicioids_richness", ]

cat("Original Poisson model:\n")
cat(sprintf("  Dispersion test:     p = %.4f (FAILED)\n", poisson_diag$Dispersion_p))
cat(sprintf("  Zero-inflation test: p = %.4f (FAILED)\n", poisson_diag$ZeroInflation_p))
cat(sprintf("  Mean richness: %.2f, Variance: %.2f (ratio: %.2f)\n", 
            mean(data$calicioids_richness), 
            var(data$calicioids_richness),
            var(data$calicioids_richness) / mean(data$calicioids_richness)))
cat(sprintf("  Zeros: %.1f%% of plots\n\n", mean(data$calicioids_richness == 0) * 100))

cat("Strategy: Negative Binomial model (handles overdispersion)\n\n")


# ==============================================================================
# SECTION 2: FIT NEGATIVE BINOMIAL MODEL
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FITTING NEGATIVE BINOMIAL MODEL                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

model_calicioids_nb <- glmmTMB(
  calicioids_richness ~ 
    elevation_scaled + volume_snags_scaled + tree_height_median_scaled + 
    canopy_cover_scaled + regeneration_scaled + decay2_scaled + 
    decay3_scaled + dbh_max_scaled + dbh_sd_scaled + 
    n_dead_50cm_scaled + ba_spruce_scaled + ba_beech_scaled,
  data = data_scaled,
  family = nbinom2
)

convergence_code <- model_calicioids_nb$fit$convergence

if (convergence_code == 0) {
  cat("✅ Model converged successfully\n\n")
} else {
  cat("⚠️  Convergence issues (code =", convergence_code, ")\n\n")
}


# ==============================================================================
# SECTION 3: MODEL SUMMARY
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("MODEL SUMMARY\n")
cat("═══════════════════════════════════════════════════════════\n\n")

summary(model_calicioids_nb)

cat("\n")
dispersion_param <- sigma(model_calicioids_nb)
cat(sprintf("Dispersion parameter: %.3f (allows variance > mean)\n", dispersion_param))

# Significant predictors
coefs <- summary(model_calicioids_nb)$coefficients$cond
sig_predictors <- rownames(coefs)[coefs[, 4] < 0.05 & rownames(coefs) != "(Intercept)"]

if (length(sig_predictors) > 0) {
  cat("\nSignificant predictors (p < 0.05):\n")
  for (pred in sig_predictors) {
    coef_est <- coefs[pred, "Estimate"]
    p_val <- coefs[pred, "Pr(>|z|)"]
    pct_change <- (exp(coef_est) - 1) * 100
    cat(sprintf("  %s: %.3f (p = %.4f) → %+.1f%% per SD\n", 
                pred, coef_est, p_val, pct_change))
  }
} else {
  cat("\nNo significant predictors (all p > 0.05)\n")
}
cat("\n\n")


# ==============================================================================
# SECTION 4: PRE-DHARMA VALIDATION
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PRE-DHARMA CHECKS                                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

coef_table_nb <- summary(model_calicioids_nb)$coefficients$cond

# Check 1: Coefficients
max_coef_nb <- max(abs(coef_table_nb[, "Estimate"]))
max_coef_name_nb <- rownames(coef_table_nb)[which.max(abs(coef_table_nb[, "Estimate"]))]

cat(sprintf("1. Coefficient magnitude: %.2f (%s)", max_coef_nb, max_coef_name_nb))
coef_flag_nb <- if (max_coef_nb > 10) {
  cat(" 🚨 HUGE\n"); "FAIL"
} else if (max_coef_nb > 5) {
  cat(" ⚠️  Large\n"); "WARNING"
} else {
  cat(" ✅\n"); "PASS"
}

# Check 2: Standard errors
max_se_nb <- max(coef_table_nb[, "Std. Error"])
max_se_name_nb <- rownames(coef_table_nb)[which.max(coef_table_nb[, "Std. Error"])]

cat(sprintf("2. Standard errors: %.2f (%s)", max_se_nb, max_se_name_nb))
se_flag_nb <- if (max_se_nb > 10) {
  cat(" 🚨 HUGE\n"); "FAIL"
} else if (max_se_nb > 5) {
  cat(" ⚠️  Large\n"); "WARNING"
} else {
  cat(" ✅\n"); "PASS"
}

# Check 3: Convergence
cat(sprintf("3. Convergence: code = %d", convergence_code))
conv_flag_nb <- if (convergence_code == 0) {
  cat(" ✅\n"); "PASS"
} else {
  cat(" 🚨\n"); "FAIL"
}

# Check 4: Dispersion
cat(sprintf("4. Dispersion: %.3f", dispersion_param))
disp_flag_nb <- if (dispersion_param < 0.1 || dispersion_param > 100) {
  cat(" ⚠️  Extreme\n"); "WARNING"
} else {
  cat(" ✅\n"); "PASS"
}

# Verdict
cat("\n")
all_flags_nb <- c(coef_flag_nb, se_flag_nb, conv_flag_nb, disp_flag_nb)

overall_status_nb <- if (all(all_flags_nb == "PASS")) {
  cat("✅ ALL PRE-DHARMA CHECKS PASSED\n\n"); "VALID"
} else if (any(all_flags_nb == "FAIL")) {
  cat("❌ CRITICAL ISSUES DETECTED\n\n"); "INVALID"
} else {
  cat("⚠️  SOME WARNINGS\n\n"); "BORDERLINE"
}


# ==============================================================================
# SECTION 5: MODEL COMPARISON (Poisson vs NB)
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON                                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Refit Poisson for comparison
model_poisson_original <- glm(
  calicioids_richness ~ elevation_scaled + volume_snags_scaled + 
    tree_height_median_scaled + canopy_cover_scaled + regeneration_scaled + 
    decay2_scaled + decay3_scaled + dbh_max_scaled + dbh_sd_scaled + 
    n_dead_50cm_scaled + ba_spruce_scaled + ba_beech_scaled,
  data = data_scaled,
  family = poisson(link = "log")
)

aic_poisson <- AIC(model_poisson_original)
aic_nb <- AIC(model_calicioids_nb)
delta_aic <- aic_poisson - aic_nb

cat(sprintf("AIC:  Poisson = %.2f, NB = %.2f, ΔAIC = %.2f", 
            aic_poisson, aic_nb, delta_aic))

if (delta_aic > 10) {
  cat(" (NB overwhelmingly better)\n")
} else if (delta_aic > 4) {
  cat(" (NB substantially better)\n")
} else if (delta_aic > 2) {
  cat(" (NB moderately better)\n")
} else {
  cat(" (Similar)\n")
}
cat("\n")


# ==============================================================================
# SECTION 6: DHARMA DIAGNOSTICS
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DHARMa DIAGNOSTICS                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

sim_resid_nb <- simulateResiduals(fittedModel = model_calicioids_nb, 
                                  n = 1000, plot = FALSE)

png("03b_DHARMa_calicioids_NB.png", width = 10, height = 6, units = "in", res = 300)
plot(sim_resid_nb, main = "DHARMa: Calicioids Richness (Negative Binomial)")
dev.off()

cat("Running tests...\n\n")

cat("1. DISPERSION TEST:\n")
disp_test_nb <- testDispersion(sim_resid_nb, plot = FALSE)
print(disp_test_nb)
cat("\n")

cat("2. ZERO-INFLATION TEST:\n")
zi_test_nb <- testZeroInflation(sim_resid_nb, plot = FALSE)
print(zi_test_nb)
cat("\n")

cat("3. OUTLIER TEST:\n")
outlier_test_nb <- testOutliers(sim_resid_nb, plot = FALSE)
print(outlier_test_nb)
cat("\n")

cat("4. SPATIAL AUTOCORRELATION TEST:\n")
spatial_test_nb <- testSpatialAutocorrelation(sim_resid_nb, 
                                              x = data$X, y = data$Y,
                                              plot = FALSE)
print(spatial_test_nb)
cat("\n")

# Compile results
diagnostics_nb <- data.frame(
  Model = "calicioids_richness",
  Dispersion_p = round(disp_test_nb$p.value, 4),
  ZeroInflation_p = round(zi_test_nb$p.value, 4),
  Outliers_p = round(outlier_test_nb$p.value, 4),
  SpatialAutocorr_p = round(spatial_test_nb$p.value, 4),
  Dispersion_OK = disp_test_nb$p.value > 0.05,
  ZeroInflation_OK = zi_test_nb$p.value > 0.05,
  Outliers_OK = outlier_test_nb$p.value > 0.05,
  Spatial_OK = spatial_test_nb$p.value > 0.05,
  stringsAsFactors = FALSE
)

all_passed <- all(diagnostics_nb[, 6:9])

cat("═══════════════════════════════════════════════════════════\n")
if (all_passed) {
  cat("🎉 ALL DIAGNOSTICS PASSED\n\n")
} else {
  failed <- names(diagnostics_nb)[6:9][!diagnostics_nb[, 6:9]]
  cat("⚠️  Failed:", paste(gsub("_OK", "", failed), collapse = ", "), "\n\n")
}


# ==============================================================================
# SECTION 7: DIAGNOSTIC COMPARISON
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  BEFORE vs AFTER COMPARISON                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

comparison <- data.frame(
  Test = c("Dispersion", "Zero-Inflation", "Outliers", "Spatial"),
  Poisson_p = c(0.022, 0.014, 1.000, 0.241),
  Poisson_OK = c(FALSE, FALSE, TRUE, TRUE),
  NB_p = c(diagnostics_nb$Dispersion_p, diagnostics_nb$ZeroInflation_p,
           diagnostics_nb$Outliers_p, diagnostics_nb$SpatialAutocorr_p),
  NB_OK = c(diagnostics_nb$Dispersion_OK, diagnostics_nb$ZeroInflation_OK,
            diagnostics_nb$Outliers_OK, diagnostics_nb$Spatial_OK)
)

print(comparison)
cat(sprintf("\n✅ Negative Binomial resolved %d issue(s)\n\n", 
            sum(comparison$NB_OK & !comparison$Poisson_OK)))


# ==============================================================================
# SECTION 8: VISUALIZE FIT
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  VISUALIZING MODEL FIT                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

data$predicted_poisson <- predict(model_poisson_original, type = "response")
data$predicted_nb <- predict(model_calicioids_nb, type = "response")

png("03b_model_comparison.png", width = 12, height = 5, units = "in", res = 300)
par(mfrow = c(1, 3))

# Poisson
plot(data$predicted_poisson, data$calicioids_richness,
     xlab = "Predicted", ylab = "Observed", main = "Poisson",
     pch = 16, col = rgb(0, 0, 1, 0.4))
abline(0, 1, col = "red", lwd = 2, lty = 2); grid()

# NB
plot(data$predicted_nb, data$calicioids_richness,
     xlab = "Predicted", ylab = "Observed", main = "Negative Binomial",
     pch = 16, col = rgb(0, 0.5, 0, 0.4))
abline(0, 1, col = "red", lwd = 2, lty = 2); grid()

# Residuals
boxplot(list(Poisson = residuals(model_poisson_original, type = "pearson"),
             NB = residuals(model_calicioids_nb, type = "pearson")),
        main = "Residuals", col = c("lightblue", "lightgreen"))
abline(h = 0, col = "red", lwd = 2, lty = 2); grid()

dev.off()
cat("✓ Plots saved\n\n")


# ==============================================================================
# SECTION 9: FINAL DECISION
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL DECISION                                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

use_nb <- all_passed && delta_aic > 2

if (use_nb) {
  cat("✅ USE NEGATIVE BINOMIAL MODEL\n\n")
  cat("Rationale:\n")
  cat(sprintf("  • All diagnostics passed\n"))
  cat(sprintf("  • Better AIC (ΔAIC = %.2f)\n", delta_aic))
  cat(sprintf("  • Converged successfully\n\n"))
  
  # Update model list
  models_glm$calicioids_richness <- model_calicioids_nb
  diagnostics_summary[diagnostics_summary$Model == "calicioids_richness", ] <- 
    diagnostics_nb[, names(diagnostics_summary)]
  
  cat("Updated:\n  6 Binomial GLMs + 1 Negative Binomial GLM ✓\n\n")
  
} else {
  cat("⚠️  Manual decision required\n\n")
}


# ==============================================================================
# SECTION 10: SAVE WORKSPACE
# ==============================================================================

if (use_nb) {
  save(data, data_scaled, models_glm, diagnostics_summary, spatial_results_raw,
       lichen_binary, lichen_count, all_lichen, predictors, predictors_scaled,
       inv_dist_weights, model_calicioids_nb,
       file = "Lichens/to_model/03b_final_models_NB.RData")
  
  cat("✅ Workspace saved: 03b_final_models_NB.RData\n\n")
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("ALL 7 MODELS FINALIZED\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  1-6. Binomial GLMs (presence/absence) ✓\n")
  cat("    7. Negative Binomial GLM (richness) ✓\n\n")
  cat("Ready for interpretation and visualization!\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
}

################################################################################
# END OF SCRIPT 03b
################################################################################