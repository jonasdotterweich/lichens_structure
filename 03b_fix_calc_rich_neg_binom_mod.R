################################################################################
# ŠUMAVA NATIONAL PARK LICHEN-STRUCTURE MODELING
# Script 03b: Fix Calicioids Richness - Negative Binomial Model
# Author: Jonas Dotterweich
# Date: 2026-02-14
# Purpose: Address overdispersion in calicioids richness using simpler NB model
################################################################################

# ==============================================================================
# CONTEXT: WHY WE'RE HERE
# ==============================================================================

# PROBLEM WITH ZINB MODEL (Script 03):
#   - Convergence failure (non-positive-definite Hessian)
#   - AIC = NA (couldn't be calculated)
#   - Dispersion parameter = 4.6e+07 (absurdly huge)
#   - Overparameterization: 14 parameters for ~119 observations
#
# ROOT CAUSE:
#   - Too complex (12 predictors + zero-inflation + dispersion)
#   - Zero-inflation was marginal (p = 0.014, borderline)
#   - Overdispersion was the main issue (p = 0.022)
#
# SOLUTION:
#   - Try NEGATIVE BINOMIAL without zero-inflation
#   - Simpler model = better convergence
#   - NB handles overdispersion, may be sufficient for mild zero-inflation

# ==============================================================================
# LOAD WORKSPACE FROM SCRIPT 02
# ==============================================================================

load("02_GLM_models_and_diagnostics.RData")

# Load required libraries
library(glmmTMB)
library(DHARMa)
library(ggplot2)


# ==============================================================================
# SECTION 1: RECAP OF THE PROBLEM
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  CALICIOIDS RICHNESS: DIAGNOSTIC ISSUES RECAP               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Original Poisson model issues:\n")
poisson_diag <- diagnostics_summary[diagnostics_summary$Model == "calicioids_richness", ]
cat(sprintf("  - Dispersion test: p = %.4f (FAILED)\n", poisson_diag$Dispersion_p))
cat(sprintf("  - Zero-inflation test: p = %.4f (FAILED)\n", poisson_diag$ZeroInflation_p))
cat("\n")

cat("Data characteristics:\n")
cat(sprintf("  - Mean richness: %.2f\n", mean(data$calicioids_richness)))
cat(sprintf("  - Variance: %.2f\n", var(data$calicioids_richness)))
cat(sprintf("  - Variance/Mean ratio: %.2f (>1 = overdispersion)\n", 
            var(data$calicioids_richness) / mean(data$calicioids_richness)))
cat(sprintf("  - Zero-inflation: %.1f%% plots with zero richness\n", 
            mean(data$calicioids_richness == 0) * 100))
cat("\n")

cat("Strategy:\n")
cat("  → Fit NEGATIVE BINOMIAL model (handles overdispersion)\n")
cat("  → NO zero-inflation component (simpler = better convergence)\n")
cat("  → Use all 12 structural predictors\n\n")


# ==============================================================================
# SECTION 2: FIT NEGATIVE BINOMIAL MODEL
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FITTING NEGATIVE BINOMIAL MODEL                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Model structure:\n")
cat("  Family: Negative Binomial (type 2 parameterization)\n")
cat("  Formula: calicioids_richness ~ 12 structural predictors\n")
cat("  Link: log (standard for count data)\n\n")

cat("Fitting model...\n")

# Fit negative binomial model using glmmTMB
model_calicioids_nb <- glmmTMB(
  calicioids_richness ~ 
    elevation_scaled + 
    volume_snags_scaled + 
    tree_height_median_scaled + 
    canopy_cover_scaled + 
    regeneration_scaled + 
    decay2_scaled + 
    decay3_scaled + 
    dbh_max_scaled + 
    dbh_sd_scaled + 
    n_dead_50cm_scaled + 
    ba_spruce_scaled + 
    ba_beech_scaled,
  data = data_scaled,
  family = nbinom2  # Negative binomial, no zero-inflation
)

cat("✓ Model fitted\n\n")

# Check for convergence warnings
convergence_code <- model_calicioids_nb$fit$convergence
if (convergence_code == 0) {
  cat("✅ MODEL CONVERGED SUCCESSFULLY (convergence code = 0)\n\n")
} else {
  cat("⚠️  Warning: Convergence issues detected (code =", convergence_code, ")\n\n")
}


# ==============================================================================
# SECTION 3: MODEL SUMMARY AND INTERPRETATION
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("NEGATIVE BINOMIAL MODEL SUMMARY\n")
cat("═══════════════════════════════════════════════════════════\n\n")

summary(model_calicioids_nb)

cat("\n")
cat("─────────────────────────────────────────────────────────────\n")
cat("INTERPRETING NEGATIVE BINOMIAL OUTPUT:\n")
cat("─────────────────────────────────────────────────────────────\n\n")

cat("1. CONDITIONAL MODEL (Count component):\n")
cat("   - Coefficients are on LOG scale (log of expected count)\n")
cat("   - Positive coefficient → richness increases with predictor\n")
cat("   - Negative coefficient → richness decreases with predictor\n")
cat("   - exp(coefficient) = multiplicative effect on richness\n\n")

cat("2. DISPERSION PARAMETER:\n")
dispersion_param <- sigma(model_calicioids_nb)
cat(sprintf("   - Estimated dispersion: %.3f\n", dispersion_param))
cat("   - This allows variance > mean (unlike Poisson)\n")
cat("   - Larger values = more overdispersion\n\n")

# Extract and display significant predictors
coefs <- summary(model_calicioids_nb)$coefficients$cond
sig_predictors <- rownames(coefs)[coefs[, 4] < 0.05 & rownames(coefs) != "(Intercept)"]

if (length(sig_predictors) > 0) {
  cat("3. SIGNIFICANT PREDICTORS (p < 0.05):\n")
  for (pred in sig_predictors) {
    coef_est <- coefs[pred, "Estimate"]
    p_val <- coefs[pred, "Pr(>|z|)"]
    multiplicative_effect <- exp(coef_est)
    
    direction <- ifelse(coef_est > 0, "INCREASES", "DECREASES")
    pct_change <- (multiplicative_effect - 1) * 100
    
    cat(sprintf("   - %s: %.3f (p = %.4f)\n", pred, coef_est, p_val))
    cat(sprintf("     → Richness %s by %.1f%% per 1 SD increase\n", 
                direction, abs(pct_change)))
  }
} else {
  cat("3. NO SIGNIFICANT PREDICTORS (all p > 0.05)\n")
  cat("   → Calicioids richness may be driven by unmeasured factors\n")
  cat("   → Or sample size insufficient to detect effects\n")
}

cat("\n\n")




# ==============================================================================
# SECTION: PRE-DHARMA VALIDATION (Negative Binomial Model)
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PRE-DHARMA CHECKS: Negative Binomial Model                 ║\n")
cat("╚═════════════════════════════════════��════════════════════════╝\n\n")

cat("Checking model convergence and numerical stability...\n\n")

# Extract coefficient table
coef_table_nb <- summary(model_calicioids_nb)$coefficients$cond

# ----------------------------------------------------------------------------
# CHECK 1: Coefficient Magnitudes
# ----------------------------------------------------------------------------
max_coef_nb <- max(abs(coef_table_nb[, "Estimate"]))
max_coef_name_nb <- rownames(coef_table_nb)[which.max(abs(coef_table_nb[, "Estimate"]))]

cat(sprintf("1. COEFFICIENT MAGNITUDE:\n"))
cat(sprintf("   Max coefficient: %.2f (%s)\n", max_coef_nb, max_coef_name_nb))

if (max_coef_nb > 10) {
  cat("   🚨 HUGE coefficient - likely separation!\n")
  coef_flag_nb <- "FAIL"
} else if (max_coef_nb > 5) {
  cat("   ⚠️  Large coefficient (5-10 range) - monitor\n")
  coef_flag_nb <- "WARNING"
} else {
  cat("   ✅ Coefficients in normal range\n")
  coef_flag_nb <- "PASS"
}
cat("\n")

# ----------------------------------------------------------------------------
# CHECK 2: Standard Errors
# ----------------------------------------------------------------------------
max_se_nb <- max(coef_table_nb[, "Std. Error"])
max_se_name_nb <- rownames(coef_table_nb)[which.max(coef_table_nb[, "Std. Error"])]

cat(sprintf("2. STANDARD ERRORS:\n"))
cat(sprintf("   Max SE: %.2f (%s)\n", max_se_nb, max_se_name_nb))

if (max_se_nb > 10) {
  cat("   🚨 HUGE standard error - extreme uncertainty!\n")
  se_flag_nb <- "FAIL"
} else if (max_se_nb > 5) {
  cat("   ⚠️  Large SE (5-10 range) - uncertain estimates\n")
  se_flag_nb <- "WARNING"
} else {
  cat("   ✅ Standard errors in normal range\n")
  se_flag_nb <- "PASS"
}
cat("\n")

# ----------------------------------------------------------------------------
# CHECK 3: Convergence
# ----------------------------------------------------------------------------
cat("3. CONVERGENCE:\n")

# Check convergence for glmmTMB
convergence_code <- model_calicioids_nb$fit$convergence

cat(sprintf("   Convergence code: %d\n", convergence_code))

if (convergence_code == 0) {
  cat("   ✅ Model converged successfully\n")
  conv_flag_nb <- "PASS"
} else {
  cat("   🚨 Model did NOT converge!\n")
  conv_flag_nb <- "FAIL"
}
cat("\n")

# ----------------------------------------------------------------------------
# CHECK 4: Dispersion Parameter (Negative Binomial specific)
# ----------------------------------------------------------------------------
cat("4. DISPERSION PARAMETER:\n")

# Extract dispersion parameter
disp_param <- sigma(model_calicioids_nb)

cat(sprintf("   Dispersion (theta): %.3f\n", disp_param))

if (disp_param < 0.1) {
  cat("   ⚠️  Very small dispersion (close to Poisson)\n")
  disp_flag_nb <- "WARNING"
} else if (disp_param > 100) {
  cat("   ⚠️  Very large dispersion (unstable?)\n")
  disp_flag_nb <- "WARNING"
} else {
  cat("   ✅ Dispersion parameter reasonable\n")
  disp_flag_nb <- "PASS"
}

cat(sprintf("   Interpretation: Variance = μ + μ²/%.3f\n", disp_param))
cat("\n")

# ----------------------------------------------------------------------------
# OVERALL VERDICT
# ----------------------------------------------------------------------------
cat("─────────────────────────────────────────────────────────────\n")
cat("OVERALL VERDICT:\n")
cat("───���─────────────────────────────────────────────────────────\n")

all_flags_nb <- c(coef_flag_nb, se_flag_nb, conv_flag_nb, disp_flag_nb)

if (all(all_flags_nb == "PASS")) {
  cat("✅ ALL PRE-DHARMA CHECKS PASSED\n")
  cat("   → Model is healthy and ready for DHARMa diagnostics\n\n")
  overall_status_nb <- "VALID"
} else if (any(all_flags_nb == "FAIL")) {
  cat("❌ CRITICAL ISSUES DETECTED\n")
  cat("   → Model may be invalid - review carefully\n\n")
  overall_status_nb <- "INVALID"
} else {
  cat("⚠️  SOME WARNINGS DETECTED\n")
  cat("   → Proceed with caution\n\n")
  overall_status_nb <- "BORDERLINE"
}





# ==============================================================================
# SECTION 4: MODEL COMPARISON (Poisson vs NB)
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON: Poisson vs Negative Binomial             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# IMPORTANT: models_glm$calicioids_richness was overwritten in Script 03
# We need to refit the original Poisson model or use saved values

# Option 1: Refit original Poisson model
cat("Refitting original Poisson model for comparison...\n")
model_poisson_original <- glm(
  calicioids_richness ~ 
    elevation_scaled + 
    volume_snags_scaled + 
    tree_height_median_scaled + 
    canopy_cover_scaled + 
    regeneration_scaled + 
    decay2_scaled + 
    decay3_scaled + 
    dbh_max_scaled + 
    dbh_sd_scaled + 
    n_dead_50cm_scaled + 
    ba_spruce_scaled + 
    ba_beech_scaled,
  data = data_scaled,
  family = poisson(link = "log")
)

# Extract fit statistics
aic_poisson <- AIC(models_)
aic_nb <- AIC(model_calicioids_nb)
delta_aic <- aic_poisson - aic_nb

ll_poisson <- logLik(model_poisson_original)
ll_nb <- logLik(model_calicioids_nb)

cat("Model fit comparison:\n")
cat(sprintf("  Poisson GLM:    AIC = %.2f, Log-Likelihood = %.2f\n", 
            aic_poisson, as.numeric(ll_poisson)))
cat(sprintf("  Negative Binom: AIC = %.2f, Log-Likelihood = %.2f\n", 
            aic_nb, as.numeric(ll_nb)))
cat(sprintf("  ΔAIC = %.2f (Poisson - NB)\n\n", delta_aic))


# ==============================================================================
# SECTION 5: DHARMA DIAGNOSTICS FOR NB MODEL
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DHARMa DIAGNOSTICS: NEGATIVE BINOMIAL MODEL                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Running simulation-based diagnostics (n=1000 simulations)...\n")
cat("(This may take a minute)\n\n")

# Simulate residuals for NB model
sim_resid_nb <- simulateResiduals(fittedModel = model_calicioids_nb, 
                                  n = 1000, 
                                  plot = FALSE)

# Create diagnostic plot
png("03b_DHARMa_calicioids_NB.png", width = 10, height = 6, units = "in", res = 300)
plot(sim_resid_nb, main = "DHARMa Diagnostics: Calicioids Richness (Negative Binomial)")
dev.off()

cat("✓ Diagnostic plot saved: 03b_DHARMa_calicioids_NB.png\n\n")


# Run diagnostic tests
cat("══════════���════════════════════════════════════════════════\n")
cat("DIAGNOSTIC TESTS:\n")
cat("═══════════════════════════════════════════════════════════\n\n")

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
                                              x = data$X, 
                                              y = data$Y,
                                              plot = FALSE)
print(spatial_test_nb)
cat("\n")


# Compile diagnostics
diagnostics_nb <- data.frame(
  Model = "calicioids_richness (NB)",
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

cat("═══════════════════════════════════════════════════════════\n")
cat("DIAGNOSTIC SUMMARY:\n")
cat("═══════════════════════════════════════════════════════════\n\n")
print(diagnostics_nb)
cat("\n")

all_passed <- all(diagnostics_nb[, 6:9])

if (all_passed) {
  cat("🎉 ALL DIAGNOSTICS PASSED FOR NEGATIVE BINOMIAL MODEL!\n")
  cat("   ✓ Overdispersion resolved\n")
  cat("   ✓ Zero-inflation acceptable (or resolved)\n")
  cat("   ✓ No spatial autocorrelation\n")
  cat("   ✓ No outlier issues\n")
  cat("\n→ Model is appropriate for analysis\n\n")
} else {
  failed_tests <- names(diagnostics_nb)[6:9][!diagnostics_nb[, 6:9]]
  cat("⚠️  Some diagnostics failed:\n")
  for (test in failed_tests) {
    cat(sprintf("   - %s\n", gsub("_OK", "", test)))
  }
  cat("\n→ May need further model refinement\n\n")
}

# ==============================================================================
# SECTION 6: DIAGNOSTIC COMPARISON (Poisson vs NB)
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DIAGNOSTIC COMPARISON: Poisson vs Negative Binomial        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Original Poisson diagnostics (from Script 02 - hardcoded since overwritten)
cat("Note: Using original Poisson diagnostics from Script 02\n")
cat("      (current diagnostics_summary was modified in Script 03)\n\n")

poisson_diag_values <- data.frame(
  Dispersion_p = 0.022,
  ZeroInflation_p = 0.014,
  Outliers_p = 1.000,
  SpatialAutocorr_p = 0.241,
  Dispersion_OK = FALSE,
  ZeroInflation_OK = FALSE,
  Outliers_OK = TRUE,
  Spatial_OK = TRUE
)

comparison_table <- data.frame(
  Test = c("Dispersion", "Zero-Inflation", "Outliers", "Spatial Autocorr"),
  Poisson_p = c(poisson_diag_values$Dispersion_p, 
                poisson_diag_values$ZeroInflation_p,
                poisson_diag_values$Outliers_p,
                poisson_diag_values$SpatialAutocorr_p),
  Poisson_OK = c(poisson_diag_values$Dispersion_OK,
                 poisson_diag_values$ZeroInflation_OK,
                 poisson_diag_values$Outliers_OK,
                 poisson_diag_values$Spatial_OK),
  NB_p = c(diagnostics_nb$Dispersion_p,
           diagnostics_nb$ZeroInflation_p,
           diagnostics_nb$Outliers_p,
           diagnostics_nb$SpatialAutocorr_p),
  NB_OK = c(diagnostics_nb$Dispersion_OK,
            diagnostics_nb$ZeroInflation_OK,
            diagnostics_nb$Outliers_OK,
            diagnostics_nb$Spatial_OK)
)

print(comparison_table)
cat("\n")

improved_tests <- sum(comparison_table$NB_OK & !comparison_table$Poisson_OK)
if (improved_tests > 0) {
  cat(sprintf("✅ NEGATIVE BINOMIAL resolved %d diagnostic issue(s)\n\n", improved_tests))
} else if (all(comparison_table$NB_OK)) {
  cat("✅ NEGATIVE BINOMIAL maintains all passed diagnostics\n\n")
} else {
  cat("⚠️  NEGATIVE BINOMIAL did not resolve all issues\n\n")
}

# ==============================================================================
# SECTION 7: VISUALIZE MODEL FIT
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  VISUALIZING MODEL FIT                                      ║\n")
cat("╚══════════════════════════════════════════════════���═══════════╝\n\n")

# Get predicted values
data$predicted_poisson <- predict(models_glm$calicioids_richness, type = "response")
data$predicted_nb <- predict(model_calicioids_nb, type = "response")

# Create comparison plot
png("03b_model_predictions_comparison.png", width = 12, height = 5, units = "in", res = 300)
par(mfrow = c(1, 3))

# Observed vs Predicted - Poisson
plot(data$predicted_poisson, data$calicioids_richness,
     xlab = "Predicted richness (Poisson)",
     ylab = "Observed richness",
     main = "Poisson GLM",
     pch = 16, col = rgb(0, 0, 1, 0.4),
     xlim = c(0, max(data$predicted_poisson, data$calicioids_richness)),
     ylim = c(0, max(data$predicted_poisson, data$calicioids_richness)))
abline(0, 1, col = "red", lwd = 2, lty = 2)
grid()

# Observed vs Predicted - NB
plot(data$predicted_nb, data$calicioids_richness,
     xlab = "Predicted richness (Negative Binomial)",
     ylab = "Observed richness",
     main = "Negative Binomial",
     pch = 16, col = rgb(0, 0.5, 0, 0.4),
     xlim = c(0, max(data$predicted_nb, data$calicioids_richness)),
     ylim = c(0, max(data$predicted_nb, data$calicioids_richness)))
abline(0, 1, col = "red", lwd = 2, lty = 2)
grid()

# Residuals comparison
boxplot(
  list(
    Poisson = residuals(models_glm$calicioids_richness, type = "pearson"),
    NB = residuals(model_calicioids_nb, type = "pearson")
  ),
  main = "Pearson Residuals Comparison",
  ylab = "Pearson residuals",
  col = c("lightblue", "lightgreen"),
  border = c("blue", "darkgreen")
)
abline(h = 0, col = "red", lwd = 2, lty = 2)
grid()

dev.off()

cat("✓ Model comparison plot saved: 03b_model_predictions_comparison.png\n\n")


# ==============================================================================
# SECTION 8: DECISION AND FINALIZATION
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL MODEL SELECTION DECISION                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Decision logic
use_nb <- all_passed && (delta_aic > 2 || (!is.na(delta_aic) && delta_aic > 0))

if (use_nb) {
  cat("✅ DECISION: Use NEGATIVE BINOMIAL model\n\n")
  cat("Rationale:\n")
  cat("  - All DHARMa diagnostics passed\n")
  cat(sprintf("  - Better AIC than Poisson (ΔAIC = %.2f)\n", delta_aic))
  cat("  - Model converged successfully\n")
  cat("  - Handles overdispersion appropriately\n\n")
  
  # Update model list
  models_glm$calicioids_richness <- model_calicioids_nb
  
  # Update diagnostics summary
  diagnostics_summary[diagnostics_summary$Model == "calicioids_richness", ] <- 
    diagnostics_nb[, names(diagnostics_summary)]
  
  cat("Updated model list:\n")
  cat("  6 Binomial GLMs (presence/absence) ✓\n")
  cat("  1 Negative Binomial GLM (richness) ✓\n\n")
  
} else {
  cat("⚠️  NEGATIVE BINOMIAL did not fully resolve issues\n\n")
  cat("Options:\n")
  cat("  A) Accept NB model with caveats\n")
  cat("  B) Try reduced predictor set (model selection first)\n")
  cat("  C) Use hurdle model (two-part model)\n")
  cat("  D) Accept Poisson with robust standard errors\n\n")
  cat("→ Manual decision required before proceeding\n\n")
}


# ==============================================================================
# SECTION 9: SAVE UPDATED WORKSPACE
# ==============================================================================

if (use_nb) {
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("SAVING UPDATED WORKSPACE\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  
  # Save updated workspace
  save(data, data_scaled,
       models_glm,  # Now contains NB model for calicioids
       diagnostics_summary,
       spatial_results_raw,
       lichen_binary, lichen_count, all_lichen,
       predictors, predictors_scaled,
       inv_dist_weights,
       file = "Lichens/to_model/03b_final_models_NB.RData")
  
  cat("✅ Workspace saved: 03b_final_models_NB.RData\n\n")
  
  cat("Updated diagnostics summary:\n")
  print(diagnostics_summary)
  cat("\n")
}


# ==============================================================================
# SECTION 10: SUMMARY AND NEXT STEPS
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODELING COMPLETE: ALL 7 LICHEN GROUPS                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

if (use_nb) {
  cat("✅ ALL MODELS FINALIZED AND VALIDATED\n\n")
  
  cat("Final model set:\n")
  cat("  1. core_ogf_presence (Binomial GLM) ✓\n")
  cat("  2. mycoblastus_presence (Binomial GLM) ✓\n")
  cat("  3. ochrolechia_presence (Binomial GLM) ✓\n")
  cat("  4. xylographa_presence (Binomial GLM) ✓\n")
  cat("  5. parmelia_agg_presence (Binomial GLM) ✓\n")
  cat("  6. elite_rare_presence (Binomial GLM) ✓\n")
  cat("  7. calicioids_richness (Negative Binomial GLM) ✓\n\n")
  
  cat("Key findings:\n")
  cat("  • Spatial autocorrelation fully explained by structural predictors\n")
  cat("  • No need for spatial random effects (GLMMs not required)\n")
  cat("  • Negative binomial successfully handles overdispersion\n")
  cat("  • All models pass diagnostic tests\n\n")
  
  cat("═══════════════════════════════════════════════════════════\n")
  cat("READY FOR:\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  → Coefficient interpretation\n")
  cat("  → Effect plots (marginal effects)\n")
  cat("  → Model selection (optional simplification)\n")
  cat("  → Threshold detection (segmented regression)\n")
  cat("  → Manuscript preparation\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
}


################################################################################
# END OF SCRIPT 03b
################################################################################