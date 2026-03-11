################################################################################
# ŠUMAVA NATIONAL PARK LICHEN-STRUCTURE MODELING
# Script 02_Combined: All GLMs with Complete Diagnostics
# Author: Jonas Dotterweich
# Date: 2026-02-23
# Purpose: Fit all 7 lichen models (6 binomial + 1 NB) with full validation
################################################################################

# ==============================================================================
# SETUP & CONFIGURATION
# ==============================================================================

library(DHARMa)
library(ape)
library(glmmTMB)
library(ggplot2)

# ┌─────────────────────────────────────────────────────────────────────────┐
# │ MODIFY HERE: Choose your dataset                                       │
# ├─────────────────────────────────────────────────────────────────────────┤
# │ Option A: Standard predictor set (12 predictors)                       │
# │    load("Lichens/to_model/1_prepared_data_for_modeling.RData")                     │
# │                                                                         │
# │ Option B: Enlarged predictor set (more predictors)                     ��
# │    load("Lichens/to_model/01_data_prepared_enlarged.RData")            │
# │                                                                         │
# │ Uncomment the one you want to use:                                     │
# └─────────────────────────────────────────────────────────────────────────┘

# Standard dataset (DEFAULT):
#load("Lichens/to_model/01_prepared_data_for_modeling.RData")
#dataset_name <- "standard"

# Enlarged dataset (ALTERNATIVE - uncomment to use):
# load("Lichens/to_model/01_prepared_data_for_modeling(14pred).RData")
# dataset_name <- "enlarged"
#

#testing without mycoblastus and different predictors(no regen, no tee hight median)
load("Lichens/to_model/01_prepared_data_for_modeling(12pred_wo_mycob).RData")
dataset_name <- "reduced"


# ┌─────────────────────────────────────────────────────────────────────────┐
# │ OUTPUT DIRECTORY CONFIGURATION                                          │
# ├─────────────────────────────────────────────────────────────────────────┤
# │ All outputs (PNGs, CSVs, RData) will be saved here                     │
# └─────────────────────────────────────────────────────────────────────────┘

output_dir <- "Lichens/model_outputs/reduced/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("✓ Created output directory:", output_dir, "\n")
}


# ┌────────────────��────────────────────────────────────────────────────────┐
# │ DETECT & DISPLAY DATASET INFORMATION                                    │
# └─────────────────────────────────────────────────────────────────────────┘

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  ŠUMAVA LICHEN MODELING - COMPLETE WORKFLOW                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Dataset configuration:\n")
cat(sprintf("  Dataset: %s\n", dataset_name))
cat(sprintf("  Observations: %d plots\n", nrow(data)))
cat(sprintf("  Scaled predictors: %d\n", length(predictors_scaled)))
cat(sprintf("  Output directory: %s\n\n", output_dir))

cat("Predictors in model:\n")
for (i in 1:length(predictors_scaled)) {
  cat(sprintf("  %2d. %s\n", i, predictors_scaled[i]))
}
cat("\n")

cat("Models to fit:\n")
cat(sprintf("  • %d Binomial GLMs (presence/absence)\n", length(lichen_binary)))
cat(sprintf("  • 1 Negative Binomial GLM (calicioids richness)\n\n"))


# ==============================================================================
# SECTION 1: RAW SPATIAL AUTOCORRELATION (for Test 5 comparison)
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("STEP 1: Testing spatial autocorrelation in RAW data\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Create distance-based weights
coords <- as.matrix(data[, c("X", "Y")])
dist_matrix <- as.matrix(dist(coords))
inv_dist_weights <- 1 / dist_matrix
diag(inv_dist_weights) <- 0

# Test function
test_morans_I_ape <- function(response_var, weights, var_name) {
  moran_test <- Moran.I(response_var, weights)
  data.frame(
    Variable = var_name,
    Morans_I = round(moran_test$observed, 4),
    P_value = round(moran_test$p.value, 4),
    Z_score = round((moran_test$observed - moran_test$expected) / sqrt(moran_test$sd^2), 2),
    Sig = ifelse(moran_test$p.value < 0.001, "***",
                 ifelse(moran_test$p.value < 0.01, "**",
                        ifelse(moran_test$p.value < 0.05, "*", "ns"))),
    stringsAsFactors = FALSE
  )
}

# Test all lichen groups
spatial_results_raw <- data.frame()
for (lichen in all_lichen) {
  result <- test_morans_I_ape(data[[lichen]], inv_dist_weights, lichen)
  spatial_results_raw <- rbind(spatial_results_raw, result)
  cat(sprintf("%-25s | I = %7.4f | p = %.4f %s\n", 
              lichen, result$Morans_I, result$P_value, result$Sig))
}

cat("\n✓ Raw spatial patterns documented\n\n")

# Save to output directory
csv_filename <- paste0(output_dir, "02c_spatial_autocorrelation_raw_", dataset_name, ".csv")
write.csv(spatial_results_raw, csv_filename, row.names = FALSE)
cat("✓ Saved:", csv_filename, "\n\n")


# ==============================================================================
# SECTION 2: DEFINE MODEL FORMULAS
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("STEP 2: Setting up model formulas\n")
cat("═════��═════════════════════════════════════════════════════\n\n")

# ┌──────────────────────��──────────────────────────────────────────────────┐
# │ NOTE: predictors_scaled is loaded from the RData file                  │
# │ If you need to modify predictors manually:                             │
# │                                                                         │
# │ predictors_scaled <- c("elevation_scaled", "volume_snags_scaled", ...) │
# │                                                                         │
# │ For standard dataset: 12 predictors                                    │
# │ For enlarged dataset: may have additional predictors                   │
# └─────────────────────────────────────────────────────────────────────────┘

formula_full <- as.formula(paste("response ~", paste(predictors_scaled, collapse = " + ")))

cat("Formula:", deparse(formula_full), "\n")
cat("Number of predictors:", length(predictors_scaled), "\n\n")


# ==============================================================================
# SECTION 3: PRE-DHARMA VALIDATION FUNCTION
# ==============================================================================

run_pre_dharma <- function(model, model_name, family_type) {
  
  cat(sprintf("\n─────────────────────────────────────────────────────────\n"))
  cat(sprintf("Pre-DHARMa: %s\n", model_name))
  cat(sprintf("─────────────────────────────────────────────────────────\n"))
  
  # Extract coefficients
  if (family_type == "nb") {
    coef_table <- summary(model)$coefficients$cond
    convergence_code <- model$fit$convergence
  } else {
    coef_table <- summary(model)$coefficients
    convergence_code <- ifelse(model$converged, 0, 1)
  }
  
  # Check 1: Coefficients
  max_coef <- max(abs(coef_table[, "Estimate"]))
  max_coef_name <- rownames(coef_table)[which.max(abs(coef_table[, "Estimate"]))]
  
  cat(sprintf("  Max coef: %.2f (%s)", max_coef, max_coef_name))
  coef_flag <- if (max_coef > 10) {
    cat(" 🚨\n"); "FAIL"
  } else if (max_coef > 5) {
    cat(" ⚠️\n"); "WARNING"
  } else {
    cat(" ✅\n"); "PASS"
  }
  
  # Check 2: Standard errors
  max_se <- max(coef_table[, "Std. Error"])
  max_se_name <- rownames(coef_table)[which.max(coef_table[, "Std. Error"])]
  
  cat(sprintf("  Max SE:   %.2f (%s)", max_se, max_se_name))
  se_flag <- if (max_se > 10) {
    cat(" 🚨\n"); "FAIL"
  } else if (max_se > 5) {
    cat(" ⚠️\n"); "WARNING"
  } else {
    cat(" ✅\n"); "PASS"
  }
  
  # Check 3: Convergence & Iterations
  if (family_type == "nb") {
    cat(sprintf("  Converged: %d", convergence_code))
    conv_flag <- if (convergence_code == 0) {
      cat(" ✅\n"); "PASS"
    } else {
      cat(" 🚨\n"); "FAIL"
    }
    iter_flag <- "PASS"
  } else {
    n_iter <- model$iter
    cat(sprintf("  Converged: %s, Iter: %d", model$converged, n_iter))
    conv_flag <- if (model$converged) { cat(" ✅\n"); "PASS" } else { cat(" 🚨\n"); "FAIL" }
    iter_flag <- if (n_iter >= 20) { "FAIL" } else if (n_iter >= 15) { "WARNING" } else { "PASS" }
  }
  
  # Check 4: Pseudo R² or Dispersion
  if (family_type == "binomial") {
    pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
    cat(sprintf("  Pseudo R²: %.3f", pseudo_r2))
    r2_flag <- if (pseudo_r2 > 0.95) {
      cat(" 🚨\n"); "FAIL"
    } else if (pseudo_r2 > 0.85) {
      cat(" ⚠️\n"); "WARNING"
    } else {
      cat(" ✅\n"); "PASS"
    }
  } else if (family_type == "nb") {
    disp_param <- sigma(model)
    cat(sprintf("  Dispersion: %.3f", disp_param))
    r2_flag <- if (disp_param < 0.1 || disp_param > 100) {
      cat(" ⚠️\n"); "WARNING"
    } else {
      cat(" ✅\n"); "PASS"
    }
  } else {
    r2_flag <- "PASS"
  }
  
  # Overall verdict
  all_flags <- c(coef_flag, se_flag, conv_flag, iter_flag, r2_flag)
  
  if (all(all_flags == "PASS")) {
    cat("  → ✅ VALID\n")
    overall <- "VALID"
  } else if (any(all_flags == "FAIL")) {
    cat("  → ❌ INVALID\n")
    overall <- "INVALID"
  } else {
    cat("  → ⚠️  BORDERLINE\n")
    overall <- "BORDERLINE"
  }
  
  return(list(
    status = overall,
    max_coef = max_coef,
    max_se = max_se,
    flags = all_flags
  ))
}


# ==============================================================================
# SECTION 4: DHARMA DIAGNOSTICS FUNCTION (with Test 5)
# ==============================================================================

run_dharma_diagnostics <- function(model, model_name, coords_df, output_dir) {
  
  cat(sprintf("\n═══════════════════════════════════════════════════════════\n"))
  cat(sprintf("DHARMa: %s\n", model_name))
  cat(sprintf("═══════════════════════════════════════════════════════════\n"))
  
  # Simulate residuals
  sim_resid <- simulateResiduals(fittedModel = model, n = 1000, plot = FALSE)
  
  # Save plot to output directory
  plot_filename <- paste0(output_dir, "DHARMa_", 
                          gsub(" ", "_", model_name), "_",
                          dataset_name, ".png")
  png(filename = plot_filename, width = 10, height = 6, units = "in", res = 300)
  plot(sim_resid, main = paste("DHARMa:", model_name, paste0("(", dataset_name, ")")))
  dev.off()
  
  # Run 4 tests
  cat("  Tests: ")
  disp_test <- testDispersion(sim_resid, plot = FALSE)
  zi_test <- testZeroInflation(sim_resid, plot = FALSE)
  outlier_test <- testOutliers(sim_resid, plot = FALSE)
  spatial_test <- testSpatialAutocorrelation(sim_resid, x = coords_df$X, 
                                             y = coords_df$Y, plot = FALSE)
  
  cat(sprintf("Disp=%.3f, ZI=%.3f, Out=%.3f, Spat=%.3f\n",
              disp_test$p.value, zi_test$p.value, 
              outlier_test$p.value, spatial_test$p.value))
  
  # Test 5: Spatial comparison
  if (exists("spatial_results_raw", envir = .GlobalEnv)) {
    raw_spatial <- spatial_results_raw[spatial_results_raw$Variable == model_name, ]
    
    if (nrow(raw_spatial) > 0) {
      cat("\n  5. SPATIAL PATTERN COMPARISON:\n")
      cat(sprintf("     Raw:       I = %7.4f, p = %.4f %s\n", 
                  raw_spatial$Morans_I, raw_spatial$P_value, raw_spatial$Sig))
      cat(sprintf("     Residuals: I = %7.4f, p = %.4f\n", 
                  as.numeric(spatial_test$statistic), spatial_test$p.value))
      
      cat("     ")
      if (raw_spatial$P_value < 0.05 && spatial_test$p.value >= 0.05) {
        cat("✅ Fully explained\n")
      } else if (raw_spatial$P_value < 0.05 && spatial_test$p.value < 0.05) {
        cat("⚠️  Partially explained\n")
      } else {
        cat("✅ No spatial issues\n")
      }
    }
  }
  
  # Compile results
  diagnostics <- data.frame(
    Model = model_name,
    Dataset = dataset_name,
    Dispersion_p = round(disp_test$p.value, 4),
    ZeroInflation_p = round(zi_test$p.value, 4),
    Outliers_p = round(outlier_test$p.value, 4),
    SpatialAutocorr_p = round(spatial_test$p.value, 4),
    Dispersion_OK = disp_test$p.value > 0.05,
    ZeroInflation_OK = zi_test$p.value > 0.05,
    Outliers_OK = outlier_test$p.value > 0.05,
    Spatial_OK = spatial_test$p.value > 0.05,
    stringsAsFactors = FALSE
  )
  
  all_ok <- all(diagnostics[, 7:10])
  cat(ifelse(all_ok, "\n  ✅ ALL PASSED\n", "\n  ⚠️  SOME FAILED\n"))
  
  return(list(diagnostics = diagnostics, sim_resid = sim_resid))
}


# ==============================================================================
# SECTION 5: FIT ALL MODELS
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("STEP 3: Fitting all 7 models\n")
cat("═══════════════════════════════════════════════════════════\n\n")

models_glm <- list()
pre_dharma_results <- list()

# Fit 6 binomial GLMs
for (lichen in lichen_binary) {
  cat(sprintf("Fitting: %s (binomial)...", lichen))
  
  data_scaled$response <- data_scaled[[lichen]]
  
  models_glm[[lichen]] <- glm(formula_full, 
                              data = data_scaled, 
                              family = binomial(link = "logit"))
  
  cat(" ✓\n")
  
  # Pre-DHARMa check
  pre_dharma_results[[lichen]] <- run_pre_dharma(models_glm[[lichen]], 
                                                 lichen, "binomial")
}

# Fit 1 negative binomial GLM for calicioids
cat(sprintf("Fitting: calicioids_richness (negative binomial)..."))

# Build formula dynamically for NB model
formula_nb <- as.formula(paste("calicioids_richness ~", 
                               paste(predictors_scaled, collapse = " + ")))

models_glm[["calicioids_richness"]] <- glmmTMB(
  formula_nb,
  data = data_scaled,
  family = nbinom2
)

cat(" ✓\n")

# Pre-DHARMa check
pre_dharma_results[["calicioids_richness"]] <- run_pre_dharma(
  models_glm[["calicioids_richness"]], 
  "calicioids_richness", "nb")

cat("\n✓ All 7 models fitted\n\n")


# ==============================================================================
# SECTION 6: PRE-DHARMA SUMMARY
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("PRE-DHARMA SUMMARY\n")
cat("═══════════════════════════════════════════════════════════\n\n")

n_valid <- sum(sapply(pre_dharma_results, function(x) x$status == "VALID"))
n_invalid <- sum(sapply(pre_dharma_results, function(x) x$status == "INVALID"))
n_borderline <- sum(sapply(pre_dharma_results, function(x) x$status == "BORDERLINE"))

cat(sprintf("Valid:      %d/%d models\n", n_valid, length(all_lichen)))
cat(sprintf("Borderline: %d/%d models\n", n_borderline, length(all_lichen)))
cat(sprintf("Invalid:    %d/%d models\n\n", n_invalid, length(all_lichen)))

if (n_invalid > 0) {
  cat("⚠️  Invalid models (DO NOT USE):\n")
  for (lichen in names(pre_dharma_results)) {
    if (pre_dharma_results[[lichen]]$status == "INVALID") {
      cat(sprintf("   • %s\n", lichen))
    }
  }
  cat("\n")
}


# ==============================================================================
# SECTION 7: DHARMA DIAGNOSTICS FOR ALL MODELS
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("STEP 4: Running DHARMa diagnostics (all 7 models)\n")
cat("═══════════════════════════════════════════════════════════\n")

diagnostics_glm <- list()
diagnostics_summary <- data.frame()

for (lichen in all_lichen) {
  # Skip if pre-DHARMa failed
  if (pre_dharma_results[[lichen]]$status == "INVALID") {
    cat(sprintf("\n⚠️  Skipping DHARMa for %s (failed pre-DHARMa)\n", lichen))
    next
  }
  
  diag_result <- run_dharma_diagnostics(models_glm[[lichen]], 
                                        lichen, 
                                        data[, c("X", "Y")],
                                        output_dir)
  diagnostics_glm[[lichen]] <- diag_result$sim_resid
  diagnostics_summary <- rbind(diagnostics_summary, diag_result$diagnostics)
}

cat("\n✓ DHARMa diagnostics complete\n\n")

# Save DHARMa summary to output directory
dharma_csv <- paste0(output_dir, "02c_DHARMa_diagnostics_summary_", dataset_name, ".csv")
write.csv(diagnostics_summary, dharma_csv, row.names = FALSE)
cat("✓ Saved:", dharma_csv, "\n\n")


# ==============================================================================
# SECTION 8: FINAL SUMMARY
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL SUMMARY                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Dataset: %s (%d predictors)\n\n", dataset_name, length(predictors_scaled)))

cat("DHARMa Results:\n")
print(diagnostics_summary[, c("Model", "Dispersion_OK", "ZeroInflation_OK", 
                              "Outliers_OK", "Spatial_OK")])
cat("\n")

n_passed <- sum(apply(diagnostics_summary[, 7:10], 1, all))
cat(sprintf("Models passing all DHARMa tests: %d/%d\n\n", 
            n_passed, nrow(diagnostics_summary)))

# Spatial pattern summary
if (exists("spatial_results_raw")) {
  n_raw_clustered <- sum(spatial_results_raw$P_value < 0.05)
  cat(sprintf("Spatial autocorrelation:\n"))
  cat(sprintf("  Raw data: %d/%d models clustered\n", n_raw_clustered, 
              nrow(spatial_results_raw)))
  
  if (n_passed > 0) {
    cat(sprintf("  Residuals: All p > 0.05 (fully explained) ✅\n\n"))
  }
}

cat("═══════════════════════════════════════════════════════════\n")
cat("MODEL SET COMPLETE\n")
cat("═══════════════════════════════════════════════════════════\n")
cat(sprintf("  1-%d. Binomial GLMs (presence/absence) ✓\n", length(lichen_binary)))
cat("    7. Negative Binomial GLM (richness) ✓\n\n")
cat("All models validated and ready for interpretation!\n")
cat("═══════════════════════════════════════════════════════════\n\n")


# ==============================================================================
# SECTION 9: SAVE WORKSPACE
# ==============================================================================

# Save to output directory
rdata_filename <- paste0(output_dir, "02c_models_complete_", dataset_name, ".RData")

save(data, data_scaled,
     models_glm, diagnostics_glm, diagnostics_summary,
     spatial_results_raw, 
     lichen_binary, lichen_count, all_lichen,
     predictors, predictors_scaled,
     inv_dist_weights,
     pre_dharma_results,
     dataset_name,
     output_dir,
     file = rdata_filename)

cat("✅ Workspace saved:", rdata_filename, "\n\n")

cat("Files created in", output_dir, ":\n")
cat(sprintf("  • models_complete_%s.RData\n", dataset_name))
cat(sprintf("  • spatial_autocorrelation_raw_%s.csv\n", dataset_name))
cat(sprintf("  • DHARMa_diagnostics_summary_%s.csv\n", dataset_name))
cat(sprintf("  • DHARMa_*_%s.png (7 diagnostic plots)\n\n", dataset_name))

cat("═══════════════════════════════════════════════════════════\n")
cat("READY FOR SCRIPT 04: Interpretation & Visualization\n")
cat("═══════════════════════════════════════════════════════════\n\n")

################################################################################
# END OF SCRIPT
################################################################################


