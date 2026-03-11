################################################################################
# ŠUMAVA NATIONAL PARK LICHEN-STRUCTURE MODELING
# Script 02: Spatial Autocorrelation Testing and GLMM Fitting (SIMPLIFIED)
# Author: Jonas Dotterweich
# Date: 2026-02-14
# Purpose: Test spatial structure and fit GLMMs for all 7 lichen groups
################################################################################

# ==============================================================================
# LOAD WORKSPACE FROM PREVIOUS SCRIPT
# ==============================================================================

load("Lichens/to_model/01_prepared_data_for_modeling.RData") ### for 12 predictors
# or 
load("Lichens/to_model/01_prepared_data_for_modeling(14pred).RData") ### for 14 predictors (with decay4 & decay5)
#or
load("Lichens/to_model/01_prepared_data_for_modeling(12pred_wo_mycob).RData") ### wo myco and wo tree median hight and regen

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

library(tidyverse)
library(ape)           # For Moran's I test (reliable method)
library(ncf)           # For spatial correlograms
library(DHARMa)        # For model diagnostics
library(performance)   # For model checking

cat("✅ All required packages loaded successfully\n\n")


# ==============================================================================
# SECTION 1: SPATIAL AUTOCORRELATION TESTING - RAW RESPONSES
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 1: TESTING FOR SPATIAL AUTOCORRELATION IN RESPONSES   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Extract coordinates
coords <- as.matrix(data[, c("X", "Y")])

cat("Creating spatial distance matrix...\n")

# Create Euclidean distance matrix
dist_matrix <- as.matrix(dist(coords))

# Inverse distance weights (closer plots have higher weights)
# Add small constant to avoid division by zero
inv_dist_weights <- 1 / (dist_matrix + 1)  # +1 meter to avoid infinity
diag(inv_dist_weights) <- 0  # No self-weighting

cat("  ✓ Distance matrix created\n")
cat("  ✓ Inverse distance weights calculated\n")
cat("  - Matrix dimensions:", dim(inv_dist_weights), "\n\n")


# ==============================================================================
# Function to test Moran's I using ape package
# ==============================================================================

test_morans_I_ape <- function(variable, weight_matrix, var_name) {
  
  # Remove any NA values (though we shouldn't have any)
  complete_idx <- which(complete.cases(variable))
  var_clean <- variable[complete_idx]
  weight_clean <- weight_matrix[complete_idx, complete_idx]
  
  # Perform Moran's I test using ape
  moran_result <- ape::Moran.I(var_clean, weight_clean)
  
  # Extract results
  morans_i <- moran_result$observed
  expected_i <- moran_result$expected
  sd_i <- moran_result$sd
  p_value <- moran_result$p.value
  
  # Calculate z-score
  z_score <- (morans_i - expected_i) / sd_i
  
  # Interpret significance
  significance <- ifelse(p_value < 0.001, "***",
                         ifelse(p_value < 0.01, "**",
                                ifelse(p_value < 0.05, "*", "ns")))
  
  interpretation <- ifelse(morans_i > expected_i & p_value < 0.05, 
                           "POSITIVE clustering (similar values cluster together)",
                           ifelse(morans_i < expected_i & p_value < 0.05,
                                  "NEGATIVE clustering (dissimilar values cluster together)",
                                  "NO spatial autocorrelation detected"))
  
  # Return as data frame
  result <- data.frame(
    Variable = var_name,
    Morans_I = round(morans_i, 4),
    Expected_I = round(expected_i, 4),
    SD = round(sd_i, 4),
    Z_score = round(z_score, 3),
    P_value = round(p_value, 4),
    Sig = significance,
    Interpretation = interpretation,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


# ==============================================================================
# Test all lichen responses for spatial autocorrelation
# ==============================================================================

cat("═══════════════════════════════════════════════════════════\n")
cat("MORAN'S I TEST FOR SPATIAL AUTOCORRELATION IN LICHEN DATA\n")
cat("(Using inverse distance weighting method)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

spatial_results_raw <- data.frame()

for (lichen in all_lichen) {
  result <- test_morans_I_ape(data[[lichen]], inv_dist_weights, lichen)
  spatial_results_raw <- rbind(spatial_results_raw, result)
  
  cat(sprintf("%-25s | Moran's I = %7.4f | Z = %6.2f | p = %.4f %-3s | %s\n",
              lichen,
              result$Morans_I,
              result$Z_score,
              result$P_value,
              result$Sig,
              result$Interpretation))
}

cat("\n")

# Save results
#write.csv(spatial_results_raw, 
#          "02_spatial_autocorrelation_raw_responses.csv", 
#          row.names = FALSE)

#cat("✅ Spatial autocorrelation results saved to CSV\n\n")


# ==============================================================================
# SECTION 2: VISUALIZE SPATIAL PATTERNS
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("CREATING SPATIAL VISUALIZATIONS\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Plot 1: Spatial distribution of core_ogf_presence
p_spatial <- ggplot(data, aes(x = X, y = Y, color = factor(core_ogf_presence))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("0" = "gray70", "1" = "forestgreen"),
                     labels = c("Absent", "Present"),
                     name = "Core OGF Lichens") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Spatial Distribution of Core Old-Growth Forest Lichens",
       subtitle = "Šumava National Park (n=120 plots)",
       x = "Easting (X)", y = "Northing (Y)") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

print(p_spatial)
ggsave("02_spatial_pattern_core_ogf.png", width = 8, height = 7, dpi = 300)
cat("  ✓ Spatial pattern plot saved\n")


# Plot 2: Spatial correlogram (Moran's I at different distance lags)
cat("\nCalculating spatial correlogram for core_ogf_presence...\n")
cat("  (This may take a minute with bootstrap resampling...)\n")

tryCatch({
  
  # Calculate correlogram using ncf package
  correlogram <- ncf::correlog(x = data$X, 
                               y = data$Y, 
                               z = data$core_ogf_presence,
                               increment = 1000,  # Distance increment in meters
                               resamp = 100,      # Bootstrap resamples
                               quiet = TRUE)
  
  # Save plot
  png("02_spatial_correlogram_core_ogf.png", width = 10, height = 6, 
      units = "in", res = 300)
  
  plot(correlogram$mean.of.class/1000,  # Convert to km
       correlogram$correlation,
       type = "b",
       pch = 16,
       xlab = "Distance (km)",
       ylab = "Moran's I",
       main = "Spatial Correlogram: Core OGF Presence",
       sub = "Red points indicate p < 0.05",
       ylim = c(min(-0.2, min(correlogram$correlation, na.rm = TRUE)), 
                max(0.3, max(correlogram$correlation, na.rm = TRUE))),
       cex = 1.2,
       lwd = 2,
       col = "black")
  
  abline(h = 0, lty = 2, col = "gray50", lwd = 2)
  
  # Add significance markers
  sig_points <- which(correlogram$p < 0.05)
  if (length(sig_points) > 0) {
    points(correlogram$mean.of.class[sig_points]/1000, 
           correlogram$correlation[sig_points],
           pch = 16, col = "red", cex = 2.5)
  }
  
  # Add grid
  grid(col = "gray90")
  
  legend("topright", 
         legend = c("Moran's I", "Significant (p<0.05)", "No autocorrelation"),
         pch = c(16, 16, NA),
         lty = c(1, NA, 2),
         col = c("black", "red", "gray50"),
         cex = 1,
         bg = "white")
  
  dev.off()
  
  cat("  ✓ Correlogram calculated and saved\n\n")
  
}, error = function(e) {
  cat("  ⚠️  Correlogram calculation failed (optional analysis):\n")
  cat("     ", e$message, "\n\n")
  correlogram <<- NULL
})


# ==============================================================================
# SECTION 3: DECIDE ON MODELING STRATEGY
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 2: DETERMINING MODELING APPROACH                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Check if ANY lichen shows significant spatial autocorrelation
has_spatial_autocorr <- any(spatial_results_raw$P_value < 0.05)

if (has_spatial_autocorr) {
  cat("⚠️  SPATIAL AUTOCORRELATION DETECTED in one or more lichen groups!\n\n")
  
  # Show which ones
  spatial_lichen <- spatial_results_raw[spatial_results_raw$P_value < 0.05, ]
  cat("Lichen groups with significant spatial structure:\n")
  print(spatial_lichen[, c("Variable", "Morans_I", "Z_score", "P_value", "Sig")])
  
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  MODELING OPTIONS:                                           ║\n")
  cat("╠══════════════════════════════════════════════════════════════╣\n")
  cat("║  A) Create spatial blocks → Use (1|block) random effect     ║\n")
  cat("║  B) Use GAMs with spatial smoothers → s(X, Y)               ║\n")
  cat("║  C) Fit GLMs first, check residuals → Add RE if needed      ║\n")
  cat("╚════════════════════════════════════════════��═════════════════╝\n\n")
  cat("RECOMMENDED: Start with Option C\n")
  cat("  1. Fit standard GLMs with all 12 predictors\n")
  cat("  2. Check DHARMa diagnostics for residual autocorrelation\n")
  cat("  3. If residuals still show autocorrelation → create spatial blocks\n\n")
  
} else {
  cat("✅ NO SPATIAL AUTOCORRELATION detected in raw lichen responses\n\n")
  cat("RECOMMENDED: Use standard GLMs (no random effects)\n")
  cat("  - Simpler models, easier interpretation\n")
  cat("  - Still check residual diagnostics to confirm!\n\n")
}


# ==============================================================================
# SECTION 4: CREATE SPATIAL BLOCKS (FOR POTENTIAL USE IN GLMMs)
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════\n")
cat("CREATING SPATIAL BLOCKS (k-means clustering)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# Use k-means to create spatial blocks
set.seed(123)  # For reproducibility
n_blocks <- 8  # Number of spatial blocks (n/15 rule: 120/15 = 8)

cat("Performing k-means clustering on coordinates...\n")

# Perform k-means clustering
kmeans_result <- kmeans(coords, centers = n_blocks, nstart = 25, iter.max = 100)

# Add to both datasets
data$spatial_block <- as.factor(kmeans_result$cluster)
data_scaled$spatial_block <- data$spatial_block

cat("  ✓ K-means clustering complete\n\n")
cat("Spatial blocks summary:\n")
cat("  - Number of blocks:", n_blocks, "\n")
cat("  - Plots per block:\n")
block_table <- table(data$spatial_block)
print(block_table)
cat("\n")

# Check balance
cat("  - Min plots per block:", min(block_table), "\n")
cat("  - Max plots per block:", max(block_table), "\n")
cat("  - Mean plots per block:", round(mean(block_table), 1), "\n\n")

# Visualize spatial blocks
p_blocks <- ggplot(data, aes(x = X, y = Y, color = spatial_block, label = spatial_block)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_brewer(palette = "Set2", name = "Spatial\nBlock") +
  coord_equal() +
  theme_minimal() +
  labs(title = "K-means Spatial Blocks (k=8)",
       subtitle = "For potential use as random effects in GLMMs",
       x = "Easting (X)", y = "Northing (Y)") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

print(p_blocks)
ggsave("02_spatial_blocks_kmeans.png", width = 9, height = 7, dpi = 300)
cat("  ✓ Spatial blocks visualization saved\n\n")


# ==============================================================================
# SECTION 5: FIT GLMs FOR ALL LICHEN GROUPS (ALL 12 PREDICTORS)
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 3: FITTING GLMs (NO RANDOM EFFECTS - BASELINE MODELS) ║\n")
cat("╚═════════════════════��════════════════════════════════════════╝\n\n")

# Create formula with all 12 scaled predictors
formula_full <- as.formula(paste("response ~", 
                                 paste(predictors_scaled, collapse = " + ")))

cat("Model formula:\n")
cat("  response ~ ", paste(predictors_scaled, collapse = " + "), "\n\n")

cat("Predictors included (n=14):\n") #11-14 predictors depending on input from above
for (i in seq_along(predictors)) {
  cat(sprintf("  %2d. %-20s (scaled: %s)\n", 
              i, predictors[i], predictors_scaled[i]))
}
cat("\n")

# Storage for models
models_glm <- list()
model_summaries <- list()


# ==============================================================================
# Function to fit GLM and extract key info
# ==============================================================================

fit_lichen_glm <- function(lichen_name, data, family_type) {
  
  cat(sprintf("\n═══════════════════════════════════════════════════════════\n"))
  cat(sprintf("Fitting GLM: %s\n", lichen_name))
  cat(sprintf("═══════════════════════════════════════════════════════════\n"))
  
  # Prepare response variable
  data$response <- data[[lichen_name]]
  
  # Check prevalence/distribution
  if (family_type == "binomial") {
    n_present <- sum(data$response)
    n_absent <- sum(data$response == 0)
    prevalence <- mean(data$response) * 100
    cat(sprintf("  Response distribution: %d present, %d absent (%.1f%% prevalence)\n", 
                n_present, n_absent, prevalence))
    
    # Warn if prevalence is very low
    if (prevalence < 10 || prevalence > 90) {
      cat("  ⚠️  Low/high prevalence may affect model stability\n")
    }
    
  } else if (family_type == "poisson") {
    cat("  Richness distribution:\n")
    richness_summary <- summary(data$response)
    print(richness_summary)
    n_zeros <- sum(data$response == 0)
    pct_zeros <- mean(data$response == 0) * 100
    cat(sprintf("  Zero-inflation: %d plots with zero richness (%.1f%%)\n", 
                n_zeros, pct_zeros))
  }
  
  # Fit model
  cat("\n  Fitting model...\n")
  
  if (family_type == "binomial") {
    model <- glm(formula_full, data = data, family = binomial(link = "logit"))
  } else if (family_type == "poisson") {
    model <- glm(formula_full, data = data, family = poisson(link = "log"))
  }
  
  cat("  ✓ Model fitted\n\n")
  
  # Model summary
  model_sum <- summary(model)
  print(model_sum)
  
  # Model fit statistics
  cat("\n─────────────────────────────────────────────────────────────\n")
  cat("MODEL FIT STATISTICS:\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat(sprintf("  AIC: %.2f\n", AIC(model)))
  cat(sprintf("  Log-Likelihood: %.2f\n", logLik(model)))
  cat(sprintf("  Null deviance: %.2f (df = %d)\n", 
              model$null.deviance, model$df.null))
  cat(sprintf("  Residual deviance: %.2f (df = %d)\n", 
              model$deviance, model$df.residual))
  
  # Pseudo R-squared for binomial
  if (family_type == "binomial") {
    null_deviance <- model$null.deviance
    residual_deviance <- model$deviance
    pseudo_r2 <- 1 - (residual_deviance / null_deviance)
    cat(sprintf("  McFadden's Pseudo R²: %.3f\n", pseudo_r2))
  }
  
  # Check overdispersion
  residual_deviance <- model$deviance
  df_residual <- model$df.residual
  dispersion_ratio <- residual_deviance / df_residual
  
  cat(sprintf("  Dispersion ratio (deviance/df): %.3f\n", dispersion_ratio))
  
  if (family_type == "poisson" && dispersion_ratio > 2) {
    cat("  ⚠️  OVERDISPERSION detected! Consider negative binomial family.\n")
  } else if (dispersion_ratio > 1.5) {
    cat("  ⚠️  Slight overdispersion present\n")
  } else {
    cat("  ✓  Dispersion acceptable\n")
  }
  
  cat("\n")
  
  return(model)
}


# ==============================================================================
# Fit GLMs for all BINARY lichen groups (6 presence/absence)
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FITTING BINOMIAL GLMs (6 presence/absence groups)          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")

for (lichen in lichen_binary) {
  models_glm[[lichen]] <- fit_lichen_glm(lichen, data_scaled, "binomial")
}


# ==============================================================================
# Fit GLM for COUNT lichen group (calicioids richness)
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FITTING POISSON GLM (calicioids richness)                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")

models_glm[["calicioids_richness"]] <- fit_lichen_glm("calicioids_richness", 
                                                      data_scaled, 
                                                      "poisson")




# ==============================================================================
# PRE-DHARMA CHECKS: Convergence and Coefficient Sanity
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PRE-DHARMA VALIDATION: Check for Model Problems           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

for (lichen in all_lichen) {
  model <- models_glm[[lichen]]
  
  cat(sprintf("%-25s: ", lichen))
  
  # Check 1: Coefficient magnitudes
  max_coef <- max(abs(coef(model)))
  max_se <- max(summary(model)$coefficients[, "Std. Error"])
  
  # Check 2: Pseudo R²
  if (inherits(model, "glm") && model$family$family == "binomial") {
    pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
  } else {
    pseudo_r2 <- NA
  }
  
  # Check 3: Convergence
  converged <- model$converged
  
  # Flags
  flags <- c()
  if (max_coef > 10) flags <- c(flags, "Huge coef")
  if (max_se > 10) flags <- c(flags, "Huge SE")
  if (!is.na(pseudo_r2) && pseudo_r2 > 0.95) flags <- c(flags, "Perfect fit")
  if (!converged) flags <- c(flags, "No convergence")
  
  if (length(flags) == 0) {
    cat("✅ OK\n")
  } else {
    cat(sprintf("🚨 WARNING: %s\n", paste(flags, collapse = ", ")))
  }
}

cat("\n")


# ==============================================================================
# ENHANCED PRE-DHARMA CHECKS
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PRE-DHARMA VALIDATION: Comprehensive Model Checks          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

for (lichen in all_lichen) {
  model <- models_glm[[lichen]]
  
  cat(sprintf("\n%-25s\n", lichen))
  cat("────────────────────────────\n")
  
  # Extract coefficient table
  coef_table <- summary(model)$coefficients
  
  # Check 1: Maximum coefficient
  max_coef <- max(abs(coef_table[, "Estimate"]))
  max_coef_name <- rownames(coef_table)[which.max(abs(coef_table[, "Estimate"]))]
  
  cat(sprintf("  Max coefficient: %.2f (%s)", max_coef, max_coef_name))
  if (max_coef > 10) {
    cat(" 🚨 HUGE!\n")
  } else if (max_coef > 5) {
    cat(" ⚠️ Large\n")
  } else {
    cat(" ✓\n")
  }
  
  # Check 2: Maximum SE
  max_se <- max(coef_table[, "Std. Error"])
  max_se_name <- rownames(coef_table)[which.max(coef_table[, "Std. Error"])]
  
  cat(sprintf("  Max SE: %.2f (%s)", max_se, max_se_name))
  if (max_se > 100) {
    cat(" 🚨 ENORMOUS!\n")
  } else if (max_se > 10) {
    cat(" 🚨 Very large\n")
  } else if (max_se > 5) {
    cat(" ⚠️ Large\n")
  } else {
    cat(" ✓\n")
  }
  
  # Check 3: Pseudo R² (binomial only)
  if (inherits(model, "glm") && model$family$family == "binomial") {
    pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
    cat(sprintf("  Pseudo R²: %.3f", pseudo_r2))
    if (pseudo_r2 > 0.95) {
      cat(" 🚨 Too perfect!\n")
    } else if (pseudo_r2 > 0.8) {
      cat(" ⚠️ Very high\n")
    } else {
      cat(" ✓\n")
    }
  }
  
  # Check 4: Convergence
  cat(sprintf("  Converged: %s", model$converged))
  if (!model$converged) {
    cat(" 🚨 FAILED!\n")
  } else {
    cat(" ✓\n")
  }
  
  # Check 5: Iterations
  n_iter <- model$iter
  cat(sprintf("  Iterations: %d", n_iter))
  if (n_iter >= 20) {
    cat(" 🚨 Too many!\n")
  } else if (n_iter >= 15) {
    cat(" ⚠️ Many\n")
  } else {
    cat(" ✓\n")
  }
  
  # Check 6: Any non-significant coefficients with huge estimates?
  nonsig_huge <- coef_table[coef_table[, "Pr(>|z|)"] > 0.5 & 
                              abs(coef_table[, "Estimate"]) > 5, ]
  if (nrow(nonsig_huge) > 0) {
    cat("  🚨 Large non-significant coefficients detected!\n")
    cat(sprintf("     %s\n", paste(rownames(nonsig_huge), collapse=", ")))
  }
  
  # Overall verdict
  has_issues <- max_coef > 10 || max_se > 10 || 
    (exists("pseudo_r2") && pseudo_r2 > 0.95) || 
    !model$converged || n_iter >= 20
  
  if (has_issues) {
    cat("\n  ⚠️ VERDICT: MODEL HAS ISSUES - CHECK CAREFULLY!\n")
  } else {
    cat("\n  ✅ VERDICT: Model appears healthy\n")
  }
}

cat("\n")



####----
# ==============================================================================
# ENHANCED PRE-DHARMA CHECKS: Comprehensive Validation
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PRE-DHARMA VALIDATION: Comprehensive Model Checks          ║\n")
cat("╚════════���═════════════════════════════════════════════════════╝\n\n")

cat("Checking for separation, convergence issues, and numerical instability...\n\n")

# Initialize results dataframe
pre_dharma_results <- data.frame(
  Model = character(),
  Max_Coef = numeric(),
  Problematic_Coef = character(),
  Max_SE = numeric(),
  Problematic_SE = character(),
  Pseudo_R2 = numeric(),
  Converged = character(),
  Iterations = integer(),
  Coef_Flag = character(),
  SE_Flag = character(),
  R2_Flag = character(),
  Conv_Flag = character(),
  Iter_Flag = character(),
  Overall_Status = character(),
  stringsAsFactors = FALSE
)

# Loop through all models
for (lichen in all_lichen) {
  model <- models_glm[[lichen]]
  
  cat(sprintf("═══════════════════════════════════════════════════════════\n"))
  cat(sprintf("Model: %s\n", lichen))
  cat(sprintf("═══════════════════════════════════════════════════════════\n"))
  
  # Extract coefficient table
  coef_table <- summary(model)$coefficients
  
  # ============================================================================
  # CHECK 1: Coefficient Magnitudes
  # ============================================================================
  max_coef <- max(abs(coef_table[, "Estimate"]))
  max_coef_name <- rownames(coef_table)[which.max(abs(coef_table[, "Estimate"]))]
  
  # Find ALL problematic coefficients (not just max)
  problem_coefs <- rownames(coef_table)[abs(coef_table[, "Estimate"]) > 10]
  
  cat(sprintf("  Max coefficient: %.2f (%s)\n", max_coef, max_coef_name))
  
  if (length(problem_coefs) > 0) {
    cat(sprintf("  🚨 Problem coefficients (|coef| > 10): %s\n", 
                paste(problem_coefs, collapse=", ")))
    coef_flag <- "FAIL"
  } else if (max_coef > 5) {
    cat("  ⚠️  Large coefficient (5-10 range) - acceptable but monitor\n")
    coef_flag <- "WARNING"
  } else {
    cat("  ✓  Coefficients in normal range\n")
    coef_flag <- "PASS"
  }
  
  # ============================================================================
  # CHECK 2: Standard Errors
  # ============================================================================
  max_se <- max(coef_table[, "Std. Error"])
  max_se_name <- rownames(coef_table)[which.max(coef_table[, "Std. Error"])]
  
  # Find ALL problematic SEs
  problem_ses <- rownames(coef_table)[coef_table[, "Std. Error"] > 10]
  
  cat(sprintf("  Max SE: %.2f (%s)\n", max_se, max_se_name))
  
  if (length(problem_ses) > 0) {
    cat(sprintf("  🚨 Problem SEs (SE > 10): %s\n", 
                paste(problem_ses, collapse=", ")))
    se_flag <- "FAIL"
  } else if (max_se > 5) {
    cat("  ⚠️  Large SE (5-10 range) - acceptable but uncertain\n")
    se_flag <- "WARNING"
  } else {
    cat("  ✓  Standard errors in normal range\n")
    se_flag <- "PASS"
  }
  
  # ============================================================================
  # CHECK 3: Pseudo R-squared (binomial only)
  # ============================================================================
  if (inherits(model, "glm") && model$family$family == "binomial") {
    pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
    cat(sprintf("  Pseudo R²: %.3f\n", pseudo_r2))
    
    if (pseudo_r2 > 0.95) {
      cat("  🚨 Perfect/near-perfect fit - likely separation!\n")
      r2_flag <- "FAIL"
    } else if (pseudo_r2 > 0.85) {
      cat("  ⚠️  Very high R² - check for overfitting\n")
      r2_flag <- "WARNING"
    } else {
      cat("  ✓  R² in acceptable range\n")
      r2_flag <- "PASS"
    }
  } else {
    pseudo_r2 <- NA
    r2_flag <- "N/A"
    cat("  Pseudo R²: N/A (not binomial model)\n")
  }
  
  # ============================================================================
  # CHECK 4: Convergence Status
  # ============================================================================
  converged <- model$converged
  cat(sprintf("  Converged: %s\n", converged))
  
  if (!converged) {
    cat("  🚨 Model did NOT converge!\n")
    conv_flag <- "FAIL"
  } else {
    cat("  ✓  Model converged\n")
    conv_flag <- "PASS"
  }
  
  # ============================================================================
  # CHECK 5: Number of Iterations (NEW!)
  # ============================================================================
  n_iter <- model$iter
  cat(sprintf("  Iterations: %d\n", n_iter))
  
  if (n_iter >= 20) {
    cat("  🚨 Many iterations (≥20) - convergence struggled!\n")
    iter_flag <- "FAIL"
  } else if (n_iter >= 15) {
    cat("  ⚠️  Moderate iterations (15-19) - borderline\n")
    iter_flag <- "WARNING"
  } else {
    cat("  ✓  Normal iteration count\n")
    iter_flag <- "PASS"
  }
  
  # ============================================================================
  # ADDITIONAL CHECK: Large non-significant coefficients
  # ============================================================================
  nonsig_huge <- coef_table[coef_table[, "Pr(>|z|)"] > 0.5 & 
                              abs(coef_table[, "Estimate"]) > 5, ]
  if (nrow(nonsig_huge) > 0) {
    cat("\n  ⚠️  WARNING: Large but non-significant coefficients detected:\n")
    cat(sprintf("     %s\n", paste(rownames(nonsig_huge), collapse=", ")))
    cat("     (May indicate numerical instability)\n")
  }
  
  # ============================================================================
  # OVERALL VERDICT
  # ============================================================================
  cat("\n")
  
  # Determine overall status
  has_failures <- any(c(coef_flag, se_flag, r2_flag, conv_flag, iter_flag) == "FAIL")
  has_warnings <- any(c(coef_flag, se_flag, r2_flag, conv_flag, iter_flag) == "WARNING")
  
  if (has_failures) {
    overall_status <- "❌ INVALID"
    cat("  ╔════════════════════════════════════════════════════════════╗\n")
    cat("  ║  ❌ VERDICT: MODEL HAS CRITICAL ISSUES                    ║\n")
    cat("  ║     → DO NOT USE for inference                            ║\n")
    cat("  ║     → Likely complete/quasi-complete separation           ��\n")
    cat("  ║     → Try: Firth's regression, simplify model, or exclude ║\n")
    cat("  ╚════════════════════════════════════════════════════════════╝\n")
  } else if (has_warnings) {
    overall_status <- "⚠️  BORDERLINE"
    cat("  ╔════════════════════════════════════════════════════════════╗\n")
    cat("  ║  ⚠️  VERDICT: MODEL ACCEPTABLE BUT USE WITH CAUTION       ║\n")
    cat("  ║     → Proceed but interpret results carefully             ║\n")
    cat("  ║     → Check DHARMa diagnostics thoroughly                 ║\n")
    cat("  ╚════════════════════════════════════════════════════════════╝\n")
  } else {
    overall_status <- "✅ VALID"
    cat("  ╔════════════════════════════════════════════════════════════╗\n")
    cat("  ║  ✅ VERDICT: MODEL HEALTHY                                ║\n")
    cat("  ║     → Safe to use for inference                           ║\n")
    cat("  ║     → Proceed to DHARMa diagnostics                       ║\n")
    cat("  ╚════════════════════════════════════════════════════════════╝\n")
  }
  
  cat("\n")
  
  # Store results in dataframe
  pre_dharma_results <- rbind(pre_dharma_results, data.frame(
    Model = lichen,
    Max_Coef = round(max_coef, 2),
    Problematic_Coef = ifelse(length(problem_coefs) > 0, 
                              paste(problem_coefs, collapse=", "), 
                              "None"),
    Max_SE = round(max_se, 2),
    Problematic_SE = ifelse(length(problem_ses) > 0, 
                            paste(problem_ses, collapse=", "), 
                            "None"),
    Pseudo_R2 = round(pseudo_r2, 3),
    Converged = ifelse(converged, "Yes", "No"),
    Iterations = n_iter,
    Coef_Flag = coef_flag,
    SE_Flag = se_flag,
    R2_Flag = r2_flag,
    Conv_Flag = conv_flag,
    Iter_Flag = iter_flag,
    Overall_Status = overall_status,
    stringsAsFactors = FALSE
  ))
}


# ==============================================================================
# SUMMARY TABLE 1: Quick Overview
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY TABLE: Pre-DHARMa Validation Results               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

print(pre_dharma_results[, c("Model", "Overall_Status", "Coef_Flag", 
                             "SE_Flag", "R2_Flag", "Conv_Flag", "Iter_Flag")])
cat("\n")


# ==============================================================================
# SUMMARY TABLE 2: Detailed Diagnostics
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DETAILED DIAGNOSTICS TABLE                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

print(pre_dharma_results[, c("Model", "Max_Coef", "Max_SE", "Pseudo_R2", 
                             "Converged", "Iterations")])
cat("\n")


# ==============================================================================
# SUMMARY TABLE 3: Problematic Parameters (if any)
# ==============================================================================

problematic_models <- pre_dharma_results[
  pre_dharma_results$Problematic_Coef != "None" | 
    pre_dharma_results$Problematic_SE != "None", 
]

if (nrow(problematic_models) > 0) {
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  ⚠️  MODELS WITH PROBLEMATIC PARAMETERS                     ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")
  
  print(problematic_models[, c("Model", "Problematic_Coef", "Problematic_SE")])
  cat("\n")
  
  cat("INTERPRETATION:\n")
  cat("• Coefficients with |β| > 10 indicate separation or numerical instability\n")
  cat("• Standard errors > 10 indicate extreme uncertainty (unreliable estimates)\n")
  cat("• These models should NOT be used for inference\n\n")
} else {
  cat("\n✅ No models with problematic parameters detected!\n\n")
}


# ==============================================================================
# SUMMARY TABLE 4: Flag Interpretation Guide
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FLAG INTERPRETATION GUIDE                                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

interpretation_table <- data.frame(
  Test = c("Coefficient Magnitude", 
           "Standard Error", 
           "Pseudo R²", 
           "Convergence", 
           "Iterations"),
  
  Pass_Threshold = c("|β| < 5",
                     "SE < 5",
                     "R² < 0.85",
                     "converged = TRUE",
                     "iterations < 15"),
  
  Warning_Threshold = c("5 ≤ |β| < 10",
                        "5 ≤ SE < 10",
                        "0.85 ≤ R² < 0.95",
                        "N/A",
                        "15 ≤ iter < 20"),
  
  Fail_Threshold = c("|β| ≥ 10",
                     "SE ≥ 10",
                     "R² ≥ 0.95",
                     "converged = FALSE",
                     "iterations ≥ 20"),
  
  Why_It_Matters = c(
    "Large coefficients suggest separation or implausible effects",
    "Large SEs indicate extreme uncertainty in estimates",
    "Perfect/near-perfect fit suggests model memorized data",
    "Non-convergence means algorithm couldn't find solution",
    "Many iterations suggest algorithm struggled (near-separation)"
  ),
  
  stringsAsFactors = FALSE
)

print(interpretation_table)
cat("\n")


# ==============================================================================
# FINAL SUMMARY: Count of Model Statuses
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("���  FINAL SUMMARY                                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

n_valid <- sum(grepl("VALID", pre_dharma_results$Overall_Status))
n_borderline <- sum(grepl("BORDERLINE", pre_dharma_results$Overall_Status))
n_invalid <- sum(grepl("INVALID", pre_dharma_results$Overall_Status))

cat(sprintf("Total models checked: %d\n\n", nrow(pre_dharma_results)))
cat(sprintf("✅ VALID models:        %d/%d (%.1f%%)\n", 
            n_valid, nrow(pre_dharma_results), 
            n_valid/nrow(pre_dharma_results)*100))
cat(sprintf("⚠️  BORDERLINE models:  %d/%d (%.1f%%)\n", 
            n_borderline, nrow(pre_dharma_results), 
            n_borderline/nrow(pre_dharma_results)*100))
cat(sprintf("❌ INVALID models:      %d/%d (%.1f%%)\n\n", 
            n_invalid, nrow(pre_dharma_results), 
            n_invalid/nrow(pre_dharma_results)*100))

if (n_invalid > 0) {
  invalid_models <- pre_dharma_results$Model[grepl("INVALID", pre_dharma_results$Overall_Status)]
  cat("⚠️  INVALID MODELS - DO NOT USE:\n")
  cat(sprintf("   • %s\n", paste(invalid_models, collapse="\n   • ")))
  cat("\nRECOMMENDATIONS:\n")
  cat("1. Check for complete separation using detectseparation package\n")
  cat("2. Try Firth's penalized logistic regression (logistf)\n")
  cat("3. Simplify model (remove predictors, especially decay4/decay5)\n")
  cat("4. Consider excluding from analysis if too problematic\n\n")
}

if (n_borderline > 0) {
  borderline_models <- pre_dharma_results$Model[grepl("BORDERLINE", pre_dharma_results$Overall_Status)]
  cat("⚠️  BORDERLINE MODELS - USE WITH CAUTION:\n")
  cat(sprintf("   • %s\n", paste(borderline_models, collapse="\n   • ")))
  cat("\nRECOMMENDATIONS:\n")
  cat("1. Proceed to DHARMa diagnostics\n")
  cat("2. Interpret results cautiously\n")
  cat("3. Report warnings in methods section\n\n")
}

if (n_valid == nrow(pre_dharma_results)) {
  cat("🎉 ALL MODELS PASSED PRE-DHARMA CHECKS!\n")
  cat("   → Proceed to DHARMa diagnostics with confidence\n\n")
}

# Save results
#write.csv(pre_dharma_results, "02_pre_dharma_validation_results.csv", row.names = FALSE)
#cat("✓ Results saved: 02_pre_dharma_validation_results.csv\n\n")





# ==============================================================================
# SECTION 6: DHARMA DIAGNOSTICS FOR ALL GLMs
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 4: MODEL DIAGNOSTICS WITH DHARMa                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Running simulation-based diagnostics for all models...\n")
cat("(This will take a few minutes)\n\n")


# Function to run comprehensive DHARMa diagnostics
run_dharma_diagnostics <- function(model, model_name, coords_df) {
  
  cat(sprintf("\n═══════════════════════════════════════════════════════════\n"))
  cat(sprintf("DHARMa Diagnostics: %s\n", model_name))
  cat(sprintf("═══════════════════════════════════════════════════════════\n"))
  
  # Simulate residuals
  cat("  Simulating residuals (n=1000)...\n")
  sim_resid <- simulateResiduals(fittedModel = model, n = 1000, plot = FALSE)
  cat("  ✓ Simulation complete\n\n")
  
  # Create diagnostic plot
  plot_filename <- paste0("02_DHARMa_", gsub(" ", "_", model_name), ".png")
  png(filename = plot_filename, width = 10, height = 6, units = "in", res = 300)
  plot(sim_resid, main = paste("DHARMa Diagnostics:", model_name))
  dev.off()
  cat(sprintf("  ✓ Diagnostic plot saved: %s\n\n", plot_filename))
  
  # Statistical tests
  cat("  Running diagnostic tests...\n\n")
  
  cat("  1. DISPERSION TEST:\n")
  disp_test <- testDispersion(sim_resid, plot = FALSE)
  print(disp_test)
  cat("\n")
  
  cat("  2. ZERO-INFLATION TEST:\n")
  zi_test <- testZeroInflation(sim_resid, plot = FALSE)
  print(zi_test)
  cat("\n")
  
  cat("  3. OUTLIER TEST:\n")
  outlier_test <- testOutliers(sim_resid, plot = FALSE)
  print(outlier_test)
  cat("\n")
  
  cat("  4. SPATIAL AUTOCORRELATION TEST (on residuals):\n")
  spatial_test <- testSpatialAutocorrelation(sim_resid, 
                                             x = coords_df$X, 
                                             y = coords_df$Y,
                                             plot = FALSE)
  print(spatial_test)
  cat("\n")
  
  # Compile results into data frame
  diagnostics <- data.frame(
    Model = model_name,
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
  
  # Summary interpretation
  all_ok <- all(diagnostics[, 6:9])
  if (all_ok) {
    cat("  ✅ ALL DIAGNOSTICS PASSED\n")
  } else {
    cat("  ⚠️  SOME DIAGNOSTICS FAILED - Review required\n")
  }
  
  cat("\n")
  
  return(list(diagnostics = diagnostics, sim_resid = sim_resid))
}


# Run diagnostics for all models
diagnostics_glm <- list()
diagnostics_summary <- data.frame()

for (lichen in all_lichen) {
  diag_result <- run_dharma_diagnostics(models_glm[[lichen]], 
                                        lichen, 
                                        data[, c("X", "Y")])
  diagnostics_glm[[lichen]] <- diag_result$sim_resid
  diagnostics_summary <- rbind(diagnostics_summary, diag_result$diagnostics)
}




# ==============================================================================
# SECTION 6: DHARMA DIAGNOSTICS FOR ALL GLMs (this section performs the same but adds a comparison with raw spatial pattern if available)
# ==============================================================================



cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STEP 4: MODEL DIAGNOSTICS WITH DHARMa                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Running simulation-based diagnostics for all models...\n")
cat("(This will take a few minutes)\n\n")


# Function to run comprehensive DHARMa diagnostics
run_dharma_diagnostics <- function(model, model_name, coords_df) {
  
  cat(sprintf("\n═══════════════════════════════════════════════════════════\n"))
  cat(sprintf("DHARMa Diagnostics: %s\n", model_name))
  cat(sprintf("═══════════════════════════════════════════════════════════\n"))
  
  # Simulate residuals
  cat("  Simulating residuals (n=1000)...\n")
  sim_resid <- simulateResiduals(fittedModel = model, n = 1000, plot = FALSE)
  cat("  ✓ Simulation complete\n\n")
  
  # Create diagnostic plot
  plot_filename <- paste0("02_DHARMa_", gsub(" ", "_", model_name), ".png")
  png(filename = plot_filename, width = 10, height = 6, units = "in", res = 300)
  plot(sim_resid, main = paste("DHARMa Diagnostics:", model_name))
  dev.off()
  cat(sprintf("  ✓ Diagnostic plot saved: %s\n\n", plot_filename))
  
  # Statistical tests
  cat("  Running diagnostic tests...\n\n")
  
  cat("  1. DISPERSION TEST:\n")
  disp_test <- testDispersion(sim_resid, plot = FALSE)
  print(disp_test)
  cat("\n")
  
  cat("  2. ZERO-INFLATION TEST:\n")
  zi_test <- testZeroInflation(sim_resid, plot = FALSE)
  print(zi_test)
  cat("\n")
  
  cat("  3. OUTLIER TEST:\n")
  outlier_test <- testOutliers(sim_resid, plot = FALSE)
  print(outlier_test)
  cat("\n")
  
  cat("  4. SPATIAL AUTOCORRELATION TEST (on residuals):\n")
  spatial_test <- testSpatialAutocorrelation(sim_resid, 
                                             x = coords_df$X, 
                                             y = coords_df$Y,
                                             plot = FALSE)
  print(spatial_test)
  cat("\n")
  
  # ============================================================================
  # TEST 5: Compare with raw spatial pattern (if available)
  # ============================================================================
  
  if (exists("spatial_results_raw", envir = .GlobalEnv)) {
    raw_spatial <- spatial_results_raw[spatial_results_raw$Variable == model_name, ]
    
    if (nrow(raw_spatial) > 0) {
      cat("  5. SPATIAL PATTERN COMPARISON (Before vs After):\n\n")
      
      cat(sprintf("     Raw data:        Moran's I = %7.4f, p = %.4f %s\n", 
                  raw_spatial$Morans_I, 
                  raw_spatial$P_value,
                  raw_spatial$Sig))
      
      cat(sprintf("     Residuals:       Moran's I = %7.4f, p = %.4f\n", 
                  as.numeric(spatial_test$statistic), 
                  spatial_test$p.value))
      
      cat("\n     ")
      if (raw_spatial$P_value < 0.05 && spatial_test$p.value >= 0.05) {
        cat("✅ Model fully explained spatial pattern!\n")
      } else if (raw_spatial$P_value < 0.05 && spatial_test$p.value < 0.05) {
        cat("⚠️  Spatial clustering remains in residuals\n")
      } else {
        cat("✅ No spatial issues\n")
      }
      cat("\n")
    }
  }
  
  # Compile results into data frame
  diagnostics <- data.frame(
    Model = model_name,
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
  
  # Summary interpretation
  all_ok <- all(diagnostics[, 6:9])
  if (all_ok) {
    cat("  ✅ ALL DIAGNOSTICS PASSED\n")
  } else {
    cat("  ⚠️  SOME DIAGNOSTICS FAILED - Review required\n")
  }
  
  cat("\n")
  
  return(list(diagnostics = diagnostics, sim_resid = sim_resid))
}


# Run diagnostics for all models
diagnostics_glm <- list()
diagnostics_summary <- data.frame()

for (lichen in all_lichen) {
  diag_result <- run_dharma_diagnostics(models_glm[[lichen]], 
                                        lichen, 
                                        data[, c("X", "Y")])
  diagnostics_glm[[lichen]] <- diag_result$sim_resid
  diagnostics_summary <- rbind(diagnostics_summary, diag_result$diagnostics)
}


# Print summary
cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DHARMA DIAGNOSTICS SUMMARY                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

print(diagnostics_summary)
cat("\n")

# Save results
write.csv(diagnostics_summary, "02_DHARMa_diagnostics_summary.csv", row.names = FALSE)
cat("✓ Results saved: 02_DHARMa_diagnostics_summary.csv\n\n")







# ==============================================================================
# SECTION 7: COMPREHENSIVE DIAGNOSTICS SUMMARY
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════���═══════╗\n")
cat("║  STEP 5: DIAGNOSTICS SUMMARY - ALL 7 MODELS                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

print(diagnostics_summary)
cat("\n")

# Save diagnostics table
write.csv(diagnostics_summary, 
          "02_GLM_diagnostics_summary.csv", 
          row.names = FALSE)

cat("✅ Diagnostics summary saved to: 02_GLM_diagnostics_summary.csv\n\n")


# Identify problems by category
cat("═══════════════════════════════════════════════════════════\n")
cat("DIAGNOSTIC ISSUES DETECTED:\n")
cat("═══════════════════════════════════════════════════════════\n\n")

problems_spatial <- diagnostics_summary[!diagnostics_summary$Spatial_OK, ]
problems_dispersion <- diagnostics_summary[!diagnostics_summary$Dispersion_OK, ]
problems_zeroinf <- diagnostics_summary[!diagnostics_summary$ZeroInflation_OK, ]
problems_outliers <- diagnostics_summary[!diagnostics_summary$Outliers_OK, ]

if (nrow(problems_spatial) > 0) {
  cat("⚠️  RESIDUAL SPATIAL AUTOCORRELATION:\n")
  for (i in 1:nrow(problems_spatial)) {
    cat(sprintf("    - %s (p = %.4f)\n", 
                problems_spatial$Model[i], 
                problems_spatial$SpatialAutocorr_p[i]))
  }
  cat("\n  → SOLUTION: Refit as GLMMs with (1|spatial_block) random effect\n\n")
} else {
  cat("✅ NO residual spatial autocorrelation detected\n\n")
}

if (nrow(problems_dispersion) > 0) {
  cat("⚠️  DISPERSION ISSUES:\n")
  for (i in 1:nrow(problems_dispersion)) {
    cat(sprintf("    - %s (p = %.4f)\n", 
                problems_dispersion$Model[i], 
                problems_dispersion$Dispersion_p[i]))
  }
  cat("\n  → SOLUTION: Use negative binomial family or quasi-models\n\n")
} else {
  cat("✅ NO dispersion issues detected\n\n")
}

if (nrow(problems_zeroinf) > 0) {
  cat("⚠️  ZERO-INFLATION:\n")
  for (i in 1:nrow(problems_zeroinf)) {
    cat(sprintf("    - %s (p = %.4f)\n", 
                problems_zeroinf$Model[i], 
                problems_zeroinf$ZeroInflation_p[i]))
  }
  cat("\n  → SOLUTION: Use zero-inflated models (glmmTMB with zi formula or glmm neg binom)\n\n") #. # turned out later that glmm neg binomial was way better.    
} else {
  cat("✅ NO zero-inflation issues detected\n\n")
}

if (nrow(problems_outliers) > 0) {
  cat("⚠️  OUTLIERS DETECTED:\n")
  for (i in 1:nrow(problems_outliers)) {
    cat(sprintf("    - %s (p = %.4f)\n", 
                problems_outliers$Model[i], 
                problems_outliers$Outliers_p[i]))
  }
  cat("\n  → Usually minor issue, but check influential points\n\n")
} else {
  cat("✅ NO outlier issues detected\n\n")
}


# Overall summary
n_problems <- sum(!diagnostics_summary$Spatial_OK) + 
  sum(!diagnostics_summary$Dispersion_OK) +
  sum(!diagnostics_summary$ZeroInflation_OK)

cat("═══════════════════════════════════════════════════════════\n")
if (n_problems == 0) {
  cat("🎉 ALL MODELS PASSED DIAGNOSTICS!\n")
  cat("   → GLMs are appropriate for your data\n")
  cat("   → Proceed to model selection and interpretation\n")
} else {
  cat(sprintf("⚠️  %d DIAGNOSTIC ISSUE(S) DETECTED\n", n_problems))
  cat("   → Review DHARMa plots (02_DHARMa_*.png files)\n")
  cat("   → Consider refitting problematic models as GLMMs\n")
}
cat("═══════════════════════════════════════════════════════════\n\n")






# ==============================================================================
# SECTION 8: SAVE WORKSPACE
# ==============================================================================

save(data, data_scaled,
     models_glm, diagnostics_glm, diagnostics_summary,
     spatial_results_raw, 
     lichen_binary, lichen_count, all_lichen,
     predictors, predictors_scaled,
     inv_dist_weights, 
     file = "Lichens/to_model/02_GLM_models_and_diagnostics.RData")

cat("\n✅ ANALYSIS COMPLETE!\n\n")
cat("Files saved:\n")
cat("  - 02_GLM_models_and_diagnostics.RData (workspace)\n")
cat("  - 02_spatial_autocorrelation_raw_responses.csv\n")
cat("  - 02_GLM_diagnostics_summary.csv\n")
cat("  - 02_spatial_pattern_core_ogf.png\n")
cat("  - 02_spatial_blocks_kmeans.png\n")
cat("  - 02_DHARMa_*.png (7 diagnostic plots)\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  NEXT STEPS:                                                 ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║  tbd                                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

################################################################################
# END OF SCRIPT 02
################################################################################