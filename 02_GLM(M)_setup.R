################################################################################
# ŠUMAVA NATIONAL PARK LICHEN-STRUCTURE MODELING
# Script 01: Setup, Data Loading, and GLMM Framework
# Author: Jonas Dotterweich
# Date: 2026-02-09
# Purpose: Prepare environment and conduct GLMM analysis for lichen indicators
################################################################################

# ==============================================================================
# SECTION 1: LIBRARY LOADING
# ==============================================================================

# Core data manipulation and visualization
library(tidyverse)      # dplyr, ggplot2, tidyr, etc.
library(here)           # For reproducible file paths

# Statistical modeling - GLMMs
library(lme4)           # Standard GLMM package (glmer function)
library(glmmTMB)        # Advanced GLMMs (better for zero-inflation, flexibility)
library(MASS)           # For negative binomial GLMs (glm.nb)

# Model diagnostics and validation
library(DHARMa)         # Residual diagnostics for GLMMs (ESSENTIAL!)
library(performance)    # Model checking (check_collinearity, check_overdispersion)
library(see)            # Visualization for performance package

# Model selection and comparison
library(MuMIn)          # Multi-model inference, AICc calculations
library(bbmle)          # AIC tables and model comparison
library(car)            # VIF rechecking, Anova() for Type II/III tests

# Spatial analysis (for autocorrelation testing)
library(spdep)          # Moran's I and spatial weights
library(sp)             # Spatial data structures
library(gstat)          # Variogram analysis

# Threshold detection (segmented regression)
library(segmented)      # Breakpoint detection in GLMs
library(chngpt)         # Change-point models

# Visualization enhancements
library(sjPlot)         # Publication-ready model tables and plots
library(effects)        # Effect plots for GLMMs
library(ggeffects)      # Marginal effects plotting

# Reporting
library(broom)          # Tidy model outputs
library(broom.mixed)    # Tidy outputs for mixed models


# ==============================================================================
# SECTION 2: DATA LOADING
# ==============================================================================

# Set working directory (modify to your project folder)
# Alternatively, use an RStudio project and here() package
# setwd("~/path/to/your/project")

# Load the three VIF-cleaned structural datasets
data_16 <- read.csv("Lichens/to_model/modeling_dataset_FULL_16predictors.csv", 
                    stringsAsFactors = FALSE)

data_14 <- read.csv("Lichens/to_model/modeling_dataset_MODERATE_14predictors.csv", 
                    stringsAsFactors = FALSE)

data_12 <- read.csv("Lichens/to_model/modeling_dataset_CONSERVATIVE_12predictors.csv", 
                    stringsAsFactors = FALSE)

# Set primary dataset for analysis (CONSERVATIVE recommended)
data <- data_16  # Modify to data_14 or data_16 if needed

# checking and showing the dataset
names(data)


#I will delete the colum tree hight median

data <- data %>% select(-tree_height_median, -regeneration)

# ==============================================================================
# SECTION 3: DATA INSPECTION
# ==============================================================================

# Check dimensions
cat("\n=== DATASET DIMENSIONS ===\n")
cat("Full (16 predictors):", dim(data_16), "\n")
cat("Moderate (14 predictors):", dim(data_14), "\n")
cat("Conservative (12 predictors):", dim(data_12), "\n")
cat("\nWorking with:", dim(data), "observations\n")

# Verify structure
str(data)

# Check for missing values (should be 0)
cat("\n=== MISSING VALUE CHECK ===\n")
cat("Total missing values:", sum(is.na(data)), "\n")

# Display first few rows
head(data)


# ==============================================================================
# SECTION 4: DEFINE RESPONSE AND PREDICTOR VARIABLES
# ==============================================================================

# Lichen response variables (7 groups)
lichen_binary <- c(
  "core_ogf_presence",
  #"mycoblastus_presence", 
  "ochrolechia_presence",
  "xylographa_presence",
  "parmelia_agg_presence",
  "elite_rare_presence"
)

lichen_count <- c(
  "calicioids_richness"
)

all_lichen <- c(lichen_binary, lichen_count)

# Structural predictors (12 in CONSERVATIVE dataset)
predictors_12 <- c(
  "elevation",
  "volume_snags",
  #"tree_height_median",
  "canopy_cover",
  #"regeneration",
  "decay2",           # Key predictor from univariate screening
  "decay3",           # Key predictor from univariate screening
  "dbh_max",
  "dbh_sd",
  "n_dead_50cm",
  "ba_spruce",
  "ba_beech"
)

# Structural predictors (14 in MODERATE dataset)
predictors_14 <- c(
  predictors_12,
  "n_living_80cm",
  "ba_late_succ"
)

# Structural predictors (16 in FULL dataset) 
predictors_16 <- c(
  predictors_12,
  "decay4",
  "decay5"
)

# Set active predictor set based on chosen dataset
predictors <- predictors_16  # Modify if using data_14 or data_16


predictors


# ==============================================================================
# SECTION 5: VERIFY LICHEN RESPONSE DISTRIBUTIONS
# ==============================================================================

cat("\n=== LICHEN RESPONSE VARIABLE SUMMARIES ===\n\n")

# Binary responses (presence/absence)
cat("BINARY RESPONSES (prevalence):\n")
for (lichen in lichen_binary) {
  prevalence <- mean(data[[lichen]], na.rm = TRUE)
  n_present <- sum(data[[lichen]], na.rm = TRUE)
  cat(sprintf("  %s: %.1f%% present (%d/%d plots)\n", 
              lichen, prevalence * 100, n_present, nrow(data)))
}

# Count response (richness)
cat("\nCOUNT RESPONSE (richness):\n")
cat("  calicioids_richness:\n")
summary(data$calicioids_richness)
cat("  Zero plots:", sum(data$calicioids_richness == 0), "\n")
cat("  Max richness:", max(data$calicioids_richness), "\n")


# ==============================================================================
# SECTION 6: VERIFY SPATIAL COORDINATES ARE AVAILABLE
# ==============================================================================

# Check if X, Y coordinates exist
if (all(c("X", "Y") %in% names(data))) {
  cat("\n=== SPATIAL COORDINATES DETECTED ===\n")
  cat("Coordinate ranges:\n")
  cat("  X (Easting):", range(data$X), "\n")
  cat("  Y (Northing):", range(data$Y), "\n")
  
  # Quick spatial plot
  ggplot(data, aes(x = X, y = Y)) +
    geom_point(size = 2, alpha = 0.6) +
    coord_equal() +
    theme_minimal() +
    labs(title = "Spatial Distribution of 120 Forest Plots",
         subtitle = "Šumava National Park")
  
} else {
  warning("No X, Y coordinates found - spatial modeling may be limited")
}


# ==============================================================================
# SECTION 7: RECHECK COLLINEARITY (VIF)
# ==============================================================================

# Double-check VIF with your chosen predictor set
cat("\n=== VIF RECALCULATION FOR SELECTED PREDICTORS ===\n")

# Create a simple linear model to calculate VIF
# (Using arbitrary response - VIF only depends on predictors)
vif_check <- lm(as.formula(paste("elevation ~", 
                                 paste(predictors[-1], collapse = " + "))), 
                data = data)

vif_values <- car::vif(vif_check)
cat("\nVIF values:\n")
print(round(vif_values, 2))
cat("\nMean VIF:", round(mean(vif_values), 2), "\n")
cat("Max VIF:", round(max(vif_values), 2), "\n")

if (max(vif_values) > 10) {
  warning("⚠️  Some predictors have VIF > 10! Reconsider model structure.")
} else if (max(vif_values) > 5) {
  message("⚠️  Some predictors have VIF between 5-10 (moderate collinearity)")
} else {
  message("✅ All VIF values < 5 (acceptable collinearity)")
}


# ==============================================================================
# SECTION 8: STANDARDIZE PREDICTORS (RECOMMENDED FOR GLMMs)
# ==============================================================================

# Create standardized versions (mean = 0, SD = 1)
# This helps with:
# - Model convergence
# - Coefficient interpretation (effect sizes)
# - Comparing relative importance

data_scaled <- data

for (pred in predictors) {
  scaled_name <- paste0(pred, "_scaled")
  data_scaled[[scaled_name]] <- scale(data[[pred]])[, 1]
}

# Create list of scaled predictor names
predictors_scaled <- paste0(predictors, "_scaled")

cat("\n=== PREDICTOR STANDARDIZATION COMPLETE ===\n")
cat("Original predictors:", paste(predictors, collapse = ", "), "\n")
cat("Scaled predictors:", paste(predictors_scaled, collapse = ", "), "\n")


# ==============================================================================
# SECTION 9: SAVE WORKSPACE FOR MODELING
# ==============================================================================

# Save prepared data for downstream analyses
save(data, data_scaled, 
     lichen_binary, lichen_count, all_lichen,
     predictors, predictors_scaled,
     file = "Lichens/to_model/01_prepared_data_for_modeling(12pred_wo_mycob).RData")  ### Modify this if only twelve predictors are used

cat("\n✅ Setup complete! Ready for GLMM modeling.\n")
cat("Workspace saved to: 01_prepared_data_for_modeling.RData\n")


################################################################################
# END OF SETUP SCRIPT
# Next: Spatial autocorrelation testing and GLMM fitting
################################################################################
