# =============================================================================
# utils.R
# Helper functions for the Šumava lichen-structure modeling project
# =============================================================================
# Provides: path helpers (here::here()), configuration loading,
#           message formatting, and generic file I/O utilities.
# =============================================================================

library(here)
library(tibble)
library(dplyr)
library(readr)

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

#' Return the default project configuration as a named list
#'
#' @return A named list with file paths, modeling thresholds, and output
#'   directory settings.
#' @examples
#' cfg <- get_project_config()
#' cfg$data$structural
get_project_config <- function() {
  list(
    data = list(
      structural   = here::here("Lichens", "Structural_data.xlsx"),
      lichen       = here::here("Lichens", "Licen_data.xlsx"),
      coordinates  = here::here("Lichens", "Biodiversity_120 plot.xlsx"),
      n_plots      = 120L
    ),
    output = list(
      root     = here::here("Lichens", "model_outputs"),
      standard = here::here("Lichens", "model_outputs", "standard"),
      reduced  = here::here("Lichens", "model_outputs", "reduced")
    ),
    modeling = list(
      prevalence_min    = 20,    # minimum species prevalence (%) for modeling
      prevalence_max    = 80,    # maximum species prevalence (%) for modeling
      missing_threshold = 2L,    # maximum missing key predictors before removing a plot
      outlier_z         = 3,     # |Z-score| threshold to flag extreme outliers
      vif_threshold     = 10,    # VIF threshold for collinearity
      cor_threshold     = 0.7,   # |r| threshold for high correlation
      n_sim_dharma      = 1000L  # DHARMa simulation iterations
    ),
    lichen_groups = list(
      binary = c(
        "parmelia_agg_presence",
        "ochrolechia_presence",
        "core_ogf_presence",
        "mycoblastus_presence",
        "xylographa_presence",
        "elite_rare_presence"
      ),
      count  = "calicioids_richness"
    )
  )
}


# -----------------------------------------------------------------------------
# PATH HELPERS
# -----------------------------------------------------------------------------

#' Build a project-relative output path using here::here()
#'
#' @param ... Path components passed to here::here().
#' @return A character string with the full path.
#' @examples
#' project_path("Lichens", "model_outputs", "results.csv")
project_path <- function(...) {
  here::here(...)
}


#' Ensure an output directory exists, creating it recursively if needed
#'
#' @param path Character. Directory path to create.
#' @return The path invisibly.
#' @examples
#' ensure_dir(project_path("Lichens", "model_outputs", "reduced"))
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    lichen_message("Created directory: ", path)
  }
  invisible(path)
}


# -----------------------------------------------------------------------------
# MESSAGING / FORMATTING
# -----------------------------------------------------------------------------

#' Print a formatted project message with a check-mark prefix
#'
#' @param ... Character strings passed to cat().
#' @return NULL invisibly.
#' @examples
#' lichen_message("Data loaded: 120 plots")
lichen_message <- function(...) {
  cat("\u2713 ", ..., "\n", sep = "")
  invisible(NULL)
}


#' Print a formatted warning with a flag prefix
#'
#' @param ... Character strings passed to cat().
#' @return NULL invisibly.
#' @examples
#' lichen_warning("3 plots removed due to missing data")
lichen_warning <- function(...) {
  cat("\u26A0\uFE0F  WARNING: ", ..., "\n", sep = "")
  invisible(NULL)
}


#' Print a section header banner
#'
#' @param title Character. Section title.
#' @return NULL invisibly.
#' @examples
#' lichen_banner("STEP 1: Data Loading")
lichen_banner <- function(title) {
  line <- paste(rep("=", 60), collapse = "")
  cat("\n", line, "\n", title, "\n", line, "\n\n", sep = "")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# FILE I/O
# -----------------------------------------------------------------------------

#' Save a data frame to CSV with a standard message
#'
#' @param data A data frame or tibble.
#' @param path Character. Full file path (use project_path()).
#' @param ... Additional arguments passed to readr::write_csv().
#' @return The data invisibly.
#' @examples
#' save_csv(results, project_path("Lichens", "model_outputs", "results.csv"))
save_csv <- function(data, path, ...) {
  ensure_dir(dirname(path))
  readr::write_csv(data, path, ...)
  lichen_message("Saved CSV: ", basename(path))
  invisible(data)
}


#' Save an R object to an .RData file with a standard message
#'
#' @param ... Objects to save (passed by name to save()).
#' @param path Character. Full .RData file path.
#' @return NULL invisibly.
#' @examples
#' save_rdata(models_glm, diagnostics_summary, path = project_path("Lichens", "models.RData"))
save_rdata <- function(..., path) {
  ensure_dir(dirname(path))
  save(..., file = path)
  lichen_message("Saved RData: ", basename(path))
  invisible(NULL)
}


#' Save a ggplot to PNG with sensible defaults
#'
#' @param plot A ggplot object.
#' @param path Character. Full PNG file path.
#' @param width Numeric. Width in inches (default 10).
#' @param height Numeric. Height in inches (default 7).
#' @param dpi Numeric. Resolution (default 300).
#' @return NULL invisibly.
#' @examples
#' save_plot(my_plot, project_path("Lichens", "model_outputs", "figure.png"))
save_plot <- function(plot, path, width = 10, height = 7, dpi = 300) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(filename = path, plot = plot,
                  width = width, height = height, dpi = dpi)
  lichen_message("Saved plot: ", basename(path))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# MISC HELPERS
# -----------------------------------------------------------------------------

#' Significance star string for a p-value
#'
#' @param p Numeric p-value.
#' @return Character: "***", "**", "*", "." or "ns".
#' @examples
#' sig_stars(0.001)
#' sig_stars(0.06)
sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE       ~ "ns"
  )
}
