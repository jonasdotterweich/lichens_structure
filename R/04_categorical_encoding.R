# =============================================================================
# 04_categorical_encoding.R
# Functions for encoding categorical predictor variables
# =============================================================================
# Provides: encode_management(), encode_exposure(),
#           encode_dominant_species(), encode_all_categoricals()
# =============================================================================

library(dplyr)
library(stringr)
library(tibble)

source(here::here("R", "utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Encode the past-management variable
#'
#' Creates three derived columns from the `management` column:
#' * `management_factor` – ordered factor (none < most biomass left < minority
#'   left < no biomass left)
#' * `logged` – binary: 0 = no logging, 1 = any logging
#' * `logging_intensity` – integer 0–3 (0 = unlogged, 3 = complete removal)
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param col Character. Name of the management column.
#'   Default "management".
#' @return The input tibble with three new columns appended.
#' @examples
#' structure_enc <- structure_clean |> encode_management()
encode_management <- function(data, col = "management") {
  stopifnot(is.data.frame(data))

  if (!col %in% colnames(data)) {
    stop("Column '", col, "' not found in data. ",
         "Available columns: ", paste(colnames(data), collapse = ", "))
  }

  management_levels <- c(
    "no logging, no stumps",
    "logging with most biomass left in place",
    "logging with minority of biomass left in place",
    "logging without leaving biomass"
  )

  data <- data |>
    dplyr::mutate(
      management_factor = factor(
        !!rlang::sym(col),
        levels  = management_levels,
        ordered = TRUE
      ),
      logged = dplyr::if_else(
        !!rlang::sym(col) == "no logging, no stumps", 0L, 1L
      ),
      logging_intensity = as.integer(management_factor) - 1L
    )

  # Report category counts
  cat_counts <- table(data$management_factor, useNA = "ifany")
  lichen_message("Management encoded. Categories:")
  print(cat_counts)

  tibble::as_tibble(data)
}


#' Encode the aspect/exposure variable
#'
#' Simplifies the eight cardinal and intercardinal directions into five groups
#' (N, S, E, W, NE_SE, NW_SW, flat) and creates:
#' * `exposure_simple` – simplified character label
#' * `exposure_factor` – unordered factor
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param col Character. Name of the exposure column. Default "exposure".
#' @return The input tibble with two new columns appended.
#' @examples
#' structure_enc <- structure_enc |> encode_exposure()
encode_exposure <- function(data, col = "exposure") {
  stopifnot(is.data.frame(data))

  if (!col %in% colnames(data)) {
    stop("Column '", col, "' not found in data.")
  }

  data <- data |>
    dplyr::mutate(
      exposure_simple = dplyr::case_when(
        stringr::str_detect(!!rlang::sym(col), "not evaluated|flat") ~ "flat",
        stringr::str_detect(!!rlang::sym(col), "^N ")                ~ "N",
        stringr::str_detect(!!rlang::sym(col), "^S ")                ~ "S",
        stringr::str_detect(!!rlang::sym(col), "^E ")                ~ "E",
        stringr::str_detect(!!rlang::sym(col), "^W ")                ~ "W",
        stringr::str_detect(!!rlang::sym(col), "NE|SE")              ~ "NE_SE",
        stringr::str_detect(!!rlang::sym(col), "NW|SW")              ~ "NW_SW",
        TRUE                                                           ~ "other"
      ),
      exposure_factor = factor(exposure_simple)
    )

  cat_counts <- table(data$exposure_simple, useNA = "ifany")
  lichen_message("Exposure encoded. Categories:")
  print(cat_counts)

  tibble::as_tibble(data)
}


#' Encode the dominant tree species variable
#'
#' Groups rare species codes into "other" and creates:
#' * `dominant_grouped` – character with values SM, BK, JD, or "other"
#' * `dominant_factor` – factor with levels c("SM", "BK", "JD", "other")
#'
#' Species codes follow the Czech/Slovak dendrological convention:
#' SM = Norway spruce (*Picea abies*), BK = European beech (*Fagus sylvatica*),
#' JD = Silver fir (*Abies alba*).
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param col Character. Name of the dominant species column.
#'   Default "dominant_species".
#' @param main_species Character vector of species codes to retain as individual
#'   levels. Default `c("SM", "BK", "JD")`.
#' @return The input tibble with two new columns appended.
#' @examples
#' structure_enc <- structure_enc |> encode_dominant_species()
encode_dominant_species <- function(data,
                                    col          = "dominant_species",
                                    main_species = c("SM", "BK", "JD")) {
  stopifnot(is.data.frame(data))

  if (!col %in% colnames(data)) {
    stop("Column '", col, "' not found in data.")
  }

  factor_levels <- c(main_species, "other")

  data <- data |>
    dplyr::mutate(
      dominant_grouped = dplyr::if_else(
        !!rlang::sym(col) %in% main_species,
        !!rlang::sym(col),
        "other"
      ),
      dominant_factor = factor(dominant_grouped, levels = factor_levels)
    )

  cat_counts <- table(data$dominant_grouped, useNA = "ifany")
  lichen_message("Dominant species encoded. Categories:")
  print(cat_counts)

  tibble::as_tibble(data)
}


#' Apply all categorical encoding steps in sequence
#'
#' Convenience wrapper that calls `encode_management()`, `encode_exposure()`,
#' and `encode_dominant_species()` with their default column names.
#'
#' @param data A data frame or tibble (pipe-friendly).
#' @param management_col Character. Column name for management. Default "management".
#' @param exposure_col Character. Column name for exposure. Default "exposure".
#' @param dominant_col Character. Column name for dominant species.
#'   Default "dominant_species".
#' @return The input tibble with all categorical encoding columns appended.
#' @examples
#' structure_enc <- structure_clean |> encode_all_categoricals()
encode_all_categoricals <- function(data,
                                    management_col = "management",
                                    exposure_col   = "exposure",
                                    dominant_col   = "dominant_species") {
  stopifnot(is.data.frame(data))

  lichen_banner("Encoding Categorical Variables")

  if (management_col %in% colnames(data)) {
    data <- encode_management(data, col = management_col)
  } else {
    lichen_warning("Column '", management_col, "' not found – skipping management encoding")
  }

  if (exposure_col %in% colnames(data)) {
    data <- encode_exposure(data, col = exposure_col)
  } else {
    lichen_warning("Column '", exposure_col, "' not found – skipping exposure encoding")
  }

  if (dominant_col %in% colnames(data)) {
    data <- encode_dominant_species(data, col = dominant_col)
  } else {
    lichen_warning("Column '", dominant_col, "' not found – skipping dominant species encoding")
  }

  lichen_message("All categorical encoding complete")
  tibble::as_tibble(data)
}
