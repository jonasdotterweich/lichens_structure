# =============================================================================
# 01_data_loading.R
# Functions for loading raw data files for lichen-structure modeling
# =============================================================================
# Provides: load_forest_structure(), load_lichen_data(), load_coordinates(),
#           validate_data_structure()
# =============================================================================

library(readxl)
library(dplyr)
library(tibble)

source(here::here("R", "utils.R"))


# -----------------------------------------------------------------------------
# INTERNAL HELPERS
# -----------------------------------------------------------------------------

.clean_colnames <- function(df) {
  colnames(df) <- trimws(colnames(df))
  colnames(df) <- gsub("[\r\n]", "", colnames(df))
  df
}

.translate_lookup <- function(x, lookup) {
  result <- lookup[x]
  # preserve untranslated values where no mapping exists
  not_found <- is.na(result)
  result[not_found] <- x[not_found]
  unname(result)
}


# -----------------------------------------------------------------------------
# TRANSLATION TABLES (Czech → English)
# -----------------------------------------------------------------------------

.exposure_translations <- c(
  "JV - jihovychodni"              = "SE - southeast",
  "J - jizni"                      = "S - south",
  "V - vychodni"                   = "E - east",
  "S - severni"                    = "N - north",
  "SV - severovychodni"            = "NE - northeast",
  "nehodnoceno (rovina do 5 st.)"  = "not evaluated (flat land up to 5 deg.)",
  "SZ - severozapadni"             = "NW - northwest",
  "Z - zapadni"                    = "W - west",
  "JZ - jihozapadni"               = "SW - southwest"
)

.coverage_translations <- c(
  "nevyskytuje se"                  = "absent",
  "bez vyskytu"                     = "absent",
  "ojedinely vyskyt"                = "isolated occurrence",
  "velmi ridky vyskyt (do 0.2%)"    = "very sparse occurrence (up to 0.2%)",
  "ridky vyskyt (od 0.2 do 1%)"     = "sparse occurrence (from 0.2 to 1%)",
  "ridky vyskyt (0.2-1%)"           = "sparse occurrence (0.2-1%)",
  "ridky vyslyt (0.2-1%)"           = "sparse occurrence (0.2-1%)",
  "malocetny vyskyt (1-5%)"         = "infrequent occurrence (1-5%)",
  "malocetny vyskyt (od 1-5%)"      = "infrequent occurrence (from 1-5%)",
  "hojny vyskyt (6-25%)"            = "abundant occurrence (6-25%)",
  "velmi hojny vyskyt (25-50%)"     = "very common occurrence (25-50%)",
  "velmi hojny vyskyt (26-50%)"     = "very common occurrence (26-50%)",
  "velkoplosny vyskyt (51-75%)"     = "large-scale occurrence (51-75%)",
  "dominantni vyskyt (76-100%)"     = "dominant occurrence (76-100%)"
)

.management_translations <- c(
  "bez tezby, zadne parezy"                   = "no logging, no stumps",
  "tezba s ponechanim vetsiny hmoty"          = "logging with most biomass left in place",
  "tezba bez ponechani hmoty"                 = "logging without leaving biomass",
  "tezba s ponechanim mensiny hmoty"          = "logging with minority of biomass left in place"
)

.surface_translations <- c(
  "zivy les"                          = "living forest",
  "zadna z predchozich moznosti"      = "none of the previous options",
  "niva potoka/reky"                  = "stream/river floodplain",
  "velkoplosne vyvraty a polomy"      = "large-scale uprooting and windthrow",
  "kurovec"                           = "bark beetle",
  "porostni mezera"                   = "canopy gap",
  "redina po pastve"                  = "pasture woodland",
  "holina po tezbe"                   = "clear-cut area"
)


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Load and translate the forest structural data
#'
#' Reads `Lichens/Structural_data.xlsx`, cleans column names, and translates
#' all Czech categorical values to English.
#'
#' @param path Character. Path to the structural data xlsx file. Defaults to
#'   the project config value via `here::here()`.
#' @return A tibble with translated structural data (one row per plot).
#' @examples
#' structure_raw <- load_forest_structure()
load_forest_structure <- function(path = get_project_config()$data$structural) {
  if (!file.exists(path)) {
    stop("Structural data file not found: ", path,
         "\nCheck that the Lichens/ directory is present and the file exists.")
  }

  lichen_banner("Loading Forest Structural Data")
  df <- readxl::read_xlsx(path) |> .clean_colnames()

  # Translate categorical columns
  df <- df |>
    dplyr::mutate(
      exposure = .translate_lookup(exposure, .exposure_translations),
      `Past management (logging history)` = .translate_lookup(
        `Past management (logging history)`, .management_translations),
      `Plot surface appearance` = .translate_lookup(
        `Plot surface appearance`, .surface_translations)
    )

  # Translate coverage columns
  coverage_cols <- grep("coverage E[12] layer", colnames(df), value = TRUE)
  for (col in coverage_cols) {
    df[[col]] <- .translate_lookup(df[[col]], .coverage_translations)
  }

  lichen_message("Structural data loaded: ", nrow(df), " plots \u00d7 ",
                 ncol(df), " variables")
  tibble::as_tibble(df)
}


#' Load and prepare the lichen species occurrence data
#'
#' Reads `Lichens/Licen_data.xlsx` and translates Czech object ID values to
#' English.
#'
#' @param path Character. Path to the lichen data xlsx file. Defaults to the
#'   project config value via `here::here()`.
#' @return A tibble with raw lichen occurrence records.
#' @examples
#' lichens_raw <- load_lichen_data()
load_lichen_data <- function(path = get_project_config()$data$lichen) {
  if (!file.exists(path)) {
    stop("Lichen data file not found: ", path,
         "\nCheck that the Lichens/ directory is present and the file exists.")
  }

  lichen_banner("Loading Lichen Occurrence Data")
  df <- readxl::read_xlsx(path) |> .clean_colnames()

  object_id_translations <- c(
    "3 kl\u00e1da"              = "3 log",
    "8 kameny"                  = "8 stones",
    "4 pah\u00fdl"              = "4 stump",
    "5 kl\u00e1das"             = "5 logs",
    "7 kameny"                  = "7 stones",
    "6 kl\u00e1da"              = "6 log",
    "7 spadl\u00e1 v\u011btev"  = "7 fallen branch",
    "5 kl\u00e1da"              = "5 log",
    "6 n\u00edzk\u00fd pah\u00fdl" = "6 low stump",
    "3 pah\u00fdl"              = "3 stump",
    "4 kl\u00e1da"              = "4 log",
    "5 pa\u0159ez"              = "5 stump",
    "6 pah\u00fdl"              = "6 stump",
    "4 Picea abies mlad\u00fd"  = "4 Picea abies young",
    "Pinus sylvestris padl\u00e1 v\u011btev" = "Pinus sylvestris fallen branch",
    "6 v\u00fdvrat"             = "6 uprooted tree",
    "7 pa\u0159ez"              = "7 stump",
    "6 pa\u0159ez"              = "6 stump",
    "5 pah\u00fdl"              = "5 stump",
    "5 n\u00edzk\u00fd pah\u00fdl" = "5 low stump",
    "10 n\u00edzk\u00fd pah\u00fdl" = "10 low stump",
    "9 pah\u00fdl"              = "9 stump",
    "8 kl\u00e1da"              = "8 log",
    "7 v\u00fdvrat"             = "7 uprooted tree",
    "7 n\u00edzk\u00fd pah\u00fdl" = "7 low stump",
    "4 n\u00edzk\u00fd pah\u00fdl" = "4 low stump",
    "1 pah\u00fdl"              = "1 stump",
    "7 pah\u00fdl"              = "7 stump",
    "7 kl\u00e1da"              = "7 log",
    "2 kl\u00e1da"              = "2 log",
    "4 pa\u0159ez"              = "4 stump",
    "8 pa\u0159ez"              = "8 stump"
  )

  if ("object ID" %in% colnames(df)) {
    df <- df |>
      dplyr::mutate(`object ID` = .translate_lookup(`object ID`, object_id_translations))
  }

  lichen_message("Lichen data loaded: ", nrow(df), " observations")
  tibble::as_tibble(df)
}


#' Load plot coordinate data
#'
#' Reads `Lichens/Biodiversity_120 plot.xlsx`, cleans column names, and
#' returns a tibble with Project_ID, X (Easting), Y (Northing) in S-JTSK.
#'
#' @param path Character. Path to the coordinates xlsx file. Defaults to the
#'   project config value via `here::here()`.
#' @return A tibble with columns: Project_ID, X, Y, coordinates.
#' @examples
#' coords <- load_coordinates()
load_coordinates <- function(path = get_project_config()$data$coordinates) {
  if (!file.exists(path)) {
    stop("Coordinate file not found: ", path,
         "\nCheck that the Lichens/ directory is present and the file exists.")
  }

  lichen_banner("Loading Plot Coordinates")
  df <- readxl::read_xlsx(path) |> .clean_colnames()

  # Rename ID and coordinate columns (Czech header "\u010d." = "No.")
  id_col   <- intersect(c("\u010d.", "c.", "No.", "Project_ID"), colnames(df))[1]
  coord_col <- intersect(c("souradnice", "coordinates"), colnames(df))[1]

  if (is.na(id_col)) {
    stop("Cannot find plot ID column in coordinate file. ",
         "Expected one of: \u010d., c., No., Project_ID\n",
         "Columns present: ", paste(colnames(df), collapse = ", "))
  }

  df <- df |>
    dplyr::rename(Project_ID = dplyr::all_of(id_col)) |>
    dplyr::select(Project_ID, X, Y,
                  dplyr::any_of(c(coord_col, "coordinates", "souradnice")))

  lichen_message("Coordinates loaded: ", nrow(df), " plots")
  tibble::as_tibble(df)
}


#' Validate that raw data tables have the expected structure
#'
#' Checks required columns, ID overlap between lichen and structural data, and
#' reports any missing plots.
#'
#' @param structure A tibble returned by `load_forest_structure()`.
#' @param lichen A tibble returned by `load_lichen_data()`.
#' @param coords A tibble returned by `load_coordinates()` (optional).
#' @return A named list with elements: `valid` (logical), `n_plots_structure`,
#'   `n_plots_lichen`, `missing_in_structure` (character vector),
#'   `missing_in_lichen` (character vector), `messages` (character vector).
#' @examples
#' v <- validate_data_structure(structure_raw, lichens_raw, coords)
#' if (!v$valid) stop("Data validation failed")
validate_data_structure <- function(structure, lichen, coords = NULL) {
  lichen_banner("Validating Data Structure")

  messages <- character()
  valid <- TRUE

  # Required columns in structure
  required_structure <- c("Project_ID")
  missing_str_cols <- setdiff(required_structure, colnames(structure))
  if (length(missing_str_cols) > 0) {
    msg <- paste("Structure data missing required columns:",
                 paste(missing_str_cols, collapse = ", "))
    lichen_warning(msg)
    messages <- c(messages, msg)
    valid <- FALSE
  }

  # Required columns in lichen
  required_lichen <- c("Plot", "species")
  missing_lich_cols <- setdiff(required_lichen, colnames(lichen))
  if (length(missing_lich_cols) > 0) {
    msg <- paste("Lichen data missing required columns:",
                 paste(missing_lich_cols, collapse = ", "))
    lichen_warning(msg)
    messages <- c(messages, msg)
    valid <- FALSE
  }

  n_plots_structure <- if ("Project_ID" %in% colnames(structure)) {
    dplyr::n_distinct(structure$Project_ID)
  } else NA_integer_

  n_plots_lichen <- if ("Plot" %in% colnames(lichen)) {
    dplyr::n_distinct(lichen$Plot)
  } else NA_integer_

  # ID overlap
  missing_in_structure <- character()
  missing_in_lichen    <- character()
  if (!is.na(n_plots_structure) && !is.na(n_plots_lichen)) {
    missing_in_structure <- setdiff(unique(lichen$Plot),
                                    unique(structure$Project_ID))
    missing_in_lichen    <- setdiff(unique(structure$Project_ID),
                                    unique(lichen$Plot))

    if (length(missing_in_structure) > 0) {
      msg <- paste(length(missing_in_structure),
                   "lichen plot(s) not found in structural data")
      lichen_warning(msg)
      messages <- c(messages, msg)
    }
  }

  # Coordinate coverage
  if (!is.null(coords)) {
    missing_coords <- setdiff(unique(structure$Project_ID),
                              unique(coords$Project_ID))
    if (length(missing_coords) > 0) {
      msg <- paste(length(missing_coords), "plot(s) have no coordinates")
      lichen_warning(msg)
      messages <- c(messages, msg)
    }
  }

  lichen_message("Structure plots: ", n_plots_structure,
                 " | Lichen plots: ", n_plots_lichen)
  if (valid && length(messages) == 0) {
    lichen_message("All data validation checks passed")
  }

  list(
    valid                = valid,
    n_plots_structure    = n_plots_structure,
    n_plots_lichen       = n_plots_lichen,
    missing_in_structure = missing_in_structure,
    missing_in_lichen    = missing_in_lichen,
    messages             = messages
  )
}
