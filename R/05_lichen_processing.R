# =============================================================================
# 05_lichen_processing.R
# Functions for preparing lichen response variables
# =============================================================================
# Provides: check_species_prevalence(), create_species_groups(),
#           validate_response_variables()
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

source(here::here("R", "utils.R"))


# -----------------------------------------------------------------------------
# PUBLIC FUNCTIONS
# -----------------------------------------------------------------------------

#' Check prevalence and modeling suitability of lichen species
#'
#' Computes per-species occurrence statistics across all plots and categorises
#' each species by its suitability for binomial GLM modeling.
#'
#' @param data A data frame or tibble of raw lichen records (pipe-friendly).
#'   Expected columns: `Plot` (plot identifier) and `species`.
#' @param n_plots Integer. Total number of surveyed plots (denominator for
#'   prevalence). Default 120.
#' @param min_prevalence Numeric. Lower bound (%) for "good" modeling range.
#'   Default 20.
#' @param max_prevalence Numeric. Upper bound (%) for "good" modeling range.
#'   Default 80.
#' @return A tibble with one row per species: species, n_plots, prevalence_pct,
#'   mean_cover, median_cover, max_cover, total_observations, model_suitability.
#' @examples
#' prevalence <- lichens_raw |> check_species_prevalence()
#' prevalence |> dplyr::filter(model_suitability == "good_balanced")
check_species_prevalence <- function(data,
                                     n_plots        = get_project_config()$data$n_plots,
                                     min_prevalence = get_project_config()$modeling$prevalence_min,
                                     max_prevalence = get_project_config()$modeling$prevalence_max) {
  stopifnot(is.data.frame(data))

  required <- c("Plot", "species")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop("check_species_prevalence() requires columns: ",
         paste(missing_cols, collapse = ", "))
  }

  cover_col <- intersect(c("cover", "Cover", "abundance"), colnames(data))[1]
  has_cover <- !is.na(cover_col)

  .cover_stat <- function(fn) {
    if (has_cover) fn(.data[[cover_col]], na.rm = TRUE) else NA_real_
  }

  summary <- data |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      n_plots            = dplyr::n_distinct(Plot),
      total_observations = dplyr::n(),
      mean_cover         = .cover_stat(mean),
      median_cover       = .cover_stat(median),
      max_cover          = .cover_stat(max),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      prevalence_pct    = round(n_plots / !!n_plots * 100, 1),
      model_suitability = dplyr::case_when(
        prevalence_pct >= min_prevalence & prevalence_pct <= max_prevalence
          ~ "good_balanced",
        prevalence_pct > max_prevalence
          ~ "ubiquitous_low_variance",
        prevalence_pct >= 10 & prevalence_pct < min_prevalence
          ~ "rare_consider_grouping",
        TRUE
          ~ "very_rare_not_suitable"
      )
    ) |>
    dplyr::arrange(dplyr::desc(n_plots))

  n_good <- sum(summary$model_suitability == "good_balanced")
  lichen_message("Species assessed: ", nrow(summary),
                 " | Suitable for modeling (", min_prevalence, "-", max_prevalence,
                 "%): ", n_good)
  summary
}


#' Create lichen species groups / response variables for modeling
#'
#' Aggregates raw lichen records into plot-level presence/absence and count
#' response variables. The groups returned match those used in the Šumava
#' lichen-structure analysis:
#'
#' * `calicioids_richness` – species richness of calicioid (pin) lichens
#' * `parmelia_agg_presence` – presence of *Parmelia* aggregate
#' * `ochrolechia_presence` – presence of *Ochrolechia* spp.
#' * `core_ogf_presence` – presence of core old-growth-forest assemblage
#' * `mycoblastus_presence` – presence of *Mycoblastus* spp.
#' * `xylographa_presence` – presence of *Xylographa vitiligo*
#' * `elite_rare_presence` – presence of elite/rare lichen assemblage
#'
#' @param data A data frame or tibble of raw lichen records (pipe-friendly).
#'   Expected columns: `Plot`, `species`, and optionally `cover`.
#' @param calicioid_species Character vector of calicioid species names to sum
#'   for richness. Uses a curated default list when `NULL`.
#' @param parmelia_species Character vector of *Parmelia* agg. species names.
#' @param ochrolechia_species Character vector of *Ochrolechia* species names.
#' @param core_ogf_species Character vector of core OGF indicator species.
#' @param mycoblastus_species Character vector of *Mycoblastus* species names.
#' @param xylographa_species Character vector of *Xylographa* species names.
#' @param elite_species Character vector of elite/rare indicator species.
#' @return A tibble with one row per plot and all response-variable columns.
#' @examples
#' lichen_groups <- lichens_raw |> create_species_groups()
create_species_groups <- function(
    data,
    calicioid_species  = NULL,
    parmelia_species   = NULL,
    ochrolechia_species = NULL,
    core_ogf_species   = NULL,
    mycoblastus_species = NULL,
    xylographa_species = NULL,
    elite_species      = NULL) {

  stopifnot(is.data.frame(data))

  required <- c("Plot", "species")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop("create_species_groups() requires columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Default species lists (based on Šumava calicioid assemblage)
  if (is.null(calicioid_species)) {
    calicioid_species <- c(
      "Calicium viride", "Chaenotheca brunneola", "Chaenotheca chrysocephala",
      "Chaenotheca ferruginea", "Chaenotheca furfuracea", "Chaenotheca hispidula",
      "Chaenotheca phaeocephala", "Chaenotheca sphaerocephala", "Chaenotheca stemonea",
      "Chaenotheca subroscida", "Chaenotheca trichialis", "Chaenothecopsis pusilla",
      "Mycocalicium subtile", "Microcalicium disseminatum"
    )
  }
  if (is.null(parmelia_species)) {
    parmelia_species <- c("Parmelia sulcata", "Parmelia saxatilis",
                          "Parmelina tiliacea", "Parmotrema perlatum")
  }
  if (is.null(ochrolechia_species)) {
    ochrolechia_species <- c("Ochrolechia arborea", "Ochrolechia pallescens",
                             "Ochrolechia turneri")
  }
  if (is.null(core_ogf_species)) {
    core_ogf_species <- c("Lobaria pulmonaria", "Lobaria scrobiculata",
                          "Sticta limbata", "Menegazzia terebrata",
                          "Hypogymnia vittata", "Alectoria sarmentosa",
                          "Bryoria nadvornikiana")
  }
  if (is.null(mycoblastus_species)) {
    mycoblastus_species <- c("Mycoblastus sanguinarius", "Mycoblastus fucatus")
  }
  if (is.null(xylographa_species)) {
    xylographa_species <- c("Xylographa vitiligo", "Xylographa parallela")
  }
  if (is.null(elite_species)) {
    elite_species <- c(core_ogf_species, "Usnea longissima", "Evernia divaricata",
                       "Cetrelia cetrarioides", "Platismatia norvegica")
  }

  .presence <- function(sp_list) {
    data |>
      dplyr::filter(species %in% sp_list) |>
      dplyr::distinct(Plot) |>
      dplyr::mutate(.present = 1L)
  }

  all_plots <- tibble::tibble(Plot = unique(data$Plot))

  # Calicioid richness (count)
  calicioid_richness <- data |>
    dplyr::filter(species %in% calicioid_species) |>
    dplyr::group_by(Plot) |>
    dplyr::summarise(calicioids_richness = dplyr::n_distinct(species),
                     .groups = "drop")

  # Binary presence variables
  make_presence_col <- function(sp_list, col_name) {
    all_plots |>
      dplyr::left_join(.presence(sp_list), by = "Plot") |>
      dplyr::mutate(!!col_name := tidyr::replace_na(.present, 0L)) |>
      dplyr::select(Plot, dplyr::all_of(col_name))
  }

  parmelia_pres    <- make_presence_col(parmelia_species,    "parmelia_agg_presence")
  ochrolechia_pres <- make_presence_col(ochrolechia_species, "ochrolechia_presence")
  core_ogf_pres    <- make_presence_col(core_ogf_species,    "core_ogf_presence")
  mycoblastus_pres <- make_presence_col(mycoblastus_species, "mycoblastus_presence")
  xylographa_pres  <- make_presence_col(xylographa_species,  "xylographa_presence")
  elite_pres       <- make_presence_col(elite_species,       "elite_rare_presence")

  groups <- all_plots |>
    dplyr::left_join(calicioid_richness, by = "Plot") |>
    dplyr::mutate(calicioids_richness = tidyr::replace_na(calicioids_richness, 0L)) |>
    dplyr::left_join(parmelia_pres,    by = "Plot") |>
    dplyr::left_join(ochrolechia_pres, by = "Plot") |>
    dplyr::left_join(core_ogf_pres,    by = "Plot") |>
    dplyr::left_join(mycoblastus_pres, by = "Plot") |>
    dplyr::left_join(xylographa_pres,  by = "Plot") |>
    dplyr::left_join(elite_pres,       by = "Plot")

  lichen_message("Species groups created for ", nrow(groups), " plots")
  tibble::as_tibble(groups)
}


#' Validate response variables in the modeling dataset
#'
#' Checks that all expected response columns are present, reports
#' prevalence/range, and flags potential issues (e.g. zero-inflation, extreme
#' imbalance).
#'
#' @param data A data frame or tibble that includes response variable columns
#'   (pipe-friendly).
#' @param binary_responses Character vector of binary (0/1) response column
#'   names. Defaults to the six binomial response variables.
#' @param count_responses Character vector of count response column names.
#'   Default `"calicioids_richness"`.
#' @return A tibble with one row per response variable: variable, type,
#'   n_present (or mean for counts), prevalence_pct (binary only), min, max,
#'   status.
#' @examples
#' response_check <- modeling_data |> validate_response_variables()
validate_response_variables <- function(
    data,
    binary_responses = get_project_config()$lichen_groups$binary,
    count_responses  = get_project_config()$lichen_groups$count) {

  stopifnot(is.data.frame(data))

  lichen_banner("Validating Response Variables")

  all_responses <- c(binary_responses, count_responses)
  missing_vars  <- setdiff(all_responses, colnames(data))

  if (length(missing_vars) > 0) {
    lichen_warning("Missing response columns: ",
                   paste(missing_vars, collapse = ", "))
  }

  present_vars <- intersect(all_responses, colnames(data))

  report <- purrr::map_dfr(present_vars, function(var) {
    x    <- data[[var]]
    type <- if (var %in% binary_responses) "binary" else "count"

    if (type == "binary") {
      n_present      <- sum(x == 1, na.rm = TRUE)
      prevalence_pct <- round(n_present / length(x) * 100, 1)
      status <- dplyr::case_when(
        prevalence_pct < 10 | prevalence_pct > 90 ~ "imbalanced",
        prevalence_pct < 20 | prevalence_pct > 80 ~ "moderate_imbalance",
        TRUE                                        ~ "ok"
      )
      tibble::tibble(variable = var, type = type,
                     n_present = n_present, prevalence_pct = prevalence_pct,
                     min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
                     status = status)
    } else {
      pct_zero <- mean(x == 0, na.rm = TRUE) * 100
      status   <- if (pct_zero > 50) "zero_inflated" else "ok"
      tibble::tibble(variable = var, type = type,
                     n_present = sum(x > 0, na.rm = TRUE),
                     prevalence_pct = NA_real_,
                     min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
                     status = status)
    }
  })

  print(report, n = 20)
  lichen_message("Response variable validation complete")
  report
}
