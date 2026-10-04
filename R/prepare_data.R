library(here)
library(tidyverse)
source(here("R", "io.R"))

#' Read MIR spectra and average replicate scans
#'
#' Loads OPUS spectra and averages the replicate scans of each sample, giving
#' one row per sample. Use it on its own for new spectra that have no soil
#' property data yet.
#'
#' @param spectra_dir Path to the directory containing OPUS `.0` files, or a
#'   vector of paths to individual OPUS files.
#'
#' @return A tibble with `smp_id` and one column per wavenumber.
#' @export
prepare_spectra <- function(spectra_dir) {
  read_opus_tbl(spectra_dir) |>
    mutate(smp_id = str_extract(spc_id, "RT_\\d+")) |>
    group_by(smp_id) |>
    summarize(across(starts_with("wn_"), mean))
}

#' Prepare Ring Trial MIR spectra and soil property data
#'
#' Loads OPUS spectra, averages replicate scans, loads soil properties,
#' sets RT_44's placeholder clay value to NA, and joins them into a single tibble.
#'
#' @param spectra_dir Path to the directory containing OPUS `.0` files.
#' @param soil_prop_file Path to the soil properties CSV file.
#'
#' @return A tibble with averaged spectra and soil properties joined by sample ID.
#' @export
prepare_rt <- function(spectra_dir, soil_prop_file) {
  soil_prop <- read_csv(soil_prop_file, show_col_types = FALSE) |>
    # RT_44's clay value is 0, but every other property of RT_44 except pH is
    # missing. A clay content of exactly 0 is implausible: treat it as missing.
    mutate(clay_perc = if_else(sample_id == "RT_44" & clay_perc == 0, NA_real_, clay_perc))

  prepare_spectra(spectra_dir) |>
    left_join(soil_prop, by = c("smp_id" = "sample_id"))
}
