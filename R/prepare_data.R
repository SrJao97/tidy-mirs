library(here)
library(tidyverse)
source(here("R", "io.R"))

#' Prepare Ring Trial MIR spectra and soil property data
#'
#' Loads OPUS spectra, averages replicate scans, loads soil properties,
#' and joins them into a single tibble.
#'
#' @param spectra_dir Path to the directory containing OPUS `.O` files.
#' @param soil_prop_file Path to the soil properties CSV file.
#'
#' @return A tibble with averaged spectra and soil properties joined by sample ID.
#' @export
prepare_rt <- function(spectra_dir, soil_prop_file) {
  spectra <- read_opus_tbl(spectra_dir)

  spectra_mean <- spectra |>
    mutate(smp_id = str_extract(spc_id, "RT_\\d+")) |>
    group_by(smp_id) |>
    summarize(across(starts_with("wn_"), mean))

  soil_prop <- read_csv(soil_prop_file, show_col_types = FALSE)

  spectra_mean |>
    left_join(soil_prop, by = c("smp_id" = "sample_id"))
}
