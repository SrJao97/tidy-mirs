library(opusreader2)
library(purrr)
library(tibble)
library(dplyr)
library(here)

#' Read OPUS files into a tidy tibble
#'
#' Reads all OPUS binary files from a folder and returns a wide tibble with
#' one row per spectrum. Wavenumbers are used as column names and absorbance
#' values as cell values. All spectra must share the same wavenumber grid.
#'
#' @param path Path to a folder containing OPUS files
#' @param wn_format Function that takes a numeric wavenumber vector and returns
#'   character column names. Defaults to `sprintf("wn_%.1f", wn)`.
#' @return A tibble with one row per spectrum, `sample_id` as the first column,
#'   followed by one column per wavenumber.
#' @export
#' @examples
#' spectra_tbl <- read_opus_tbl(here("data/RT_spectra"))
read_opus_tbl <- function(path, wn_format = function(wn) sprintf("wn_%.1f", wn)) {
  spectra <- read_opus(dsn = path)
  wns <- wn_format(spectra[[1]]$ab$wavenumbers)

  map_dfr(names(spectra), function(id) {
    ab_vec <- as.vector(spectra[[id]]$ab$data)
    names(ab_vec) <- wns
    as_tibble_row(ab_vec) |>
      mutate(spc_id = tools::file_path_sans_ext(id), .before = 1)
  })
}
