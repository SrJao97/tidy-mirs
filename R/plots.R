library(tidyverse)
theme_set(theme_bw())

# ── Plotting spectra: Overview and highlighting ─────────────────────────────────────────────────
plot_spectra <- function(
    data,
    cols = tidyselect::starts_with("wn_"),
    names_prefix = "wn_",
    id_col = "smp_id",
    highlight = NULL,
    alpha = 0.3
) {
  spectra <- data |>
    tidyr::pivot_longer(
      cols = {{ cols }},
      names_to = "wavenumber",
      names_prefix = names_prefix,
      names_transform = list(wavenumber = as.numeric),
      values_to = "absorbance"
    )

  p <- ggplot2::ggplot(
    spectra,
    ggplot2::aes(
      x = wavenumber,
      y = absorbance,
      group = .data[[id_col]]
    )
  )

  if (is.null(highlight)) {
    p <- p +
      ggplot2::geom_line(
        alpha = alpha,
        colour = "steelblue"
      )
  } else {
    background <- spectra |>
      dplyr::filter(!.data[[id_col]] %in% highlight)

    selected <- spectra |>
      dplyr::filter(.data[[id_col]] %in% highlight)

    p <- p +
      ggplot2::geom_line(
        data = background,
        alpha = alpha,
        colour = "grey75"
      ) +
      ggplot2::geom_line(
        data = selected,
        ggplot2::aes(colour = .data[[id_col]]),
        alpha = 1
      ) +
      ggplot2::labs(colour = "Sample")
  }

  p +
    ggplot2::scale_x_reverse(
      breaks = seq(500, 4000, by = 250),
      minor_breaks = seq(500, 4000, by = 125)
    ) +
    ggplot2::labs(
      x = "Wavenumber (cm-1)",
      y = "Absorbance"
    ) +
    ggplot2::theme_minimal()
}