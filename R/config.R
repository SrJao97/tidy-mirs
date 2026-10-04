# Shared settings for the course pages.
# Each page loads it with source(..., echo = TRUE), which prints every line,
# so the rendered page shows each value it uses.
# Change a value here, re-render, and every page uses the new value.

# Samples left out of modelling.
# Outlier detection flags their spectra as the two most unusual.
EXCLUDED_SAMPLES <- c("RT_44", "RT_42")

# Spectral preprocessing
RESAMPLING_GRID <- seq(3995, 600, by = -5)  # target wavenumbers, in cm-1
SG_W <- 15L  # Savitzky-Golay window width, in grid points
SG_P <- 3L   # Savitzky-Golay polynomial order
SG_M <- 1L   # Savitzky-Golay derivative order (0 = smoothing only)

# Cutoff on the scaled Mahalanobis distance of a spectrum.
# Chosen by eye on the Outlier detection page, reused at inference.
OUTLIER_CUTOFF <- 1.9
