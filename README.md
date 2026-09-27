# Soil mid-infrared spectroscopy with tidymodels

A one-week course on soil mid-infrared (MIR) spectroscopy in R, built as a [Quarto](https://quarto.org) website. It is for soil scientists and spectroscopists who want to build spectral models with the [tidyverse](https://tidyverse.org/) and [tidymodels](https://www.tidymodels.org). You load raw spectrometer output, preprocess it, find outliers, tune and test a model, and then run the same pipeline on the [Open Soil Spectral Library (OSSL)](https://docs.soilspectroscopy.org/). Each page shows a complete step first, then explains how it works. You need little R experience to start.

## Set up

You need R 4.4 or later. The course was built with R 4.5.3 and [Positron](https://positron.posit.co). To build the website, you also need [Quarto](https://quarto.org/docs/get-started/).

Get the code:

```bash
git clone https://github.com/franckalbinet/tidy-mirs.git
```

Open the folder in Positron, then install the pinned package versions from `renv.lock` in the R console:

```r
renv::restore()
```

`renv` installs the packages into a library inside the project. It does not change your main R library. If `renv::restore()` fails, `0-gettting-setup.qmd` lists the packages to install by hand.

## Course pages

| Page | Content |
|---|---|
| `index.qmd` | Course home page and teaching approach |
| `0-gettting-setup.qmd` | Install R, Positron and the packages |
| `1-full-pipeline.qmd` | The whole pipeline once, end to end |
| `2-loading-data.qmd` | Read OPUS files and a wet chemistry CSV into one tibble |
| `3-eda.qmd` | Explore the spectra and the soil properties |
| `4-feature-engineering.qmd` | Resampling, SNV and Savitzky-Golay as `recipes` steps |
| `5-outlier-detection.qmd` | PCA and Mahalanobis distance to find outliers |
| `6-evaluating-performance.qmd` | Repeated cross-validation, tuning and the one-standard-error rule |
| `7-preparing-ossl.qmd` | Download OSSL and save a working sample for exchangeable potassium |
| `8-kex-ossl-pipeline.qmd` | Predict exchangeable potassium from OSSL with PLS and Cubist |

## Project layout

- `R/` holds the helper scripts that the pages load with `source()`. `R/steps.R` defines the spectral `recipes` steps `step_resample()`, `step_snv()` and `step_sg()`, which wrap functions from `prospectr`. `R/io.R` reads Bruker OPUS files. `R/prepare_data.R` builds the RT dataset. `R/split.R` has the Kennard-Stone split. `R/plots.R` has `plot_spectra()`.
- `data/RT_spectra/` holds the OPUS files for the 70 RT soil samples, with 4 scans per sample. `data/RT_wetchem_soildata.csv` holds their wet chemistry.
- `data/ossl_mir_k_sample.rds` is the OSSL working sample that `7-preparing-ossl.qmd` saves: 4,548 samples with exchangeable potassium and a MIR spectrum.
- `_freeze/` stores the results of each page's code. Commit it with the pages.

## Build the website

Preview the whole site:

```bash
quarto preview
```

Render it to `_site/`:

```bash
quarto render
```

The project uses `freeze: auto`. Quarto runs a page's code again only when that page's source changes. Some pages take minutes to run:

- `7-preparing-ossl.qmd` downloads and reads the whole OSSL file. `data/ossl_mir_k_sample.rds` is in the repository, so you need to run this page only to change the sample.
- `8-kex-ossl-pipeline.qmd` tunes PLS and Cubist models in about 3 minutes.
- `6-evaluating-performance.qmd` fits 1,900 PLS models for repeated cross-validation.

To work on one page, preview only that page, for example `quarto preview 5-outlier-detection.qmd`.

## References

- Safanelli, J. L. et al. (2025). [Open Soil Spectral Library (OSSL): Building reproducible soil calibration models through open development and community engagement](https://doi.org/10.1371/journal.pone.0296545). *PLOS ONE*, 20(1).
- Albinet, F., Peng, Y., Eguchi, T., Smolders, E. and Dercon, G. (2022). [Prediction of exchangeable potassium in soil through mid-infrared spectroscopy and deep learning: From prediction to explainability](https://doi.org/10.1016/j.aiia.2022.10.001). *Artificial Intelligence in Agriculture*, 6, 230-241.
- Kuhn, M. and Silge, J. [Tidy Modeling with R](https://www.tmwr.org/).
