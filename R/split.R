pca_scores <- function(data, cols, max_expl_var, use_snv = TRUE) {
  wn_matrix <- as.matrix(dplyr::select(data, {{ cols }}))

  xr <- if (use_snv) prospectr::standardNormalVariate(wn_matrix) else wn_matrix

  resemble::ortho_projection(
    Xr     = xr,
    ncomp  = resemble::ncomp_by_cumvar(max_expl_var),
    center = TRUE,
    scale  = FALSE
  )$scores
}

kennard_stone_split <- function(data,
                                 prop         = 0.75,
                                 cols,
                                 max_expl_var = 0.99,
                                 seed         = 123L) {
  n <- nrow(data)

  set.seed(seed)
  ks <- prospectr::kenStone(
    X = pca_scores(data, {{ cols }}, max_expl_var),
    k = floor(n * prop), metric = "mahal",
    .center = TRUE, .scale = FALSE
  )
  pool_idx <- as.integer(sort(ks$model))
  test_idx <- as.integer(sort(ks$test))

  message(glue::glue(
    "KS split: {length(pool_idx)} train | {length(test_idx)} test"
  ))

  rsample::make_splits(
    x    = list(analysis = pool_idx, assessment = test_idx),
    data = data
  )
}
