library(prospectr)

################################
# SNV
################################
# ── Constructor (user-facing) ─────────────────────────────────────────────────
step_snv <- function(
  recipe,
  ...,
  role    = NA,
  trained = FALSE,
  skip    = FALSE,
  id      = rand_id("snv")
) {
  terms <- enquos(...)
  if (is_empty(terms)) terms <- quos(all_predictors())  # <-- default
  
  add_step(
    recipe,
    step_snv_new(
      terms   = terms,
      cols    = NULL,
      role    = role,
      trained = trained,
      skip    = skip,
      id      = id
    )
  )
}

# ── Internal constructor ───────────────────────────────────────────────────────
step_snv_new <- function(terms, cols, role, trained, skip, id) {
  step(
    subclass = "snv",
    terms    = terms,
    cols     = cols,
    role     = role,
    trained  = trained,
    skip     = skip,
    id       = id
  )
}

# ── prep: resolve selectors → column names (nothing to estimate) ──────────────
prep.step_snv <- function(x, training, info = NULL, ...) {
  cols <- recipes_eval_select(x$terms, training, info)
  
  step_snv_new(
    terms   = x$terms,
    cols    = cols,
    role    = x$role,
    trained = TRUE,
    skip    = x$skip,
    id      = x$id
  )
}

# ── bake: apply SNV row-wise ───────────────────────────────────────────────────
bake.step_snv <- function(object, new_data, ...) {
  check_new_data(object$cols, object, new_data)
  
  spc     <- as.matrix(new_data[, object$cols])
  spc_snv <- prospectr::standardNormalVariate(spc)
  colnames(spc_snv) <- object$cols
  
  new_data |>
    select(-all_of(object$cols)) |>
    bind_cols(as_tibble(spc_snv))
}

# ── print ──────────────────────────────────────────────────────────────────────
print.step_snv <- function(x, width = max(20, options()$width - 35), ...) {
  print_step(
    tr_obj   = x$cols,
    untr_obj = x$terms,
    trained  = x$trained,
    title    = "Standard Normal Variate (SNV) on ",
    width    = width
  )
  invisible(x)
}

# ── tidy: one row per selected column ─────────────────────────────────────────
tidy.step_snv <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble::tibble(terms = x$cols)
  } else {
    res <- tibble::tibble(terms = recipes::sel2char(x$terms))
  }
  res$id <- x$id
  res
}


################################
# SAVITSKY-GOLAY
################################
# ── Constructor (user-facing) ─────────────────────────────────────────────────
step_sg <- function(
  recipe,
  ...,
  w       = 11L,
  p       = 2L,
  m       = 1L,
  role    = NA,
  trained = FALSE,
  skip    = FALSE,
  id      = rand_id("sg")
) {
  terms <- enquos(...)
  if (is_empty(terms)) terms <- quos(all_predictors())
  
  add_step(
    recipe,
    step_sg_new(
      terms   = terms,
      w       = w,
      p       = p,
      m       = m,
      cols    = NULL,
      role    = role,
      trained = trained,
      skip    = skip,
      id      = id
    )
  )
}

# ── Internal constructor ───────────────────────────────────────────────────────
step_sg_new <- function(terms, w, p, m, cols, role, trained, skip, id) {
  step(
    subclass = "sg",
    terms    = terms,
    w        = w,
    p        = p,
    m        = m,
    cols     = cols,
    role     = role,
    trained  = trained,
    skip     = skip,
    id       = id
  )
}

# ── prep: resolve selectors → column names ────────────────────────────────────
prep.step_sg <- function(x, training, info = NULL, ...) {
  cols <- recipes_eval_select(x$terms, training, info)
  
  step_sg_new(
    terms   = x$terms,
    w       = x$w,
    p       = x$p,
    m       = x$m,
    cols    = cols,
    role    = x$role,
    trained = TRUE,
    skip    = x$skip,
    id      = x$id
  )
}

# ── bake: apply Savitzky-Golay, preserve trimmed wavenumber names ─────────────
bake.step_sg <- function(object, new_data, ...) {
  check_new_data(object$cols, object, new_data)
  
  spc    <- as.matrix(new_data[, object$cols])
  spc_sg <- prospectr::savitzkyGolay(
    spc,
    w = object$w,
    p = object$p,
    m = object$m
  )
  
  # savitzkyGolay() trims (w-1)/2 points from each end
  trim      <- (object$w - 1L) / 2L
  kept_cols <- object$cols[seq(trim + 1L, length(object$cols) - trim)]
  colnames(spc_sg) <- kept_cols
  
  new_data |>
    select(-all_of(object$cols)) |>
    bind_cols(as_tibble(spc_sg))
}

# ── print ─────────────────────────────────────────────────────────────────────
print.step_sg <- function(x, width = max(20, options()$width - 35), ...) {
  title <- glue::glue("Savitzky-Golay (w={x$w}, p={x$p}, m={x$m}) on ")
  print_step(
    tr_obj   = x$cols,
    untr_obj = x$terms,
    trained  = x$trained,
    title    = title,
    width    = width
  )
  invisible(x)
}

# ── tidy: one row per selected column ─────────────────────────────────────────
tidy.step_sg <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble::tibble(
      terms = x$cols,
      w     = x$w,
      p     = x$p,
      m     = x$m
    )
  } else {
    res <- tibble::tibble(
      terms = recipes::sel2char(x$terms),
      w     = x$w,
      p     = x$p,
      m     = x$m
    )
  }
  res$id <- x$id
  res
}

# ── tunable: w, p, and m are all tunable ──────────────────────────────────────
tunable.step_sg <- function(x, ...) {
  tibble::tibble(
    name      = c("w", "p", "m"),
    call_info = list(
      list(pkg = NULL, fun = "sg_window"),
      list(pkg = NULL, fun = "sg_degree"),
      list(pkg = NULL, fun = "sg_diff_order")
    ),
    source       = "recipe",
    component    = "step_sg",
    component_id = x$id
  )
}


#################################
# RESAMPLE
#################################

# # ── Constructor (user-facing) ─────────────────────────────────────────────────
# step_resample <- function(
#   recipe,
#   ...,
#   by      = 5,
#   role    = NA,
#   trained = FALSE,
#   skip    = FALSE,
#   id      = rand_id("resample")
# ) {
#   terms <- enquos(...)
#   if (is_empty(terms)) terms <- quos(all_predictors())
  
#   add_step(
#     recipe,
#     step_resample_new(
#       terms    = terms,
#       by       = by,
#       cols     = NULL,
#       wns      = NULL,   # original wavenumbers, learned in prep
#       new_wns  = NULL,   # target wavenumber grid, learned in prep
#       role     = role,
#       trained  = trained,
#       skip     = skip,
#       id       = id
#     )
#   )
# }

# # ── Internal constructor ───────────────────────────────────────────────────────
# step_resample_new <- function(terms, by, cols, wns, new_wns, role, trained, skip, id) {
#   step(
#     subclass = "resample",
#     terms    = terms,
#     by       = by,
#     cols     = cols,
#     wns      = wns,
#     new_wns  = new_wns,
#     role     = role,
#     trained  = trained,
#     skip     = skip,
#     id       = id
#   )
# }

# # ── prep: resolve selectors and fix wavenumber grids ─────────────────────────
# prep.step_resample <- function(x, training, info = NULL, ...) {
#   cols    <- recipes_eval_select(x$terms, training, info)
#   wns     <- get_wns(training)
#   new_wns <- seq(min(wns), max(wns), by = x$by)
  
#   # Inherit the role from the original spectral columns
#   # so downstream steps can find them via all_predictors()
#   col_info <- info[info$variable %in% cols, ]
#   inherited_role <- unique(col_info$role)  # should be "predictor"
  
#   step_resample_new(
#     terms   = x$terms,
#     by      = x$by,
#     cols    = cols,
#     wns     = wns,
#     new_wns = new_wns,
#     role    = inherited_role,   # <-- was NA, now "predictor"
#     trained = TRUE,
#     skip    = x$skip,
#     id      = x$id
#   )
# }

# # ── bake: resample onto the fixed grid learned in prep ────────────────────────
# bake.step_resample <- function(object, new_data, ...) {
#   check_new_data(object$cols, object, new_data)
  
#   spc    <- as.matrix(new_data[, object$cols])
#   spc_rs <- prospectr::resample(
#     spc,
#     wav      = object$wns,
#     new.wav  = object$new_wns,
#     interpol = "linear"
#   )
#   new_col_names <- paste0("X", round(object$new_wns))
#   colnames(spc_rs) <- new_col_names
  
#   new_data |>
#     select(-all_of(object$cols)) |>
#     bind_cols(as_tibble(spc_rs))
# }

# # ── print ─────────────────────────────────────────────────────────────────────
# print.step_resample <- function(x, width = max(20, options()$width - 35), ...) {
#   title <- glue::glue("Spectral resampling (by={x$by}) on ")
#   print_step(
#     tr_obj   = x$cols,
#     untr_obj = x$terms,
#     trained  = x$trained,
#     title    = title,
#     width    = width
#   )
#   invisible(x)
# }

# # ── tidy: one row per selected column ─────────────────────────────────────────
# tidy.step_resample <- function(x, ...) {
#   if (is_trained(x)) {
#     res <- tibble::tibble(
#       terms   = x$cols,
#       by      = x$by,
#       wns     = x$wns,
#       new_wns = x$new_wns
#     )
#   } else {
#     res <- tibble::tibble(
#       terms   = recipes::sel2char(x$terms),
#       by      = x$by,
#       wns     = NA_real_,
#       new_wns = NA_real_
#     )
#   }
#   res$id <- x$id
#   res
# }



###################################
# RESAMPLE RE-FACTORED
###################################
# ── Helpers ────────────────────────────────────────────────────────────────────

infer_wn_spec <- function(cols) {
  # Matches decimal, integer, signed, and scientific-notation numbers.
  pattern <- "[-+]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?"

  matches <- lapply(cols, function(x) {
    loc <- gregexpr(pattern, x, perl = TRUE)[[1]]

    if (length(loc) != 1L || loc[[1]] == -1L) {
      rlang::abort(
        "Each spectral predictor name must contain exactly one numeric wavenumber."
      )
    }

    match_length <- attr(loc, "match.length")

    if (length(match_length) != 1L) {
      rlang::abort(
        "Each spectral predictor name must contain exactly one numeric wavenumber."
      )
    }

    token <- substr(
      x,
      loc[[1]],
      loc[[1]] + match_length[[1]] - 1L
    )

    list(
      value  = as.numeric(token),
      token  = token,
      prefix = substr(x, 1L, loc[[1]] - 1L),
      suffix = substr(
        x,
        loc[[1]] + match_length[[1]],
        nchar(x)
      )
    )
  })

  values <- vapply(matches, `[[`, numeric(1), "value")
  tokens <- vapply(matches, `[[`, character(1), "token")
  prefix <- vapply(matches, `[[`, character(1), "prefix")
  suffix <- vapply(matches, `[[`, character(1), "suffix")

  if (anyNA(values) || any(!is.finite(values))) {
    rlang::abort(
      "Could not extract finite numeric wavenumbers from predictor names."
    )
  }

  if (length(unique(prefix)) != 1L ||
      length(unique(suffix)) != 1L) {
    rlang::abort(
      "Spectral predictor names must share a common prefix and suffix."
    )
  }

  decimal_places <- function(x) {
    mantissa <- sub("[eE].*$", "", x)

    if (!grepl("\\.", mantissa)) {
      return(0L)
    }

    nchar(sub("^[^.]*\\.", "", mantissa))
  }

  list(
    values = values,
    prefix = prefix[[1]],
    suffix = suffix[[1]],
    digits = max(vapply(tokens, decimal_places, integer(1)))
  )
}


digits_required <- function(x) {
  x <- format(
    x,
    scientific = FALSE,
    trim = TRUE,
    digits = 15
  )

  if (!grepl("\\.", x)) {
    return(0L)
  }

  nchar(sub("^[^.]*\\.", "", x))
}


format_wn_names <- function(wns, spec, digits) {
  values <- formatC(
    wns,
    format = "f",
    digits = digits
  )

  has_decimal <- grepl("\\.", values)

  values[has_decimal] <- sub(
    "0+$",
    "",
    values[has_decimal]
  )

  values[has_decimal] <- sub(
    "\\.$",
    "",
    values[has_decimal]
  )

  paste0(
    spec$prefix,
    values,
    spec$suffix
  )
}



# ── Constructor ────────────────────────────────────────────────────────────────

step_resample <- function(
  recipe,
  ...,
  grid = NULL,
  by = NULL,
  interpol = c("linear", "spline"),
  name_fun = NULL,
  role = NA,
  trained = FALSE,
  skip = FALSE,
  id = recipes::rand_id("resample")
) {
  interpol <- match.arg(interpol)
  terms <- rlang::enquos(...)

  if (rlang::is_empty(terms)) {
    terms <- rlang::quos(all_predictors())
  }

  recipes::add_step(
    recipe,
    step_resample_new(
      terms = terms,
      grid = grid,
      by = by,
      interpol = interpol,
      name_fun = name_fun,
      spec = NULL,
      cols = NULL,
      wns = NULL,
      new_wns = NULL,
      new_cols = NULL,
      role = role,
      trained = trained,
      skip = skip,
      id = id
    )
  )
}


# ── Internal constructor ───────────────────────────────────────────────────────

step_resample_new <- function(
  terms,
  grid,
  by,
  interpol,
  name_fun,
  spec,
  cols,
  wns,
  new_wns,
  new_cols,
  role,
  trained,
  skip,
  id
) {
  recipes::step(
    subclass = "resample",
    terms = terms,
    grid = grid,
    by = by,
    interpol = interpol,
    name_fun = name_fun,
    spec = spec,
    cols = cols,
    wns = wns,
    new_wns = new_wns,
    new_cols = new_cols,
    role = role,
    trained = trained,
    skip = skip,
    id = id
  )
}


# ── Prep ───────────────────────────────────────────────────────────────────────

prep.step_resample <- function(x, training, info = NULL, ...) {
  if (!is.null(x$grid) && !is.null(x$by)) {
    rlang::abort("Supply either `grid` or `by`, not both.")
  }

  if (is.null(x$grid) && is.null(x$by)) {
    rlang::abort("Supply either `grid` or `by`.")
  }

  if (!is.null(x$by)) {
    if (length(x$by) != 1L ||
        !is.numeric(x$by) ||
        !is.finite(x$by) ||
        x$by <= 0) {
      rlang::abort("`by` must be one positive, finite number.")
    }
  }

  cols <- recipes::recipes_eval_select(
    x$terms,
    training,
    info
  )

  if (length(cols) < 2L) {
    rlang::abort(
      "step_resample() requires at least two spectral predictor columns."
    )
  }

  spec <- infer_wn_spec(cols)
  wns <- spec$values

  if (anyDuplicated(wns)) {
    rlang::abort(
      "Spectral predictor names must represent unique wavenumbers."
    )
  }

  axis_direction <- sign(wns[[length(wns)]] - wns[[1]])

  if (axis_direction == 0L ||
      !all(sign(diff(wns)) == axis_direction)) {
    rlang::abort(
      "Spectral predictor names must form a strictly monotonic wavenumber grid."
    )
  }

  observed_min <- min(wns)
  observed_max <- max(wns)

  # Explicit target grid.
  if (!is.null(x$grid)) {
    new_wns <- as.numeric(x$grid)

    if (length(new_wns) < 2L ||
        anyNA(new_wns) ||
        any(!is.finite(new_wns))) {
      rlang::abort(
        "`grid` must contain at least two finite numeric wavenumbers."
      )
    }

  # Convenient inferred grid.
  } else {
    if (axis_direction < 0) {
      start <- observed_max
      end <- observed_min
      step <- -x$by
    } else {
      start <- observed_min
      end <- observed_max
      step <- x$by
    }

    new_wns <- seq(
      from = start,
      to = end,
      by = step
    )
  }

  if (length(new_wns) < 2L ||
      anyNA(new_wns) ||
      any(!is.finite(new_wns))) {
    rlang::abort(
      "The resampling grid must contain at least two finite values."
    )
  }

  grid_diffs <- diff(new_wns)

  if (any(grid_diffs == 0) ||
      !all(grid_diffs > 0) && !all(grid_diffs < 0)) {
    rlang::abort(
      "The resampling grid must be strictly monotonic."
    )
  }

  # Prevent extrapolation by default.
  if (any(new_wns < observed_min) ||
      any(new_wns > observed_max)) {
    rlang::abort(
      paste0(
        "`grid` extends beyond the observed wavenumber range ",
        "[", observed_min, ", ", observed_max, "]. ",
        "Resampling would require extrapolation."
      )
    )
  }

  # Preserve the inferred decimal precision, while also allowing enough
  # precision for an explicitly supplied `by`.
  digits <- spec$digits

  if (!is.null(x$by)) {
    digits <- max(digits, digits_required(x$by))
  }

  new_cols <- if (is.null(x$name_fun)) {
    format_wn_names(
      wns = new_wns,
      spec = spec,
      digits = digits
    )
  } else {
    x$name_fun(new_wns)
  }

  new_cols <- as.character(new_cols)

  if (length(new_cols) != length(new_wns) ||
      anyNA(new_cols) ||
      anyDuplicated(new_cols)) {
    rlang::abort(
      "`name_fun` must return one unique, non-missing name per wavenumber."
    )
  }

  other_cols <- setdiff(names(training), cols)

  if (any(new_cols %in% other_cols)) {
    rlang::abort(
      "Generated resampled names collide with existing non-spectral columns."
    )
  }

  selected_roles <- unique(
    info$role[match(cols, info$variable)]
  )

  inherited_role <- if (length(selected_roles) == 1L &&
                        !is.na(selected_roles)) {
    selected_roles
  } else {
    "predictor"
  }

  step_resample_new(
    terms = x$terms,
    grid = x$grid,
    by = x$by,
    interpol = x$interpol,
    name_fun = x$name_fun,
    spec = spec,
    cols = cols,
    wns = wns,
    new_wns = new_wns,
    new_cols = new_cols,
    role = inherited_role,
    trained = TRUE,
    skip = x$skip,
    id = x$id
  )
}


# ── Bake ───────────────────────────────────────────────────────────────────────

bake.step_resample <- function(object, new_data, ...) {
  recipes::check_new_data(
    object$cols,
    object,
    new_data
  )

  spc <- as.matrix(
    new_data[, object$cols, drop = FALSE]
  )

  spc_resampled <- prospectr::resample(
    spc,
    wav = object$wns,
    new.wav = object$new_wns,
    interpol = object$interpol
  )

  colnames(spc_resampled) <- object$new_cols

  new_data |>
    dplyr::select(-dplyr::all_of(object$cols)) |>
    dplyr::bind_cols(
      tibble::as_tibble(spc_resampled)
    )
}


# ── Print ──────────────────────────────────────────────────────────────────────

print.step_resample <- function(
  x,
  width = max(20, options()$width - 35),
  ...
) {
  title <- if (is.null(x$grid)) {
    glue::glue(
      "Spectral resampling (by = {x$by}, interpol = {x$interpol}) on "
    )
  } else {
    glue::glue(
      "Spectral resampling (explicit grid, interpol = {x$interpol}) on "
    )
  }

  recipes::print_step(
    tr_obj = x$cols,
    untr_obj = x$terms,
    trained = x$trained,
    title = title,
    width = width
  )

  invisible(x)
}


# ── Tidy ───────────────────────────────────────────────────────────────────────

tidy.step_resample <- function(x, ...) {
  if (recipes::is_trained(x)) {
    tibble::tibble(
      terms = x$cols,
      interpol = x$interpol,
      wns = list(x$wns),
      new_wns = list(x$new_wns),
      new_cols = list(x$new_cols),
      id = x$id
    )
  } else {
    tibble::tibble(
      terms = recipes::sel2char(x$terms),
      interpol = x$interpol,
      wns = list(NA_real_),
      new_wns = list(NA_real_),
      new_cols = list(NA_character_),
      id = x$id
    )
  }
}


####################################
# ── tunable ───────────────────────────────────────────────────────────────────
###################################
tunable.step_resample <- function(x, ...) {
  tibble::tibble(
    name = "by",
    call_info = list(
      list(pkg = "your_pkg", fun = "resample_by")
    ),
    source       = "recipe",
    component    = "step_resample",
    component_id = x$id
  )
}

sg_window <- function(range = c(5L, 21L), trans = NULL) {
  dials::new_quant_param(
    type      = "integer",
    range     = range,
    inclusive = c(TRUE, TRUE),
    trans     = trans,
    label     = c(sg_window = "SG Window Size")
  )
}

sg_degree <- function(range = c(1L, 3L), trans = NULL) {
  dials::new_quant_param(
    type      = "integer",
    range     = range,
    inclusive = c(TRUE, TRUE),
    trans     = trans,
    label     = c(sg_degree = "SG Polynomial Degree")
  )
}

sg_diff_order <- function(range = c(0L, 2L), trans = NULL) {
  dials::new_quant_param(
    type      = "integer",
    range     = range,
    inclusive = c(TRUE, TRUE),
    trans     = trans,
    label     = c(sg_diff_order = "SG Derivative Order")
  )
}

resample_by <- function(range = c(1, 20), trans = NULL) {
  dials::new_quant_param(
    type      = "double",
    range     = range,
    inclusive = c(TRUE, TRUE),
    trans     = trans,
    label     = c(resample_by = "Resampling Interval")
  )
}


# ── Constructor (user-facing) ─────────────────────────────────────────────────
step_rm_co2 <- function(
  recipe,
  ...,
  min_wn  = 2269,
  max_wn  = 2389,
  role    = NA,
  trained = FALSE,
  skip    = FALSE,
  id      = rand_id("rm_co2")
) {
  terms <- enquos(...)
  if (is_empty(terms)) terms <- quos(all_predictors())

  add_step(
    recipe,
    step_rm_co2_new(
      terms   = terms,
      min_wn  = min_wn,
      max_wn  = max_wn,
      cols    = NULL,
      role    = role,
      trained = trained,
      skip    = skip,
      id      = id
    )
  )
}

# ── Internal constructor ───────────────────────────────────────────────────────
step_rm_co2_new <- function(terms, min_wn, max_wn, cols, role, trained, skip, id) {
  step(
    subclass = "rm_co2",
    terms    = terms,
    min_wn   = min_wn,
    max_wn   = max_wn,
    cols     = cols,
    role     = role,
    trained  = trained,
    skip     = skip,
    id       = id
  )
}

# ── prep: identify wavenumber columns within the CO2 band ─────────────────────
prep.step_rm_co2 <- function(x, training, info = NULL, ...) {
  all_cols <- recipes_eval_select(x$terms, training, info)
  wns      <- as.numeric(sub("^X", "", all_cols))
  cols     <- all_cols[!is.na(wns) & wns >= x$min_wn & wns <= x$max_wn]

  step_rm_co2_new(
    terms   = x$terms,
    min_wn  = x$min_wn,
    max_wn  = x$max_wn,
    cols    = cols,
    role    = x$role,
    trained = TRUE,
    skip    = x$skip,
    id      = x$id
  )
}

# ── bake: drop the CO2 band columns ───────────────────────────────────────────
bake.step_rm_co2 <- function(object, new_data, ...) {
  check_new_data(object$cols, object, new_data)
  new_data |> select(-all_of(object$cols))
}

# ── print ──────────────────────────────────────────────────────────────────────
print.step_rm_co2 <- function(x, width = max(20, options()$width - 35), ...) {
  title <- glue::glue("CO2 band removal ({x$min_wn}–{x$max_wn} cm⁻¹) on ")
  print_step(
    tr_obj   = x$cols,
    untr_obj = x$terms,
    trained  = x$trained,
    title    = title,
    width    = width
  )
  invisible(x)
}

# ── tidy: one row per removed column ──────────────────────────────────────────
tidy.step_rm_co2 <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble::tibble(terms = x$cols)
  } else {
    res <- tibble::tibble(terms = recipes::sel2char(x$terms))
  }
  res$id <- x$id
  res
}
