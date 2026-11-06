# Spatial cross-validation ------------------------------------------------------

#' Assign points to spatial cross-validation folds
#'
#' Random cross-validation of species distribution models is over-optimistic
#' because training and test points lie next to each other and share the same
#' environment (spatial autocorrelation). Spatial blocking splits the study
#' area into square blocks and keeps each block entirely in one fold, so
#' every model is tested on areas it has not seen.
#'
#' Points are projected to an equal-area projection centred on the data (when
#' they are in longitude/latitude), cut into square blocks of `block_size` km,
#' and the blocks are assigned to `k` folds at random. Out of `n_iter` random
#' assignments the one with the most even numbers of presences and of
#' absences/background points per fold is kept.
#'
#' **Choosing `block_size`.** Blocks should be larger than the range of
#' spatial autocorrelation of the predictors, which can be estimated with
#' `blockCV::cv_spatial_autocor()`. As a rule of thumb, values between 50 and
#' 500 km work for regional to continental studies. With few presences, use
#' fewer folds (e.g. `k = 3`) or smaller blocks so that every fold still
#' contains presences.
#'
#' @param points Points to split: a data.frame with coordinates, a CSV path,
#'   an `sf` object or a [terra::SpatVector].
#' @param k Number of folds.
#' @param block_size Block side length in kilometres.
#' @param occ Optional 0/1 vector (1 = presence) of the same length as
#'   `points`. When given, presences and absences are balanced separately.
#' @param n_iter Number of random block-to-fold assignments tried.
#' @param seed Optional random seed.
#' @inheritParams extract_predictions
#' @return An integer vector of fold numbers (1 to `k`), one per point, with
#'   attributes `block` (block id of each point) and `block_size_km`.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))
#' pts$fold <- spatial_folds(pts, k = 5, block_size = 100, occ = pts$occ,
#'                           seed = 1)
#' table(pts$fold, pts$occ)
spatial_folds <- function(points, k = 5, block_size, occ = NULL,
                          n_iter = 100, seed = NULL, coords = NULL,
                          crs = "EPSG:4326") {
  if (missing(block_size) || !is.numeric(block_size) || block_size <= 0) {
    stop("Give `block_size` in km (e.g. 100). See ?spatial_folds for how ",
         "to choose it.", call. = FALSE)
  }
  pts <- as_points(points, coords, crs)
  n <- nrow(pts)
  if (is.null(occ)) occ <- rep(1L, n)
  if (length(occ) != n) {
    stop("`occ` must have one value per point (", n, ").", call. = FALSE)
  }
  if (terra::is.lonlat(pts)) {
    xy <- terra::crds(pts)
    laea <- sprintf(paste0("+proj=laea +lat_0=%.4f +lon_0=%.4f +x_0=0 ",
                           "+y_0=0 +datum=WGS84 +units=m +no_defs"),
                    mean(xy[, 2]), mean(xy[, 1]))
    pts <- terra::project(pts, laea)
  }
  xy <- terra::crds(pts)
  size <- block_size * 1000
  block <- paste(floor(xy[, 1] / size), floor(xy[, 2] / size))
  ub <- unique(block)
  bi <- match(block, ub)
  if (length(ub) < k) {
    stop("Only ", length(ub), " blocks contain points, fewer than k = ", k,
         ". Use a smaller `block_size` or fewer folds.", call. = FALSE)
  }
  np <- tabulate(bi[occ == 1], length(ub))
  na <- tabulate(bi[occ != 1], length(ub))
  cv <- function(v) if (mean(v) == 0) 0 else stats::sd(v) / mean(v)

  best <- NULL
  best_score <- Inf
  with_seed(seed, {
    for (i in seq_len(max(1, n_iter))) {
      f <- sample(rep_len(seq_len(k), length(ub)))
      sp <- tabulate(rep(f, np), k)
      sa <- tabulate(rep(f, na), k)
      # a fold without presences cannot be evaluated: heavy penalty
      score <- cv(sp) + cv(sa) + 100 * sum(sp == 0 & sum(np) >= k)
      if (score < best_score) {
        best <- f
        best_score <- score
      }
    }
  })
  fold <- best[bi]
  per_fold <- tabulate(fold[occ == 1], k)
  if (any(per_fold == 0) && any(occ != 1)) {
    warning("Fold(s) ", paste(which(per_fold == 0), collapse = ", "),
            " contain no presences; AUC cannot be computed for them. Use ",
            "fewer folds or smaller blocks.", call. = FALSE)
  } else if (any(per_fold < 5) && any(occ != 1)) {
    warning("Some folds have fewer than 5 presences (", paste(per_fold,
            collapse = ", "), "); their metrics will be noisy.",
            call. = FALSE)
  }
  attr(fold, "block") <- bi
  attr(fold, "block_size_km") <- block_size
  fold
}

#' Evaluate cross-validated models on their held-out folds
#'
#' Summarises spatial (or any) cross-validation results per model or
#' algorithm. You supply, for every cross-validation run, the predictions of
#' the model *on the points it was not trained on*; the function computes
#' AUC, MAE, bias, TSS and the Boyce index for each run and averages them.
#'
#' Runs can also be repeated over several background (pseudo-absence) sets via
#' `set`. The output then shows whether the background sets or the spatial
#' folds cause more variation (`AUC_sd_sets` vs `AUC_sd_folds`): if the sets
#' vary much less than the folds, adding more background sets gains little.
#'
#' @param data A data.frame in long format with one row per held-out point and
#'   model run, containing the columns named by `model`, `fold`, `obs` and
#'   `pred` (and `set`, if used). Predictions from biomod2, ENMeval, Maxent,
#'   tidymodels or your own code can all be brought into this form (see the
#'   README).
#' @param model Column with the model or algorithm name (results are
#'   summarised per value of this column).
#' @param fold Column with the fold that was held out.
#' @param set Optional column with the background set / replicate.
#' @param obs Column with observations (1 = presence, 0 = absence or
#'   background).
#' @param pred Column with predicted suitability (0-1).
#' @param threshold Optional column with a fixed threshold per run (e.g. the
#'   TSS-optimal or 10th-percentile threshold found on that run's
#'   *calibration* data, as biomod2 does for its validation TSS). Without it,
#'   TSS uses the best threshold on the held-out data, which is slightly
#'   optimistic.
#' @return An object of class `enm_cv`: a data.frame with one row per model,
#'   the mean of each metric over runs, their SDs (`AUC_sd`, ...), `n_runs`
#'   and, when available, `AUC_sd_folds` and `AUC_sd_sets`. Per-run metrics
#'   are in `attr(x, "runs")` and mean AUC per fold in `attr(x, "by_fold")`
#'   (a fold that is low for every model points to a region that is hard to
#'   predict, not to a bad model).
#' @export
#' @examples
#' # see the README for a full example with spatial folds and two GLMs
#' set.seed(1)
#' d <- data.frame(model = rep(c("good", "poor"), each = 200),
#'                 fold = rep(1:4, 100), obs = rep(rbinom(200, 1, 0.3), 2))
#' d$pred <- ifelse(d$model == "good", plogis(3 * d$obs - 1.5 +
#'                  rnorm(400)), runif(400))
#' evaluate_cv(d)
evaluate_cv <- function(data, model = "model", fold = "fold", set = NULL,
                        obs = "obs", pred = "pred", threshold = NULL) {
  cols <- c(model, fold, set, obs, pred, threshold)
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  if (!all(data[[obs]] %in% c(0, 1))) {
    stop("`", obs, "` must contain only 0 and 1.", call. = FALSE)
  }
  check_prediction_scale(data[pred])
  set_v <- if (is.null(set)) rep(1L, nrow(data)) else data[[set]]
  key <- interaction(data[[model]], set_v, data[[fold]], drop = TRUE,
                     lex.order = TRUE, sep = "\r")
  idx <- split(seq_len(nrow(data)), key)
  runs <- do.call(rbind, lapply(idx, function(i) {
    m <- metric_row(data[[obs]][i], data[[pred]][i])
    if (!is.null(threshold)) {
      thr <- unique(data[[threshold]][i])
      if (length(thr) != 1) {
        stop("`", threshold, "` must be constant within each run.",
             call. = FALSE)
      }
      ok <- !is.na(data[[pred]][i])
      m[["TSS"]] <- tss_score(data[[obs]][i][ok], data[[pred]][i][ok],
                              threshold = thr)[1]
    }
    data.frame(model = as.character(data[[model]][i[1]]), set = set_v[i[1]],
               fold = data[[fold]][i[1]], t(m), check.names = FALSE)
  }))
  rownames(runs) <- NULL
  runs$n_presence <- as.integer(runs$n_presence)
  runs$n_absence <- as.integer(runs$n_absence)
  bad <- runs$n_presence == 0 | runs$n_absence == 0
  if (any(bad)) {
    warning(sum(bad), " run(s) have no presences or no absences in the ",
            "held-out fold; their AUC/TSS/Boyce are NA and are ignored in ",
            "the means.", call. = FALSE)
  }

  mets <- c("AUC", "MAE", "bias", "TSS", "Boyce")
  out <- do.call(rbind, lapply(split(runs, factor(runs$model,
                                                  unique(runs$model))),
                               function(r) {
    mn <- vapply(r[mets], mean, numeric(1), na.rm = TRUE)
    sdv <- vapply(r[mets], stats::sd, numeric(1), na.rm = TRUE)
    d <- data.frame(model = r$model[1], n_runs = sum(!is.na(r$AUC)),
                    t(mn), t(stats::setNames(sdv, paste0(mets, "_sd"))),
                    check.names = FALSE)
    by_set <- split(r$AUC, r$set)
    d$AUC_sd_folds <- mean(vapply(by_set, stats::sd, numeric(1),
                                  na.rm = TRUE), na.rm = TRUE)
    if (!is.null(set)) {
      d$AUC_sd_sets <- stats::sd(vapply(by_set, mean, numeric(1),
                                        na.rm = TRUE))
    }
    d
  }))
  rownames(out) <- NULL
  out[] <- lapply(out, function(v) if (is.numeric(v)) {
    v[is.nan(v)] <- NA
    v
  } else v)

  by_fold <- stats::aggregate(runs$AUC, list(model = runs$model,
                                             fold = runs$fold),
                              mean, na.rm = TRUE)
  by_fold <- stats::reshape(by_fold, idvar = "model", timevar = "fold",
                            direction = "wide")
  names(by_fold) <- sub("^x\\.", "fold_", names(by_fold))
  rownames(by_fold) <- NULL

  attr(out, "runs") <- runs
  attr(out, "by_fold") <- by_fold
  attr(out, "has_sets") <- !is.null(set)
  class(out) <- c("enm_cv", "data.frame")
  out
}

#' @export
print.enm_cv <- function(x, digits = 3, ...) {
  if (is.null(attr(x, "runs"))) return(print_plain(x, digits, ...))
  runs <- attr(x, "runs")
  cat("Cross-validated evaluation: ", nrow(x), " model(s), ",
      length(unique(runs$fold)), " folds",
      if (isTRUE(attr(x, "has_sets"))) {
        paste0(", ", length(unique(runs$set)), " background sets")
      }, "\n\n", sep = "")
  df <- as.data.frame(x)
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE, ...)
  bf <- attr(x, "by_fold")
  if (!is.null(bf)) {
    cat("\nMean hold-out AUC per fold:\n")
    bf[-1] <- lapply(bf[-1], round, digits = digits)
    print(bf, row.names = FALSE)
  }
  if (isTRUE(attr(x, "has_sets"))) {
    cat("\nMedian SD of AUC between background sets: ",
        round(stats::median(x$AUC_sd_sets, na.rm = TRUE), 4),
        "; between folds: ",
        round(stats::median(x$AUC_sd_folds, na.rm = TRUE), 4), "\n", sep = "")
  }
  invisible(x)
}

# Evaluate `expr` with a temporary random seed, restoring the user's RNG state
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(invisible(force(expr)))
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env)
  on.exit(if (had) assign(".Random.seed", old, envir = env)
          else if (exists(".Random.seed", envir = env, inherits = FALSE))
            rm(".Random.seed", envir = env))
  set.seed(seed)
  invisible(force(expr))
}

# Print a (subsetted) result as a plain rounded data.frame (internal)
print_plain <- function(x, digits = 3, ...) {
  df <- as.data.frame(unclass(x), check.names = FALSE)
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE, ...)
  invisible(x)
}
