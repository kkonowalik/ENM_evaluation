# Model evaluation and bootstrap ------------------------------------------------

#' Evaluate one or many models and rank them by Pareto optimality
#'
#' The main entry point of the package. Extracts predictions at presence and
#' absence/background points, computes AUC, MAE, bias, TSS and the Boyce
#' index (plus any metrics you add), and ranks models by Pareto optimality on
#' the chosen `criteria`. The default, AUC (maximised) and MAE (minimised), is
#' the combination recommended by Konowalik & Nosol (2021).
#'
#' **Choosing the metrics to compare.** `criteria` can be any of the computed
#' metrics, each with a direction: `"max"` (higher is better), `"min"`
#' (lower is better) or `"zero"` (closer to 0 is better). For example
#' `c(AUC = "max", TSS = "max")`, `c(Boyce = "max", MAE = "min")` or three
#' criteria `c(AUC = "max", MAE = "min", bias = "zero")`. Your own metrics are
#' added with `metrics`, e.g.
#' `metrics = list(Brier = function(obs, pred) mean((obs - pred)^2))`, and can
#' then be used in `criteria` (`c(AUC = "max", Brier = "min")`).
#'
#' When `n_boot > 0`, evaluation points are resampled with replacement
#' (separately for presences and absences, so their numbers stay fixed) and all
#' models are re-evaluated on each resample. This gives:
#'   * percentile confidence intervals for every criterion, and
#'   * `p_pareto`: the proportion of resamples in which the model lies on the
#'     Pareto front. It shows how robust the selection is: a model with
#'     `p_pareto = 0.95` is almost always among the best trade-offs, while
#'     `0.20` means it reached the front partly by chance.
#'
#' @inheritParams extract_predictions
#' @param n_boot Number of bootstrap resamples (0 = no bootstrap). 500 to 1000
#'   is a good choice for final results.
#' @param conf Confidence level for bootstrap intervals (default 0.95).
#' @param criteria Named vector of the metrics to compare and their
#'   direction (`"max"`, `"min"` or `"zero"`); see Details. Default
#'   `c(AUC = "max", MAE = "min")`.
#' @param metrics Optional named list of extra metric functions
#'   `function(obs, pred)` returning a single number; `obs` is 0/1, `pred`
#'   the predicted suitability.
#' @param weights Optional named vector of the relative importance of the
#'   criteria in the compromise score (default: equal); see [pareto_rank()].
#' @param tol Tolerance for Pareto ties; see [pareto_rank()].
#' @param seed Optional random seed for reproducible bootstraps.
#' @param data Optionally, a table returned by [extract_predictions()]. When
#'   supplied, `models`, `presence` and `absence` are ignored. This is also
#'   useful when predictions were extracted with other software: supply a
#'   data.frame with a 0/1 column `obs` and one numeric column per model.
#'
#' @return An object of class `enm_evaluation` (a data.frame) with one row per
#'   model and columns `model`, `AUC`, `MAE`, `bias`, `TSS`, `Boyce`, any
#'   extra metrics, `n_presence`, `n_absence`, `pareto_front` (1 =
#'   non-dominated on the criteria), `distance` (normalised distance to the
#'   ideal point, lower is better) and `rank`. With bootstrapping there are
#'   also `<criterion>_lower`, `<criterion>_upper` and `<criterion>_sd` for
#'   each criterion, and `p_pareto`. The bootstrap replicates are stored in
#'   `attr(x, "bootstrap")` and the criteria in `attr(x, "criteria")`.
#' @export
#' @references Konowalik, K., Nosol, A. (2021) Evaluation metrics and
#'   validation of presence-only species distribution models based on
#'   distributional maps with varying coverage. *Scientific Reports* 11, 1482.
#'   \doi{10.1038/s41598-020-80062-1}
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' res <- evaluate_models(sim$models, sim$presence, sim$absence,
#'                        n_boot = 100, seed = 1)
#' res
#' plot_pareto(res)
#'
#' # compare other metrics, e.g. TSS and the Boyce index
#' res2 <- evaluate_models(sim$models, sim$presence, sim$absence,
#'                         criteria = c(TSS = "max", Boyce = "max"))
#' plot_pareto(res2)
#'
#' # add your own metric and use it as a criterion
#' brier <- function(obs, pred) mean((obs - pred)^2)
#' res3 <- evaluate_models(sim$models, sim$presence, sim$absence,
#'                         metrics = list(Brier = brier),
#'                         criteria = c(AUC = "max", Brier = "min"))
evaluate_models <- function(models = NULL, presence = NULL, absence = NULL,
                            criteria = c(AUC = "max", MAE = "min"),
                            metrics = NULL,
                            n_boot = 0, conf = 0.95,
                            weights = NULL, tol = 0,
                            coords = NULL, crs = "EPSG:4326",
                            na_action = c("common", "per_model"),
                            seed = NULL, data = NULL) {
  na_action <- match.arg(na_action)
  if (is.null(data)) {
    if (is.null(models) || is.null(presence) || is.null(absence)) {
      stop("Supply `models`, `presence` and `absence`, or `data`.",
           call. = FALSE)
    }
    data <- extract_predictions(models, presence, absence, coords = coords,
                                crs = crs, na_action = na_action)
  }
  validate_table(data)
  obs <- data$obs
  preds <- data[setdiff(names(data), "obs")]

  extra <- check_extra_metrics(metrics)
  res <- do.call(rbind, lapply(preds, function(p) metric_row(obs, p, extra)))
  res <- data.frame(model = names(preds), res, row.names = NULL,
                    check.names = FALSE)
  res$n_presence <- as.integer(res$n_presence)
  res$n_absence <- as.integer(res$n_absence)
  criteria <- check_criteria(criteria, c(names(builtin_metrics()),
                                         names(extra)))

  pc <- pareto_core(to_minimise(res, criteria), weights = weights, tol = tol)
  res$pareto_front <- pc$front
  res$distance <- pc$distance

  boot <- NULL
  if (n_boot > 0) {
    funs <- c(builtin_metrics(), extra)[names(criteria)]
    boot <- with_seed(seed, bootstrap_metrics(obs, preds, n_boot, funs,
                                              criteria, tol))
    alpha <- (1 - conf) / 2
    q <- function(x) stats::quantile(x, c(alpha, 1 - alpha), na.rm = TRUE,
                                     names = FALSE)
    stats_by_model <- lapply(split(boot, boot$model), function(b) {
      d <- data.frame(model = b$model[1])
      for (k in names(criteria)) {
        qq <- q(b[[k]])
        d[[paste0(k, "_lower")]] <- qq[1]
        d[[paste0(k, "_upper")]] <- qq[2]
        d[[paste0(k, "_sd")]] <- stats::sd(b[[k]], na.rm = TRUE)
      }
      d$p_pareto <- mean(b$pareto_front == 1)
      d
    })
    bs <- do.call(rbind, stats_by_model)
    res <- cbind(res, bs[match(res$model, bs$model), -1, drop = FALSE])
  }

  # Final ranking: Pareto front first, then bootstrap robustness (if any),
  # then compromise distance to the ideal point.
  ord <- if (!is.null(boot)) {
    order(res$pareto_front, -res$p_pareto, res$distance)
  } else {
    order(res$pareto_front, res$distance)
  }
  res$rank <- NA_integer_
  res$rank[ord] <- seq_along(ord)
  res <- res[ord, , drop = FALSE]
  rownames(res) <- NULL

  attr(res, "bootstrap") <- boot
  attr(res, "conf") <- conf
  attr(res, "n_boot") <- n_boot
  attr(res, "criteria") <- criteria
  attr(res, "tol") <- tol
  class(res) <- c("enm_evaluation", "data.frame")
  res
}

validate_table <- function(data) {
  if (!is.data.frame(data) || !"obs" %in% names(data)) {
    stop("`data` must be a data.frame with an `obs` column (1/0).",
         call. = FALSE)
  }
  if (!all(data$obs %in% c(0, 1))) {
    stop("`obs` must contain only 0 (absence/background) and 1 (presence).",
         call. = FALSE)
  }
  if (ncol(data) < 2) stop("`data` has no prediction columns.", call. = FALSE)
  if (sum(data$obs == 1) < 2 || sum(data$obs == 0) < 2) {
    stop("At least two presence and two absence points are required.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Stratified bootstrap of evaluation points (internal)
bootstrap_metrics <- function(obs, preds, n_boot, funs, criteria, tol = 0) {
  i1 <- which(obs == 1)
  i0 <- which(obs == 0)
  out <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    idx <- c(i1[sample.int(length(i1), replace = TRUE)],
             i0[sample.int(length(i0), replace = TRUE)])
    o <- obs[idx]
    m <- do.call(rbind, lapply(preds, function(p) {
      compute_metrics(o, p[idx], funs)
    }))
    d <- data.frame(replicate = b, model = names(preds), m, row.names = NULL,
                    check.names = FALSE)
    d$pareto_front <- pareto_fronts(to_minimise(d, criteria),
                                    tol = expand_tol(tol, names(criteria)))
    out[[b]] <- d
  }
  do.call(rbind, out)
}

#' @export
print.enm_evaluation <- function(x, digits = 3, ...) {
  if (!all(c("model", "pareto_front", "n_presence") %in% names(x))) {
    return(print_plain(x, digits, ...))
  }
  nb <- attr(x, "n_boot")
  cat("ENM evaluation: ", nrow(x), " model(s), ", x$n_presence[1],
      " presence and ", x$n_absence[1], " absence/background points\n",
      sep = "")
  if (!is.null(nb) && nb > 0) {
    cat("Bootstrap: ", nb, " resamples, ", 100 * attr(x, "conf"),
        "% percentile intervals\n", sep = "")
  }
  crit <- attr(x, "criteria")
  if (is.null(crit)) crit <- c(AUC = "max", MAE = "min")
  cat("Criteria: ", paste0(names(crit), " (", vapply(crit, direction_label,
                                                     character(1)), ")",
                           collapse = ", "), "\n", sep = "")
  cat("Pareto front (best trade-offs): ",
      paste(x$model[x$pareto_front == 1], collapse = ", "), "\n\n", sep = "")
  crit_cols <- unlist(lapply(names(crit), function(k) {
    c(k, paste0(k, c("_lower", "_upper")))
  }))
  other <- setdiff(names(x), c("rank", "model", crit_cols, "pareto_front",
                               "p_pareto", "distance", "n_presence",
                               "n_absence", grep("_sd$", names(x),
                                                 value = TRUE)))
  cols <- intersect(c("rank", "model", crit_cols, other, "pareto_front",
                      "p_pareto", "distance"), names(x))
  df <- as.data.frame(x)[cols]
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE, ...)
  invisible(x)
}
