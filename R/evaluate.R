# Model evaluation and bootstrap ------------------------------------------------

#' Evaluate one or many models with AUC, MAE and bias
#'
#' The main entry point of the package. Extracts predictions at presence and
#' absence/background points, computes the metrics recommended by Konowalik &
#' Nosol (2021), and ranks models by Pareto optimality on AUC (maximised) and
#' MAE (minimised).
#'
#' When `n_boot > 0`, evaluation points are resampled with replacement
#' (separately for presences and absences, so their numbers stay fixed) and all
#' models are re-evaluated on each resample. This gives:
#'   * percentile confidence intervals for AUC and MAE, and
#'   * `p_pareto`: the proportion of resamples in which the model lies on the
#'     Pareto front. It shows how robust the selection is: a model with
#'     `p_pareto = 0.95` is almost always among the best trade-offs, while
#'     `0.20` means it reached the front partly by chance.
#'
#' @inheritParams extract_predictions
#' @param n_boot Number of bootstrap resamples (0 = no bootstrap). 500 to 1000
#'   is a good choice for final results.
#' @param conf Confidence level for bootstrap intervals (default 0.95).
#' @param weights Relative importance of AUC and MAE in the compromise score;
#'   see [pareto_rank()].
#' @param seed Optional random seed for reproducible bootstraps.
#' @param data Optionally, a table returned by [extract_predictions()]. When
#'   supplied, `models`, `presence` and `absence` are ignored. This is also
#'   useful when predictions were extracted with other software: supply a
#'   data.frame with a 0/1 column `obs` and one numeric column per model.
#'
#' @return An object of class `enm_evaluation` (a data.frame) with one row per
#'   model and columns `model`, `AUC`, `MAE`, `bias`, `n_presence`,
#'   `n_absence`, `pareto_front` (1 = non-dominated), `distance` (normalised
#'   distance to the ideal point, lower is better) and `rank`. With
#'   bootstrapping there are also `AUC_lower`, `AUC_upper`, `MAE_lower`,
#'   `MAE_upper`, `AUC_sd`, `MAE_sd` and `p_pareto`. The bootstrap replicates
#'   are stored in `attr(x, "bootstrap")`.
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
evaluate_models <- function(models = NULL, presence = NULL, absence = NULL,
                            n_boot = 0, conf = 0.95,
                            weights = c(AUC = 1, MAE = 1),
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

  res <- do.call(rbind, lapply(preds, function(p) metric_row(obs, p)))
  res <- data.frame(model = names(preds), res, row.names = NULL,
                    check.names = FALSE)
  res$n_presence <- as.integer(res$n_presence)
  res$n_absence <- as.integer(res$n_absence)

  ranked <- pareto_rank(res, maximize = "AUC", minimize = "MAE",
                        weights = weights)
  idx <- match(res$model, ranked$model)  # pareto_rank() returns sorted rows
  res$pareto_front <- ranked$pareto_front[idx]
  res$distance <- ranked$distance[idx]

  boot <- NULL
  if (n_boot > 0) {
    if (!is.null(seed)) {
      old <- if (exists(".Random.seed", envir = globalenv())) {
        get(".Random.seed", envir = globalenv())
      }
      on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
              else assign(".Random.seed", old, envir = globalenv()),
              add = TRUE)
      set.seed(seed)
    }
    boot <- bootstrap_metrics(obs, preds, n_boot)
    alpha <- (1 - conf) / 2
    q <- function(x) stats::quantile(x, c(alpha, 1 - alpha), na.rm = TRUE,
                                     names = FALSE)
    stats_by_model <- lapply(split(boot, boot$model), function(b) {
      data.frame(model = b$model[1],
                 AUC_lower = q(b$AUC)[1], AUC_upper = q(b$AUC)[2],
                 MAE_lower = q(b$MAE)[1], MAE_upper = q(b$MAE)[2],
                 AUC_sd = stats::sd(b$AUC, na.rm = TRUE),
                 MAE_sd = stats::sd(b$MAE, na.rm = TRUE),
                 p_pareto = mean(b$pareto_front == 1))
    })
    bs <- do.call(rbind, stats_by_model)
    res <- cbind(res, bs[match(res$model, bs$model), -1])
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
bootstrap_metrics <- function(obs, preds, n_boot) {
  i1 <- which(obs == 1)
  i0 <- which(obs == 0)
  out <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    idx <- c(i1[sample.int(length(i1), replace = TRUE)],
             i0[sample.int(length(i0), replace = TRUE)])
    o <- obs[idx]
    m <- do.call(rbind, lapply(preds, function(p) {
      p <- p[idx]
      ok <- !is.na(p)
      c(AUC = auc_score(o[ok], p[ok]), MAE = mae_score(o[ok], p[ok]))
    }))
    d <- data.frame(replicate = b, model = names(preds), m, row.names = NULL)
    d$pareto_front <- pareto_fronts(cbind(-d$AUC, d$MAE))
    out[[b]] <- d
  }
  do.call(rbind, out)
}

#' @export
print.enm_evaluation <- function(x, digits = 3, ...) {
  nb <- attr(x, "n_boot")
  cat("ENM evaluation: ", nrow(x), " model(s), ", x$n_presence[1],
      " presence and ", x$n_absence[1], " absence/background points\n",
      sep = "")
  if (!is.null(nb) && nb > 0) {
    cat("Bootstrap: ", nb, " resamples, ", 100 * attr(x, "conf"),
        "% percentile intervals\n", sep = "")
  }
  cat("Pareto front (best AUC/MAE trade-offs): ",
      paste(x$model[x$pareto_front == 1], collapse = ", "), "\n\n", sep = "")
  cols <- intersect(c("rank", "model", "AUC", "AUC_lower", "AUC_upper", "MAE",
                      "MAE_lower", "MAE_upper", "bias", "pareto_front",
                      "p_pareto", "distance"), names(x))
  df <- as.data.frame(x)[cols]
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE, ...)
  invisible(x)
}
