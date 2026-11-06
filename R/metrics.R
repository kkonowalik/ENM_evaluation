# Evaluation metrics -----------------------------------------------------------
# Implemented in base R (no dependency on the 'Metrics' package) and numerically
# identical to Metrics::auc(), Metrics::mae() and Metrics::bias().

#' Area under the ROC curve
#'
#' Computed with the Mann-Whitney U statistic, so ties are handled correctly
#' and no threshold is needed.
#'
#' @param obs Numeric vector of observations (1 = presence, 0 = absence).
#' @param pred Numeric vector of predicted suitability.
#' @return A single number between 0 and 1.
#' @export
#' @examples
#' auc_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
auc_score <- function(obs, pred) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  n1 <- sum(obs == 1)
  n0 <- sum(obs == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(pred)
  (sum(r[obs == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Mean absolute error
#'
#' Mean absolute difference between the observations (0/1) and the predicted
#' suitability. Lower is better.
#'
#' @inheritParams auc_score
#' @export
#' @examples
#' mae_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
mae_score <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

#' Bias
#'
#' Mean of `obs - pred`. Positive values mean the model under-predicts
#' suitability on average and negative values mean it over-predicts. Values
#' closer to 0 are better.
#'
#' @inheritParams auc_score
#' @export
#' @examples
#' bias_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
bias_score <- function(obs, pred) mean(obs - pred, na.rm = TRUE)

#' True skill statistic (TSS)
#'
#' TSS = sensitivity + specificity - 1, ranging from -1 to 1 (0 = random).
#' By default the threshold that maximises TSS is used (the same as the
#' maximum of Youden's J). This is the usual way TSS is reported, but it is
#' slightly optimistic on hold-out data. Supply a `threshold` fixed beforehand
#' (e.g. from calibration data, see [threshold_p10()]) for a strict hold-out
#' value.
#'
#' @inheritParams auc_score
#' @param threshold Optional fixed threshold; predictions `>= threshold` count
#'   as presences. `NULL` (default) uses the TSS-maximising threshold.
#' @return A single number. The threshold used is stored in
#'   `attr(, "threshold")`.
#' @export
#' @examples
#' tss_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
#' tss_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1), threshold = 0.5)
tss_score <- function(obs, pred, threshold = NULL) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  n1 <- sum(obs == 1)
  n0 <- sum(obs == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  if (!is.null(threshold)) {
    out <- mean(pred[obs == 1] >= threshold) +
      mean(pred[obs == 0] < threshold) - 1
    attr(out, "threshold") <- threshold
    return(out)
  }
  o <- order(pred, decreasing = TRUE)
  p <- pred[o]
  tp <- cumsum(obs[o] == 1)
  fp <- cumsum(obs[o] == 0)
  last <- c(p[-1] != p[-length(p)], TRUE)   # last index of each tied value
  tss <- tp[last] / n1 - fp[last] / n0
  best <- which.max(tss)
  out <- max(0, tss[best])
  attr(out, "threshold") <- if (tss[best] > 0) p[last][best] else Inf
  out
}

#' Continuous Boyce index
#'
#' Measures how well predictions match the *frequency* of presences, which
#' is especially suited to presence-only data (Hirzel et al. 2006). A moving
#' window slides along the range of predicted values; in each window the
#' ratio of the share of presences to the share of background points
#' (predicted/expected ratio) is computed. The index is the Spearman
#' correlation between this ratio and the window position: values close to 1
#' mean presences become steadily more frequent as suitability increases,
#' 0 means no better than random and negative values mean a wrong model. The
#' implementation reproduces `ecospat::ecospat.boyce()` (moving window,
#' Spearman correlation), which biomod2 also uses.
#'
#' Use background points that represent the whole study area (or all cells);
#' absences sampled away from presences make the index unreliable.
#'
#' @param presence Predicted values at presence points.
#' @param background Predicted values at background (or absence) points, or
#'   at all cells of the study area.
#' @param n_windows Number of window positions (default 100).
#' @param window_width Window width as a proportion of the range of the
#'   background predictions (default 0.1).
#' @return A single number between -1 and 1, or `NA` when it cannot be
#'   computed (e.g. too few distinct window ratios).
#' @references Hirzel, A.H., Le Lay, G., Helfer, V., Randin, C., Guisan, A.
#'   (2006) Evaluating the ability of habitat suitability models to predict
#'   species presences. *Ecological Modelling* 199, 142-152.
#' @export
#' @examples
#' set.seed(1)
#' bg <- runif(1000)
#' pres <- rbeta(100, 4, 1)          # presences mostly at high suitability
#' boyce_index(pres, bg)
boyce_index <- function(presence, background, n_windows = 100,
                        window_width = 0.1) {
  presence <- presence[!is.na(presence)]
  background <- background[!is.na(background)]
  if (length(presence) < 2 || length(background) < 2) return(NA_real_)
  lo <- min(background, presence)
  hi <- max(background, presence)
  w <- (max(background) - min(background)) * window_width
  if (hi == lo || w <= 0 || hi - lo - w <= 0) return(NA_real_)
  start <- seq(from = lo, to = hi - w, by = (hi - lo - w) / n_windows)
  # reproduces ecospat exactly: the last window start is shifted by +1
  # (an ecospat "trick" for closed intervals that in effect drops it)
  start[length(start)] <- start[length(start)] + 1
  start_end <- start + w
  ratio <- vapply(seq_along(start), function(i) {
    p <- mean(presence >= start[i] & presence <= start_end[i])
    e <- mean(background >= start[i] & background <= start_end[i])
    round(p / e, 10)          # p > 0, e = 0 gives Inf (kept); 0/0 = NaN
  }, numeric(1))
  keep <- which(!is.nan(ratio))
  ratio <- ratio[keep]
  start <- start[keep]
  if (length(ratio) < 2) return(NA_real_)
  # drop successive duplicated ratios (flat stretches), exactly as ecospat
  r <- which(ratio != c(ratio[-1], TRUE))
  if (length(r) < 2) return(NA_real_)
  out <- suppressWarnings(stats::cor(ratio[r], start[r], method = "spearman"))
  if (is.na(out)) NA_real_ else out
}

#' Threshold at a given percentile of presence predictions
#'
#' The "10th percentile training presence" threshold: the suitability value
#' exceeded by 90% of presence points. Usually computed on the calibration
#' (training) presences and then applied to independent data, e.g. to compute
#' sensitivity of a hindcast in [temporal_validation()].
#'
#' @param pred_presence Predicted values at presence points.
#' @param p Proportion of presences allowed below the threshold (0.1 = 10%).
#' @export
#' @examples
#' threshold_p10(c(0.2, 0.5, 0.6, 0.8, 0.9))
threshold_p10 <- function(pred_presence, p = 0.1) {
  unname(stats::quantile(pred_presence, p, na.rm = TRUE))
}

# Compute all metrics for one model (internal)
metric_row <- function(obs, pred) {
  ok <- !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  c(AUC = auc_score(obs, pred),
    MAE = mae_score(obs, pred),
    bias = bias_score(obs, pred),
    TSS = tss_score(obs, pred)[1],
    Boyce = boyce_index(pred[obs == 1], pred[obs == 0]),
    n_presence = sum(obs == 1),
    n_absence = sum(obs == 0))
}
