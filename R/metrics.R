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

# Compute all metrics for one model (internal)
metric_row <- function(obs, pred) {
  ok <- !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  c(AUC = auc_score(obs, pred),
    MAE = mae_score(obs, pred),
    bias = bias_score(obs, pred),
    n_presence = sum(obs == 1),
    n_absence = sum(obs == 0))
}
