# Pareto (non-dominated) ranking ------------------------------------------------

#' Pareto ranking and compromise score for multi-criteria model selection
#'
#' A model **dominates** another if it is at least as good on every criterion
#' and strictly better on at least one. Models that no other model dominates
#' form the *Pareto front* (front 1). These are the defensible choices: each
#' is a different trade-off between the criteria. Removing front 1 and
#' repeating gives front 2, and so on (non-dominated sorting, as in NSGA-II).
#'
#' Within and across fronts, models are also scored by their Euclidean
#' `distance` to the ideal point, the hypothetical model that is best on every
#' criterion. Criteria are first min-max normalised to the range 0 to 1, so
#' AUC and MAE contribute on a comparable scale. `weights` let you make one
#' criterion count more than another. This is the "compromise programming"
#' (TOPSIS-like) approach to choose a single model from the front.
#'
#' The function is generic: any number of criteria (e.g. AUC, MAE, TSS,
#' Boyce index, AICc, omission rate) can be used.
#'
#' **Tolerance.** Differences of a few thousandths in AUC are usually within
#' sampling noise (compare them with the SD of AUC across folds). With
#' `tol > 0`, a model is dominated only if another model is better by more
#' than `tol` on at least one criterion *and* not worse by more than `tol` on
#' any, so near-ties no longer remove models from the front.
#'
#' @param x A data.frame with one row per model.
#' @param maximize Character vector of columns where larger is better.
#' @param minimize Character vector of columns where smaller is better.
#' @param zero Character vector of columns where values closer to 0 are
#'   better (e.g. `bias`).
#' @param weights Optional named numeric vector of criterion weights (names
#'   must match `maximize`/`minimize`); missing names get weight 1.
#' @param tol Differences smaller than this count as ties: one number for all
#'   criteria or a named vector, e.g. `c(AUC = 0.01, MAE = 0.005)`.
#'   Default 0 (strict Pareto dominance).
#' @return `x` with added columns `pareto_front` (1 = non-dominated) and
#'   `distance` (0 = ideal), sorted by front then distance.
#' @export
#' @examples
#' df <- data.frame(model = c("A", "B", "C", "D"),
#'                  AUC = c(0.90, 0.85, 0.80, 0.88),
#'                  MAE = c(0.30, 0.20, 0.25, 0.32))
#' pareto_rank(df, maximize = "AUC", minimize = "MAE")
pareto_rank <- function(x, maximize = "AUC", minimize = "MAE", zero = NULL,
                        weights = NULL, tol = 0) {
  crit <- c(maximize, minimize, zero)
  if (length(crit) < 1) stop("Specify at least one criterion.", call. = FALSE)
  miss <- setdiff(crit, names(x))
  if (length(miss)) {
    stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  criteria <- stats::setNames(
    rep(c("max", "min", "zero"),
        c(length(maximize), length(minimize), length(zero))), crit)
  pc <- pareto_core(to_minimise(x, criteria), weights = weights, tol = tol)
  x$pareto_front <- pc$front
  x$distance <- pc$distance

  x[order(x$pareto_front, x$distance), , drop = FALSE]
}

#' Non-dominated sorting (internal)
#'
#' @param m Numeric matrix, rows = solutions, columns = criteria to MINIMISE.
#' @return Integer vector with the front number of each row (1 = best).
#'   Rows with missing values get `NA`.
#' @keywords internal
#' @noRd
pareto_fronts <- function(m, tol = 0) {
  m <- as.matrix(m)
  tol <- rep_len(tol, ncol(m))
  n <- nrow(m)
  front <- rep(NA_integer_, n)
  remaining <- which(stats::complete.cases(m))
  k <- 1L
  while (length(remaining) > 0) {
    sub <- m[remaining, , drop = FALSE]
    dominated <- vapply(seq_along(remaining), function(i) {
      # j dominates i: not worse by more than tol anywhere, better by more
      # than tol somewhere (tol = 0: classic Pareto dominance)
      better_eq <- sweep(sub, 2, sub[i, ] + tol, `<=`)
      strictly <- sweep(sub, 2, sub[i, ] - tol, `<`)
      any(rowSums(better_eq) == ncol(sub) & rowSums(strictly) > 0)
    }, logical(1))
    front[remaining[!dominated]] <- k
    remaining <- remaining[dominated]
    k <- k + 1L
  }
  front
}

# Tolerance as one value per criterion, in the order of `crit` (internal)
expand_tol <- function(tol, crit) {
  if (is.null(names(tol))) return(rep_len(tol, length(crit)))
  out <- stats::setNames(rep(0, length(crit)), crit)
  hit <- intersect(names(tol), crit)
  out[hit] <- tol[hit]
  unname(out)
}
