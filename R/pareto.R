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
#' @param x A data.frame with one row per model.
#' @param maximize Character vector of columns where larger is better.
#' @param minimize Character vector of columns where smaller is better.
#' @param weights Optional named numeric vector of criterion weights (names
#'   must match `maximize`/`minimize`); missing names get weight 1.
#' @return `x` with added columns `pareto_front` (1 = non-dominated) and
#'   `distance` (0 = ideal), sorted by front then distance.
#' @export
#' @examples
#' df <- data.frame(model = c("A", "B", "C", "D"),
#'                  AUC = c(0.90, 0.85, 0.80, 0.88),
#'                  MAE = c(0.30, 0.20, 0.25, 0.32))
#' pareto_rank(df, maximize = "AUC", minimize = "MAE")
pareto_rank <- function(x, maximize = "AUC", minimize = "MAE",
                        weights = NULL) {
  crit <- c(maximize, minimize)
  if (length(crit) < 1) stop("Specify at least one criterion.", call. = FALSE)
  miss <- setdiff(crit, names(x))
  if (length(miss)) {
    stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  # Convert everything to "smaller is better"
  m <- as.matrix(x[crit])
  storage.mode(m) <- "double"
  m[, seq_along(maximize)] <- -m[, seq_along(maximize)]

  x$pareto_front <- pareto_fronts(m)

  w <- stats::setNames(rep(1, length(crit)), crit)
  if (!is.null(weights)) {
    w[intersect(names(weights), crit)] <- weights[intersect(names(weights),
                                                           crit)]
  }
  norm <- apply(m, 2, function(v) {
    r <- range(v, na.rm = TRUE)
    if (diff(r) == 0) rep(0, length(v)) else (v - r[1]) / diff(r)
  })
  norm <- matrix(norm, nrow = nrow(m))
  x$distance <- sqrt(rowSums(sweep(norm^2, 2, w, `*`))) / sqrt(sum(w))

  x[order(x$pareto_front, x$distance), , drop = FALSE]
}

#' Non-dominated sorting (internal)
#'
#' @param m Numeric matrix, rows = solutions, columns = criteria to MINIMISE.
#' @return Integer vector with the front number of each row (1 = best).
#'   Rows with missing values get `NA`.
#' @keywords internal
#' @noRd
pareto_fronts <- function(m) {
  m <- as.matrix(m)
  n <- nrow(m)
  front <- rep(NA_integer_, n)
  remaining <- which(stats::complete.cases(m))
  k <- 1L
  while (length(remaining) > 0) {
    sub <- m[remaining, , drop = FALSE]
    dominated <- vapply(seq_along(remaining), function(i) {
      better_eq <- sweep(sub, 2, sub[i, ], `<=`)
      strictly <- sweep(sub, 2, sub[i, ], `<`)
      any(rowSums(better_eq) == ncol(sub) & rowSums(strictly) > 0)
    }, logical(1))
    front[remaining[!dominated]] <- k
    remaining <- remaining[dominated]
    k <- k + 1L
  }
  front
}
