# Model / algorithm selection ---------------------------------------------------

#' Select models or algorithms by minimum values or Pareto optimality
#'
#' Combines the two selection rules used in practice:
#'
#' * **`"threshold"`** (default): keep every model reaching conventional
#'   minimum values on all `min_metrics`. This drops only clearly poor models
#'   and keeps several good ones for an ensemble. Typical values: AUC >= 0.8
#'   ("excellent" discrimination; 0.7-0.8 "acceptable", Hosmer & Lemeshow
#'   2000), TSS >= 0.6 ("substantial" agreement on the Landis & Koch 1977
#'   scale), Boyce >= 0.5.
#' * **`"pareto"`**: keep only models on the Pareto front of `criteria`: no
#'   other model is better on one criterion without being worse on another.
#'   Stricter; use `tol` so that differences within sampling noise count as
#'   ties.
#'
#' Both rules are always computed and reported (`pass`, `pareto`), so you can
#' compare them; `rule` only decides which one fills `selected`. If no model
#' passes all minimum values, the Pareto front is used instead, with a
#' warning. `rank_mean` (mean rank over `criteria`, 1 = best) orders the
#' output.
#'
#' Any table with a `model` column and metric columns works: the output of
#' [evaluate_cv()], [evaluate_models()] or [temporal_validation()], or
#' several of them joined with [merge()], e.g. hold-out metrics together with
#' fossil validation.
#'
#' @param x A data.frame with one row per model and a `model` column.
#' @param rule `"threshold"` or `"pareto"` (see Details).
#' @param min_metrics Named vector of minimum values (larger = better).
#'   Metrics that are not columns of `x` are ignored with a message. To set a
#'   maximum for a "lower is better" metric, add a negated column, e.g.
#'   `x$neg_MAE <- -x$MAE` and `min_metrics = c(neg_MAE = -0.3)`.
#' @param criteria Named character vector giving, for each Pareto criterion,
#'   `"max"` (larger is better), `"min"` (smaller is better) or `"zero"`
#'   (closer to 0 is better, e.g. `bias`).
#' @param tol Tolerance for Pareto ties; see [pareto_rank()].
#' @param override Optional character vector of model names to select by
#'   hand. The rules are still computed and reported.
#' @return An object of class `enm_selection` (a data.frame): `x` with
#'   columns `pass`, `pareto`, `rank_mean` and `selected`, sorted with the
#'   selected models first.
#' @export
#' @examples
#' df <- data.frame(model = c("GLM", "GBM", "RF", "SRE"),
#'                  AUC = c(0.86, 0.88, 0.87, 0.71),
#'                  TSS = c(0.62, 0.65, 0.66, 0.38),
#'                  MAE = c(0.30, 0.27, 0.24, 0.41))
#' select_models(df)                                   # minimum values
#' select_models(df, rule = "pareto", tol = 0.01)      # Pareto front
select_models <- function(x, rule = c("threshold", "pareto"),
                          min_metrics = c(AUC = 0.8, TSS = 0.6),
                          criteria = c(AUC = "max", MAE = "min"),
                          tol = 0, override = NULL) {
  rule <- match.arg(rule)
  x <- as.data.frame(x)
  if (!"model" %in% names(x)) {
    stop("`x` needs a `model` column.", call. = FALSE)
  }
  check_criteria(criteria, names(criteria))

  # rule 1: minimum values
  miss <- setdiff(names(min_metrics), names(x))
  if (length(miss)) {
    message("Minimum values ignored (no such column): ",
            paste(miss, collapse = ", "))
  }
  mins <- min_metrics[names(min_metrics) %in% names(x)]
  x$pass <- if (length(mins)) {
    Reduce(`&`, lapply(names(mins), function(k) {
      !is.na(x[[k]]) & x[[k]] >= mins[[k]]
    }))
  } else {
    rep(TRUE, nrow(x))
  }

  # rule 2: Pareto front (missing values count as worst)
  miss_c <- setdiff(names(criteria), names(x))
  if (length(miss_c)) {
    message("Pareto criteria ignored (no such column): ",
            paste(miss_c, collapse = ", "))
  }
  crit <- criteria[names(criteria) %in% names(x)]
  if (!length(crit)) stop("None of the `criteria` are columns of `x`.",
                          call. = FALSE)
  m <- to_minimise(x, crit)
  x$pareto <- pareto_fronts(m, tol = expand_tol(tol, names(crit))) == 1
  x$rank_mean <- rowMeans(apply(m, 2, rank, ties.method = "average"))

  sel <- if (!is.null(override)) {
    unknown <- setdiff(override, x$model)
    if (length(unknown)) {
      stop("Unknown model(s) in `override`: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    x$model %in% override
  } else if (rule == "pareto") {
    x$pareto
  } else {
    x$pass
  }
  if (!any(sel)) {
    warning("No model reaches all minimum values (",
            paste(names(mins), mins, sep = " >= ", collapse = ", "),
            "); the Pareto front is selected instead.", call. = FALSE)
    sel <- x$pareto
  }
  x$selected <- sel
  x <- x[order(!x$selected, x$rank_mean), , drop = FALSE]
  rownames(x) <- NULL
  attr(x, "rule") <- if (!is.null(override)) "manual" else rule
  attr(x, "min_metrics") <- mins
  attr(x, "criteria") <- crit
  attr(x, "tol") <- tol
  class(x) <- c("enm_selection", "data.frame")
  x
}

#' @export
print.enm_selection <- function(x, digits = 3, ...) {
  if (is.null(attr(x, "criteria")) ||
      !all(c("model", "pass", "pareto", "selected") %in% names(x))) {
    return(print_plain(x, digits, ...))
  }
  mins <- attr(x, "min_metrics")
  crit <- attr(x, "criteria")
  cat("Model selection (rule: ", attr(x, "rule"), ")\n", sep = "")
  if (length(mins)) {
    cat("  Minimum values: ", paste(names(mins), mins, sep = " >= ",
                                    collapse = ", "),
        " | pass: ", paste(x$model[x$pass], collapse = ", "), "\n", sep = "")
  }
  cat("  Pareto front (", paste(names(crit), crit, sep = ":", collapse = ", "),
      "; tol ", paste(attr(x, "tol"), collapse = "/"), "): ",
      paste(x$model[x$pareto], collapse = ", "), "\n", sep = "")
  cat("  SELECTED: ", paste(x$model[x$selected], collapse = ", "), "\n\n",
      sep = "")
  df <- as.data.frame(x)
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE, ...)
  invisible(x)
}
