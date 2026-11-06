# Metrics registry and selection criteria (internal) ----------------------------

# Built-in metrics: each is function(obs, pred) -> single number
builtin_metrics <- function() {
  list(
    AUC = auc_score,
    MAE = mae_score,
    bias = bias_score,
    TSS = function(obs, pred) as.numeric(tss_score(obs, pred)),
    Boyce = function(obs, pred) boyce_index(pred[obs == 1], pred[obs == 0])
  )
}

# Which direction is "better" for the built-in metrics
default_directions <- c(AUC = "max", MAE = "min", bias = "zero", TSS = "max",
                        Boyce = "max")

# Validate user-supplied metric functions
check_extra_metrics <- function(metrics) {
  if (is.null(metrics)) return(list())
  if (is.function(metrics)) {
    stop("`metrics` must be a *named list* of functions, e.g. ",
         "list(Brier = function(obs, pred) mean((obs - pred)^2)).",
         call. = FALSE)
  }
  if (!is.list(metrics) || is.null(names(metrics)) ||
      any(names(metrics) == "") ||
      !all(vapply(metrics, is.function, logical(1)))) {
    stop("`metrics` must be a named list of functions(obs, pred).",
         call. = FALSE)
  }
  clash <- intersect(names(metrics), c(names(builtin_metrics()), "model",
                                       "obs", "pred", "n_presence",
                                       "n_absence"))
  if (length(clash)) {
    stop("Metric name(s) already used by the package: ",
         paste(clash, collapse = ", "), ". Please rename them.",
         call. = FALSE)
  }
  metrics
}

# Compute a set of metrics for one prediction vector
compute_metrics <- function(obs, pred, funs) {
  ok <- !is.na(pred) & !is.na(obs)
  obs <- obs[ok]
  pred <- pred[ok]
  vapply(names(funs), function(nm) {
    v <- tryCatch(funs[[nm]](obs, pred), error = function(e) {
      stop("Metric '", nm, "' failed: ", conditionMessage(e), call. = FALSE)
    })
    if (length(v) != 1 || !is.numeric(v)) {
      stop("Metric '", nm, "' must return a single number.", call. = FALSE)
    }
    as.numeric(v)
  }, numeric(1))
}

# Validate selection criteria against the available columns
check_criteria <- function(criteria, available) {
  if (is.null(names(criteria)) || any(names(criteria) == "") ||
      !all(criteria %in% c("max", "min", "zero"))) {
    stop("`criteria` must be a named vector with values \"max\" (higher is ",
         "better), \"min\" (lower is better) or \"zero\" (closer to 0 is ",
         "better), e.g. c(AUC = \"max\", MAE = \"min\").", call. = FALSE)
  }
  miss <- setdiff(names(criteria), available)
  if (length(miss)) {
    stop("Unknown metric(s) in `criteria`: ", paste(miss, collapse = ", "),
         ". Available: ", paste(available, collapse = ", "),
         ". Add your own with `metrics = list(name = function(obs, pred) ...)`.",
         call. = FALSE)
  }
  criteria
}

# Matrix of criteria transformed so that smaller is always better; missing
# values count as worst
to_minimise <- function(x, criteria) {
  m <- vapply(names(criteria), function(k) {
    v <- as.numeric(x[[k]])
    v <- switch(criteria[[k]], max = -v, min = v, zero = abs(v))
    v[is.na(v)] <- Inf
    v
  }, numeric(nrow(x)))
  matrix(m, nrow = nrow(x), dimnames = list(NULL, names(criteria)))
}

# Pareto front + compromise distance on a "smaller is better" matrix
pareto_core <- function(m, weights = NULL, tol = 0) {
  crit <- colnames(m)
  front <- pareto_fronts(m, tol = expand_tol(tol, crit))
  w <- stats::setNames(rep(1, length(crit)), crit)
  if (!is.null(weights)) {
    hit <- intersect(names(weights), crit)
    w[hit] <- weights[hit]
  }
  norm <- apply(m, 2, function(v) {
    fin <- is.finite(v)
    if (!any(fin)) return(rep(1, length(v)))
    r <- range(v[fin])
    out <- if (diff(r) == 0) rep(0, length(v)) else (v - r[1]) / diff(r)
    out[!fin] <- 1
    out
  })
  norm <- matrix(norm, nrow = nrow(m))
  dist <- sqrt(rowSums(sweep(norm^2, 2, w, `*`))) / sqrt(sum(w))
  list(front = front, distance = dist)
}

# Human-readable direction of a criterion
direction_label <- function(dir) {
  switch(dir, max = "higher is better", min = "lower is better",
         zero = "closer to 0 is better", "")
}
