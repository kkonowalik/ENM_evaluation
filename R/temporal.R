# Temporal (e.g. fossil) validation --------------------------------------------

#' Validate hindcasts against dated records (e.g. fossils)
#'
#' Tests whether models projected to past climates predict independent,
#' dated occurrences, e.g. fossil or pollen records. Records have no true
#' absences and each is matched to the climate of its own time, so instead of
#' a classical AUC each record is compared with **background cells of the
#' same time slice**:
#'
#' 1. For every record and every time slice it falls into, the score is the
#'    proportion of background cells of that slice with lower predicted
#'    suitability (ties count half). 0.5 = random, 1 = the record is more
#'    suitable than every background cell.
#' 2. A record that spans several slices (e.g. a 1000-year bin covering ten
#'    100-year climate slices) gets the mean of its slice scores.
#' 3. The temporal AUC of a model is the mean over records. It equals the
#'    probability that a record is more suitable than a random background
#'    cell of its time, which is the definition of AUC.
#'
#' Optionally the score is also given per period (e.g. Holocene,
#' Lateglacial, Glacial) and for groups of periods (`combine`). Older periods
#' are usually what hindcasts are about, while the far more numerous recent
#' records dominate the overall mean. Further options: an "analog" score that
#' ignores record-slices with climate outside the calibration range (`novel`,
#' see [outside_range()]), and sensitivity, the share of records predicted
#' above a present-day threshold (`thresholds`, see [threshold_p10()]).
#'
#' @param sites Long data.frame of predictions at the records: one row per
#'   record x time slice x model (x run), with the columns named by `record`,
#'   `slice`, `model`, `pred` and optionally `run`, `period` and `novel`.
#' @param background Long data.frame of predictions at background cells: one
#'   row per cell x time slice x model (x run), with columns `slice`, `model`,
#'   `pred` (and `run`). Use the same cells in every slice and several hundred
#'   or more cells per slice.
#' @param model,run,record,slice,pred,period,novel Column names. `run` (e.g.
#'   the cross-validation model or background set) is optional: when given,
#'   scores are computed per run and then averaged per model, with an SD.
#'   `period` must be constant within a record; `novel` is logical (TRUE =
#'   non-analog climate).
#' @param thresholds Optional thresholds for sensitivity: a named numeric
#'   vector (names = models) or a data.frame with columns `model`, `run`
#'   (if used) and `threshold`.
#' @param combine Optional named list of period groups, e.g.
#'   `list(deep = c("Lateglacial", "Glacial"))`. Records of these periods are
#'   pooled (so periods with more records weigh more).
#' @param prefix Prefix of the output columns, e.g. `"fossil"` gives
#'   `fossil_AUC`.
#' @return A data.frame with one row per model: `<prefix>_AUC`, per-period and
#'   combined columns `<prefix>_AUC_<period>`, `<prefix>_AUC_analog`,
#'   `<prefix>_sens`, `<prefix>_n_records`, and with `run` the SD of each
#'   metric across runs (`_sd`) and `<prefix>_n_runs`. Per-run results are
#'   in `attr(x, "runs")`.
#' @export
#' @examples
#' set.seed(1)
#' # 3 time slices, 500 background cells each, 40 records
#' bg <- expand.grid(cell = 1:500, slice = 1:3, model = c("good", "random"))
#' bg$pred <- runif(nrow(bg))
#' st <- expand.grid(record = 1:40, model = c("good", "random"))
#' st$slice <- sample(1:3, nrow(st), replace = TRUE)
#' st$pred <- ifelse(st$model == "good", rbeta(nrow(st), 5, 1),
#'                   runif(nrow(st)))
#' st$period <- ifelse(st$record <= 30, "Holocene", "Glacial")
#' temporal_validation(st, bg, period = "period", prefix = "fossil")
temporal_validation <- function(sites, background, model = "model",
                                run = NULL, record = "record",
                                slice = "slice", pred = "pred",
                                period = NULL, novel = NULL,
                                thresholds = NULL, combine = NULL,
                                prefix = "temporal") {
  need_s <- c(model, run, record, slice, pred, period, novel)
  need_b <- c(model, run, slice, pred)
  miss <- c(setdiff(need_s, names(sites)), setdiff(need_b, names(background)))
  if (length(miss)) {
    stop("Columns not found: ", paste(unique(miss), collapse = ", "),
         call. = FALSE)
  }
  run_s <- if (is.null(run)) rep("1", nrow(sites)) else as.character(sites[[run]])
  run_b <- if (is.null(run)) rep("1", nrow(background)) else
    as.character(background[[run]])
  gs <- paste(sites[[model]], run_s, sep = "\r")
  gb <- paste(background[[model]], run_b, sep = "\r")
  bg_by_group <- split(seq_len(nrow(background)), gb)
  if (!is.null(period)) {
    np <- tapply(sites[[period]], sites[[record]],
                 function(v) length(unique(v)))
    if (any(np > 1)) stop("`", period, "` must be constant within a record.",
                          call. = FALSE)
  }
  thr_lookup <- function(mod, rn) {
    if (is.null(thresholds)) return(NA_real_)
    if (is.numeric(thresholds)) {
      v <- thresholds[as.character(mod)]
      return(if (is.na(v)) NA_real_ else unname(v))
    }
    hit <- thresholds$model == mod
    if (!is.null(run) && "run" %in% names(thresholds)) {
      hit <- hit & as.character(thresholds$run) == rn
    }
    if (any(hit)) thresholds$threshold[which(hit)[1]] else NA_real_
  }
  p <- prefix
  empty_slices <- character(0)

  runs <- lapply(split(seq_len(nrow(sites)), gs), function(si) {
    g <- gs[si[1]]
    bi <- bg_by_group[[g]]
    if (is.null(bi)) {
      stop("No background predictions for model/run '",
           sub("\r", " / ", g), "'.", call. = FALSE)
    }
    bvals <- split(background[[pred]][bi], background[[slice]][bi])
    bvals <- lapply(bvals, function(v) sort(v[!is.na(v)]))
    x <- sites[[pred]][si]
    sl <- as.character(sites[[slice]][si])
    score <- rep(NA_real_, length(si))
    for (s in unique(sl)) {
      rows <- which(sl == s)
      v <- bvals[[s]]
      if (is.null(v) || !length(v)) {
        empty_slices <<- c(empty_slices, s)
        next
      }
      score[rows] <- (findInterval(x[rows], v, left.open = TRUE) +
                        findInterval(x[rows], v)) / (2 * length(v))
    }
    score[is.na(x)] <- NA
    rec <- sites[[record]][si]
    rec_score <- tapply(score, rec, mean, na.rm = TRUE)
    rec_score[is.nan(rec_score)] <- NA
    out <- list(model = as.character(sites[[model]][si[1]]),
                run = run_s[si[1]])
    out[[paste0(p, "_AUC")]] <- mean(rec_score, na.rm = TRUE)
    if (!is.null(period)) {
      rec_per <- tapply(as.character(sites[[period]][si]), rec, `[`, 1)
      for (pp in sort(unique(stats::na.omit(rec_per)))) {
        out[[paste0(p, "_AUC_", pp)]] <- mean(rec_score[rec_per == pp],
                                              na.rm = TRUE)
      }
      for (nm in names(combine)) {
        out[[paste0(p, "_AUC_", nm)]] <-
          mean(rec_score[rec_per %in% combine[[nm]]], na.rm = TRUE)
      }
    }
    if (!is.null(novel)) {
      an <- !as.logical(sites[[novel]][si])
      ra <- tapply(score[an], rec[an], mean, na.rm = TRUE)
      out[[paste0(p, "_AUC_analog")]] <- mean(ra, na.rm = TRUE)
    }
    if (!is.null(thresholds)) {
      thr <- thr_lookup(out$model, out$run)
      rec_pred <- tapply(x, rec, mean, na.rm = TRUE)
      out[[paste0(p, "_sens")]] <- if (is.na(thr)) NA_real_ else
        mean(rec_pred >= thr, na.rm = TRUE)
    }
    out[[paste0(p, "_n_records")]] <- sum(!is.na(rec_score))
    out
  })
  if (length(empty_slices)) {
    warning("No background cells for slice(s): ",
            paste(unique(empty_slices), collapse = ", "),
            ". Records in them were ignored.", call. = FALSE)
  }
  cols <- unique(unlist(lapply(runs, names)))
  runs <- do.call(rbind, lapply(runs, function(r) {
    r[setdiff(cols, names(r))] <- NA
    as.data.frame(r[cols], check.names = FALSE)
  }))
  rownames(runs) <- NULL
  nrec <- paste0(p, "_n_records")
  mets <- setdiff(cols, c("model", "run", nrec))
  runs[mets] <- lapply(runs[mets], function(v) {
    v[is.nan(v)] <- NA
    v
  })

  if (is.null(run)) {
    out <- runs[setdiff(names(runs), "run")]
  } else {
    out <- do.call(rbind, lapply(split(runs, factor(runs$model,
                                                    unique(runs$model))),
                                 function(r) {
      mn <- vapply(r[mets], mean, numeric(1), na.rm = TRUE)
      sdv <- vapply(r[mets], stats::sd, numeric(1), na.rm = TRUE)
      d <- data.frame(model = r$model[1], t(mn),
                      t(stats::setNames(sdv, paste0(mets, "_sd"))),
                      check.names = FALSE)
      d[[nrec]] <- max(r[[nrec]])
      d[[paste0(p, "_n_runs")]] <- nrow(r)
      d
    }))
  }
  rownames(out) <- NULL
  attr(out, "runs") <- runs
  out
}

#' Assign ages to named periods
#'
#' @param age Numeric ages (e.g. calibrated years BP).
#' @param upper Named vector of the upper (oldest) age limit of each period,
#'   from youngest to oldest. The default follows the INTIMATE event
#'   stratigraphy: Holocene up to 11,700, Lateglacial up to 14,700, Glacial
#'   older.
#' @return Character vector of period names.
#' @export
#' @examples
#' age_to_period(c(500, 11700, 12000, 20000))
age_to_period <- function(age, upper = c(Holocene = 11700,
                                         Lateglacial = 14700,
                                         Glacial = Inf)) {
  if (is.null(names(upper)) || is.unsorted(upper)) {
    stop("`upper` must be a named, increasing vector of age limits.",
         call. = FALSE)
  }
  names(upper)[findInterval(age, c(-Inf, upper), left.open = TRUE)]
}

#' Past climate from present climate and a climate-model anomaly
#'
#' The "delta" or anomaly method: the change simulated by a (coarse) climate
#' model between a past period and its own present is added to a
#' high-resolution present-day climatology (e.g. CHELSA or WorldClim). This
#' removes much of the climate model's bias. Temperature-like variables use
#' the difference (`current + past - present`); precipitation amounts use the
#' ratio (`current * (past + 1) / (present + 1)`, the +1 avoids division by
#' zero in dry cells), limited to `ratio_limits` to avoid extreme values.
#'
#' Works on matrices or data.frames of values at points (fast; no maps
#' needed) and on [terra::SpatRaster]s with layers of the same names.
#'
#' @param current Present-day high-resolution values.
#' @param past Climate model values for the past period.
#' @param present Climate model values for the present (the reference period
#'   of the same model).
#' @param ratio_vars Names of variables that use the ratio (e.g. the
#'   precipitation variables). Tip: `grep("^bio1[2-9]$", names(x),
#'   value = TRUE)` selects bioclim precipitation variables 12-19.
#' @param ratio_limits Lower and upper limits of the ratio.
#' @return An object like `current` with past values.
#' @export
#' @examples
#' cur  <- data.frame(bio1 = c(5, 10), bio12 = c(600, 900))
#' past <- data.frame(bio1 = c(-3, 2), bio12 = c(300, 500))
#' pres <- data.frame(bio1 = c(4, 9), bio12 = c(650, 850))
#' anomaly_climate(cur, past, pres, ratio_vars = "bio12")
anomaly_climate <- function(current, past, present, ratio_vars = NULL,
                            ratio_limits = c(0.1, 10)) {
  vars <- names(current)
  if (!identical(sort(vars), sort(names(past))) ||
      !identical(sort(vars), sort(names(present)))) {
    stop("`current`, `past` and `present` must have the same variable names.",
         call. = FALSE)
  }
  if (inherits(current, "SpatRaster")) {
    layers <- lapply(vars, function(v) {
      if (v %in% ratio_vars) {
        r <- (past[[v]] + 1) / (present[[v]] + 1)
        current[[v]] * terra::clamp(r, ratio_limits[1], ratio_limits[2],
                                    values = TRUE)
      } else {
        current[[v]] + past[[v]] - present[[v]]
      }
    })
    out <- terra::rast(layers)
    names(out) <- vars
    return(out)
  }
  out <- current
  for (v in vars) {
    out[, v] <- if (v %in% ratio_vars) {
      r <- (past[, v] + 1) / (present[, v] + 1)
      current[, v] * pmin(pmax(r, ratio_limits[1]), ratio_limits[2])
    } else {
      current[, v] + past[, v] - present[, v]
    }
  }
  out
}

#' Flag values outside the calibration range (non-analog climate)
#'
#' Models extrapolate when projected to climates absent from the calibration
#' data. This flags rows where at least one variable lies outside the range
#' found in `reference` (a simple univariate check, like the "clamping mask"
#' of Maxent and biomod2).
#'
#' @param x Matrix or data.frame of values to check.
#' @param reference Matrix or data.frame of calibration values (columns
#'   matched by name).
#' @return Logical vector, `TRUE` = outside the range of at least one
#'   variable.
#' @export
#' @examples
#' ref <- data.frame(t = c(0, 10), p = c(100, 900))
#' outside_range(data.frame(t = c(5, -2), p = c(500, 500)), ref)
outside_range <- function(x, reference) {
  vars <- intersect(colnames(x), colnames(reference))
  if (!length(vars)) stop("No common columns in `x` and `reference`.",
                          call. = FALSE)
  x <- as.matrix(x[, vars, drop = FALSE])
  rng <- apply(as.matrix(reference[, vars, drop = FALSE]), 2, range,
               na.rm = TRUE)
  rowSums(sweep(x, 2, rng[1, ], `<`) | sweep(x, 2, rng[2, ], `>`),
          na.rm = TRUE) > 0
}
