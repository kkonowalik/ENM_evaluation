# Reading and harmonising inputs ------------------------------------------------

#' Read model predictions as a SpatRaster
#'
#' Accepts the many ways a user may hold model output and returns a single
#' multi-layer [terra::SpatRaster] with one named layer per model.
#'
#' @param models One of:
#'   * a [terra::SpatRaster] (each layer is treated as one model),
#'   * a character vector of file paths (e.g. GeoTIFF / ASCII grids exported
#'     from Maxent, biomod2, ENMeval, QGIS, ...),
#'   * a (named) list mixing `SpatRaster` objects and file paths,
#'   * a legacy `raster::RasterLayer`/`RasterStack` (converted with
#'     [terra::rast()]).
#' @return A named list of single-layer `SpatRaster`s.
#' @keywords internal
#' @noRd
as_model_list <- function(models) {
  if (inherits(models, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    models <- terra::rast(models)
  }
  if (inherits(models, "SpatRaster")) {
    out <- lapply(seq_len(terra::nlyr(models)), function(i) models[[i]])
    names(out) <- names(models)
  } else if (is.character(models)) {
    nms <- names(models)
    if (is.null(nms)) nms <- sub("\\.[^.]*$", "", basename(models))
    out <- lapply(models, terra::rast)
    out <- lapply(out, function(r) {
      if (terra::nlyr(r) > 1) {
        warning("A file has several layers; only the first one is used.",
                call. = FALSE)
        r <- r[[1]]
      }
      r
    })
    names(out) <- nms
  } else if (is.list(models)) {
    out <- lapply(models, function(m) {
      if (is.character(m)) m <- terra::rast(m)
      if (inherits(m, c("RasterLayer", "RasterStack", "RasterBrick"))) {
        m <- terra::rast(m)
      }
      if (!inherits(m, "SpatRaster")) {
        stop("Every element of `models` must be a SpatRaster or a file path.",
             call. = FALSE)
      }
      m[[1]]
    })
    nms <- names(models)
    if (is.null(nms) || any(nms == "")) {
      auto <- vapply(out, function(r) names(r)[1], character(1))
      nms <- if (is.null(nms)) auto else ifelse(nms == "", auto, nms)
    }
    names(out) <- nms
  } else {
    stop("`models` must be a SpatRaster, file path(s) or a list of them.",
         call. = FALSE)
  }
  if (anyDuplicated(names(out))) names(out) <- make.unique(names(out))
  out
}

#' Read occurrence points as a SpatVector
#'
#' @param x A data.frame / matrix with coordinates, a path to a CSV file,
#'   an `sf` object or a [terra::SpatVector] of points.
#' @param coords Optional character vector of length 2 giving the names of
#'   the x (longitude) and y (latitude) columns. If `NULL`, common names such
#'   as `lon/lat`, `longitude/latitude`, `x/y` or `decimalLongitude/
#'   decimalLatitude` (GBIF) are detected automatically.
#' @param crs Coordinate reference system of tabular coordinates. Defaults to
#'   WGS84 longitude/latitude (`"EPSG:4326"`).
#' @keywords internal
#' @noRd
as_points <- function(x, coords = NULL, crs = "EPSG:4326") {
  if (inherits(x, "SpatVector")) return(x)
  if (inherits(x, "sf")) return(terra::vect(x))
  if (inherits(x, "SpatialPointsDataFrame") || inherits(x, "SpatialPoints")) {
    return(terra::vect(x))
  }
  if (is.character(x) && length(x) == 1) {
    x <- utils::read.csv(x, check.names = FALSE)
  }
  if (is.matrix(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) {
    stop("Points must be a data.frame, matrix, CSV path, sf or SpatVector.",
         call. = FALSE)
  }
  if (is.null(coords)) coords <- guess_coords(names(x))
  if (!all(coords %in% names(x))) {
    stop("Coordinate columns not found: ",
         paste(setdiff(coords, names(x)), collapse = ", "), call. = FALSE)
  }
  xy <- data.frame(x = as.numeric(x[[coords[1]]]),
                   y = as.numeric(x[[coords[2]]]))
  xy <- xy[stats::complete.cases(xy), , drop = FALSE]
  terra::vect(xy, geom = c("x", "y"), crs = crs)
}

guess_coords <- function(nms) {
  low <- tolower(nms)
  candidates <- list(
    c("decimallongitude", "decimallatitude"),
    c("longitude", "latitude"),
    c("lon", "lat"),
    c("long", "lat"),
    c("lng", "lat"),
    c("x", "y")
  )
  for (cand in candidates) {
    if (all(cand %in% low)) return(nms[match(cand, low)])
  }
  stop("Could not detect coordinate columns. Please supply ",
       "`coords = c(\"<x/longitude column>\", \"<y/latitude column>\")`.\n",
       "Available columns: ", paste(nms, collapse = ", "), call. = FALSE)
}

#' Extract model predictions at presence and absence/background points
#'
#' Builds the observation/prediction table that all evaluation metrics are
#' computed from. Points are re-projected to each raster's CRS when needed.
#'
#' @param models Model predictions: a [terra::SpatRaster] (one layer per
#'   model), file path(s) to rasters, or a named list of these.
#' @param presence Presence points (data.frame, CSV path, `sf` or
#'   `SpatVector`).
#' @param absence Absence points or, for presence-only data, background /
#'   pseudo-absence points (same formats as `presence`).
#' @param coords Optional names of the x/longitude and y/latitude columns
#'   when points are tabular. Detected automatically if `NULL`.
#' @param crs CRS of tabular coordinates (default WGS84, `"EPSG:4326"`).
#' @param na_action How to handle points that fall on `NA` cells.
#'   `"common"` (default) drops a point for *all* models if it is `NA` in any
#'   model, so every model is evaluated on exactly the same points (a fair
#'   comparison). `"per_model"` keeps points independently for each model
#'   (this reproduces the original 2021 script).
#'
#' @return A data.frame with columns `obs` (1 = presence, 0 = absence) and one
#'   column of predicted values per model. Attribute `"n_dropped"` holds the
#'   number of points removed because of `NA` predictions.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models, sim$presence, sim$absence)
#' head(tab)
extract_predictions <- function(models, presence, absence, coords = NULL,
                                crs = "EPSG:4326",
                                na_action = c("common", "per_model")) {
  na_action <- match.arg(na_action)
  mods <- as_model_list(models)
  pres <- as_points(presence, coords, crs)
  abs <- as_points(absence, coords, crs)
  pts <- rbind(pres, abs)
  obs <- c(rep(1L, nrow(pres)), rep(0L, nrow(abs)))

  vals <- lapply(mods, function(r) {
    p <- pts
    if (!is.na(terra::crs(r)) && terra::crs(r) != "" &&
        !is.na(terra::crs(p)) && terra::crs(p) != "" &&
        !terra::same.crs(r, p)) {
      p <- terra::project(p, r)
    }
    terra::extract(r, p, ID = FALSE)[[1]]
  })
  tab <- data.frame(obs = obs, vals, check.names = FALSE)
  names(tab) <- c("obs", names(mods))

  n_before <- nrow(tab)
  if (na_action == "common") {
    tab <- tab[stats::complete.cases(tab), , drop = FALSE]
  }
  dropped <- n_before - nrow(tab)
  if (dropped > 0) {
    message(dropped, " point(s) fell on NA cells (e.g. outside the raster) ",
            "and were removed.")
  }
  rownames(tab) <- NULL
  attr(tab, "n_dropped") <- dropped
  check_prediction_scale(tab[-1])
  tab
}

check_prediction_scale <- function(preds) {
  rng <- range(unlist(preds), na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    warning("Some predictions lie outside [0, 1] (range: ",
            signif(rng[1], 3), " to ", signif(rng[2], 3), "). MAE and bias ",
            "assume probability-like suitability values; consider rescaling ",
            "(e.g. Maxent 'cloglog'/'logistic' output).", call. = FALSE)
  }
  invisible(NULL)
}
