#' Simulate example data for trying the package
#'
#' Creates a virtual species on a small grid, samples presence and absence
#' points, and produces several "model" rasters of varying quality. This lets
#' you test the workflow without any of your own data.
#'
#' @param n_presence,n_absence Number of points to sample.
#' @param seed Optional random seed.
#' @return A list with `models` (a multi-layer `SpatRaster`, values 0-1),
#'   `presence` and `absence` (data.frames with `lon` and `lat` columns) and
#'   `truth` (the true suitability raster).
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 42)
#' terra::plot(sim$models)
simulate_enm_data <- function(n_presence = 100, n_absence = 300, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  r <- terra::rast(nrows = 60, ncols = 80, xmin = 14, xmax = 24, ymin = 49,
                   ymax = 55, crs = "EPSG:4326")
  xy <- terra::crds(r)
  env1 <- (xy[, 1] - 14) / 10
  env2 <- (xy[, 2] - 49) / 6
  noise <- function(sd) stats::rnorm(nrow(xy), 0, sd)
  logistic <- function(z) 1 / (1 + exp(-z))
  truth <- logistic(6 - 30 * (env1 - 0.55)^2 - 25 * (env2 - 0.45)^2)

  make <- function(v) {
    out <- terra::rast(r)
    terra::values(out) <- pmin(pmax(v, 0), 1)
    out
  }
  models <- c(
    make(logistic(stats::qlogis(pmin(pmax(truth, 1e-4), 1 - 1e-4)) + noise(0.6))),
    make(logistic(stats::qlogis(pmin(pmax(truth, 1e-4), 1 - 1e-4)) * 0.5 + noise(0.4))),
    make(truth^0.4 + noise(0.05)),
    make(truth * 0.6 + noise(0.08)),
    make(logistic(5 - 20 * (env1 - 0.5)^2) + noise(0.05)),
    make(stats::runif(nrow(xy))),
    make(0.7 * truth + 0.3 * logistic(4 - 25 * (env2 - 0.3)^2) + noise(0.1))
  )
  names(models) <- c("Maxent_LQ", "Maxent_H_rm4", "GLM", "RandomForest",
                     "BIOCLIM", "Random", "GBM")

  tr <- make(truth)
  names(tr) <- "truth"
  pv <- terra::values(tr)[, 1]
  pres_cells <- sample(seq_along(pv), n_presence, prob = pv^2)
  abs_cells <- sample(seq_along(pv), n_absence, prob = (1 - pv)^2)
  jitter <- function(cells) {
    p <- terra::xyFromCell(r, cells)
    data.frame(lon = p[, 1] + stats::runif(length(cells), -0.06, 0.06),
               lat = p[, 2] + stats::runif(length(cells), -0.05, 0.05))
  }
  list(models = models, presence = jitter(pres_cells),
       absence = jitter(abs_cells), truth = tr)
}
