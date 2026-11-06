sim <- simulate_enm_data(seed = 7, n_presence = 60, n_absence = 150)

test_that("evaluate_models runs end to end on rasters and data.frames", {
  res <- evaluate_models(sim$models, sim$presence, sim$absence)
  expect_s3_class(res, "enm_evaluation")
  expect_equal(nrow(res), terra::nlyr(sim$models))
  expect_true(all(c("AUC", "MAE", "bias", "pareto_front", "rank") %in%
                    names(res)))
  expect_true(any(res$pareto_front == 1))
  # the random model must be clearly worse than the others
  expect_lt(res$AUC[res$model == "Random"], 0.65)
  expect_equal(sort(res$rank), seq_len(nrow(res)))
})

test_that("results equal the original 2021 script logic", {
  tab <- extract_predictions(sim$models, sim$presence, sim$absence)
  res <- evaluate_models(data = tab)
  g <- res[res$model == "GLM", ]
  expect_equal(g$MAE, mean(abs(tab$obs - tab$GLM)))
  expect_equal(g$bias, mean(tab$obs - tab$GLM))
})

test_that("bootstrap adds intervals and Pareto probabilities", {
  res <- evaluate_models(sim$models, sim$presence, sim$absence,
                         n_boot = 30, seed = 1)
  expect_true(all(c("AUC_lower", "AUC_upper", "p_pareto") %in% names(res)))
  expect_true(all(res$AUC_lower <= res$AUC_upper))
  expect_true(all(res$p_pareto >= 0 & res$p_pareto <= 1))
  expect_equal(nrow(attr(res, "bootstrap")), 30 * nrow(res))
  res2 <- evaluate_models(sim$models, sim$presence, sim$absence,
                          n_boot = 30, seed = 1)
  expect_equal(res$p_pareto, res2$p_pareto)
})

test_that("inputs: file paths, CSV, sf, projected CRS", {
  td <- tempfile("enm"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE))
  f <- file.path(td, paste0(names(sim$models), ".tif"))
  for (i in seq_along(f)) terra::writeRaster(sim$models[[i]], f[i])
  csv_p <- file.path(td, "pres.csv"); csv_a <- file.path(td, "abs.csv")
  utils::write.csv(data.frame(id = 1, Longitude = sim$presence$lon,
                              Latitude = sim$presence$lat), csv_p)
  utils::write.csv(data.frame(id = 1, Longitude = sim$absence$lon,
                              Latitude = sim$absence$lat), csv_a)
  ref <- evaluate_models(sim$models, sim$presence, sim$absence)
  res <- evaluate_models(f, csv_p, csv_a)
  expect_equal(res$AUC[match(ref$model, res$model)], ref$AUC)

  # rasters in a projected CRS, points in WGS84 -> reprojected automatically
  proj <- terra::project(sim$models, "EPSG:3035", method = "near")
  res_p <- evaluate_models(proj, sim$presence, sim$absence)
  expect_equal(nrow(res_p), nrow(ref))

  skip_if_not_installed("sf")
  p_sf <- sf::st_as_sf(sim$presence, coords = c("lon", "lat"), crs = 4326)
  a_sf <- sf::st_as_sf(sim$absence, coords = c("lon", "lat"), crs = 4326)
  res_sf <- evaluate_models(sim$models, p_sf, a_sf)
  expect_equal(res_sf$AUC[match(ref$model, res_sf$model)], ref$AUC)
})

test_that("informative errors", {
  expect_error(evaluate_models(), "Supply")
  expect_error(evaluate_models(sim$models, data.frame(a = 1, b = 2),
                               sim$absence), "coordinate columns")
  expect_error(evaluate_models(data = data.frame(obs = c(1, 2), m = 1:2)),
               "only 0")
})

test_that("plot_pareto returns a ggplot", {
  res <- evaluate_models(sim$models, sim$presence, sim$absence, n_boot = 10,
                         seed = 2)
  expect_s3_class(plot_pareto(res), "ggplot")
  tab <- extract_predictions(sim$models, sim$presence, sim$absence)
  expect_s3_class(plot_pareto(evaluate_models(data = tab)), "ggplot")
  expect_error(plot_pareto(data.frame()), "evaluate_models")
})

test_that("fronts and distances stay attached to the right model", {
  res <- evaluate_models(sim$models, sim$presence, sim$absence)
  # brute force: a model is on front 1 iff no other model dominates it
  dominated <- vapply(seq_len(nrow(res)), function(i) {
    any(res$AUC >= res$AUC[i] & res$MAE <= res$MAE[i] &
          (res$AUC > res$AUC[i] | res$MAE < res$MAE[i]))
  }, logical(1))
  expect_equal(res$pareto_front == 1, !dominated)
  ref <- pareto_rank(as.data.frame(res)[c("model", "AUC", "MAE")])
  expect_equal(res$distance, ref$distance[match(res$model, ref$model)])
})
