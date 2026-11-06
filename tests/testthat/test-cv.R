sim <- simulate_enm_data(seed = 3, n_presence = 120, n_absence = 400)
pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))

test_that("spatial_folds gives k balanced folds and keeps blocks intact", {
  f <- spatial_folds(pts, k = 4, block_size = 100, occ = pts$occ, seed = 1)
  expect_equal(sort(unique(as.integer(f))), 1:4)
  expect_true(all(tabulate(f[pts$occ == 1], 4) > 10))
  blk <- attr(f, "block")
  expect_true(all(tapply(f, blk, function(v) length(unique(v))) == 1))
  expect_equal(as.integer(f), as.integer(spatial_folds(
    pts, k = 4, block_size = 100, occ = pts$occ, seed = 1)))
})

test_that("spatial_folds does not change the user's random numbers", {
  set.seed(99); a <- runif(1)
  set.seed(99); spatial_folds(pts, 4, 100, occ = pts$occ, seed = 1)
  expect_equal(runif(1), a)
})

test_that("spatial_folds: errors and small-data warnings", {
  expect_error(spatial_folds(pts, k = 5), "block_size")
  expect_error(spatial_folds(pts, k = 5, block_size = 5000), "fewer than k")
  small <- rbind(cbind(sim$presence[1:6, ], occ = 1),
                 cbind(sim$absence[1:60, ], occ = 0))
  expect_warning(spatial_folds(small, k = 5, block_size = 60,
                               occ = small$occ, seed = 1), "presences")
})

test_that("spatial_folds accepts projected points", {
  v <- terra::project(terra::vect(pts, geom = c("lon", "lat"),
                                  crs = "EPSG:4326"), "EPSG:3035")
  f <- spatial_folds(v, k = 4, block_size = 100, occ = pts$occ, seed = 1)
  expect_length(f, nrow(pts))
})

# hold-out predictions of two GLMs on spatial folds, two background sets
make_holdout <- function() {
  pts$fold <- spatial_folds(pts, k = 4, block_size = 100, occ = pts$occ,
                            seed = 1)
  pts <- cbind(pts, terra::extract(sim$env, pts[c("lon", "lat")], ID = FALSE))
  pts$set <- ifelse(pts$occ == 1, 0, rep(1:2, length.out = nrow(pts)))
  out <- list()
  for (s in 1:2) for (f in 1:4) {
    d <- pts[pts$set %in% c(0, s), ]
    tr <- d[d$fold != f, ]
    te <- d[d$fold == f, ]
    good <- glm(occ ~ env1 + I(env1^2) + env2 + I(env2^2), binomial, tr)
    poor <- glm(occ ~ env1, binomial, tr)
    out[[length(out) + 1]] <- rbind(
      data.frame(model = "quadratic", set = s, fold = f, obs = te$occ,
                 pred = predict(good, te, type = "response")),
      data.frame(model = "linear", set = s, fold = f, obs = te$occ,
                 pred = predict(poor, te, type = "response")))
  }
  do.call(rbind, out)
}
ho <- suppressWarnings(make_holdout())  # glm: perfect separation is fine here

test_that("evaluate_cv averages per-run metrics correctly", {
  cv <- evaluate_cv(ho, set = "set")
  expect_s3_class(cv, "enm_cv")
  runs <- attr(cv, "runs")
  expect_equal(nrow(runs), 2 * 2 * 4)
  r1 <- ho[ho$model == "quadratic" & ho$set == 1 & ho$fold == 2, ]
  expect_equal(runs$AUC[runs$model == "quadratic" & runs$set == 1 &
                          runs$fold == 2], auc_score(r1$obs, r1$pred))
  q <- cv[cv$model == "quadratic", ]
  expect_equal(q$AUC, mean(runs$AUC[runs$model == "quadratic"]))
  expect_equal(q$AUC_sd, sd(runs$AUC[runs$model == "quadratic"]))
  expect_gt(q$AUC, cv$AUC[cv$model == "linear"])
  expect_true(all(c("AUC_sd_sets", "AUC_sd_folds") %in% names(cv)))
  expect_equal(ncol(attr(cv, "by_fold")), 1 + 4)
  expect_output(print(cv), "background sets")
})

test_that("evaluate_cv: fixed thresholds, factors, empty folds", {
  ho$thr <- 0.5
  a <- evaluate_cv(ho, set = "set")
  b <- evaluate_cv(ho, set = "set", threshold = "thr")
  expect_true(all(b$TSS <= a$TSS + 1e-12))
  ho$model <- factor(ho$model)
  expect_equal(sort(evaluate_cv(ho)$model), c("linear", "quadratic"))
  bad <- ho
  bad$obs[bad$fold == 1] <- 0
  expect_warning(evaluate_cv(bad), "no presences")
  expect_error(evaluate_cv(ho, fold = "nope"), "not found")
})

test_that("full chain: cv -> select -> plot", {
  cv <- evaluate_cv(ho, set = "set")
  s <- select_models(cv, min_metrics = c(AUC = 0.8))
  expect_equal(s$model[s$selected], "quadratic")
  expect_s3_class(plot_selection(s), "ggplot")
})
