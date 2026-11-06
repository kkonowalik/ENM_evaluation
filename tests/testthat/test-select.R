df <- data.frame(model = c("GLM", "GBM", "RF", "SRE"),
                 AUC = c(0.86, 0.88, 0.87, 0.71),
                 TSS = c(0.62, 0.65, 0.66, 0.38),
                 MAE = c(0.30, 0.27, 0.24, 0.41))

test_that("threshold rule keeps all models above the minimum values", {
  s <- select_models(df)
  expect_s3_class(s, "enm_selection")
  expect_setequal(s$model[s$selected], c("GLM", "GBM", "RF"))
  expect_equal(s$model[s$selected], s$model[seq_len(3)])   # selected first
  expect_false(s$pass[s$model == "SRE"])
})

test_that("Pareto rule and tolerance", {
  s <- select_models(df, rule = "pareto")
  expect_setequal(s$model[s$selected], c("GBM", "RF"))
  s2 <- select_models(df, rule = "pareto", criteria = c(AUC = "max"),
                      tol = 0.025)
  expect_setequal(s2$model[s2$selected], c("GBM", "RF", "GLM"))
})

test_that("fallback, override, missing metrics", {
  expect_warning(s <- select_models(df, min_metrics = c(AUC = 0.95)),
                 "Pareto front is selected")
  expect_setequal(s$model[s$selected], c("GBM", "RF"))
  s <- select_models(df, override = "SRE")
  expect_equal(s$model[s$selected], "SRE")
  expect_error(select_models(df, override = "XX"), "Unknown")
  expect_message(select_models(df, min_metrics = c(AUC = 0.8, Boyce = 0.5)),
                 "Boyce")
  df$AUC[2] <- NA            # missing metric: never passes, worst in Pareto
  s <- select_models(df)
  expect_false(s$selected[s$model == "GBM"])
})

test_that("plot_selection returns a ggplot", {
  s <- select_models(df)
  expect_s3_class(plot_selection(s), "ggplot")
  expect_s3_class(plot_selection(s, x_metric = "TSS"), "ggplot")
  expect_error(plot_selection(df), "select_models")
  expect_error(plot_selection(s, x_metric = "nope"), "not found")
})
