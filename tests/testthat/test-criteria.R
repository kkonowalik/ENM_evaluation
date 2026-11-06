sim <- simulate_enm_data(seed = 8, n_presence = 80, n_absence = 250)
tab <- extract_predictions(sim$models, sim$presence, sim$absence)

test_that("default criteria are AUC (max) and MAE (min)", {
  r <- evaluate_models(data = tab)
  expect_equal(attr(r, "criteria"), c(AUC = "max", MAE = "min"))
  expect_true(all(c("TSS", "Boyce") %in% names(r)))
})

test_that("front follows the chosen criteria (brute force)", {
  check_front <- function(r, crit) {
    m <- sapply(names(crit), function(k) {
      v <- r[[k]]
      switch(crit[[k]], max = -v, min = v, zero = abs(v))
    })
    dominated <- sapply(seq_len(nrow(m)), function(i) any(
      apply(m, 1, function(o) all(o <= m[i, ]) && any(o < m[i, ]))))
    expect_equal(r$pareto_front == 1, !dominated)
  }
  for (crit in list(c(TSS = "max", Boyce = "max"),
                    c(AUC = "max", bias = "zero"),
                    c(AUC = "max", MAE = "min", Boyce = "max"))) {
    r <- evaluate_models(data = tab, criteria = crit)
    expect_equal(attr(r, "criteria"), crit)
    check_front(r, crit)
  }
})

test_that("custom metrics can be computed and used as criteria", {
  brier <- function(obs, pred) mean((obs - pred)^2)
  r <- evaluate_models(data = tab, metrics = list(Brier = brier),
                       criteria = c(AUC = "max", Brier = "min"))
  expect_equal(r$Brier[r$model == "GLM"], brier(tab$obs, tab$GLM))
  expect_s3_class(plot_pareto(r), "ggplot")
  cv <- evaluate_cv(data.frame(model = "GLM", fold = rep(1:2, length.out =
                                 nrow(tab)), obs = tab$obs, pred = tab$GLM),
                    metrics = list(Brier = brier))
  expect_true(all(c("Brier", "Brier_sd") %in% names(cv)))
})

test_that("bootstrap gives intervals for every chosen criterion", {
  r <- evaluate_models(data = tab, criteria = c(TSS = "max", Boyce = "max"),
                       n_boot = 20, seed = 1)
  expect_true(all(c("TSS_lower", "TSS_upper", "TSS_sd", "Boyce_lower",
                    "Boyce_upper", "Boyce_sd", "p_pareto") %in% names(r)))
  expect_false("AUC_lower" %in% names(r))
  expect_true(all(r$TSS_lower <= r$TSS_upper))
})

test_that("helpful errors for wrong criteria or metrics", {
  expect_error(evaluate_models(data = tab, criteria = c(AUCC = "max")),
               "Available: AUC, MAE")
  expect_error(evaluate_models(data = tab, criteria = c(AUC = "high")),
               "\"max\"")
  expect_error(evaluate_models(data = tab, criteria = c("AUC", "MAE")),
               "named vector")
  expect_error(evaluate_models(data = tab, metrics = function(o, p) 1),
               "named list")
  expect_error(evaluate_models(data = tab, metrics = list(AUC = mean)),
               "already used")
  expect_error(evaluate_models(data = tab, metrics = list(B = function(o, p)
    c(1, 2))), "single number")
})

test_that("plots work for any metric pair and direction", {
  r <- evaluate_models(data = tab, criteria = c(AUC = "max", bias = "zero"),
                       n_boot = 10, seed = 1)
  expect_s3_class(plot_pareto(r), "ggplot")
  expect_s3_class(plot_pareto(r, x_metric = "TSS", zoom = TRUE), "ggplot")
  expect_error(plot_pareto(r, x_metric = "nope"), "not found")
  r3 <- evaluate_models(data = tab,
                        criteria = c(AUC = "max", MAE = "min", TSS = "max"))
  p <- plot_pareto(r3)
  expect_match(p$labels$subtitle, "3 criteria")
  s <- select_models(r3, criteria = c(AUC = "max", bias = "zero"))
  expect_s3_class(plot_selection(s, "bias", "AUC"), "ggplot")
})

test_that("select_models and pareto_rank support 'zero'", {
  df <- data.frame(model = c("a", "b", "c"), AUC = c(0.9, 0.9, 0.8),
                   bias = c(-0.3, 0.1, 0))
  s <- select_models(df, criteria = c(AUC = "max", bias = "zero"),
                     rule = "pareto")
  expect_setequal(s$model[s$pareto], c("b", "c"))
  r <- pareto_rank(df, maximize = "AUC", minimize = NULL, zero = "bias")
  expect_setequal(r$model[r$pareto_front == 1], c("b", "c"))
})
