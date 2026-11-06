test_that("AUC matches the Mann-Whitney definition, incl. ties", {
  obs <- c(1, 1, 1, 0, 0, 0, 0)
  pred <- c(0.9, 0.5, 0.5, 0.5, 0.2, 0.1, 0.6)
  # brute-force: P(pred_pres > pred_abs) + 0.5 * P(tie)
  pp <- pred[obs == 1]; pa <- pred[obs == 0]
  cmp <- outer(pp, pa, ">") + 0.5 * outer(pp, pa, "==")
  expect_equal(auc_score(obs, pred), mean(cmp))
  expect_equal(auc_score(c(1, 0), c(1, 0)), 1)
  expect_true(is.na(auc_score(c(1, 1), c(0.2, 0.3))))
})

test_that("MAE and bias match the original Metrics-package formulas", {
  obs <- c(1, 0, 1, 0)
  pred <- c(0.8, 0.3, 0.4, 0.1)
  expect_equal(mae_score(obs, pred), mean(abs(obs - pred)))
  expect_equal(bias_score(obs, pred), mean(obs - pred))
  if (requireNamespace("Metrics", quietly = TRUE)) {
    expect_equal(mae_score(obs, pred), Metrics::mae(obs, pred))
    expect_equal(bias_score(obs, pred), Metrics::bias(obs, pred))
    expect_equal(auc_score(obs, pred), Metrics::auc(obs, pred))
  }
})

test_that("TSS equals the brute-force maximum over thresholds", {
  set.seed(4)
  obs <- rbinom(300, 1, 0.3)
  pred <- round(plogis(obs * 1.5 + rnorm(300)), 2)   # rounded: many ties
  brute <- max(vapply(c(sort(unique(pred)), Inf), function(t) {
    mean(pred[obs == 1] >= t) + mean(pred[obs == 0] < t) - 1
  }, numeric(1)))
  expect_equal(as.numeric(tss_score(obs, pred)), brute)
  expect_equal(as.numeric(tss_score(c(1, 1, 0, 0), c(0.9, 0.8, 0.2, 0.1))), 1)
  expect_equal(as.numeric(tss_score(c(1, 0, 1, 0), c(0.9, 0.2, 0.3, 0.6),
                                    threshold = 0.5)), 0)
  expect_true(is.na(tss_score(c(1, 1), c(0.5, 0.6))))
})

test_that("Boyce index reproduces ecospat::ecospat.boyce()", {
  # reference values computed with ecospat 4.1.4 (rounded to 3 decimals)
  set.seed(10); bg <- runif(500); pr <- rbeta(60, 4, 1)
  expect_equal(round(boyce_index(pr, bg), 3), 0.961)
  set.seed(11); bg <- rbeta(300, 1, 3); pr <- runif(25)
  expect_equal(round(boyce_index(pr, bg), 3), 0.743)
  set.seed(12); bg <- round(runif(200), 2); pr <- round(rbeta(40, 2, 2), 2)
  expect_equal(round(boyce_index(pr, bg), 3), 0.268)
})

test_that("Boyce index handles degenerate input", {
  expect_true(is.na(boyce_index(0.5, runif(10))))           # 1 presence
  expect_true(is.na(boyce_index(rep(0.5, 5), rep(0.5, 50)))) # no range
  expect_true(is.na(boyce_index(c(NA, NA), runif(10))))
})

test_that("threshold_p10 is the 10th percentile", {
  expect_equal(threshold_p10(1:11 / 11), unname(quantile(1:11 / 11, 0.1)))
})
