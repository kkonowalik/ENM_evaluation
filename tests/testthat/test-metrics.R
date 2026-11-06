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
