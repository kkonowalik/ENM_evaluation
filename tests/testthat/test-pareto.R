test_that("Pareto fronts are identified correctly", {
  df <- data.frame(model = c("A", "B", "C", "D", "E"),
                   AUC = c(0.90, 0.85, 0.80, 0.88, 0.90),
                   MAE = c(0.30, 0.20, 0.25, 0.32, 0.30))
  r <- pareto_rank(df)
  f <- setNames(r$pareto_front, r$model)
  expect_equal(unname(f[c("A", "B", "E")]), c(1L, 1L, 1L))  # ties both stay
  expect_equal(unname(f["C"]), 2L)  # dominated by B
  expect_equal(unname(f["D"]), 2L)  # dominated by A
  expect_true(all(r$distance >= 0 & r$distance <= 1))
})

test_that("a model best on all criteria has distance 0 and rank 1", {
  df <- data.frame(model = c("best", "x", "y"),
                   AUC = c(0.95, 0.7, 0.8), MAE = c(0.1, 0.3, 0.2))
  r <- pareto_rank(df)
  expect_equal(r$model[1], "best")
  expect_equal(r$distance[1], 0)
  expect_equal(sum(r$pareto_front == 1), 1)
})

test_that("pareto_rank works with more than two criteria", {
  df <- data.frame(m = 1:3, a = c(1, 2, 3), b = c(3, 2, 1), c = c(1, 1, 0))
  r <- pareto_rank(df, maximize = c("a", "c"), minimize = "b")
  # 3 is best on a and b; 1 and 2 tie on c (1 > 0), 1 is dominated by 2
  expect_setequal(r$m[r$pareto_front == 1], c(2, 3))
  expect_equal(r$pareto_front[r$m == 1], 2L)
  expect_error(pareto_rank(df, maximize = "zzz"), "not found")
})

test_that("tolerance turns near-ties into ties", {
  df <- data.frame(model = c("A", "B"), AUC = c(0.900, 0.895),
                   MAE = c(0.200, 0.201))
  expect_equal(sum(pareto_rank(df)$pareto_front == 1), 1)        # A wins
  expect_equal(sum(pareto_rank(df, tol = 0.01)$pareto_front == 1), 2)
  # named tolerance: only AUC differences are forgiven
  df2 <- data.frame(model = c("A", "B"), AUC = c(0.900, 0.895),
                    MAE = c(0.20, 0.25))
  r <- pareto_rank(df2, tol = c(AUC = 0.01))
  expect_equal(r$model[r$pareto_front == 1], "A")
})
