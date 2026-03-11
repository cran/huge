test_that("huge.inference Gaussian produces valid p-values", {
  set.seed(70)
  L <- huge.generator(n = 100, d = 15, graph = "hub", g = 3, verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)
  T_hat <- tail(fit$icov, 1)[[1]]
  inf <- huge.inference(L$data, T_hat, L$theta)

  expect_true(!is.null(inf$p))
  expect_equal(dim(inf$p), c(15, 15))
  # p-values should be in [0, 1]
  expect_true(all(inf$p >= 0 & inf$p <= 1))
  # error rate should be in [0, 1]
  expect_true(inf$error >= 0 && inf$error <= 1)
})
