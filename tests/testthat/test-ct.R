test_that("ct returns valid structure", {
  set.seed(40)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "ct", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "ct")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
})

test_that("ct sparsity is non-decreasing", {
  set.seed(41)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "ct", verbose = FALSE)
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("ct path matrices are symmetric and binary", {
  set.seed(42)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "ct", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p))
    expect_true(all(p %in% c(0, 1)))
    expect_true(all(diag(p) == 0))
  }
})
