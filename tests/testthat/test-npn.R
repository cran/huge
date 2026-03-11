test_that("huge.npn shrinkage transformation works", {
  set.seed(60)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  Q <- huge.npn(L$data, npn.func = "shrinkage", verbose = FALSE)

  expect_equal(dim(Q), dim(L$data))
  expect_true(all(is.finite(Q)))
  # columns should have unit variance (approximately)
  sds <- apply(Q, 2, sd)
  expect_true(all(abs(sds - 1) < 0.1))
})

test_that("huge.npn truncation transformation works", {
  set.seed(61)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  Q <- huge.npn(L$data, npn.func = "truncation", verbose = FALSE)

  expect_equal(dim(Q), dim(L$data))
  expect_true(all(is.finite(Q)))
})

test_that("huge.npn skeptic returns a correlation matrix", {
  set.seed(62)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  Q <- huge.npn(L$data, npn.func = "skeptic", verbose = FALSE)

  expect_equal(dim(Q), c(20, 20))
  expect_equal(Q, t(Q), tolerance = 1e-10)
  expect_true(all(abs(diag(Q) - 1) < 1e-10))
})

test_that("npn-transformed data works with estimation methods", {
  set.seed(63)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  Q <- huge.npn(L$data^3, verbose = FALSE)  # nonlinear distortion
  fit <- huge(Q, method = "mb", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})
