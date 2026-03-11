test_that("glasso returns valid structure", {
  set.seed(10)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "glasso")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$icov), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
  expect_equal(length(fit$loglik), length(fit$lambda))
  expect_equal(length(fit$df), length(fit$lambda))
})

test_that("glasso lambda path is decreasing and sparsity is non-decreasing", {
  set.seed(11)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  expect_true(all(diff(fit$lambda) < 0))
  expect_true(all(diff(fit$sparsity) >= -1e-10))
  expect_true(all(diff(fit$df) >= 0))
})

test_that("glasso loglik values are finite", {
  set.seed(12)
  L <- huge.generator(n = 80, d = 30, graph = "band", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)
  expect_true(all(is.finite(fit$loglik)))
})

test_that("glasso path matrices are symmetric", {
  set.seed(13)
  L <- huge.generator(n = 80, d = 30, graph = "cluster", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p), info = paste("path asymmetric at k =", k))
    expect_true(all(diag(p) == 0))
  }
})

test_that("glasso icov matrices are symmetric", {
  set.seed(14)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)

  for (k in seq_along(fit$icov)) {
    ic <- fit$icov[[k]]
    expect_equal(ic, t(ic), tolerance = 1e-4,
                 info = paste("icov asymmetric at k =", k))
  }
})

test_that("glasso with scr=TRUE produces valid results", {
  set.seed(15)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", scr = TRUE, verbose = FALSE)

  expect_true(all(diff(fit$sparsity) >= -1e-10))
  expect_true(all(is.finite(fit$loglik)))
})

test_that("glasso cov.output works", {
  set.seed(16)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", cov.output = TRUE, verbose = FALSE)

  expect_true(!is.null(fit$cov))
  expect_equal(length(fit$cov), length(fit$lambda))
  for (k in seq_along(fit$cov)) {
    co <- fit$cov[[k]]
    expect_equal(dim(co), c(30, 30))
    expect_equal(co, t(co), tolerance = 1e-10)
  }
})

test_that("glasso accepts covariance matrix input", {
  set.seed(17)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  S <- cor(L$data)
  fit <- huge(S, method = "glasso", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_true(fit$cov.input)
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("glasso works across graph types", {
  set.seed(18)
  for (g in c("hub", "band", "cluster", "random")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    fit <- huge(L$data, method = "glasso", verbose = FALSE)
    expect_true(all(diff(fit$sparsity) >= -1e-10),
                info = paste("non-monotone sparsity for graph =", g))
  }
})
