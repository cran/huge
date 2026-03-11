test_that("tiger returns valid structure", {
  set.seed(30)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "tiger", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "tiger")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$icov), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
})

test_that("tiger sparsity is non-decreasing", {
  set.seed(31)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "tiger", verbose = FALSE)

  expect_true(all(diff(fit$lambda) < 0))
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("tiger path matrices are symmetric and binary", {
  set.seed(32)
  L <- huge.generator(n = 80, d = 30, graph = "band", verbose = FALSE)
  fit <- huge(L$data, method = "tiger", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p), info = paste("path asymmetric at k =", k))
    expect_true(all(p %in% c(0, 1)))
    expect_true(all(diag(p) == 0))
  }
})

test_that("tiger works across graph types", {
  set.seed(33)
  for (g in c("hub", "band", "cluster")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    fit <- huge(L$data, method = "tiger", verbose = FALSE)
    expect_true(all(diff(fit$sparsity) >= -1e-10),
                info = paste("non-monotone sparsity for graph =", g))
  }
})
