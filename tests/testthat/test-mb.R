test_that("mb returns valid structure", {
  set.seed(20)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_equal(fit$method, "mb")
  expect_equal(length(fit$path), length(fit$lambda))
  expect_equal(length(fit$sparsity), length(fit$lambda))
})

test_that("mb sparsity is non-decreasing", {
  set.seed(21)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)

  expect_true(all(diff(fit$lambda) < 0))
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("mb path matrices are symmetric and binary", {
  set.seed(22)
  L <- huge.generator(n = 80, d = 30, graph = "band", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)

  for (k in seq_along(fit$path)) {
    p <- as.matrix(fit$path[[k]])
    expect_equal(p, t(p), info = paste("path asymmetric at k =", k))
    expect_true(all(p %in% c(0, 1)))
    expect_true(all(diag(p) == 0))
  }
})

test_that("mb with scr=TRUE works", {
  set.seed(23)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", scr = TRUE, verbose = FALSE)

  expect_s3_class(fit, "huge")
  expect_true(all(diff(fit$sparsity) >= -1e-10))
})

test_that("mb sym='and' gives sparser graphs than sym='or'", {
  set.seed(24)
  L <- huge.generator(n = 80, d = 30, graph = "hub", verbose = FALSE)
  fit_or  <- huge(L$data, method = "mb", sym = "or",  verbose = FALSE)
  fit_and <- huge(L$data, method = "mb", sym = "and", verbose = FALSE)

  # "and" should be at least as sparse as "or" at each lambda
  for (k in seq_along(fit_or$path)) {
    expect_true(sum(fit_and$path[[k]]) <= sum(fit_or$path[[k]]) + 1e-10,
                info = paste("and not sparser than or at k =", k))
  }
})

test_that("mb works across graph types", {
  set.seed(25)
  for (g in c("hub", "band", "cluster", "random")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    fit <- huge(L$data, method = "mb", verbose = FALSE)
    expect_true(all(diff(fit$sparsity) >= -1e-10),
                info = paste("non-monotone sparsity for graph =", g))
  }
})
