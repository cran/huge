test_that("huge.select with ebic works for glasso", {
  set.seed(50)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "glasso", verbose = FALSE)
  sel <- huge.select(fit, criterion = "ebic", verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_true(!is.null(sel$refit))
  expect_true(!is.null(sel$opt.lambda))
  expect_true(!is.null(sel$opt.index))
  expect_true(sel$opt.lambda > 0)
  expect_true(sel$opt.index >= 1 && sel$opt.index <= length(sel$lambda))
  expect_equal(dim(sel$refit), c(20, 20))
})

test_that("huge.select with ric works for mb", {
  set.seed(51)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)
  sel <- huge.select(fit, criterion = "ric", verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_true(!is.null(sel$refit))
  expect_true(sel$opt.lambda > 0)
})

test_that("huge.select with stars works for mb", {
  set.seed(52)
  L <- huge.generator(n = 80, d = 20, graph = "hub", verbose = FALSE)
  fit <- huge(L$data, method = "mb", verbose = FALSE)
  sel <- huge.select(fit, criterion = "stars", rep.num = 5, verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_true(!is.null(sel$refit))
  expect_true(!is.null(sel$variability))
  expect_equal(length(sel$variability), length(sel$lambda))
  expect_true(all(sel$variability >= 0))
})
