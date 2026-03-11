test_that("huge.generator produces valid output for all graph types", {
  set.seed(1)
  for (g in c("hub", "band", "cluster", "random", "scale-free")) {
    L <- huge.generator(n = 60, d = 20, graph = g, verbose = FALSE)
    expect_s3_class(L, "sim")
    expect_equal(dim(L$data), c(60, 20))
    expect_equal(dim(L$sigma), c(20, 20))
    expect_equal(dim(L$omega), c(20, 20))

    # sigma should be positive definite
    eig <- eigen(L$sigma, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig > 0))

    # theta should be symmetric binary
    theta <- as.matrix(L$theta)
    expect_equal(theta, t(theta))
    expect_true(all(theta %in% c(0, 1)))
    expect_true(all(diag(theta) == 0))

    # sparsity should be in [0, 1]
    expect_true(L$sparsity >= 0 && L$sparsity <= 1)
  }
})

test_that("huge.generator respects dimension parameters", {
  set.seed(2)
  L <- huge.generator(n = 100, d = 30, graph = "hub", verbose = FALSE)
  expect_equal(nrow(L$data), 100)
  expect_equal(ncol(L$data), 30)
  expect_equal(nrow(L$sigma), 30)
})
