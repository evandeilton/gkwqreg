## V13 -- the estimator recovers the truth, and the fitted quantile really is a
## quantile.

test_that("V13: coefficients recover the data-generating truth", {
  ## The design study's pilot: 200 replicates at n = 400 gave bias below 0.004
  ## with 100% convergence. Here a smaller version, with tolerances scaled to
  ## the Monte Carlo error of the replicate count.
  skip_on_cran()
  R <- 40L; n <- 400L; tau <- 0.5
  est <- t(vapply(seq_len(R), function(r) {
    d <- sim_kw(n = n, tau = tau, seed = 1000 + r)
    f <- try(gkwqreg(y ~ x, data = d, tau = tau, family = "kw"), silent = TRUE)
    if (inherits(f, "try-error")) return(c(NA, NA, NA))
    c(f$coef_list$mu[[1]], f$coef_list$mu[[2]], exp(f$coef_list$alpha[[1]]))
  }, numeric(3)))
  conv <- mean(stats::complete.cases(est))
  expect_gt(conv, 0.95)
  m <- colMeans(est, na.rm = TRUE)
  s <- apply(est, 2, stats::sd, na.rm = TRUE) / sqrt(sum(stats::complete.cases(est)))
  expect_lt(abs(m[1] - 0.4), 4 * s[1])
  expect_lt(abs(m[2] - 1.1), 4 * s[2])
  expect_lt(abs(m[3] - 2.0), 4 * s[3])
})

test_that("the fitted values really are conditional tau-quantiles", {
  ## Empirical coverage must sit at tau, by construction of the parametrization.
  d <- sim_kw(n = 1000, seed = 7)
  for (tau in c(0.1, 0.5, 0.9)) {
    f <- gkwqreg(y ~ x, data = d, tau = tau, family = "kw")
    cov_emp <- mean(f$y <= fitted(f))
    expect_lt(abs(cov_emp - tau), 4 * sqrt(tau * (1 - tau) / 1000))
  }
})

test_that("fitted() returns quantiles, not means", {
  ## The single most likely mistake a gkwreg user can make on this class.
  d <- sim_kw(n = 300)
  f <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
  expect_false(isTRUE(all.equal(fitted(f), fitted(f, type = "mean"))))
  expect_true(all(fitted(f) > fitted(f, type = "mean")))   # 90th pct above mean
})
