## V4 -- the reparametrized log-likelihood equals the original.

test_that("V4: the tape's log-likelihood equals gkwdist's at matching parameters", {
  d <- sim_kw(n = 200)
  fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  pv <- fit$parameter_vectors
  direct <- sum(gkwdist::dkw(fit$y, pv$alpha, pv$beta, log = TRUE))
  expect_equal(fit$loglik, direct, tolerance = 1e-10)
})

test_that("V4: the R cascade equals the tape for every family", {
  d <- sim_kw(n = 150)
  logdens <- gq(".gkwq_logdens")
  for (fam in ALL_FAMILIES) {
    fit <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = 0.5, family = fam))
    pv <- fit$parameter_vectors
    ll_r <- sum(logdens(fit$y, pv$alpha, pv$beta, pv$gamma, pv$delta, pv$lambda,
                        delta_is_zero = fit$spec$delta_is_zero == 1L))
    ## Relative, and deliberately looser than the tape-vs-gkwdist check above.
    ## Those two agree to 1e-10 because they evaluate the same compiled path;
    ## here the R cascade and the tape are separate implementations, so the sum
    ## over n observations picks up whatever the platform's compiler does with
    ## operation order and FMA contraction. macOS arm64 lands ~2e-9 away from
    ## x86_64 Linux. A genuine cascade/tape divergence is orders of magnitude
    ## larger than that, so 1e-7 still fails loudly for a real defect.
    expect_equal(ll_r, fit$loglik, tolerance = 1e-7,
                 info = sprintf("family %s", fam))
  }
})

test_that("the log-domain cascade survives responses that break gkwdist::dgkw", {
  ## Responses within 1e-10 of 1 are not exotic in bounded-response work, and
  ## dgkw throws those observations away. The expected value is derived
  ## independently: for the Kumaraswamy, log f = log(alpha*beta) +
  ## (alpha-1) log y + (beta-1) log(1 - y^alpha), and 1 - y^2 = (1-y)(1+y) is
  ## computed without cancellation.
  logdens <- gq(".gkwq_logdens")
  alpha <- 2; beta <- 3
  for (eps in c(1e-8, 1e-10, 1e-12)) {
    y <- 1 - eps
    gap <- 1 - y          # the exactly representable distance from 1
    want <- log(alpha * beta) + (alpha - 1) * log(y) +
      (beta - 1) * (log(gap) + log(1 + y))
    got <- logdens(y, alpha, beta, 1, 0, 1, delta_is_zero = TRUE)
    expect_true(is.finite(got))
    expect_lt(abs(got - want), 1e-6)
  }
  ## And the comparison that motivates the cascade at all.
  expect_true(is.infinite(gkwdist::dgkw(1 - 1e-10, 2, 3, 1, 0, 1, log = TRUE)))
})
