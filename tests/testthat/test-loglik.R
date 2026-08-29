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
    ## The cascade and the tape are separate implementations, so the sum over n
    ## observations picks up whatever each compiler does with operation order
    ## and FMA contraction; macOS arm64 lands ~2e-9 from x86_64 Linux.
    ##
    ## There is also a hard floor neither side can beat. The normalization
    ## carries -log B(gamma, delta+1), and both sides build it as
    ## lgamma(g) + lgamma(d+1) - lgamma(g+d+1) -- a difference of numbers of
    ## magnitude lgamma(g+d+1) whose true value is a few hundred. The absolute
    ## precision available is therefore eps * lgamma(g+d+1) per observation,
    ## and a comparison tighter than that tests the arithmetic rather than the
    ## two implementations. It bites in exactly one place: a gkw fit that does
    ## not converge runs off to gamma = 6.3e9, where lgamma is 1.4e11 and the
    ## floor is 4.6e-3 over 150 observations -- which is the gap seen there, to
    ## within one per cent.
    ##
    ## The factor of four covers the four lgamma evaluations across the two
    ## sides. A genuine cascade/tape divergence is orders of magnitude larger
    ## than either bound and still fails loudly.
    fp_floor <- 4 * length(fit$y) * .Machine$double.eps *
      max(abs(lgamma(pv$gamma + pv$delta + 1)), 1)
    expect_lt(abs(ll_r - fit$loglik),
              max(fp_floor, 1e-7 * abs(fit$loglik)),
              label = sprintf("cascade vs tape, family %s", fam))
  }
})

test_that("the log-domain cascade holds at responses within 1e-12 of 1", {
  ## Responses within 1e-10 of 1 are not exotic in bounded-response work. The
  ## expected value is derived independently: for the Kumaraswamy, log f =
  ## log(alpha*beta) + (alpha-1) log y + (beta-1) log(1 - y^alpha), and
  ## 1 - y^2 = (1-y)(1+y) is computed without cancellation.
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
  ## This block used to assert that gkwdist::dgkw returned -Inf here, which was
  ## the reason the cascade was written. Whether it still does depends on the
  ## version installed: 1.1.5 returns -Inf, 1.1.6 handles the boundary. Neither
  ## is something this package's tests should pin, so the comparison is made
  ## only where gkwdist itself produces a finite value -- where both can compute
  ## it, they must agree, and that is a check on the cascade rather than a claim
  ## about someone else's release.
  skip_if_not_installed("gkwdist")
  ## Mid-range first, so the block cannot pass by comparing nothing: every
  ## version computes this one, and agreement there is the anchor.
  expect_equal(logdens(0.6, alpha, beta, 1, 0, 1, delta_is_zero = TRUE),
               gkwdist::dgkw(0.6, alpha, beta, 1, 0, 1, log = TRUE),
               tolerance = 1e-10)
  ## Then the boundary, wherever the installed version can still reach it.
  for (eps in c(1e-8, 1e-10, 1e-12, 1e-14)) {
    y <- 1 - eps
    ref <- gkwdist::dgkw(y, alpha, beta, 1, 0, 1, log = TRUE)
    if (!is.finite(ref)) next
    expect_equal(logdens(y, alpha, beta, 1, 0, 1, delta_is_zero = TRUE), ref,
                 tolerance = 1e-10, info = sprintf("1 - y = %.0e", eps))
  }
})

test_that("the log-CDF keeps the upper tail that pbeta's argument destroys", {
  ## log F and log S were both taken from pbeta(exp(lambda * lv), ...). The
  ## argument rounds to 1 as soon as lambda*lv falls below about -1e-17, so the
  ## upper tail collapsed to -Inf while its true value stayed finite: for
  ## kw with alpha = 2, beta = 60 the function returned -Inf at y = 0.8, where
  ## log S is -201.93, and at y = 1 - 1e-10, where it is -1339.96.
  ##
  ## For the Kumaraswamy sub-family the survival function is elementary,
  ## log S(y) = beta * log(1 - y^alpha), which is what the result must equal.
  lcdf <- gq(".gkwq_logcdf"); l1me <- gq(".log1mexp")
  alpha <- 2; beta <- 60
  for (y in c(0.8, 0.95, 1 - 1e-2, 1 - 1e-6, 1 - 1e-10, 1 - 1e-14)) {
    exact <- beta * l1me(alpha * log(y))
    expect_equal(lcdf(y, alpha, beta, 1, 0, 1, lower.tail = FALSE), exact,
                 tolerance = 1e-10, info = sprintf("y = %.17g", y))
  }

  ## The two tails still sum to one wherever both are representable, and the
  ## result is unchanged from the previous route in that regime.
  set.seed(3)
  for (i in seq_len(200)) {
    y <- stats::runif(1, 0.05, 0.95)
    a <- stats::runif(1, 0.5, 3); b <- stats::runif(1, 0.5, 5)
    g <- stats::runif(1, 0.5, 3); d <- stats::runif(1, 0, 3)
    L <- stats::runif(1, 0.5, 2)
    lF <- lcdf(y, a, b, g, d, L, lower.tail = TRUE)
    lS <- lcdf(y, a, b, g, d, L, lower.tail = FALSE)
    expect_equal(exp(lF) + exp(lS), 1, tolerance = 1e-9)
  }
})
