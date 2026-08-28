## SPEC 1.6 -- the three identifiability facts, each as an enforced guard.

test_that("delta = 0 confounds gamma and lambda", {
  ## With delta = 0, I_u(gamma, 1) = u^gamma, so F(y) = v^(gamma*lambda) and
  ## only the PRODUCT is identified.
  yv <- c(.2, .5, .8)
  a <- gkwdist::pgkw(yv, 1.5, 2, 3, 0, 0.7)
  b <- gkwdist::pgkw(yv, 1.5, 2, 0.7, 0, 3)   # same product 2.1
  expect_lt(max(abs(a - b)), 1e-15)
})

test_that("the lambda anchor is refused when delta is fixed at 0 with gamma free", {
  ## No shipped family reaches this combination, because every family that fixes
  ## delta at 0 also fixes gamma at 1. The guard is enforced at the registry so
  ## it stays correct if a `fixed` argument is ever added.
  check <- gq(".gkwq_check_identifiability")
  spec <- list(family = "synthetic", anchor = "lambda",
               full = c("alpha", "beta", "gamma", "lambda"),
               fixed = list(delta = 0))
  expect_error(check(spec, warn = FALSE), "not identified")
})

test_that("family = 'gkw' warns that it is weakly identified", {
  d <- sim_kw(n = 150)
  expect_warning(gkwqreg(y ~ x, data = d, tau = 0.5, family = "gkw"),
                 "weakly identified")
})

test_that("tau is validated and never estimated", {
  d <- sim_kw(n = 100)
  expect_error(gkwqreg(y ~ x, data = d, tau = 0, family = "kw"), "inside \\(0,1\\)")
  expect_error(gkwqreg(y ~ x, data = d, tau = 1, family = "kw"), "inside \\(0,1\\)")
  expect_error(gkwqreg(y ~ x, data = d, tau = NA, family = "kw"), "inside \\(0,1\\)")
  fit <- gkwqreg(y ~ x, data = d, tau = 0.3, family = "kw")
  expect_false(any(grepl("tau", names(coef(fit)))))
})

test_that("responses at the boundary are refused rather than silently clamped", {
  d <- sim_kw(n = 100)
  d$y[1] <- 0
  expect_error(gkwqreg(y ~ x, data = d, family = "kw"), "strictly in \\(0,1\\)")
  d$y[1] <- 1
  expect_error(gkwqreg(y ~ x, data = d, family = "kw"), "strictly in \\(0,1\\)")
})

test_that("an empty sample is refused with a clear message", {
  d <- sim_kw(n = 50)
  expect_error(gkwqreg(y ~ x, data = d, subset = x > 1e6, family = "kw"),
               "no observations remain")
})

test_that("the condition number separates well- from ill-conditioned fits", {
  ## summary() reports the information matrix's condition number precisely so
  ## that a weakly identified fit is visible rather than silent. `mc` with the
  ## lambda anchor rides a gamma-lambda ridge even when delta > 0: only the
  ## product is well determined, and the condition number says so.
  set.seed(3)
  n <- 400
  x <- runif(n, -1.5, 1.5)
  mu <- plogis(0.3 + 0.9 * x)
  spec <- gq(".gkwq_family_info")("mc", NULL)
  P <- gq(".gkwq_reconstruct")(mu, 0.5, spec,
                               list(gamma = rep(1.4, n), delta = rep(0.8, n)))
  d <- data.frame(y = gkwdist::rgkw(n, P[, 1], P[, 2], P[, 3], P[, 4], P[, 5]),
                  x = x)

  f_kw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  f_mc <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = 0.5, family = "mc"))

  expect_lt(f_kw$cond_number, 1e4)
  expect_gt(f_mc$cond_number, 1e6)
  expect_output(print(summary(f_mc)), "ill-conditioned")

  ## The quantile coefficients survive the ridge even though the nuisance
  ## parameters do not: that is the whole point of anchoring on the quantile.
  expect_lt(abs(coef(f_mc)[["mu:x"]] - 0.9), 0.25)
})
