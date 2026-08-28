## V5, V6, V7 -- the three invariance properties. Together they encode the
## design study's central finding: the anchor is a modeling choice, and it is a
## pure reparametrization exactly when the nuisance parameters are saturated.

test_that("V5: tau does not enter an intercept-only likelihood", {
  ## The profile likelihood in tau is exactly flat: for any tau' some parameter
  ## value reproduces the same distribution. tau indexes the question, not the
  ## model, which is why it is never estimated.
  d <- sim_kw(n = 300)
  f25 <- gkwqreg(y ~ 1, data = d, tau = 0.25, family = "kw")
  f75 <- gkwqreg(y ~ 1, data = d, tau = 0.75, family = "kw")
  expect_equal(f25$loglik, f75$loglik, tolerance = 1e-6)
})

test_that("V6: under saturation every anchor gives the same likelihood", {
  ## mu intercept-only, nuisance intercept-only: the nuisance is saturated
  ## relative to mu's predictor, so the anchors are a bijection.
  d <- sim_kw(n = 400)
  lls <- vapply(c("beta", "alpha"), function(a) {
    gkwqreg(y ~ 1, data = d, tau = 0.5, family = "kw", anchor = a)$loglik
  }, numeric(1))
  expect_equal(unname(diff(lls)), 0, tolerance = 1e-6)
})

test_that("V7: with a covariate and homogeneous nuisance the anchors DIVERGE", {
  ## Not a defect: the two specify genuinely different models. Anchoring on beta
  ## asserts "alpha constant, beta_i = f(mu_i, alpha)"; anchoring on alpha
  ## asserts the reverse. This test exists to stop anyone "fixing" the
  ## divergence away.
  d <- sim_kw(n = 500)
  fb <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", anchor = "beta")
  fa <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", anchor = "alpha")
  expect_gt(abs(fb$loglik - fa$loglik), 1)
  expect_equal(fb$npar, fa$npar)   # same dimension, so NOT nested
})

test_that("saturating the nuisance restores agreement between anchors", {
  ## The practical recommendation, as a test: regress the nuisance on the same
  ## covariates as mu and the coefficient of interest stops depending on the
  ## anchor.
  d <- sim_kw(n = 500)
  fb <- gkwqreg(y ~ x | x, data = d, tau = 0.5, family = "kw", anchor = "beta")
  fa <- gkwqreg(y ~ x | x, data = d, tau = 0.5, family = "kw", anchor = "alpha")
  expect_equal(unname(fb$coef_list$mu[2]), unname(fa$coef_list$mu[2]),
               tolerance = 0.05)
})

test_that("anova() refuses to compare different taus or anchors", {
  d <- sim_kw(n = 200)
  f1 <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  f2 <- gkwqreg(y ~ x, data = d, tau = 0.6, family = "kw")
  f3 <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", anchor = "alpha")
  expect_error(anova(f1, f2), "separate likelihood")
  expect_error(anova(f1, f3), "non-nested")
})
