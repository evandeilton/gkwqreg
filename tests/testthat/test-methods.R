## V14 -- every S3 method runs on every implemented family.

test_that("V14: the full S3 surface runs on all seven families", {
  d <- sim_kw(n = 150)
  for (fam in ALL_FAMILIES) {
    f <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = 0.5, family = fam))
    lab <- sprintf("family %s", fam)
    expect_type(coef(f), "double")
    expect_true(is.matrix(vcov(f)), info = lab)
    expect_s3_class(logLik(f), "logLik")
    expect_equal(nobs(f), 150L, info = lab)
    expect_type(AIC(f), "double")
    expect_type(BIC(f), "double")
    expect_length(fitted(f), 150L)
    expect_s3_class(family(f), "gkwq_family")
    expect_s3_class(terms(f), "terms")
    expect_true(is.matrix(model.matrix(f)), info = lab)
    expect_true(is.matrix(confint(f)), info = lab)
    expect_s3_class(summary(f), "summary.gkwqreg")
    expect_output(print(f), "Quantile Regression")
    expect_output(print(summary(f)), "anchor")
    expect_length(predict(f), 150L)
    expect_length(residuals(f), 150L)
    expect_type(getCall(f), "language")
    expect_silent(invisible(capture.output(print(family(f)))))
  }
})

test_that("every residual type is finite on every family", {
  d <- sim_kw(n = 150)
  types <- c("quantile", "cox-snell", "pearson", "deviance", "response",
             "check", "tau-sign")
  for (fam in ALL_FAMILIES) {
    f <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = 0.5, family = fam))
    for (ty in types) {
      r <- residuals(f, type = ty)
      expect_length(r, 150L)
      expect_true(mean(is.finite(r)) > 0.99,
                  info = sprintf("%s / %s: %d non-finite", fam, ty,
                                 sum(!is.finite(r))))
    }
  }
})

test_that("every predict type returns the documented shape", {
  d <- sim_kw(n = 100)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  expect_length(predict(f, type = "quantile"), 100L)
  expect_length(predict(f, type = "mu"), 100L)
  expect_s3_class(predict(f, type = "link"), "data.frame")
  expect_equal(dim(predict(f, type = "parameter")), c(100L, 5L))
  expect_length(predict(f, type = "mean"), 100L)
  expect_length(predict(f, type = "variance"), 100L)
  expect_equal(dim(predict(f, type = "quantile", tau = c(.1, .5, .9))),
               c(100L, 3L))
  expect_type(predict(f, type = "terms"), "list")
  nd <- data.frame(x = c(-1, 0, 1))
  expect_length(predict(f, newdata = nd), 3L)
})

test_that("predicting at a new tau reads the SAME fitted distribution", {
  ## Distinct from refitting at that level: this is one conditional
  ## distribution read at other probabilities, and it cannot cross.
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  Q <- predict(f, type = "quantile", tau = c(.1, .3, .5, .7, .9))
  expect_true(all(apply(Q, 1, function(r) all(diff(r) > 0))))
  expect_equal(unname(Q[, 3]), unname(fitted(f)), tolerance = 1e-8)
})

test_that("the sandwich covariance and estfun agree with the score identity", {
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  S <- estfun.gkwqreg(f)
  expect_equal(dim(S), c(200L, 3L))
  ## The per-observation scores must add up to the tape's own gradient. That is
  ## the real identity; "close to zero at the optimum" would only be testing how
  ## tightly nlminb happened to converge.
  ## Absolute, not relative: at the optimum both sides are near zero, where a
  ## relative tolerance measures nothing but how tightly nlminb converged.
  expect_lt(max(abs(colSums(S) + as.numeric(f$obj$gr(f$coefficients)))), 1e-6)
  Vs <- vcov(f, type = "sandwich")
  expect_equal(dim(Vs), c(3L, 3L))
  expect_true(all(diag(Vs) > 0))
})

test_that("marginal effects report effects on the quantile, not the odds", {
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  me <- marginal_effects(f)
  expect_s3_class(me, "gkwq_meff")
  ## d Q / d x = b * mean(mu (1 - mu)) under the logit link.
  mu <- fitted(f)
  expect_equal(me$table$effect[1],
               unname(coef(f)["mu:x"] * mean(mu * (1 - mu))), tolerance = 1e-8)
  expect_output(print(me), "Marginal effects")
})

test_that("update() and model.frame() round-trip", {
  d <- sim_kw(n = 120)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  f2 <- update(f, tau = 0.7)
  expect_equal(f2$tau, 0.7)
  expect_s3_class(model.frame(f), "data.frame")
})

test_that("weights and offsets reach the likelihood", {
  d <- sim_kw(n = 200)
  w <- rep(c(1, 2), length.out = 200)
  f0 <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  fw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", weights = w)
  expect_false(isTRUE(all.equal(coef(f0), coef(fw))))
  ## Doubling every weight doubles the log-likelihood exactly.
  f2 <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw",
                weights = rep(2, 200))
  expect_equal(f2$loglik, 2 * f0$loglik, tolerance = 1e-5)
})

test_that("simulate() draws samples whose tau-quantile is the fitted one", {
  d <- sim_kw(n = 300)
  f <- gkwqreg(y ~ x, data = d, tau = 0.75, family = "kw")
  s <- simulate(f, nsim = 20, seed = 7)
  expect_equal(dim(s), c(300L, 20L))
  expect_true(all(s > 0 & s < 1))
  ## Across replicates, the fraction below the fitted quantile must sit at tau.
  cov_emp <- mean(vapply(s, function(z) mean(z <= fitted(f)), numeric(1)))
  expect_lt(abs(cov_emp - 0.75), 4 * sqrt(0.75 * 0.25 / (300 * 20)))
  expect_equal(s, simulate(f, nsim = 20, seed = 7))   # seed is honoured
})

test_that("the bootstrap runs and agrees with the model covariance", {
  skip_on_cran()
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  b <- gkwq_boot(f, R = 60L, seed = 3)
  expect_s3_class(b, "gkwq_boot")
  expect_gt(b$n_ok, 50L)
  ## Bootstrap and model-based standard errors should be the same order of
  ## magnitude when the family is correctly specified.
  ratio <- b$se[c("mu:(Intercept)", "mu:x")] / f$se[c("mu:(Intercept)", "mu:x")]
  expect_true(all(ratio > 0.5 & ratio < 2))
  expect_output(print(b), "bootstrap")
  ci <- confint(f, method = "boot", R = 60L)
  expect_equal(dim(ci), c(3L, 2L))
})

test_that("lrtest() is anova() under the name lmtest users expect", {
  d <- sim_kw(n = 200)
  a <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  b <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "ekw")
  expect_equal(lrtest(a, b), anova(a, b))
})

test_that("profile intervals bracket the estimate, or report that they cannot", {
  skip_on_cran()
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  ci <- suppressWarnings(confint(f, parm = "mu:x", method = "profile"))
  expect_equal(dim(ci), c(1L, 2L))
  if (all(is.finite(ci))) {
    expect_lt(ci[1, 1], coef(f)[["mu:x"]])
    expect_gt(ci[1, 2], coef(f)[["mu:x"]])
    ## A profile interval should be in the same ballpark as the Wald one.
    w <- confint(f, parm = "mu:x", method = "wald")
    expect_lt(abs(diff(ci[1, ]) / diff(w[1, ]) - 1), 0.5)
  } else {
    ## Failure must be loud, not a silent NA.
    expect_warning(confint(f, parm = "mu:x", method = "profile"),
                   "cut-off")
  }
})

test_that("summary() still works when the response was not retained", {
  ## y = FALSE is a documented option, so summary() must report what it can
  ## rather than refusing to print a coefficient table it already holds.
  d <- sim_kw(n = 120)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", y = FALSE)
  s <- summary(f)
  expect_s3_class(s, "summary.gkwqreg")
  expect_true(is.na(s$coverage))
  expect_output(print(s), "unavailable")
  expect_error(residuals(f), "refit with y = TRUE")
})
