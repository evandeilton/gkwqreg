## V14 -- every S3 method runs on every implemented family.

test_that("V14: the full S3 surface runs on all seven families", {
  d <- sim_kw(n = 150)
  for (fam in ALL_FAMILIES) {
    ## Warnings are collected rather than suppressed. suppressWarnings() here
    ## hid that three of the seven families do not converge on this fixture,
    ## while expect_true(is.matrix(vcov(f))) passed on a covariance matrix that
    ## was not positive definite and a confint() containing NA. The surface has
    ## to run on all seven -- that is what V14 claims -- but the optimizer's
    ## verdict has to be visible, and where it did converge the covariance has
    ## to be usable. Which families struggle depends on the starting values, so
    ## the invariant is conditional rather than a list of names.
    warns <- character(0)
    f <- withCallingHandlers(
      gkwqreg(y ~ x, data = d, tau = 0.5, family = fam),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      })
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

    if (isTRUE(f$convergence == 0)) {
      ## A converged fit must carry a covariance one can actually use.
      ev <- eigen(vcov(f), symmetric = TRUE, only.values = TRUE)$values
      expect_gt(min(ev), 0)
      expect_false(anyNA(suppressWarnings(confint(f))))
      expect_true(all(is.finite(coef(f))), info = lab)
    } else {
      ## And a fit that did not converge must say so rather than pass quietly.
      expect_true(any(grepl("converge", warns, fixed = TRUE)),
                  label = paste("non-convergence is reported for", lab))
    }
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
      ## Was mean(is.finite(r)) > 0.99, which with n = 150 tolerates one
      ## non-finite residual: pure slack, and the only thing it could ever hide
      ## is the first regression. None of the seven families produces one.
      expect_true(all(is.finite(r)),
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

test_that("profile intervals are finite and bracket the estimate", {
  ## This used to return NA for every coefficient. TMB::tmbprofile() hands back
  ## a frame whose columns are named "value" and NA, with the PARAMETER grid
  ## under "value"; TMB's own confint method reads them the other way round.
  ## The interval is now computed here, and the earlier version of this test was
  ## permissive enough to pass while the feature was entirely broken.
  skip_on_cran()
  d <- sim_kw(n = 250)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")

  ci <- confint(f, method = "profile")
  expect_equal(dim(ci), c(length(coef(f)), 2L))
  expect_true(all(is.finite(ci)))                       # not NA -- the whole point
  expect_true(all(ci[, 1] < coef(f) & coef(f) < ci[, 2]))

  ## For a well-behaved model the profile and Wald intervals should agree
  ## closely; a large gap would mean one of them is wrong.
  w <- confint(f, method = "wald")
  expect_true(all(abs((ci[, 2] - ci[, 1]) / (w[, 2] - w[, 1]) - 1) < 0.25))
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

test_that("predict() with a factor in newdata is silent and respects levels", {
  ## model.frame was being handed the levels for EVERY part, so it warned
  ## "variable 'g' is not a factor" for parts that do not use g -- noise on
  ## essentially every predict(newdata=) call with a factor in the model.
  d <- sim_kw(n = 200)
  d$g <- factor(rep(c("a", "b"), length.out = 200))
  f <- gkwqreg(y ~ x + g | x, data = d, tau = 0.5, family = "kw")

  nd <- data.frame(x = c(-1, 0, 1),
                   g = factor("a", levels = c("a", "b")))
  expect_silent(p <- predict(f, newdata = nd))
  expect_length(p, 3L)
  expect_true(all(p > 0 & p < 1))

  ## Levels must still be enforced: the two groups differ, and an unseen level
  ## is still an error rather than a silent recode.
  pa <- predict(f, newdata = data.frame(x = 0, g = factor("a", levels = c("a", "b"))))
  pb <- predict(f, newdata = data.frame(x = 0, g = factor("b", levels = c("a", "b"))))
  expect_false(isTRUE(all.equal(pa, pb)))
  expect_error(predict(f, newdata = data.frame(x = 0, g = factor("z"))))
})

test_that("the documented sandwich integration is exactly what works", {
  skip_if_not_installed("sandwich")
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")

  ## sandwich::sandwich() must reproduce vcov(type = "sandwich") exactly.
  expect_equal(sandwich::sandwich(f), vcov(f, type = "sandwich"),
               ignore_attr = TRUE, tolerance = 1e-8)
  expect_true(is.matrix(sandwich::vcovCL(f, cluster = seq_len(nobs(f)))))

  ## vcovHC is documented as NOT working; pin that so the claim stays true.
  expect_error(sandwich::vcovHC(f), "model.matrix|estfun")
})

test_that("deviance residuals use a genuine supremum, not a relocation", {
  ## The reference must be sup_m l_i(m), not l_i(y_i). The density is generally
  ## maximised at some m != y, so relocating the quantile onto the observation
  ## can give a "saturated" value BELOW the fitted one, which the square root
  ## then floors to zero. That silently zeroed 93 of 400 residuals.
  d <- sim_kw(n = 400)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  r <- residuals(f, type = "deviance")

  expect_length(r, 400L)
  expect_true(all(is.finite(r)))
  expect_equal(sum(r == 0), 0L)                       # nothing floored
  ## Sign must track which side of the fitted quantile the observation is on.
  expect_true(all(sign(r) == sign(f$y - fitted(f))))
})

test_that("a non-positive variance gives NA, never a standard error of zero", {
  ## When the information matrix is not positive definite some coefficient has
  ## no standard error. Flooring the variance at zero instead -- which
  ## summary() did through sqrt(pmax(diag(V), 0)) -- reports se = 0, hence
  ## z = +-Inf and p = 0: a fit too degenerate to carry an interval reads as
  ## overwhelmingly significant. Seen on a gkw fit whose whole diagonal went
  ## negative, where the fit reported NA for every coefficient and summary()
  ## reported p = 0 for every one of them.
  ##
  ## The variance is injected rather than fished out of a badly conditioned
  ## fit, so the case is identical on every platform and does not depend on
  ## where an optimizer happened to stop.
  f <- gkwqreg(y ~ x, data = sim_kw(n = 150), tau = 0.5, family = "kw")

  ## The two paths through the package must agree in the ordinary case.
  expect_equal(unname(summary(f)$coefficients[, "Std. Error"]),
               unname(f$se), tolerance = 1e-12)

  for (bad_var in c(-1e-8, 0)) {
    bad <- f
    bad$vcov[1L, 1L] <- bad_var
    tab <- summary(bad)$coefficients
    lab <- sprintf("variance = %g", bad_var)
    expect_true(is.na(tab[1L, "Std. Error"]), label = paste("se is NA,", lab))
    expect_true(is.na(tab[1L, "z value"]), label = paste("z is NA,", lab))
    expect_true(is.na(tab[1L, "Pr(>|z|)"]), label = paste("p is NA,", lab))
    ## and the coefficients with a usable variance are left alone
    expect_true(all(is.finite(tab[-1L, "Std. Error"])),
                label = paste("others unaffected,", lab))
  }
})

test_that("summary() weights the loss it compares against the weighted pinball", {
  ## object$pinball comes from the tape, where it is weighted and divided by
  ## sum(w). summary() built its null loss and its coverage with plain means,
  ## so Pseudo-R1 = 1 - pinball/pin0 had a weighted numerator over an
  ## unweighted denominator whenever weights were supplied.
  set.seed(17); n <- 200; x <- stats::runif(n)
  d <- data.frame(y = stats::qbeta(stats::runif(n), 2 + x, 3), x = x,
                  w = rep(c(1, 5), length.out = n))

  fw <- gkwqreg(y ~ x, d, tau = 0.7, family = "kw", weights = w)
  sw <- summary(fw)
  expect_equal(sw$coverage,
               sum(d$w * (fw$y <= fitted(fw))) / sum(d$w), tolerance = 1e-12)

  null_q <- stats::quantile(fw$y, probs = 0.7, names = FALSE, type = 7)
  e0 <- fw$y - null_q
  pin0 <- sum(d$w * e0 * (0.7 - (e0 < 0))) / sum(d$w)
  expect_equal(sw$pseudo_r1, 1 - fw$pinball / pin0, tolerance = 1e-12)

  ## Unit weights must reproduce the plain means exactly.
  f0 <- gkwqreg(y ~ x, d, tau = 0.7, family = "kw")
  s0 <- summary(f0)
  expect_equal(s0$coverage, mean(f0$y <= fitted(f0)), tolerance = 1e-12)
})

test_that("a frequency weight counts as the replication it is", {
  ## Weights multiply the per-observation log-density, so w = 3 throughout is
  ## the same likelihood as replicating every row three times -- the suite
  ## already asserts the log-likelihood scales exactly. The BIC penalty was not
  ## following: it used the row count while the log-likelihood was on the
  ## sum(w) scale, so a weighted fit was penalised as though it had seen a
  ## third of its evidence. Replication is the definition, so it is the oracle.
  set.seed(4); n <- 300; x <- stats::runif(n)
  d <- data.frame(y = stats::qbeta(stats::runif(n), 2 + x, 3), x = x)

  fw <- gkwqreg(y ~ x, d, tau = 0.5, family = "kw", weights = rep(3, n))
  fr <- gkwqreg(y ~ x, d[rep(seq_len(n), each = 3), ], tau = 0.5, family = "kw")

  expect_equal(fw$loglik, fr$loglik, tolerance = 1e-5)
  expect_equal(unname(coef(fw)), unname(coef(fr)), tolerance = 1e-6)
  expect_equal(BIC(fw), BIC(fr), tolerance = 1e-4)
  expect_equal(fw$nobs_eff, fr$nobs_eff)
  expect_equal(unname(attr(logLik(fw), "nobs")), 3 * n)

  ## nobs() keeps counting rows: it is what sizes the residuals and the scores.
  expect_equal(nobs(fw), n)
  expect_length(residuals(fw), n)

  ## An unweighted fit is untouched, and there the two scales coincide.
  f1 <- gkwqreg(y ~ x, d, tau = 0.5, family = "kw")
  expect_equal(f1$nobs_eff, nobs(f1))
  expect_equal(BIC(f1), -2 * f1$loglik + log(n) * f1$npar, tolerance = 1e-8)
})
