## Formula handling: transformed terms and offsets. Both of these were broken.

test_that("transformed terms in a formula work", {
  ## The per-part offset was read by rebuilding a model frame from that part's
  ## terms, which re-evaluates predvars against a frame whose columns are
  ## already named "log(z)" rather than "z". Every transformed term therefore
  ## died with "object 'z' not found" -- log(), I(), poly(), splines, all of it.
  d <- sim_kw(n = 250)
  d$z <- runif(250, 1, 5)
  d$g <- factor(rep(c("a", "b"), length.out = 250))

  for (fml in list(y ~ log(z), y ~ I(x^2), y ~ poly(x, 2),
                   y ~ x + log(z) + g, y ~ x | log(z))) {
    f <- gkwqreg(fml, data = d, tau = 0.5, family = "kw")
    expect_true(is.finite(f$loglik), info = deparse(fml))
    expect_true(all(fitted(f) > 0 & fitted(f) < 1), info = deparse(fml))
  }

  ## and prediction must survive the transformation too
  f <- gkwqreg(y ~ log(z), data = d, tau = 0.5, family = "kw")
  p <- predict(f, newdata = data.frame(z = c(2, 4)))
  expect_length(p, 2L)
  expect_true(all(p > 0 & p < 1))
})

test_that("an offset enters the likelihood exactly once", {
  ## offset() terms were collected per part AND again from the model frame, so
  ## the conditional quantile received the offset twice.
  d <- sim_kw(n = 250)
  d$w <- runif(250, -0.5, 0.5)          # must VARY: a constant offset is
                                        # absorbed by the intercept
  a <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  b <- gkwqreg(y ~ x + offset(w), data = d, tau = 0.5, family = "kw")
  expect_false(isTRUE(all.equal(a$loglik, b$loglik)))

  ## The stored linear predictor must equal X b + offset, exactly once.
  eta <- as.numeric(cbind(1, d$x) %*% b$coef_list$mu) + d$w
  expect_equal(eta, b$linear.predictors$mu, tolerance = 1e-12)

  ## The `offset =` argument and an offset() term must agree.
  cc <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", offset = w)
  expect_equal(cc$loglik, b$loglik, tolerance = 1e-6)

  ## An offset on a non-quantile part is honoured as well.
  e <- gkwqreg(y ~ x | offset(w), data = d, tau = 0.5, family = "kw")
  expect_true(is.finite(e$loglik))
  expect_false(isTRUE(all.equal(e$loglik, a$loglik)))

  ## Length mismatches are refused rather than recycled.
  expect_error(gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw",
                       offset = rep(0, 3)), "length")
})
