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

test_that("an offset() term is resampled with the rows it belongs to", {
  ## The model frame stores an offset() term under its deparsed call --
  ## "offset(zz)", "offset(log(v))" -- which no formula can resolve on a refit.
  ## The bootstrap therefore used to re-evaluate the expression against the
  ## caller's environment, reaching past the resampled rows to the ORIGINAL
  ## variable, so replicate row i carried observation i's offset whatever row it
  ## had drawn. That produced replicates centred far from the estimate while
  ## every one of them reported convergence; when the variable happened not to
  ## be in scope it failed instead, under a message blaming the fit's stability.
  ##
  ## The invariant is that the two ways of supplying the same offset -- the
  ## argument and the term -- must bootstrap alike, since they fit alike.
  set.seed(77); n <- 150
  x <- stats::runif(n); zz <- stats::rnorm(n, 0, 0.6)
  mu <- stats::plogis(0.2 + 0.8 * x + zz)
  y <- pmin(pmax(stats::qbeta(stats::runif(n), 2, 3) * .15 + mu * .85,
                 1e-6), 1 - 1e-6)
  d <- data.frame(y = y, x = x, zz = zz, v = exp(zz))

  f_fml <- gkwqreg(y ~ x + offset(zz), data = d, tau = 0.5, family = "kw")
  f_arg <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", offset = zz)
  expect_equal(f_fml$loglik, f_arg$loglik, tolerance = 1e-8)

  b_fml <- suppressWarnings(gkwq_boot(f_fml, R = 40, seed = 5))
  b_arg <- suppressWarnings(gkwq_boot(f_arg, R = 40, seed = 5))
  expect_equal(b_fml$n_ok, b_arg$n_ok)
  expect_equal(unname(colMeans(b_fml$replicates, na.rm = TRUE)),
               unname(colMeans(b_arg$replicates, na.rm = TRUE)),
               tolerance = 1e-6)
  ## And it centres on the estimate. The broken version missed by 1.67 here.
  expect_lt(max(abs(colMeans(b_fml$replicates, na.rm = TRUE) - coef(f_fml))), 0.2)

  ## An offset() holding an expression rather than a bare name: the model frame
  ## column is "offset(log(v))", and log(v) must never be evaluated again.
  f_exp <- gkwqreg(y ~ x + offset(log(v)), data = d, tau = 0.5, family = "kw")
  expect_equal(f_exp$loglik, f_fml$loglik, tolerance = 1e-8)
  b_exp <- suppressWarnings(gkwq_boot(f_exp, R = 40, seed = 5))
  expect_equal(b_exp$n_ok, b_fml$n_ok)
  expect_equal(unname(colMeans(b_exp$replicates, na.rm = TRUE)),
               unname(colMeans(b_fml$replicates, na.rm = TRUE)),
               tolerance = 1e-6)

  ## An offset() on a nuisance part is reached too: the argument only ever
  ## attaches to mu, so this route has no substitute.
  f_nui <- gkwqreg(y ~ x | offset(zz), data = d, tau = 0.5, family = "kw")
  b_nui <- suppressWarnings(gkwq_boot(f_nui, R = 40, seed = 5))
  expect_gt(b_nui$n_ok, 30L)
  expect_lt(max(abs(colMeans(b_nui$replicates, na.rm = TRUE) - coef(f_nui))), 0.3)
})
