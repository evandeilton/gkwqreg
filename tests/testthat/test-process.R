## The quantile-regression-specific surface.

test_that("a vector tau returns a container of independent fits", {
  d <- sim_kw(n = 200)
  fits <- gkwqreg(y ~ x, data = d, tau = c(0.25, 0.5, 0.75), family = "kw")
  expect_s3_class(fits, "gkwqregs")
  expect_length(fits$fits, 3L)
  expect_equal(dim(coef(fits)), c(3L, 3L))
  expect_equal(dim(fitted(fits)), c(200L, 3L))
  ## Each level is its own likelihood, so a pooled logLik would be a lie.
  expect_error(logLik(fits), "one likelihood per tau")
})

test_that("quantiles from ONE fit cannot cross", {
  d <- sim_kw(n = 300)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  cr <- check_crossing(f, taus = seq(0.05, 0.95, by = 0.05))
  expect_s3_class(cr, "gkwq_crossing")
  expect_equal(cr$n_crossing, 0L)
  expect_equal(cr$mode, "implied")
  expect_output(print(cr), "not luck")
})

test_that("crossing across separately fitted levels is detected and fixable", {
  d <- sim_kw(n = 200)
  fits <- gkwqreg(y ~ x, data = d, tau = c(0.1, 0.3, 0.5, 0.7, 0.9),
                  family = "kw")
  cr <- check_crossing(fits)
  expect_equal(cr$mode, "separate")
  expect_true(cr$n_crossing >= 0L)
  R <- rearrange(fits)
  expect_true(all(apply(R, 1, function(r) all(diff(r) >= -1e-12))))
})

test_that("the quantile process collects coefficient paths", {
  d <- sim_kw(n = 200)
  qp <- quantile_process(gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw"),
                         taus = c(0.2, 0.4, 0.6, 0.8))
  expect_s3_class(qp, "gkwq_process")
  expect_equal(dim(qp$coef), c(3L, 4L))
  expect_true(all(qp$lower <= qp$coef & qp$coef <= qp$upper, na.rm = TRUE))
  expect_output(print(qp), "Quantile process")
})

test_that("pinball loss matches the check-loss definition", {
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.7, family = "kw")
  e <- f$y - fitted(f)
  expect_equal(pinball(f), mean(e * (0.7 - (e < 0))), tolerance = 1e-10)
  expect_equal(pinball(f), f$pinball, tolerance = 1e-10)
  expect_equal(pinball(f), mean(residuals(f, "check")), tolerance = 1e-10)
})

test_that("compare_families ranks families and Vuong compares anchors", {
  d <- sim_kw(n = 200)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  cmp <- suppressWarnings(compare_families(f, families = c("kw", "ekw", "beta")))
  expect_s3_class(cmp, "data.frame")
  expect_true(all(c("family", "AIC", "pinball") %in% names(cmp)))
  expect_false(is.unsorted(cmp$AIC, na.rm = TRUE))

  fa <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw", anchor = "alpha")
  vt <- vuong_test(f, fa)
  expect_s3_class(vt, "gkwq_vuong")
  expect_true(is.finite(vt$statistic))
  expect_output(print(vt), "Vuong")
})

test_that("plots run without error", {
  skip_if_not(capabilities("png"))
  d <- sim_kw(n = 150)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  pf <- tempfile(fileext = ".png")
  grDevices::png(pf); on.exit(unlink(pf), add = TRUE)
  expect_silent(plot(f, which = 1:6, nsim = 10))
  qp <- quantile_process(f, taus = c(0.25, 0.5, 0.75))
  expect_silent(plot(qp))
  grDevices::dev.off()
})
