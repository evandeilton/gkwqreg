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

test_that("check_crossing() is right at exactly two tau levels (regression)", {
  ## Bug: when Q has exactly two columns, diff() applied row-wise degenerates
  ## to a scalar per row, apply(Q, 1, diff) simplifies to a plain vector (not
  ## a matrix), and t() of that vector yields a 1 x n matrix instead of the
  ## intended n x 1 -- the wrong orientation. It slips past
  ## `if (is.null(dim(D))) ...` because a 1 x n matrix already has a dim().
  ## Two levels is not a hypothetical edge case: it is the minimum this
  ## function itself accepts (see the `length(tv) < 2L` guard above).
  set.seed(1)
  n <- 120
  x <- runif(n, -1, 1)
  mu <- plogis(0.2 + 1.2 * x)
  y <- rbeta(n, 4 * mu, 4 * (1 - mu))
  d <- data.frame(y = y, x = x)

  fits2 <- gkwqreg(y ~ x | x, data = d, tau = c(0.90, 0.95), family = "kw")
  ## Extrapolated grid: the same device the check_crossing() examples use to
  ## surface real crossings between independently fitted levels.
  grid <- data.frame(x = seq(-3, 3, length.out = 201))

  cr <- check_crossing(fits2, newdata = grid)
  expect_equal(cr$mode, "separate")
  expect_equal(ncol(cr$Q), 2L)

  ## Ground truth computed by hand from the (bug-free) Q matrix, exactly as
  ## the audit prescribes: Q[, 2] - Q[, 1] < -tol.
  Dm <- cr$Q[, 2] - cr$Q[, 1]
  violm <- Dm < -0  # tol = 0, the default
  rows_manual <- which(violm)

  ## Not a vacuous scenario: this grid is known to cross at these levels.
  expect_true(length(rows_manual) > 0L)

  expect_equal(cr$n_crossing, length(rows_manual))
  expect_equal(cr$frac, length(rows_manual) / nrow(cr$Q))
  expect_equal(cr$which, rows_manual)
  expect_equal(cr$worst, max(-Dm[violm]))
  expect_equal(nrow(cr$pairs), 1L)
  expect_equal(cr$pairs$tau_lo, 0.90)
  expect_equal(cr$pairs$tau_hi, 0.95)
  expect_equal(cr$pairs$n, length(rows_manual))

  ## attr(rearrange(...), "crossing") calls check_crossing() internally and
  ## must inherit the same, now-correct, count -- not the pre-fix one.
  R <- rearrange(fits2, newdata = grid)
  crR <- attr(R, "crossing")
  expect_equal(crR$n_crossing, cr$n_crossing)
  expect_true(crR$n_crossing > 0L)
  expect_true(all(apply(R, 1, function(r) all(diff(r) >= -1e-12))))
})

test_that("check_crossing() at three-plus levels still matches the manual count", {
  ## Same scenario, full 19-level grid: guards against the two-level fix
  ## having disturbed the case that already worked.
  set.seed(1)
  n <- 120
  x <- runif(n, -1, 1)
  mu <- plogis(0.2 + 1.2 * x)
  y <- rbeta(n, 4 * mu, 4 * (1 - mu))
  d <- data.frame(y = y, x = x)

  fits <- gkwqreg(y ~ x | x, data = d, tau = seq(0.05, 0.95, by = 0.05),
                  family = "kw")
  grid <- data.frame(x = seq(-3, 3, length.out = 201))
  cr <- check_crossing(fits, newdata = grid)

  m <- ncol(cr$Q)
  expect_true(m >= 3L)
  Dm <- cr$Q[, -1L, drop = FALSE] - cr$Q[, -m, drop = FALSE]
  violm <- Dm < 0
  rows_manual <- which(apply(violm, 1L, any))

  expect_true(length(rows_manual) > 0L)
  expect_equal(cr$n_crossing, length(rows_manual))
  expect_equal(cr$which, rows_manual)
  expect_equal(nrow(cr$pairs), m - 1L)
  ## data.frame() repurposes colSums(violm)'s names as rownames(cr$pairs) (see
  ## the tau-labelled row names in the check_crossing() examples), so the `n`
  ## column itself comes back unnamed -- compare values only.
  expect_equal(cr$pairs$n, unname(colSums(violm)))
})

test_that("check_crossing() validates tol", {
  d <- sim_kw(n = 150)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  expect_error(check_crossing(f, tol = -1), "non-negative")
  expect_error(check_crossing(f, tol = c(0, 1)), "non-negative")
  expect_error(check_crossing(f, tol = "a"), "non-negative")
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

test_that("quantile_process() refuses a single quantile level", {
  ## A process is a curve over tau; one level does not define one. Previously
  ## this failed deep inside with an opaque
  ## "length of 'dimnames' [2] not equal to array extent" error.
  d <- sim_kw(n = 150)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  expect_error(quantile_process(f, taus = 0.5), "at least two")
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

test_that("compare_families refuses a family override through ...", {
  ## cl$family <- fm was set per row, then silently clobbered by
  ## `for (nm in names(extra)) cl[[nm]] <- extra[[nm]]` whenever `family` was
  ## also passed in `...` -- every row was refit under the override while the
  ## `family` column kept reporting the one that was intended.
  d <- sim_kw(n = 100)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  expect_error(compare_families(f, families = c("kw", "ekw"), family = "beta"),
               "cannot be passed through")
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

test_that("plot.gkwqregs() returns its argument invisibly", {
  ## Its documented @return is shared with plot.gkwq_process: "x, invisibly."
  ## The previous body, `plot(quantile_process(x), ...)`, returned whatever
  ## plot.gkwq_process() handed back -- the derived "gkwq_process" object, not
  ## the "gkwqregs" container `x` the doc promises.
  skip_if_not(capabilities("png"))
  d <- sim_kw(n = 150)
  fits2 <- gkwqreg(y ~ x, data = d, tau = c(0.25, 0.75), family = "kw")
  pf <- tempfile(fileext = ".png")
  grDevices::png(pf); on.exit(unlink(pf), add = TRUE)
  out <- plot(fits2)
  grDevices::dev.off()
  expect_identical(out, fits2)
})

test_that("panel 6 degrades gracefully, rather than crashing, without a hessian", {
  ## vcov(x) errors when the fit skipped optimHess(), so D in panel 6 used to
  ## come back entirely NA. plot(..., type = "h") on all-NA data raised "need
  ## finite 'ylim' values" instead of drawing the "unavailable" panel the docs
  ## promise for a fit obtained with gkwq_control(hessian = FALSE).
  skip_if_not(capabilities("png"))
  d <- sim_kw(n = 150)
  f <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw",
               control = gkwq_control(hessian = FALSE))
  expect_true(is.null(f$vcov))
  pf <- tempfile(fileext = ".png")
  grDevices::png(pf); on.exit(unlink(pf), add = TRUE)
  expect_silent(plot(f, which = 6))
  expect_silent(plot(f, which = 1:6, nsim = 10))
  grDevices::dev.off()
})
