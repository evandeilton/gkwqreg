## V12 -- cross-validation against the prior art. The kw slice of this package
## IS the Kumaraswamy quantile regression of Mitnik and Baek (2013), which
## unitquantreg implements as its `kum` family. Agreement there is the strongest
## external check available.

test_that("V12: family 'kw' agrees with unitquantreg's 'kum'", {
  skip_if_not_installed("unitquantreg")
  d <- sim_kw(n = 400)
  ours <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  theirs <- unitquantreg::unitquantreg(y ~ x, data = d, tau = 0.5,
                                       family = "kum")
  expect_equal(unname(coef(ours)[1:2]), unname(coef(theirs)[1:2]),
               tolerance = 1e-3)
})

test_that("the fit is comparable to quantreg on the check loss it targets", {
  skip_if_not_installed("quantreg")
  d <- sim_kw(n = 400)
  ours <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  rq <- quantreg::rq(y ~ x, data = d, tau = 0.5)
  e <- d$y - stats::fitted(rq)
  loss_rq <- mean(e * (0.5 - (e < 0)))
  ## Under a correctly specified family the parametric fit should not be worse
  ## in-sample by any meaningful margin.
  expect_lt(pinball(ours), loss_rq * 1.05)
})

test_that("it runs on a real bounded data set", {
  skip_if_not_installed("betareg")
  data("GasolineYield", package = "betareg", envir = environment())
  f <- gkwqreg(yield ~ temp, data = GasolineYield, tau = 0.9, family = "kw")
  expect_true(all(fitted(f) > 0 & fitted(f) < 1))
  expect_true(is.finite(logLik(f)))
})
