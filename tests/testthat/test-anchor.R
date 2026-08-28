## V1, V2, V3 -- the reparametrization itself, checked against gkwdist rather
## than against our own TMB template. An oracle that shares an implementation
## with the thing under test passes by construction.

test_that("V1: the closed-form quantile matches gkwdist over a parameter grid", {
  g <- expand.grid(t = c(.01, .1, .25, .5, .75, .9, .99),
                   a = c(.5, 1, 2.5), b = c(.7, 1.5, 3),
                   gm = c(.8, 1, 2), d = c(0, .5, 2), L = c(.6, 1, 2))
  got <- gkwq_quantile(g$t, g$a, g$b, g$gm, g$d, g$L)
  want <- gkwdist::qgkw(g$t, g$a, g$b, g$gm, g$d, g$L)
  expect_lt(max(abs(got - want)), 1e-12)
})

test_that("V3: every anchor is finite and strictly positive, all families", {
  info <- gq(".gkwq_family_info")
  recon <- gq(".gkwq_reconstruct")
  for (fam in ALL_FAMILIES) {
    for (anc in gkwq_anchors(fam)) {
      spec <- info(fam, anc)
      free <- setdiff(spec$full, anc)
      g <- expand.grid(mu = c(.02, .1, .3, .5, .7, .9, .98),
                       tau = c(.05, .25, .5, .75, .95),
                       v1 = c(.5, 1.5, 2.5), v2 = c(.7, 2))
      pars <- list()
      if (length(free) >= 1) pars[[free[1]]] <- g$v1
      if (length(free) >= 2) pars[[free[2]]] <- g$v2
      if (length(free) >= 3) pars[[free[3]]] <- 1.3
      if (length(free) >= 4) pars[[free[4]]] <- 0.6
      P <- recon(g$mu, g$tau, spec, pars)
      col <- match(anc, c("alpha", "beta", "gamma", "delta", "lambda"))
      expect_true(all(is.finite(P[, col])),
                  info = sprintf("%s / %s: non-finite anchor", fam, anc))
      expect_true(all(P[, col] > 0),
                  info = sprintf("%s / %s: non-positive anchor", fam, anc))
    }
  }
})

test_that("V2: the anchored parameter round-trips through the quantile function", {
  info <- gq(".gkwq_family_info")
  recon <- gq(".gkwq_reconstruct")
  for (fam in ALL_FAMILIES) {
    for (anc in gkwq_anchors(fam)) {
      spec <- info(fam, anc)
      free <- setdiff(spec$full, anc)
      g <- expand.grid(mu = c(.05, .3, .5, .7, .95),
                       tau = c(.1, .5, .9),
                       v1 = c(.8, 2), v2 = c(1.2, 2.5))
      pars <- list()
      if (length(free) >= 1) pars[[free[1]]] <- g$v1
      if (length(free) >= 2) pars[[free[2]]] <- g$v2
      if (length(free) >= 3) pars[[free[3]]] <- 1.3
      if (length(free) >= 4) pars[[free[4]]] <- 0.6
      P <- recon(g$mu, g$tau, spec, pars)
      back <- gkwdist::qgkw(g$tau, P[, 1], P[, 2], P[, 3], P[, 4], P[, 5])
      expect_lt(max(abs(back - g$mu)), 1e-10)
    }
  }
})

test_that("the kw beta-solve reduces exactly to Mitnik and Baek (2013)", {
  ## The general formula must CONTAIN the established prior art as its
  ## two-parameter slice; that is the strongest validation the construction has.
  anchor <- gq(".gkwq_anchor_value")
  mu <- c(.15, .35, .6, .85); tau <- c(.1, .25, .5, .9); a <- c(1.2, 1.8, 2.5, 3)
  general <- anchor("beta", mu, tau, alpha = a, gamma = 1, delta = 0,
                    lambda = 1, z_mode = 1L)
  mitnik_baek <- log1p(-tau) / log1p(-mu^a)
  expect_equal(general, mitnik_baek, tolerance = 1e-14)
})

test_that("the anchors overflow in OPPOSITE tails, and neither ever fails", {
  ## SPEC N6. This is why `anchor` cannot have a universally right default:
  ## the beta solve degrades as mu -> 0 and the alpha solve as mu -> 1, so the
  ## better choice depends on where the response mass lies. Both stay finite and
  ## strictly positive throughout, which is what keeps the optimizer safe.
  anchor <- gq(".gkwq_anchor_value")
  logdens <- gq(".gkwq_logdens")
  mus <- c(1e-8, 1e-6, 1e-3, .01, .5, .99, 1 - 1e-3, 1 - 1e-6, 1 - 1e-8)
  b <- anchor("beta", mus, 0.5, alpha = 2, gamma = 1, delta = 0, lambda = 1,
              z_mode = 1L)
  a <- anchor("alpha", mus, 0.5, beta = 2, gamma = 1, delta = 0, lambda = 1,
              z_mode = 1L)

  expect_true(all(is.finite(b) & b > 0))
  expect_true(all(is.finite(a) & a > 0))

  ## Opposite tails: beta is largest at the smallest mu, alpha at the largest.
  expect_equal(which.max(b), 1L)
  expect_equal(which.max(a), length(mus))

  ## And the log-density stays finite at both extremes.
  expect_true(all(is.finite(logdens(0.5, 2, b, 1, 0, 1, delta_is_zero = TRUE))))
  expect_true(all(is.finite(logdens(0.5, a, 2, 1, 0, 1, delta_is_zero = TRUE))))
})
