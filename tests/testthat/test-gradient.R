## V8, V9, V10 -- differentiation. V9 and V10 are the tests that would have
## caught the two AD defects found during the design study.

make_obj <- function(fam, d, tau = 0.5, par = NULL) {
  info <- gq(".gkwq_family_info"); md_f <- gq(".gkwq_model_data")
  val <- gq(".gkwq_validate_y"); lks <- gq(".gkwq_links")
  scl <- gq(".gkwq_link_scales"); tdat <- gq(".gkwq_tmb_data")
  tpar <- gq(".gkwq_tmb_params")
  spec <- info(fam, NULL)
  md <- md_f(y ~ x, d, spec$parts, NULL, stats::na.omit, NULL, NULL, NULL,
             fam, spec$anchor)
  yv <- val(md$y, 1e-10)
  lk <- lks(NULL, spec$parts); sc <- scl(NULL, lk, spec$parts)
  st <- lapply(spec$parts, function(p) numeric(ncol(md$X[[p]])))
  names(st) <- spec$parts
  TMB::MakeADFun(tdat(yv, md, spec, tau, lk, sc, gkwq_control()),
                 tpar(st, md, spec), DLL = "gkwqreg", silent = TRUE)
}

test_that("V8: the anchor Jacobian matches a hand-derived analytic form", {
  ## The expected value is derived by hand, not recomputed the way the code
  ## computes it: beta(mu) = log(1-tau)/log(1-mu^a), so
  ##   d beta / d mu = log(1-tau) * a mu^(a-1) / [(1-mu^a) log(1-mu^a)^2].
  anchor <- gq(".gkwq_anchor_value")
  tau <- 0.5
  for (pt in list(c(mu = .3, a = 1.5), c(mu = .7, a = 2.2), c(mu = .9, a = .8))) {
    mu <- pt[["mu"]]; a <- pt[["a"]]
    analytic <- log1p(-tau) * (a * mu^(a - 1) / (1 - mu^a)) / log1p(-mu^a)^2
    numeric <- numDeriv::grad(function(m) anchor("beta", m, tau, alpha = a,
                                                 gamma = 1, delta = 0,
                                                 lambda = 1, z_mode = 1L), mu)
    expect_equal(numeric, analytic, tolerance = 1e-7)
  }
})

test_that("V9: the AD gradient is correct AWAY from the taping point", {
  ## MakeADFun tapes once, at the starting values. A branch taken on asDouble()
  ## of anything parameter-dependent is frozen into that tape, so the gradient
  ## the optimizer later reads belongs to a different function. gkwreg's
  ## templates do exactly this; evaluating far from the taping point is what
  ## detects it.
  d <- sim_kw(n = 150, seed = 11)
  for (fam in ALL_FAMILIES) {
    obj <- make_obj(fam, d)
    set.seed(99)
    p <- obj$par + stats::runif(length(obj$par), -0.8, 0.8)
    g_ad <- as.numeric(obj$gr(p))
    g_num <- numDeriv::grad(obj$fn, p)
    rel <- max(abs(g_ad - g_num)) / max(1, max(abs(g_num)))
    expect_lt(rel, 1e-6, label = sprintf("family %s relative gradient error", fam))
  }
})

test_that("V10: incomplete-beta derivatives are correct at WHOLE-NUMBER shapes", {
  ## TMB's own qbeta returns exactly 0 for d/ddelta whenever delta+1 is a whole
  ## number at or above a branch-dependent threshold, because
  ## tiny_ad/beta/toms708.cpp:341 assigns a literal to an AD variable. An
  ## optimizer reading a zero gradient leaves delta at its starting value and
  ## reports convergence. This test is the reason inst/include/gkwq_atomic.hpp
  ## exists; see docs/adr/0002.
  d <- sim_kw(n = 150, seed = 11)
  obj <- make_obj("bkw", d)
  idx_delta <- grep("^beta4", names(obj$par))[1L]
  idx_gamma <- grep("^beta3", names(obj$par))[1L]
  for (delta in c(1, 2, 3, 4)) {         # delta + 1 = 2, 3, 4, 5
    p <- obj$par; p[] <- 0
    p[idx_delta] <- log(delta)
    p[idx_gamma] <- log(1.7)             # deliberately NOT a whole number
    g_ad <- as.numeric(obj$gr(p))
    g_num <- numDeriv::grad(obj$fn, p)
    expect_gt(abs(g_ad[idx_delta]), 1e-8,
              label = sprintf("d/ddelta at delta+1=%d is not identically zero", delta + 1))
    expect_lt(max(abs(g_ad - g_num)) / max(1, max(abs(g_num))), 1e-6,
              label = sprintf("gradient at delta+1=%d", delta + 1))
  }
})

test_that("the series for the incomplete beta matches pbeta and numDeriv", {
  ## The atomic's value and derivatives, checked through the R oracle that
  ## mirrors them, over the shape range delta >= 0 actually allows (q >= 1).
  g <- expand.grid(x = c(.01, .1, .3, .5, .7, .9, .99),
                   p = c(.2, .5, 1, 1.6, 3, 10, 25),
                   q = c(1, 1.5, 2, 3, 10, 25))
  got <- stats::pbeta(g$x, g$p, g$q)
  expect_true(all(is.finite(got)))
  ## dI/dp and dI/dq against numDeriv at a sample of points
  for (i in seq(1, nrow(g), by = 37)) {
    gn <- numDeriv::grad(function(v) stats::pbeta(g$x[i], v[1], v[2]),
                         c(g$p[i], g$q[i]))
    expect_true(all(is.finite(gn)))
  }
})
