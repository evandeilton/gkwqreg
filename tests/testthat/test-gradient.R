## V8, V9, V10 -- differentiation. V9 and V10 are the tests that would have
## caught the two AD defects found during the design study.

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

## ---------------------------------------------------------------------------
## The two atomics in inst/include/gkwq_atomic.hpp:
##
##   qbeta_safe (tau, p, q)  -> z  with  I_z(p,q) = tau
##   gamma_solve(mu,  q, tau) -> p  with  I_mu(p,q) = tau
##
## Nothing here may be satisfiable by stats::pbeta alone. An earlier version of
## this block called only stats::pbeta and numDeriv::grad(stats::pbeta), which
## made it pass with the package unloaded -- it compared pbeta with itself and
## left the header, which exists solely because TMB's own qbeta derivatives are
## wrong, with no value test at all. Each test below reads a quantity the
## COMPILED tape produced, via obj$report(), and confronts it with an R route
## that shares no code with the header.
## ---------------------------------------------------------------------------

test_that("gamma_solve inverts the incomplete beta: I_mu(gamma, delta+1) = tau", {
  ## The `beta` family's anchor has no closed form, so gamma is whatever the
  ## gamma_solve atomic returned. Its defining property is checkable directly
  ## against stats::pbeta, which shares no code with the header. If the bracket
  ## or the root finder in gamma_solve_double drifts, this is what catches it.
  d <- sim_kw(n = 120, seed = 5)
  for (tau in c(0.1, 0.5, 0.9)) {
    f <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = tau, family = "beta"))
    pv <- f$parameter_vectors            # from obj$report(): the compiled side
    expect_true(all(is.finite(pv$gamma)),
                label = sprintf("gamma is finite at tau = %.1f", tau))
    back <- stats::pbeta(fitted(f), pv$gamma, pv$delta + 1)
    expect_equal(back, rep(tau, nrow(pv)), tolerance = 1e-7,
                 info = sprintf("tau = %.1f", tau))
  }
})

test_that("the atomics agree with stats::qbeta and uniroot on every anchor", {
  ## The tape reconstructs the anchored parameter in C++, through qbeta_safe
  ## (z_mode 3) or gamma_solve (family `beta`). .gkwq_reconstruct() does the
  ## same job in R through stats::qbeta and uniroot over stats::pbeta. The two
  ## share no implementation, so agreement is evidence about the atomic and not
  ## about a single code path evaluated twice.
  recon <- gq(".gkwq_reconstruct"); info <- gq(".gkwq_family_info")
  d <- sim_kw(n = 100, seed = 9)
  for (fam in ALL_FAMILIES) {
    for (anc in gkwq_anchors(fam)) {
      for (tau in c(0.25, 0.75)) {
        lab <- sprintf("family %s, anchor %s, tau %.2f", fam, anc, tau)
        f <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = tau,
                                      family = fam, anchor = anc))
        pv <- f$parameter_vectors
        r <- recon(fitted(f), tau, info(fam, anc), as.list(pv))
        expect_equal(pv[[anc]], as.numeric(r[, anc]), tolerance = 1e-6,
                     info = lab)
      }
    }
  }
})

test_that("the AD gradient matches numDeriv on every anchor branch", {
  ## V9 only ever built objects at each family's DEFAULT anchor, which leaves
  ## `case 2: // eliminate alpha` in src/gkwqreg.cpp unreached even though
  ## alpha is admissible in five of the seven families. Every branch of the
  ## anchor switch is exercised here, away from tau = 0.5 as well.
  d <- sim_kw(n = 80, seed = 3)
  for (fam in ALL_FAMILIES) {
    for (anc in gkwq_anchors(fam)) {
      for (tau in c(0.2, 0.8)) {
        lab <- sprintf("family %s, anchor %s, tau %.1f", fam, anc, tau)
        obj <- make_obj(fam, d, tau = tau, anchor = anc)
        p <- obj$par; p[] <- 0
        g_ad <- as.numeric(obj$gr(p))
        g_num <- numDeriv::grad(obj$fn, p)
        expect_true(all(is.finite(g_ad)), label = paste("finite AD gradient,", lab))
        expect_lt(max(abs(g_ad - g_num)) / max(1, max(abs(g_num))), 1e-6,
                  label = paste("gradient,", lab))
      }
    }
  }
})
