## ---------------------------------------------------------------------------
## The reparametrization, in pure R.
##
## This file is deliberately a SECOND implementation of what src/gkwqreg.cpp
## does. It serves three roles:
##
##   1. the independent oracle the tests check the C++ against (tests never
##      compare TMB to TMB, which would pass by construction);
##   2. the engine behind predict(), residuals() and marginal_effects(), none
##      of which need a tape;
##   3. the feasibility sweep over starting values, which must run before the
##      tape exists.
##
## It is written from the mathematics in SPEC 1.2-1.3, not transcribed from the
## template, so a mistake in one is unlikely to be mirrored in the other.
## ---------------------------------------------------------------------------

## log(1 - exp(x)) for x < 0, Maechler (2012). The two-branch form is what
## keeps the whole cascade exact: gkwreg's `if (log_x < -30) return 0` shortcut
## loses up to 6.5e4 of log-likelihood once beta is large, which the beta anchor
## makes routine (SPEC N3).
.log1mexp <- function(x) {
  x <- pmin(x, -.Machine$double.xmin)
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

## log(z_tau). The elementary cases keep the incomplete beta out of the
## computation entirely for kw, ekw and kkw (SPEC 1.2).
.gkwq_log_z <- function(tau, gamma = 1, delta = 0, z_mode = 3L) {
  if (z_mode == 1L) {
    return(log(tau))
  }
  if (z_mode == 2L) {
    return(.log1mexp(log1p(-tau) / (delta + 1)))
  }
  stats::qbeta(tau, gamma, delta + 1, log.p = FALSE, lower.tail = TRUE) |>
    log()
}

#' The Generalized Kumaraswamy quantile function, in closed form
#'
#' `Q(tau) = {1 - [1 - z^(1/lambda)]^(1/beta)}^(1/alpha)` with
#' `z = qbeta(tau, gamma, delta + 1)`. Computed entirely in the log domain.
#'
#' @param tau Quantile level(s) in `(0,1)`.
#' @param alpha,beta,gamma,delta,lambda Distribution parameters.
#' @return A numeric vector of quantiles.
#' @examples
#' gkwq_quantile(c(0.1, 0.5, 0.9), alpha = 2, beta = 3)
#' @export
gkwq_quantile <- function(tau, alpha = 1, beta = 1, gamma = 1, delta = 0,
                          lambda = 1) {
  lz <- log(stats::qbeta(tau, gamma, delta + 1))
  lw <- .log1mexp(lz / lambda)            # log(1 - z^(1/lambda))
  lm <- .log1mexp(lw / beta)              # log(mu^alpha)
  exp(lm / alpha)
}

## ---------------------------------------------------------------------------
## The anchor solves (SPEC 1.3).
##
## Each is a ratio of two logarithms of quantities in (0,1), hence provably
## positive and finite for every mu, tau in (0,1) and every admissible value of
## the remaining parameters. That is what makes an unconstrained optimizer safe:
## no box constraint is ever needed.
## ---------------------------------------------------------------------------

.gkwq_anchor_value <- function(anchor, mu, tau, alpha = 1, beta = 1, gamma = 1,
                               delta = 0, lambda = 1, z_mode = 3L) {
  if (anchor == "gamma") {
    ## No closed form. I_mu(gamma, delta+1) is strictly decreasing in gamma and
    ## sweeps all of (0,1), so a root always exists.
    return(.gkwq_gamma_solve(mu, delta + 1, tau))
  }
  lz <- .gkwq_log_z(tau, gamma, delta, z_mode)
  lmu <- log(mu)
  lw <- .log1mexp(lz / lambda)            # log(1 - z^(1/lambda))
  switch(anchor,
    beta   = lw / .log1mexp(alpha * lmu),
    alpha  = .log1mexp(lw / beta) / lmu,
    lambda = lz / .log1mexp(beta * .log1mexp(alpha * lmu)),
    stop("unknown anchor: ", sQuote(anchor), call. = FALSE)
  )
}

## Vectorised bracketed root solve for the `beta` family's implicit anchor.
.gkwq_gamma_solve <- function(mu, q, tau) {
  n <- max(length(mu), length(q), length(tau))
  mu <- rep_len(mu, n); q <- rep_len(q, n); tau <- rep_len(tau, n)
  vapply(seq_len(n), function(i) {
    f <- function(lg) stats::pbeta(mu[i], exp(lg), q[i]) - tau[i]
    out <- try(stats::uniroot(f, lower = log(1e-8), upper = log(1e8),
                              tol = .Machine$double.eps^0.75)$root,
               silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else exp(out)
  }, numeric(1))
}

## Reconstruct the full (alpha, beta, gamma, delta, lambda) vector from the
## conditional quantile and the non-anchored parameters. Vectorised in every
## argument; returns an n x 5 matrix in gkwdist's parameter order.
.gkwq_reconstruct <- function(mu, tau, spec, pars) {
  n <- length(mu)
  get <- function(nm) {
    if (!is.null(pars[[nm]])) rep_len(pars[[nm]], n) else rep_len(spec$fixed_val[[
      match(nm, c("alpha", "beta", "gamma", "delta", "lambda"))]], n)
  }
  a <- get("alpha"); b <- get("beta"); g <- get("gamma")
  d <- get("delta"); L <- get("lambda")

  anchored <- .gkwq_anchor_value(spec$anchor, mu, tau, alpha = a, beta = b,
                                 gamma = g, delta = d, lambda = L,
                                 z_mode = spec$z_mode)
  switch(spec$anchor,
    beta   = b <- anchored,
    alpha  = a <- anchored,
    lambda = L <- anchored,
    gamma  = g <- anchored
  )
  cbind(alpha = a, beta = b, gamma = g, delta = d, lambda = L)
}

## ---------------------------------------------------------------------------
## The log-density cascade (SPEC 3.5, hazards H1-H6).
##
## Never routed through gkwdist::dgkw: that returns -Inf in reachable corners
## where the true log-density is around -88.3, and responses within 1e-10 of 1
## are not exotic in bounded data (SPEC N2).
## ---------------------------------------------------------------------------

.gkwq_logdens <- function(y, alpha, beta, gamma, delta, lambda,
                          delta_is_zero = FALSE) {
  ly <- log(y)
  l1mya <- .log1mexp(alpha * ly)                    # log(1 - y^alpha)
  lv <- .log1mexp(beta * l1mya)                     # log v
  d1 <- delta + 1

  out <- log(lambda) + log(alpha) + log(beta) -
    (lgamma(gamma) + lgamma(d1) - lgamma(gamma + d1)) +
    (alpha - 1) * ly + (beta - 1) * l1mya + (gamma * lambda - 1) * lv

  ## The delta term is DROPPED, never multiplied by zero: as v -> 1,
  ## log u -> -Inf and 0 * (-Inf) is NaN (hazard H4).
  if (!delta_is_zero) {
    lu <- .log1mexp(lambda * lv)                    # log(1 - v^lambda)
    out <- out + delta * lu
  }
  out
}

## Distribution function, in the log domain throughout.
.gkwq_logcdf <- function(y, alpha, beta, gamma, delta, lambda,
                         lower.tail = TRUE) {
  lv <- .log1mexp(beta * .log1mexp(alpha * log(y)))
  stats::pbeta(exp(lambda * lv), gamma, delta + 1,
               lower.tail = lower.tail, log.p = TRUE)
}
