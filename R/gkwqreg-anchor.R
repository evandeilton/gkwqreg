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
#' Evaluates \eqn{Q(\tau)}, the \eqn{\tau}-quantile of the five-parameter
#' Generalized Kumaraswamy distribution on the unit interval. The distribution
#' function inverts analytically, so no root finding is involved. That closed
#' form is what makes the fixed-level reparametrization behind [gkwqreg()]
#' possible at all: because \eqn{Q(\tau)} can be written down, one parameter can
#' be solved for so that the conditional quantile becomes a parameter itself.
#'
#' @details
#' The Generalized Kumaraswamy distribution function at \eqn{y \in (0,1)} is
#'
#' \deqn{F(y) = I_{v^\lambda}(\gamma, \delta + 1), \qquad
#'       v = 1 - (1 - y^\alpha)^\beta,}
#'
#' where \eqn{I} denotes the regularized incomplete beta function, that is
#' `pbeta`. Inverting this expression one layer at a time, and writing
#' \eqn{z_\tau = I^{-1}_\tau(\gamma, \delta + 1)}, which in R is
#' `qbeta(tau, gamma, delta + 1)`, gives
#'
#' \deqn{Q(\tau) = \left[\, 1 -
#'   \left\{ 1 - z_\tau^{1/\lambda} \right\}^{1/\beta} \,\right]^{1/\alpha}.}
#'
#' Two special cases of \eqn{z_\tau} are worth knowing, because the
#' sub-families that enjoy them never touch the incomplete beta function at
#' all. When \eqn{\gamma = 1} and \eqn{\delta = 0} the generator is the identity
#' and \eqn{z_\tau = \tau}; this covers the `"kw"` and `"ekw"` families. When
#' \eqn{\gamma = 1} with \eqn{\delta} free,
#' \eqn{z_\tau = 1 - (1 - \tau)^{1/(\delta + 1)}}; this covers `"kkw"`.
#'
#' @section Numerical behaviour:
#' The expression above is a cascade of two nested complements, each of the
#' form \eqn{1 - w} with \eqn{w} close to one in the regions that matter. It is
#' therefore evaluated entirely in the log domain, using a two-branch stable
#' evaluation of \eqn{\log(1 - e^{x})} for \eqn{x < 0} (Maechler 2012), which
#' switches between `log(-expm1(x))` and `log1p(-exp(x))` according to the
#' magnitude of \eqn{x}. Naive evaluation loses precision precisely where a
#' bounded-response model needs it most: when \eqn{\beta} is large the inner
#' complement approaches one, a regime the default `beta` anchor visits
#' routinely.
#'
#' Checked against `gkwdist::qgkw` over a 1701-point parameter grid -- seven
#' levels of \eqn{\tau} spanning 0.01 to 0.99, with each of the five parameters
#' taking three values -- the largest absolute discrepancy is `2.7e-14`.
#'
#' @param tau Numeric vector of quantile levels, strictly inside \eqn{(0,1)}.
#'   Values outside \eqn{[0,1]} produce `NaN`.
#' @param alpha,beta Strictly positive shape parameters of the Kumaraswamy
#'   kernel. \eqn{\alpha} acts on the left tail through \eqn{y^\alpha} and
#'   \eqn{\beta} on the right tail through \eqn{(1 - y^\alpha)^\beta}.
#' @param gamma,delta Shape parameters of the beta generator, which enter only
#'   through \eqn{z_\tau}: \eqn{\gamma} strictly positive, \eqn{\delta}
#'   non-negative. Their defaults, \eqn{\gamma = 1} and \eqn{\delta = 0}, remove
#'   the generator.
#' @param lambda Strictly positive exponentiation parameter, acting on
#'   \eqn{v^\lambda}. Its default \eqn{\lambda = 1} removes it.
#'
#' @return A numeric vector of quantiles in \eqn{(0,1)}, as long as the longest
#'   argument; all arguments are recycled to that length under the usual R
#'   rules. Inadmissible or missing inputs propagate as `NaN` or `NA`.
#'
#' @seealso [gkwq_anchors()] for the solves that invert this identity in each
#'   parameter in turn; [gkwq_parts()] for the family parameter sets;
#'   [gkwqreg()] for the regression model built on it.
#'
#' @examples
#' ## The median of a Kumaraswamy(2, 3) variate. With gamma = 1, delta = 0 and
#' ## lambda = 1 the generator drops out and Q collapses to the Kumaraswamy
#' ## quantile [1 - (1 - tau)^(1/beta)]^(1/alpha).
#' gkwq_quantile(0.5, alpha = 2, beta = 3)
#' ## [1] 0.454202
#' all.equal(gkwq_quantile(0.4, alpha = 2, beta = 3),
#'           (1 - (1 - 0.4)^(1 / 3))^(1 / 2))
#' ## [1] TRUE
#'
#' ## Q is increasing in tau, as any quantile function must be.
#' round(gkwq_quantile(c(0.05, 0.25, 0.5, 0.75, 0.95), alpha = 2, beta = 3), 4)
#' ## [1] 0.1302 0.3024 0.4542 0.6083 0.7947
#'
#' ## Arguments recycle: one level, three shapes.
#' round(gkwq_quantile(0.5, alpha = c(1, 2, 4), beta = 3), 4)
#' ## [1] 0.2063 0.4542 0.6739
#'
#' ## Q genuinely inverts F, to machine precision.
#' y <- c(0.05, 0.30, 0.62, 0.95)
#' p <- gkwdist::pgkw(y, 2, 3, 1.5, 0.5, 0.8)
#' max(abs(gkwq_quantile(p, 2, 3, 1.5, 0.5, 0.8) - y))
#' ## [1] 2.831069e-14
#'
#' ## Agreement with the reference implementation across the parameter space:
#' ## 1701 combinations of level and the five parameters.
#' g <- expand.grid(t = c(.01, .1, .25, .5, .75, .9, .99), a = c(.5, 1, 2.5),
#'                  b = c(.7, 1.5, 3), gm = c(.8, 1, 2), d = c(0, .5, 2),
#'                  L = c(.6, 1, 2))
#' max(abs(gkwq_quantile(g$t, g$a, g$b, g$gm, g$d, g$L) -
#'         gkwdist::qgkw(g$t, g$a, g$b, g$gm, g$d, g$L)))
#' ## [1] 2.738261e-14
#'
#' ## The identity that anchoring inverts. Fix a level and a target quantile,
#' ## then choose beta so that Q(tau) hits the target exactly. This is the
#' ## "beta" anchor of the "kw" family, the default used by gkwqreg().
#' tau <- 0.5; mu <- 0.7; alpha <- 2
#' beta <- log1p(-tau) / log1p(-mu^alpha)
#' beta
#' ## [1] 1.029409
#' gkwq_quantile(tau, alpha = alpha, beta = beta)
#' ## [1] 0.7   -- the median is the target mu, exactly
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
