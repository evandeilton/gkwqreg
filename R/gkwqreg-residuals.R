#' Residuals from a fitted quantile regression
#'
#' Extracts one of seven residual types from a `"gkwqreg"` fit. Because the
#' model targets a conditional quantile rather than a conditional mean, some of
#' these residuals have no counterpart in mean regression, and one of them --
#' `"response"` -- carries a familiar name with a different meaning. The table
#' in Details gives the definition, the reference distribution and the intended
#' diagnostic use of each.
#'
#' @details
#' Throughout, \eqn{Q_i = Q(\tau \mid x_i)} is the fitted conditional
#' \eqn{\tau}-quantile of observation \eqn{i}, the quantity returned by
#' [fitted()]; \eqn{\theta_i} is the corresponding fitted parameter vector
#' \eqn{(\alpha_i, \beta_i, \gamma_i, \delta_i, \lambda_i)} after the anchor has
#' been solved; \eqn{F(\cdot; \theta_i)} is the fitted conditional distribution
#' function; \eqn{\Phi} is the standard normal distribution function; and
#' \eqn{\rho_\tau(u) = u\{\tau - 1(u < 0)\}} is the check, or pinball, loss.
#' "Under correct specification" means that the fitted model is the
#' data-generating one.
#'
#' | `type` | definition | reference distribution | what it diagnoses |
#' |:--|:--|:--|:--|
#' | `"quantile"` | \eqn{\Phi^{-1}\{F(y_i; \theta_i)\}} | \eqn{N(0,1)}, exactly | the conditional distribution as a whole: shape, tails, links |
#' | `"cox-snell"` | \eqn{-\log\{1 - F(y_i; \theta_i)\}} | unit exponential, exactly | the same information, on a cumulative-hazard scale |
#' | `"pearson"` | \eqn{\{y_i - E(Y \mid x_i)\} / \mathrm{sd}(Y \mid x_i)} | mean 0, variance 1 | location and scale only; comparability with mean regression |
#' | `"deviance"` | \eqn{\mathrm{sign}(y_i - Q_i)\sqrt{2\{\sup_m \ell_i(m) - \ell_i(Q_i)\}}} | none exact | which observations the fitted quantile explains worst |
#' | `"response"` | \eqn{y_i - Q_i} | \eqn{\tau}-quantile equal to 0 | fit on the response scale, in the response's own units |
#' | `"check"` | \eqn{\rho_\tau(y_i - Q_i)} | non-negative; sample mean is the pinball loss | where in covariate space the loss accumulates |
#' | `"tau-sign"` | \eqn{1\{y_i \le Q_i\} - \tau} | mean 0 exactly, variance \eqn{\tau(1-\tau)} | calibration of the quantile level itself |
#'
#' @section Notes on individual types:
#' \describe{
#'   \item{`"quantile"`}{Dunn and Smyth (1996) proposed these as *randomized*
#'     quantile residuals, the randomization being needed only because
#'     \eqn{F(y)} is not uniform for a discrete family. Every family here is
#'     absolutely continuous, so \eqn{F(Y; \theta)} is exactly
#'     \eqn{\mathrm{Uniform}(0,1)} and no randomization is used: the residuals
#'     are a deterministic function of the data and reproduce exactly from call
#'     to call. They are computed from whichever tail of \eqn{F} is better
#'     conditioned, on the log scale, so that extreme observations do not
#'     collapse onto \eqn{\pm\infty}. This is the default and the right first
#'     look: a normal quantile-quantile plot of them tests the whole
#'     specification at once.}
#'   \item{`"cox-snell"`}{A monotone transformation of the same probability
#'     integral transform, so it carries exactly the same information. Prefer
#'     `"quantile"` for a normal quantile-quantile plot; prefer `"cox-snell"`
#'     when an exponential probability plot or a cumulative-hazard plot is
#'     wanted, as is conventional in survival work.}
#'   \item{`"pearson"`}{Neither the conditional mean nor the conditional
#'     standard deviation has a closed form under this reparametrization. Both
#'     are obtained by Gauss-Legendre quadrature of the quantile function, using
#'     \eqn{E(Y^k) = \int_0^1 Q(u)^k \, du}. The type exists for comparability
#'     with mean-regression packages; it is not the natural residual for a
#'     quantile model, since it standardises about a moment the model never
#'     targeted.}
#'   \item{`"deviance"`}{For a continuous density the usual saturated
#'     log-likelihood is unbounded above -- letting *every* parameter float
#'     freely drives it to infinity -- so the standard construction has no
#'     finite reference. The reference used here is the largest log density the
#'     model can give this observation by moving the conditional quantile alone.
#'     Writing \eqn{\ell_i(m)} for the log density of \eqn{y_i} when the
#'     conditional \eqn{\tau}-quantile is set to \eqn{m} and every nuisance
#'     parameter is held at its fitted value, the residual is
#'     \eqn{d_i = \mathrm{sign}(y_i - Q_i)[2\{\sup_m \ell_i(m) -
#'     \ell_i(Q_i)\}]^{1/2}}. The supremum is taken numerically, one bounded
#'     one-dimensional maximisation per observation, and it is genuinely a
#'     maximum, so the difference is non-negative by construction.
#'     It is tempting to use \eqn{\ell_i(y_i)} instead -- to relocate the fitted
#'     quantile onto the observation -- and that is what this package did before
#'     version 0.1.0. It is not the same thing: the density is generally
#'     maximised at some \eqn{m \neq y_i}, so \eqn{\ell_i(y_i)} can fall *below*
#'     the fitted value and the residual is then floored to zero. On a
#'     400-observation `kw` fit that silently zeroed 93 residuals; taking the
#'     supremum leaves none. There is still no exact null distribution: use
#'     `"quantile"` for distributional checking and `"deviance"` for ranking
#'     observations by how badly the fitted quantile explains them.}
#'   \item{`"response"`}{**The reference point is the fitted quantile, not the
#'     fitted mean.** In mean-regression packages `type = "response"` returns
#'     \eqn{y_i - E(Y \mid x_i)}; here it returns \eqn{y_i - Q(\tau \mid x_i)}.
#'     Under correct specification a proportion \eqn{\tau} of these values falls
#'     below zero, so their \eqn{\tau}-quantile is zero -- their mean is not
#'     zero, and there is no reason it should be.}
#'   \item{`"check"`}{The loss the model is fitted to minimise, evaluated
#'     observation by observation. Its sample mean is exactly the pinball loss
#'     returned by [pinball()], so plotting it against a covariate decomposes
#'     that single number and shows where the fit deteriorates. It is
#'     non-negative and strongly right-skewed; comparing its mean between two
#'     fits is comparing their pinball losses.}
#'   \item{`"tau-sign"`}{The sharpest check the model admits, and the one with
#'     no analogue in mean regression. If the conditional quantile is correctly
#'     specified then
#'     \deqn{E\left[\,1\{Y \le Q_\tau(X)\} \mid X\,\right] = \tau}
#'     holds exactly, for every \eqn{X}, with no asymptotics and no assumption
#'     about the rest of the distribution. Each residual takes only two values,
#'     \eqn{1 - \tau} with probability \eqn{\tau} and \eqn{-\tau} with
#'     probability \eqn{1 - \tau}, so conditionally on \eqn{X} it has mean zero
#'     and variance \eqn{\tau(1 - \tau)}. Averaging within a bin of covariate
#'     values therefore gives a direct, distribution-free calibration test: a
#'     bin mean far from zero is evidence of miscalibration in that region of
#'     covariate space. Under the null and independence, a bin of \eqn{m}
#'     observations has a mean that is approximately normal with standard error
#'     \eqn{\sqrt{\tau(1 - \tau)/m}}, which supplies an immediate yardstick.}
#' }
#'
#' @section Which type to use:
#' Start with `"quantile"`: it is exact, has a single unambiguous reference
#' distribution, and a departure in its quantile-quantile plot implicates the
#' specification as a whole. Follow it with `"tau-sign"` averaged within
#' covariate bins, which localises any failure and does so without assuming the
#' rest of the distribution is right. Use `"check"` to see where the fitted loss
#' is concentrated, `"deviance"` to identify individual poorly explained
#' observations, and `"response"` when the residual must be read in the units of
#' the response. Reserve `"pearson"` for comparisons against mean-regression
#' output, and `"cox-snell"` for plots conventionally drawn on the exponential
#' scale.
#'
#' @param object A fitted model object of class `"gkwqreg"`, as returned by
#'   [gkwqreg()]. The response must have been retained, which is the default
#'   (`y = TRUE`); otherwise an error asks for a refit.
#' @param type Character; the residual to compute, one of `"quantile"`,
#'   `"cox-snell"`, `"pearson"`, `"deviance"`, `"response"`, `"check"` or
#'   `"tau-sign"`. Partially matched, and defaulting to `"quantile"`. See the
#'   table in Details.
#' @param ... Not used; present for compatibility with the [residuals()]
#'   generic.
#'
#' @return A numeric vector with one element per observation used in the fit,
#'   in the order of the model frame after `subset` and `na.action` have been
#'   applied, so that its length is `nobs(object)` rather than the number of
#'   rows of the original data. Fitting at several levels at once produces a
#'   `"gkwqregs"` container instead, whose `residuals()` method returns a matrix
#'   with one column per level, the columns labelled by the level.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' *Journal of Computational and Graphical Statistics* **5**, 236-244.
#'
#' @seealso [pinball()] for the loss whose mean the `"check"` residual
#'   decomposes; [plot.gkwqreg()] for the standard diagnostic panels;
#'   [fitted()] and [predict()] for the fitted quantiles these residuals are
#'   measured against.
#'
#' @examples
#' ## A Kumaraswamy response whose true conditional median is plogis(0.3 + 0.9x).
#' ## The beta anchor at tau = 0.5 makes that exact: choosing
#' ## beta_i = log(1 - tau) / log(1 - mu_i^alpha) puts the median at mu_i.
#' set.seed(2024)
#' n  <- 400
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.3 + 0.9 * x)
#' y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
#'
#' fit <- gkwqreg(y ~ x, tau = 0.5, family = "kw")
#'
#' ## 1. Quantile residuals are exactly standard normal when the model is right.
#' rq <- residuals(fit)                    # type = "quantile" is the default
#' round(c(mean = mean(rq), sd = sd(rq)), 3)
#' ##   mean     sd
#' ## -0.001  1.005          -- indistinguishable from N(0, 1)
#'
#' round(quantile(rq, c(0.025, 0.25, 0.5, 0.75, 0.975)), 2)
#' ##  2.5%   25%   50%   75% 97.5%
#' ## -2.06 -0.65  0.05  0.63  2.01
#' round(qnorm(c(0.025, 0.25, 0.5, 0.75, 0.975)), 2)
#' ## [1] -1.96 -0.67  0.00  0.67  1.96      -- the theoretical counterparts
#'
#' round(shapiro.test(rq)$p.value, 3)
#' ## [1] 0.484                              -- no evidence against normality
#'
#' ## 2. Cox-Snell residuals are unit exponential: mean and sd both 1.
#' round(c(mean = mean(residuals(fit, "cox-snell")),
#'         sd   = sd(residuals(fit, "cox-snell"))), 3)
#' ##  mean    sd
#' ## 0.999 1.010
#'
#' ## 3. tau-sign residuals average to zero, with sd sqrt(tau(1 - tau)) = 0.5.
#' rts <- residuals(fit, "tau-sign")
#' round(c(mean = mean(rts), sd = sd(rts)), 4)
#' ##    mean      sd
#' ## -0.0125  0.5005
#' ## The standard error of that mean under the null is sqrt(.25 / 400) = 0.025,
#' ## so -0.0125 is half a standard error from zero: the fit is calibrated.
#'
#' ## 4. What miscalibration looks like. Here the true median is curved in x,
#' ## but one of the two fitted models is linear on the logit scale.
#' set.seed(7)
#' n  <- 500
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.2 + 0.8 * x - 0.9 * x^2)
#' y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
#' wrong <- gkwqreg(y ~ x,          tau = 0.5, family = "kw")   # omits x^2
#' right <- gkwqreg(y ~ x + I(x^2), tau = 0.5, family = "kw")
#'
#' bin <- cut(x, quantile(x, 0:4 / 4), include.lowest = TRUE,
#'            labels = paste0("Q", 1:4))
#' round(tapply(residuals(wrong, "tau-sign"), bin, mean), 3)
#' ##     Q1     Q2     Q3     Q4
#' ##  0.252 -0.324 -0.188  0.348
#' ## Far from zero, and sign-alternating: the fitted median sits above the data
#' ## in the outer quartiles of x and below it in the middle two. With about 125
#' ## observations per bin the null standard error is sqrt(.25 / 125) = 0.045,
#' ## so these bin means are four to eight standard errors from zero.
#'
#' round(tapply(residuals(right, "tau-sign"), bin, mean), 3)
#' ##     Q1     Q2     Q3     Q4
#' ## -0.028  0.004  0.004  0.028      -- all well inside one standard error
#'
#' ## The pinball loss agrees, being the mean of the check residual.
#' round(c(wrong = mean(residuals(wrong, "check")),
#'         right = mean(residuals(right, "check"))), 4)
#' ##  wrong  right
#' ## 0.1042 0.0643
#' @export
residuals.gkwqreg <- function(object,
                              type = c("quantile", "cox-snell", "pearson",
                                       "deviance", "response", "check",
                                       "tau-sign"),
                              ...) {
  type <- match.arg(type)
  y <- object$y
  if (is.null(y)) {
    stop("the response was not retained; refit with y = TRUE.", call. = FALSE)
  }
  q <- object$fitted.values
  pv <- object$parameter_vectors
  tau <- object$tau

  switch(type,
    response = y - q,
    check = (y - q) * (tau - (y < q)),
    `tau-sign` = as.numeric(y <= q) - tau,
    quantile = {
      lp <- .gkwq_logcdf(y, pv$alpha, pv$beta, pv$gamma, pv$delta, pv$lambda)
      lq <- .gkwq_logcdf(y, pv$alpha, pv$beta, pv$gamma, pv$delta, pv$lambda,
                         lower.tail = FALSE)
      ## Take whichever tail is better conditioned, so extreme observations do
      ## not collapse onto +/-Inf.
      ifelse(lp < lq,
             stats::qnorm(lp, log.p = TRUE),
             stats::qnorm(lq, log.p = TRUE, lower.tail = FALSE))
    },
    `cox-snell` = -.gkwq_logcdf(y, pv$alpha, pv$beta, pv$gamma, pv$delta,
                                pv$lambda, lower.tail = FALSE),
    pearson = {
      m <- .gkwq_moment(pv, 1L)
      s <- sqrt(.gkwq_variance(pv))
      (y - m) / s
    },
    deviance = {
      dz <- object$spec$delta_is_zero == 1L
      l_hat <- .gkwq_logdens(y, pv$alpha, pv$beta, pv$gamma, pv$delta,
                             pv$lambda, delta_is_zero = dz)
      ## The reference is the SUPREMUM of the observation's own log-density
      ## over the conditional quantile, holding the nuisance parameters at
      ## their fitted values. Relocating the quantile onto the observation --
      ## setting mu_i = y_i -- is the obvious shortcut and is NOT that
      ## supremum: the density is generally maximised at some mu != y, so the
      ## "saturated" value can fall BELOW the fitted one and the deviance
      ## residual gets floored to zero. On a 400-observation kw fit that
      ## happened for 93 observations. Maximising properly leaves none.
      l_sat <- .gkwq_saturated_loglik(object, y, dz)
      sign(y - q) * sqrt(pmax(2 * (l_sat - l_hat), 0))
    }
  )
}

## The non-anchored parameters as fitted, in the shape .gkwq_reconstruct wants.
.gkwq_theta_from_fit <- function(object) {
  nm <- object$parts[-1L]
  stats::setNames(lapply(nm, function(p) object$parameter_vectors[[p]]), nm)
}

#' @rdname gkwqregs-methods
#' @export
residuals.gkwqregs <- function(object, ...) {
  out <- sapply(object$fits, residuals, ...)
  colnames(out) <- format(object$taus, trim = TRUE)
  out
}

## Per-observation supremum of the log-density over the conditional quantile,
## with the nuisance parameters held at their fitted values. One bounded 1-D
## maximisation per observation; the objective is smooth and unimodal in mu on
## (0,1), so optimize() is adequate and no starting value is needed.
.gkwq_saturated_loglik <- function(object, y, delta_is_zero) {
  theta <- .gkwq_theta_from_fit(object)
  eps <- object$control$eps_mu
  tau <- object$tau
  spec <- object$spec
  vapply(seq_along(y), function(i) {
    th_i <- lapply(theta, function(v) v[i])
    f <- function(m) {
      P <- .gkwq_reconstruct(m, tau, spec, th_i)
      .gkwq_logdens(y[i], P[1], P[2], P[3], P[4], P[5],
                    delta_is_zero = delta_is_zero)
    }
    out <- tryCatch(
      stats::optimize(f, c(eps, 1 - eps), maximum = TRUE)$objective,
      error = function(e) NA_real_)
    if (is.finite(out)) out else f(min(max(y[i], eps), 1 - eps))
  }, numeric(1))
}
