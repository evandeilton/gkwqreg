#' Tuning parameters for a Generalized Kumaraswamy quantile regression
#'
#' @description
#' Collects the optimizer settings, the starting-value strategy, the numerical
#' guards and the default covariance estimator into one object, passed to
#' [gkwqreg()] through its `control` argument. Every element has a default that
#' works for the well-behaved families without adjustment; this page explains
#' what each one does and the circumstances in which a user would change it.
#'
#' The same settings can be given directly to [gkwqreg()] through `...`, which
#' is convenient for a one-off fit. Note that doing so **replaces** any `control`
#' object passed alongside it, so use one route or the other.
#'
#' @section The optimizer:
#' `method` selects the maximizer of the log-likelihood. `"nlminb"` (the default)
#' is the quasi-Newton routine of [stats::nlminb()], driven by the exact
#' automatic-differentiation gradient; the remaining choices are the
#' corresponding methods of [stats::optim()]. No method needs bounds: the anchor
#' solve is a ratio of logarithms of numbers in \eqn{(0,1)} and is therefore
#' positive by construction, so the parameter space is genuinely unconstrained.
#'
#' `maxit` and `reltol` are translated to whichever optimizer is in use
#' (`iter.max` and `eval.max = 2 * maxit` with `rel.tol` for `nlminb`; `maxit`
#' and `reltol` for `optim`). Raise `maxit` when the optimizer reports a non-zero
#' convergence code, which [gkwqreg()] turns into a warning and `summary()`
#' repeats. If more iterations do not help, `"BFGS"` is the usual second attempt
#' and `"Nelder-Mead"` a slow, derivative-free fallback for a badly scaled
#' problem.
#'
#' @section Starting values:
#' `start_method` chooses how the initial coefficient vector is built. Unlike
#' mean regression, starting every coefficient at zero is not safe here: it puts
#' the conditional quantile at \eqn{g_\mu^{-1}(0)=0.5}{g_mu^-1(0) = 0.5} whatever `tau` is, and
#' because each anchor is a ratio of logarithms, a quantile starting near 0 or 1
#' makes the derived parameter explode. The default construction therefore
#'
#' 1. takes marginal estimates for the family from
#'    `gkwdist::gkwgetstartvalues()` on the response, and uses them as the
#'    nuisance intercepts with every nuisance slope at zero;
#' 2. fits a *quantile-shifted* least-squares regression for the quantile part:
#'    ordinary least squares of \eqn{g_\mu(y)}{g_mu(y)} on the quantile design matrix,
#'    then the intercept shifted by the `tau`-quantile of its residuals. With an
#'    intercept-only design this returns exactly
#'    \eqn{g_\mu(\mathrm{quantile}(y,\tau))}{g_mu(quantile(y, tau))}, so the two steps agree by
#'    construction rather than by luck;
#' 3. runs a feasibility sweep, evaluating the anchor in R at the proposed start
#'    and falling back to a neutral start if it is not finite and positive
#'    everywhere -- much cheaper than discovering the problem as a tape full of
#'    `NaN`s.
#'
#' The options then differ only in whether an intercept-only pre-fit is offered
#' as a *second candidate*. It is a candidate, never a replacement: the families
#' that benefit from it also have flat ridges (`kkw` nests `ekw` nests `kw`), so
#' a pre-fit can run off along one and land at an absurd corner. The objective is
#' evaluated at every candidate and the best is kept.
#'
#' \describe{
#'   \item{`"auto"`}{the default. The construction above, plus the
#'     intercept-only pre-fit for the harder families `kkw`, `bkw`, `gkw`, `mc`
#'     and `beta`, where weak identifiability turns a large share of imperfect
#'     starts into failures. The pre-fit costs a few milliseconds and is skipped
#'     for `kw` and `ekw`, which do not need it.}
#'   \item{`"ols"`, `"gkwdist"`}{the same construction with the pre-fit
#'     suppressed for every family. Use these when fitting a great many models in
#'     a simulation, or when the pre-fit is itself the source of trouble.}
#'   \item{`"intercept"`}{always add the pre-fit, `kw` and `ekw` included. This
#'     is the first thing to try when a fit fails to converge or settles at an
#'     implausible optimum.}
#'   \item{`"user"`}{take `start` verbatim and compute nothing. `start` is a
#'     named list with one numeric vector per part, in [gkwq_parts()] order, each
#'     as long as that part's design matrix has columns, and **on the link
#'     scale** -- with the default log link a starting \eqn{\alpha} of 2 is
#'     supplied as `log(2)`. A missing part is an error. The `coef_list`
#'     component of an earlier fit has exactly this shape, which makes it easy to
#'     restart a model or to reproduce a fit exactly.}
#' }
#'
#' @section Numerical guards:
#' `eps_y` clamps the response into \eqn{[\epsilon_y, 1-\epsilon_y]}. Values
#' exactly at 0 or 1 lie outside the family's support and are rejected with an
#' error rather than clamped; this guard exists for values that are merely
#' indistinguishable from the boundary in double precision.
#'
#' `eps_mu` is the one worth understanding. It is a smooth clamp confining the
#' fitted conditional quantile to \eqn{[\epsilon_\mu, 1-\epsilon_\mu]} inside the
#' automatic-differentiation tape, and it is what keeps the anchor usable near
#' the boundary. The two closed-form solves overflow in **opposite tails**. With
#' `family = "kw"` at `tau = 0.5`, the \eqn{\beta}-solve
#' \eqn{\log(1-\tau)/\log(1-\mu^{\alpha})}{log(1-tau) / log(1-mu^alpha)} grows without bound as
#' \eqn{\mu\to 0}{mu -> 0}, while the \eqn{\alpha}-solve
#' \eqn{\log(1-(1-\tau)^{1/\beta})/\log\mu}{log(1-(1-tau)^(1/beta)) / log(mu)} grows without bound as
#' \eqn{\mu\to 1}{mu -> 1}:
#'
#' ```
#'      mu      beta-solve   alpha-solve
#'  1e-06        6.93e+11      8.89e-02
#'  1e-03        6.93e+05      1.78e-01
#'  0.5          2.41e+00      1.77e+00
#'  1 - 1e-03    1.12e-01      1.23e+03
#'  1 - 1e-06    5.28e-02      1.23e+06
#' ```
#'
#' (the examples below reproduce this table in three lines of arithmetic).
#' Neither anchor dominates, so two things follow. First, `eps_mu` should be
#' *raised*, not lowered, when a fit at an extreme `tau` produces `NaN`s or
#' astronomical shape parameters: `1e-6` is a reasonable next value. Second,
#' where the response mass lies is a legitimate reason to prefer one anchor over
#' the other -- see the anchor section of [gkwqreg()].
#'
#' `tiny` is a strict-positivity floor applied inside the tape to each *modelled*
#' shape parameter (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\lambda}) after
#' its inverse link, so that a linear predictor wandering far negative during the
#' search cannot hand the density a zero or negative shape. \eqn{\delta} is
#' exempt, needing only \eqn{\delta\ge0}. It is not the floor used for log-scale
#' quantities, which is a fixed constant chosen so as not to distort
#' \eqn{\log v}{log v} exactly where the \eqn{\beta}-anchor operates. Changing `tiny` is
#' rarely useful; raise it if the optimizer stalls with one shape parameter
#' pinned at zero.
#'
#' @section The covariance matrix:
#' `vcov_type` sets the estimator that `summary()` reports by default;
#' [vcov.gkwqreg()] can always be asked for another one after the fact.
#' `"expected"` and `"observed"` both return the inverse observed information,
#' computed by [stats::optimHess()] on the automatic-differentiation gradient --
#' the same matrix under two names, so that scripts written against either
#' convention run unchanged. `"sandwich"` returns
#' \eqn{H^{-1}\bigl(\sum_i s_i s_i^{\top}\bigr)H^{-1}}{H^-1 (sum_i s_i s_i') H^-1} from the per-observation
#' scores, which is the estimator to use when the family is a working model
#' rather than a belief; it costs one numeric Jacobian of the per-observation
#' log-likelihood. Because the conditional quantile is not orthogonal to the
#' remaining parameters in the Cox-Reid sense, the off-diagonal blocks matter and
#' the full matrix should be reported rather than per-part standard errors alone.
#'
#' `hessian = FALSE` skips the [stats::optimHess()] call altogether. The fit then
#' carries point estimates, the log-likelihood and fitted quantiles, but no
#' `vcov`, no standard errors and no condition number: `summary()` prints `NA`
#' in those columns and [vcov.gkwqreg()] and [confint.gkwqreg()] raise an error.
#' Worth it inside a simulation or bootstrap loop where only the point estimates
#' are used.
#'
#' @param method Optimizer passed to [stats::nlminb()] or [stats::optim()].
#'   `"nlminb"` (the default) is a quasi-Newton method driven by the exact
#'   automatic-differentiation gradient.
#' @param maxit Maximum optimizer iterations. Raise it when the fit reports a
#'   non-zero convergence code.
#' @param reltol Relative convergence tolerance on the objective.
#' @param start_method How starting values are built: `"auto"` (the default) adds
#'   an intercept-only pre-fit as a second candidate for the harder families,
#'   `"ols"` and `"gkwdist"` suppress that pre-fit, `"intercept"` always runs it,
#'   and `"user"` takes `start` verbatim. See the section above.
#' @param start Named list of starting coefficients, one numeric vector per part
#'   in [gkwq_parts()] order and on the link scale. Required when
#'   `start_method = "user"`, ignored otherwise. A previous fit's `coef_list` has
#'   exactly this shape.
#' @param eps_y Response values are clamped into `[eps_y, 1 - eps_y]`. Values at
#'   exactly 0 or 1 lie outside the family's support and are an error, not a
#'   clamping target.
#' @param eps_mu Smooth clamp applied to the fitted conditional quantile inside
#'   the tape. The anchor solves overflow in opposite tails (the beta solve as
#'   `mu -> 0`, the alpha solve as `mu -> 1`), so this bound is what keeps either
#'   usable near the boundary; raise it, do not lower it, when an extreme `tau`
#'   misbehaves.
#' @param tiny Strict-positivity floor applied inside the tape to each modelled
#'   shape parameter after its inverse link. Not the floor used for log-scale
#'   quantities, which is a fixed constant.
#' @param vcov_type Default covariance reported by `summary()`: `"expected"` and
#'   `"observed"` both invert the observed information, `"sandwich"` builds the
#'   empirical sandwich from per-observation scores.
#' @param hessian Whether to compute the covariance matrix at all. `FALSE` gives
#'   a fit with point estimates but no standard errors.
#' @param warm_start Only relevant when `tau` is a vector. `TRUE` fits the level
#'   nearest the median first and warm-starts outward in both directions, each
#'   sweep carrying its own neighbour forward. `FALSE` gives every level an
#'   independent cold start, which is the check that the warm start is a
#'   convenience rather than an assumption.
#' @param silent Suppress TMB's tracing output while the tape is built and
#'   evaluated.
#' @param ... Ignored, for forward compatibility.
#'
#' @return An object of class `"gkwq_control"`: a list holding the validated
#'   elements `method`, `maxit`, `reltol`, `start_method`, `start`, `eps_y`,
#'   `eps_mu`, `tiny`, `vcov_type`, `hessian`, `warm_start` and `silent`. It has
#'   a `print()` method summarising the optimizer, the starting-value strategy,
#'   the covariance type and the two guards, and it is stored in the `control`
#'   component of every fit.
#'
#' @examples
#' ## -------------------------------------------------------------------------
#' ## 1. The defaults, and how to change them
#' ## -------------------------------------------------------------------------
#' gkwq_control()
#' gkwq_control(method = "BFGS", maxit = 500, vcov_type = "sandwich")
#'
#' ## -------------------------------------------------------------------------
#' ## 2. Why eps_mu exists: the two anchors overflow in OPPOSITE tails
#' ## -------------------------------------------------------------------------
#' ## Kumaraswamy family at tau = 0.5. The beta-solve blows up as mu -> 0, the
#' ## alpha-solve as mu -> 1, so neither anchor is safe at both boundaries and
#' ## eps_mu is what keeps the fitted quantile away from whichever one bites.
#' mu <- c(1e-6, 1e-3, 0.5, 1 - 1e-3, 1 - 1e-6)
#' signif(cbind(mu = mu,
#'              beta_solve  = log(1 - 0.5) / log(1 - mu^2),        # alpha = 2
#'              alpha_solve = log(1 - (1 - 0.5)^(1 / 2)) / log(mu) # beta  = 2
#'              ), 3)
#' ## If a fit at an extreme tau returns NaN or shape parameters of order 1e11,
#' ## raise eps_mu (say to 1e-6) or switch to the anchor that fails in the other
#' ## tail -- see the anchor section of ?gkwqreg.
#'
#' ## -------------------------------------------------------------------------
#' ## 3. Robust standard errors
#' ## -------------------------------------------------------------------------
#' set.seed(2024)
#' n  <- 250
#' x1 <- runif(n, -2, 2)
#' mu_true <- plogis(0.4 + 1.1 * x1)
#' y  <- gkwdist::rkw(n, alpha = 2,
#'                    beta = log(1 - 0.5) / log(1 - mu_true^2))
#' dat <- data.frame(y = y, x1 = x1)
#'
#' fit <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw")
#' round(cbind(model    = sqrt(diag(vcov(fit, type = "expected"))),
#'             sandwich = sqrt(diag(vcov(fit, type = "sandwich")))), 4)
#' ## Close agreement here is expected: the data really were generated from the
#' ## fitted family. A large gap is evidence of misspecification, and then the
#' ## sandwich column is the one to report.
#'
#' ## Make the sandwich the default for a fit, so that summary() uses it:
#' fit_s <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw",
#'                  control = gkwq_control(vcov_type = "sandwich"))
#' summary(fit_s)$vcov_type
#'
#' \donttest{
#' ## -------------------------------------------------------------------------
#' ## 4. Supplying starting values
#' ## -------------------------------------------------------------------------
#' ## One vector per part, in gkwq_parts() order, on the LINK scale: the quantile
#' ## part is logit, the alpha part is log, so alpha = 2 is entered as log(2).
#' gkwq_parts("kw")
#' st <- list(mu = c(0, 0), alpha = log(2))
#' fit_u <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw",
#'                  control = gkwq_control(start_method = "user", start = st))
#' round(coef(fit_u), 4)
#'
#' ## A previous fit's coef_list already has that shape, so a model can be
#' ## restarted from where it stopped:
#' fit_r <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw",
#'                  control = gkwq_control(start_method = "user",
#'                                         start = fit$coef_list))
#' all.equal(coef(fit), coef(fit_r), tolerance = 1e-6)
#'
#' ## -------------------------------------------------------------------------
#' ## 5. Warm starts across a grid of quantile levels
#' ## -------------------------------------------------------------------------
#' g_warm <- gkwqreg(y ~ x1, data = dat, tau = c(0.1, 0.5, 0.9), family = "kw")
#' g_cold <- gkwqreg(y ~ x1, data = dat, tau = c(0.1, 0.5, 0.9), family = "kw",
#'                   control = gkwq_control(warm_start = FALSE))
#' max(abs(coef(g_warm) - coef(g_cold)))
#' ## A negligible difference is the check that warm starting only saves time.
#'
#' ## -------------------------------------------------------------------------
#' ## 6. Skipping the covariance matrix
#' ## -------------------------------------------------------------------------
#' ## Worth it inside a simulation loop: point estimates only, no optimHess.
#' fit_fast <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw",
#'                     control = gkwq_control(hessian = FALSE))
#' coef(fit_fast)
#' is.null(fit_fast$vcov)
#' }
#'
#' @seealso [gkwqreg()] for the model these settings control, [gkwq_parts()] for
#'   the part order that `start` must follow, [vcov.gkwqreg()] and
#'   [confint.gkwqreg()] for the inference the covariance feeds, and
#'   [gkwq_boot()] when neither the model-based nor the sandwich covariance is
#'   trusted.
#'
#' @family model fitting
#' @export
gkwq_control <- function(method = c("nlminb", "BFGS", "L-BFGS-B", "CG", "Nelder-Mead"),
                         maxit = 1000L, reltol = 1e-10,
                         start_method = c("auto", "ols", "intercept", "gkwdist", "user"),
                         start = NULL,
                         eps_y = 1e-10, eps_mu = 1e-8, tiny = 1e-12,
                         vcov_type = c("expected", "observed", "sandwich"),
                         hessian = TRUE, warm_start = TRUE, silent = TRUE, ...) {
  method <- match.arg(method)
  start_method <- match.arg(start_method)
  vcov_type <- match.arg(vcov_type)
  if (start_method == "user" && is.null(start)) {
    stop("start_method = \"user\" requires `start`.", call. = FALSE)
  }
  structure(
    list(method = method, maxit = as.integer(maxit), reltol = reltol,
         start_method = start_method, start = start,
         eps_y = eps_y, eps_mu = eps_mu, tiny = tiny,
         vcov_type = vcov_type, hessian = isTRUE(hessian),
         warm_start = isTRUE(warm_start), silent = isTRUE(silent)),
    class = "gkwq_control"
  )
}

#' @export
print.gkwq_control <- function(x, ...) {
  cat("gkwqreg control\n")
  cat(sprintf("  optimizer    : %s (maxit %d, reltol %.1e)\n",
              x$method, x$maxit, x$reltol))
  cat(sprintf("  start method : %s\n", x$start_method))
  cat(sprintf("  covariance   : %s\n", x$vcov_type))
  cat(sprintf("  guards       : eps_y %.1e, eps_mu %.1e\n", x$eps_y, x$eps_mu))
  invisible(x)
}
