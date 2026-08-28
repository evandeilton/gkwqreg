## ---------------------------------------------------------------------------
## S3 surface for class "gkwqreg".
##
## The class deliberately does NOT inherit from "gkwreg". Inheritance would make
## every missing method silently dispatch mean semantics on a quantile object,
## and the most dangerous case is fitted(): here it returns conditional
## quantiles, there conditional means. A missing method must error, not guess.
## See docs/adr/0003.
## ---------------------------------------------------------------------------

#' @export
coef.gkwqreg <- function(object, part = NULL, ...) {
  if (is.null(part)) return(object$coefficients)
  part <- match.arg(part, object$parts, several.ok = TRUE)
  if (length(part) == 1L) object$coef_list[[part]] else object$coef_list[part]
}

#' Covariance matrix of a quantile regression fit
#'
#' @param object A `"gkwqreg"` fit.
#' @param type `"expected"` and `"observed"` both invert the observed
#'   information from [stats::optimHess()] applied to the automatic-
#'   differentiation gradient; `"sandwich"` returns
#'   `H^-1 (sum_i s_i s_i') H^-1` with per-observation scores.
#' @param ... Unused.
#'
#' @details
#' The observed information is never taken from `obj$he()` (that would need
#' second-order differentiation through the incomplete-beta atomic) and never
#' from the naive `J' H J` sandwich of the unreparametrized Hessian, which omits
#' the curvature term `sum_k g_k grad^2 theta_k` and is wrong in a regression.
#'
#' The conditional quantile is *not* orthogonal to the remaining parameters in
#' the Cox-Reid sense, so the off-diagonal blocks matter: report the full matrix
#' rather than per-part standard errors alone.
#'
#' @return A covariance matrix.
#' @export
vcov.gkwqreg <- function(object, type = c("expected", "observed", "sandwich"),
                         ...) {
  type <- match.arg(type)
  if (type != "sandwich") {
    if (is.null(object$vcov)) {
      stop("no covariance matrix available; refit with control = gkwq_control(hessian = TRUE).",
           call. = FALSE)
    }
    return(object$vcov)
  }
  bread <- object$vcov
  if (is.null(bread)) {
    stop("the sandwich estimator needs the observed information; refit with hessian = TRUE.",
         call. = FALSE)
  }
  S <- estfun.gkwqreg(object)
  meat <- crossprod(S)
  out <- bread %*% meat %*% bread
  dimnames(out) <- dimnames(bread)
  out
}

#' Per-observation score contributions
#'
#' Supplies the empirical estimating functions so that `sandwich::vcovHC()`,
#' `sandwich::vcovCL()` and `lmtest::coeftest()` work on `"gkwqreg"` fits.
#'
#' @param x A `"gkwqreg"` fit.
#' @param ... Unused.
#' @return An `n` by `p` matrix of score contributions.
#' @export
estfun.gkwqreg <- function(x, ...) {
  obj <- x$obj
  if (is.null(obj)) {
    stop("the fitted object does not carry its TMB object; refit in this session.",
         call. = FALSE)
  }
  f <- function(par) as.numeric(obj$report(par)$loglik_i)
  S <- numDeriv::jacobian(f, x$coefficients)
  colnames(S) <- names(x$coefficients)
  S
}

#' Bread matrix for sandwich covariance
#'
#' The inverse observed information scaled by the sample size, as
#' `sandwich::sandwich()` expects. Together with [estfun.gkwqreg()] this makes
#' `sandwich::vcovHC()`, `sandwich::vcovCL()` and `lmtest::coeftest()` work on
#' `"gkwqreg"` fits without any further glue.
#'
#' @param x A `"gkwqreg"` fit.
#' @param ... Unused.
#' @return A `p` by `p` matrix.
#' @export
bread.gkwqreg <- function(x, ...) {
  if (is.null(x$vcov)) stop("no observed information available.", call. = FALSE)
  x$vcov * x$nobs
}

#' @export
logLik.gkwqreg <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$nobs,
            class = "logLik")
}

#' @export
nobs.gkwqreg <- function(object, ...) object$nobs

#' @export
AIC.gkwqreg <- function(object, ..., k = 2) {
  -2 * object$loglik + k * object$npar
}

#' @export
BIC.gkwqreg <- function(object, ...) object$bic

#' Fitted conditional quantiles
#'
#' **Returns the fitted conditional `tau`-quantile, not the conditional mean.**
#'
#' @param object A `"gkwqreg"` fit.
#' @param type `"quantile"` (the default) gives `Q(tau | x_i)`; `"mean"` gives
#'   `E[Y | x_i]` by quadrature; `"parameter"` gives the reconstructed
#'   five-parameter vector per observation.
#' @param ... Unused.
#' @return A numeric vector, or a data frame for `type = "parameter"`.
#' @export
fitted.gkwqreg <- function(object, type = c("quantile", "mean", "parameter"),
                           ...) {
  type <- match.arg(type)
  switch(type,
    quantile = object$fitted.values,
    mean = .gkwq_moment(object$parameter_vectors, order = 1L),
    parameter = object$parameter_vectors
  )
}

#' @export
family.gkwqreg <- function(object, ...) {
  structure(list(family = object$family, anchor = object$anchor,
                 tau = object$tau, parts = object$parts,
                 link = object$link),
            class = "gkwq_family")
}

#' @export
print.gkwq_family <- function(x, ...) {
  cat(sprintf("Generalized Kumaraswamy quantile family\n  family %s | tau %s | anchor %s\n  parts  %s\n",
              x$family, format(x$tau), x$anchor, paste(x$parts, collapse = " | ")))
  invisible(x)
}

#' @export
formula.gkwqreg <- function(x, ...) x$formula

#' @export
terms.gkwqreg <- function(x, part = "mu", ...) {
  part <- match.arg(part, x$parts)
  x$terms[[part]]
}

#' @export
model.frame.gkwqreg <- function(formula, ...) {
  if (is.null(formula$model)) {
    stop("the model frame was not retained; refit with model = TRUE.",
         call. = FALSE)
  }
  formula$model
}

#' @export
model.matrix.gkwqreg <- function(object, part = "mu", ...) {
  part <- match.arg(part, object$parts)
  if (!is.null(object$x)) return(object$x[[part]])
  if (is.null(object$model)) {
    stop("neither design matrices nor the model frame were retained; refit with x = TRUE.",
         call. = FALSE)
  }
  stats::model.matrix(object$terms[[part]], object$model,
                      contrasts.arg = object$contrasts[[part]])
}

#' @export
getCall.gkwqreg <- function(x, ...) x$call

#' @export
update.gkwqreg <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  if (!missing(formula.)) {
    call$formula <- stats::update(Formula::as.Formula(object$formula), formula.)
  }
  extras <- list(...)
  if (length(extras)) {
    for (nm in names(extras)) call[[nm]] <- extras[[nm]]
  }
  if (evaluate) eval(call, parent.frame()) else call
}

#' Confidence intervals for a quantile regression fit
#'
#' @param object A `"gkwqreg"` fit.
#' @param parm Coefficients to report; defaults to all.
#' @param level Confidence level.
#' @param method `"wald"` from the covariance matrix, `"profile"` from the TMB
#'   objective one coefficient at a time, or `"boot"` percentile intervals.
#' @param R Bootstrap replicates when `method = "boot"`.
#' @param ... Passed to [gkwq_boot()] when bootstrapping.
#'
#' @details
#' The conditional quantile is not orthogonal to the remaining parameters in the
#' Cox-Reid sense, so Wald intervals can be noticeably off in small samples.
#' `"profile"` is the cheap remedy and costs one TMB profile per coefficient.
#'
#' @return A matrix of lower and upper limits.
#' @export
confint.gkwqreg <- function(object, parm, level = 0.95,
                            method = c("wald", "profile", "boot"), R = 200L,
                            ...) {
  method <- match.arg(method)
  cf <- object$coefficients
  if (missing(parm)) parm <- names(cf) else if (is.numeric(parm)) parm <- names(cf)[parm]
  a <- (1 - level) / 2
  if (method == "wald") {
    se <- object$se[parm]
    z <- stats::qnorm(1 - a)
    ci <- cbind(cf[parm] - z * se, cf[parm] + z * se)
  } else if (method == "boot") {
    b <- gkwq_boot(object, R = R, ...)
    ci <- b$percentile[parm, , drop = FALSE]
  } else {
    ## Profile the TMB objective one coefficient at a time. Worth the cost here:
    ## the conditional quantile is not orthogonal to the nuisance parameters, so
    ## Wald intervals can be noticeably off in small samples.
    idx <- match(parm, names(cf))
    ## The profile itself can succeed and the interpolation onto the likelihood
    ## cut-off still fail, when the profile does not reach the cut-off on one
    ## side. Both stages must be guarded, and a failure has to be visible rather
    ## than silently NA.
    ci <- t(vapply(idx, function(i) {
      out <- tryCatch({
        pr <- TMB::tmbprofile(object$obj, i, trace = FALSE)
        as.numeric(stats::confint(pr, level = level))
      }, error = function(e) c(NA_real_, NA_real_))
      if (length(out) != 2L) out <- c(NA_real_, NA_real_)
      out
    }, numeric(2)))
    if (anyNA(ci)) {
      warning("the profile did not reach the likelihood cut-off for ",
              sum(!stats::complete.cases(ci)), " coefficient(s); those limits ",
              "are NA. Widen the profile or use method = \"wald\".",
              call. = FALSE)
    }
  }
  dimnames(ci) <- list(parm, paste0(format(100 * c(a, 1 - a), trim = TRUE), " %"))
  ci
}

## ---------------------------------------------------------------------------
## Moments by Gauss-Legendre quadrature of the quantile function.
##
## E[Y] = int_0^1 Q(u) du.  Q is smooth on (0,1) and available in closed form,
## so integrating the quantile function is both simpler and better conditioned
## than integrating the density -- no boundary singularity to work around.
## ---------------------------------------------------------------------------

.gkwq_gauss_legendre <- function(n = 64L) {
  ## Newton iteration on the Legendre polynomial roots.
  i <- seq_len(n)
  x <- cos(pi * (i - 0.25) / (n + 0.5))
  for (it in 1:100) {
    p0 <- rep(1, n); p1 <- x
    for (k in 2:n) {
      p2 <- p1
      p1 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
      p0 <- p2
    }
    dp <- n * (x * p1 - p0) / (x^2 - 1)
    dx <- -p1 / dp
    x <- x + dx
    if (max(abs(dx)) < 1e-15) break
  }
  p0 <- rep(1, n); p1 <- x
  for (k in 2:n) { p2 <- p1; p1 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k; p0 <- p2 }
  dp <- n * (x * p1 - p0) / (x^2 - 1)
  w <- 2 / ((1 - x^2) * dp^2)
  list(nodes = (x + 1) / 2, weights = w / 2)   # mapped to (0,1)
}

.gkwq_moment <- function(pv, order = 1L, n_nodes = 64L) {
  gl <- .gkwq_gauss_legendre(n_nodes)
  n <- nrow(pv)
  out <- numeric(n)
  for (i in seq_len(n)) {
    q <- gkwq_quantile(gl$nodes, pv$alpha[i], pv$beta[i], pv$gamma[i],
                       pv$delta[i], pv$lambda[i])
    out[i] <- sum(gl$weights * q^order)
  }
  out
}

.gkwq_variance <- function(pv, n_nodes = 64L) {
  m1 <- .gkwq_moment(pv, 1L, n_nodes)
  m2 <- .gkwq_moment(pv, 2L, n_nodes)
  pmax(m2 - m1^2, 0)
}

## ---------------------------------------------------------------------------
## summary / print
## ---------------------------------------------------------------------------

#' @export
summary.gkwqreg <- function(object, level = 0.95,
                            vcov_type = NULL, ...) {
  vt <- vcov_type %||% object$control$vcov_type
  V <- tryCatch(vcov(object, type = vt), error = function(e) object$vcov)
  cf <- object$coefficients
  se <- if (is.null(V)) object$se else sqrt(pmax(diag(V), 0))
  names(se) <- names(cf)
  z <- cf / se
  p <- 2 * stats::pnorm(-abs(z))
  tab <- cbind(Estimate = cf, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = p)

  part_of <- sub(":.*$", "", names(cf))

  ## Everything from here needs the response. It is optional (y = FALSE), so
  ## summary() reports what it can rather than refusing to print a coefficient
  ## table it already has.
  y <- object$y
  if (is.null(y)) {
    rq <- NA_real_
    pseudo_r1 <- NA_real_
    coverage <- NA_real_
    resid_summary <- rep(NA_real_, 5L)
  } else {
    rq <- residuals(object, type = "quantile")
    ## Koenker-Machado style goodness of fit on the scale the quantile actually
    ## targets: check loss against the best constant quantile.
    null_q <- stats::quantile(y, probs = object$tau, names = FALSE, type = 7)
    e0 <- y - null_q
    pin0 <- mean(e0 * (object$tau - (e0 < 0)))
    pseudo_r1 <- if (pin0 > 0) 1 - object$pinball / pin0 else NA_real_
    coverage <- mean(y <= object$fitted.values)
    resid_summary <- stats::quantile(rq[is.finite(rq)], c(0, .25, .5, .75, 1),
                                     names = FALSE)
  }

  structure(list(
    call = object$call, family = object$family, tau = object$tau,
    anchor = object$anchor, parts = object$parts, link = object$link,
    coefficients = tab, part = part_of, level = level, vcov_type = vt,
    loglik = object$loglik, npar = object$npar, nobs = object$nobs,
    aic = object$aic, bic = object$bic, pinball = object$pinball,
    pseudo_r1 = pseudo_r1, cond_number = object$cond_number,
    coverage = coverage,
    residual_summary = resid_summary,
    convergence = object$convergence,
    parameter_summary = vapply(object$parameter_vectors, mean, numeric(1))
  ), class = "summary.gkwqreg")
}

#' @export
print.summary.gkwqreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nGeneralized Kumaraswamy quantile regression\n")
  cat(sprintf("family: %s   tau: %s   anchor: %s\n",
              x$family, format(x$tau), x$anchor))
  cat("\nCall:\n"); print(x$call)

  if (!all(is.na(x$residual_summary))) {
    cat("\nQuantile residuals:\n")
    print(stats::setNames(round(x$residual_summary, digits),
                          c("Min", "1Q", "Median", "3Q", "Max")))
  }

  for (p in x$parts) {
    idx <- which(x$part == p)
    if (!length(idx)) next
    lab <- if (p == "mu") {
      sprintf("\nConditional %s-quantile (link %s) -- coefficients are effects on the LOG QUANTILE ODDS log(mu/(1-mu)):\n",
              format(x$tau), x$link[[p]])
    } else {
      sprintf("\n%s (link %s):\n", p, x$link[[p]])
    }
    cat(lab)
    stats::printCoefmat(x$coefficients[idx, , drop = FALSE], digits = digits,
                        signif.stars = getOption("show.signif.stars"))
  }

  cat(sprintf("\nAnchored parameter: %s (computed from the quantile, not estimated)\n",
              x$anchor))
  cat(sprintf("Log-likelihood: %s on %d Df   AIC: %s   BIC: %s\n",
              format(x$loglik, digits = digits), x$npar,
              format(x$aic, digits = digits), format(x$bic, digits = digits)))
  cat(sprintf("Pinball loss: %s   Pseudo-R1: %s\n",
              format(x$pinball, digits = digits),
              format(x$pseudo_r1, digits = digits)))
  if (is.na(x$coverage)) {
    cat("Empirical coverage: unavailable (refit with y = TRUE)\n")
  } else {
    cat(sprintf("Empirical coverage: %s (target %s)\n",
                format(x$coverage, digits = digits), format(x$tau)))
  }
  cat(sprintf("Covariance: %s   Information condition number: %s\n",
              x$vcov_type, format(x$cond_number, digits = 4)))

  if (is.finite(x$cond_number) && x$cond_number > 1e8) {
    cat("\nNote: the information matrix is ill-conditioned (condition number above 1e8).\n")
    cat("      Standard errors are unreliable; consider a smaller sub-family.\n")
  }
  if (!identical(x$convergence, 0L) && !identical(x$convergence, 0)) {
    cat("\nWarning: the optimizer did not report convergence.\n")
  }
  invisible(x)
}

#' @export
print.gkwqreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("\nGeneralized Kumaraswamy Quantile Regression  (family: %s, tau = %s, anchor: %s)\n",
              x$family, format(x$tau), x$anchor))
  cat("\nCall:\n"); print(x$call)
  cat("\nCoefficients:\n")
  for (p in x$parts) {
    cat(sprintf("  %s:\n", p))
    print(round(x$coef_list[[p]], digits))
  }
  cat(sprintf("\nLog-likelihood: %s   AIC: %s   n: %d\n",
              format(x$loglik, digits = digits), format(x$aic, digits = digits),
              x$nobs))
  cat("fitted() returns conditional quantiles, not means.\n")
  invisible(x)
}

#' @export
print.gkwqregs <- function(x, ...) {
  cat(sprintf("\n%d Generalized Kumaraswamy quantile regressions (family: %s, anchor: %s)\n",
              length(x$taus), x$family, x$anchor))
  cat("\nCall:\n"); print(x$call)
  cf <- sapply(x$fits, function(f) f$coefficients)
  colnames(cf) <- format(x$taus, trim = TRUE)
  cat("\nCoefficients by tau:\n")
  print(round(cf, 4))
  invisible(x)
}

#' @export
coef.gkwqregs <- function(object, ...) {
  cf <- sapply(object$fits, function(f) f$coefficients)
  colnames(cf) <- format(object$taus, trim = TRUE)
  cf
}

#' @export
fitted.gkwqregs <- function(object, ...) {
  out <- sapply(object$fits, function(f) f$fitted.values)
  colnames(out) <- format(object$taus, trim = TRUE)
  out
}

#' @export
logLik.gkwqregs <- function(object, ...) {
  stop("a `gkwqregs` container holds one likelihood per tau, and they are not ",
       "comparable or additive. Take logLik() of an individual fit, e.g. ",
       "logLik(object$fits[[1]]).", call. = FALSE)
}

#' @export
summary.gkwqregs <- function(object, ...) {
  lapply(object$fits, summary, ...)
}
