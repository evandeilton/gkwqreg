#' Residuals from a quantile regression fit
#'
#' @param object A `"gkwqreg"` fit.
#' @param type The residual to compute; see Details.
#' @param ... Unused.
#'
#' @details
#' \describe{
#'   \item{`"quantile"`}{Dunn and Smyth (1996), `qnorm(F(y_i; theta_i))`. Exact,
#'     with no randomization: the family is continuous. Standard normal under
#'     correct specification, and the default.}
#'   \item{`"cox-snell"`}{`-log(1 - F(y_i))`; unit exponential under correct
#'     specification.}
#'   \item{`"pearson"`}{`(y_i - E[Y|x_i]) / sd(Y|x_i)`, both by quadrature.
#'     Provided for comparability with mean-regression packages.}
#'   \item{`"deviance"`}{`sign(y_i - Q_i) sqrt(2 (l_i(y_i) - l_i(Q_i)))`, where
#'     `l_i(y_i)` is the log-likelihood with the fitted quantile moved onto the
#'     observation. For a continuous density the usual saturated likelihood is
#'     unbounded, so the reference here is the quantile that would place `y_i`
#'     exactly at level `tau` -- which is the quantity this model is about.}
#'   \item{`"response"`}{`y_i - Q(tau|x_i)`. **The reference point is the
#'     quantile, not the mean**, unlike the same name in mean regression.}
#'   \item{`"check"`}{`(y_i - Q_i)(tau - 1{y_i < Q_i})`. Its mean is the pinball
#'     loss; plotted against a covariate it shows where the fit deteriorates.}
#'   \item{`"tau-sign"`}{`1{y_i <= Q_i} - tau`. Under correct specification
#'     `E[1{Y <= Q_tau(X)} | X] = tau` exactly, so a non-zero mean within a
#'     covariate bin is direct, distribution-free evidence of quantile
#'     miscalibration -- a check that mean regression has no analogue for.}
#' }
#'
#' @return A numeric vector of residuals.
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' *Journal of Computational and Graphical Statistics* **5**, 236-244.
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
      ## The quantile relocated onto the observation.
      Py <- .gkwq_reconstruct(pmin(pmax(y, object$control$eps_mu),
                                   1 - object$control$eps_mu),
                              tau, object$spec,
                              .gkwq_theta_from_fit(object))
      l_sat <- .gkwq_logdens(y, Py[, 1], Py[, 2], Py[, 3], Py[, 4], Py[, 5],
                             delta_is_zero = dz)
      sign(y - q) * sqrt(pmax(2 * (l_sat - l_hat), 0))
    }
  )
}

## The non-anchored parameters as fitted, in the shape .gkwq_reconstruct wants.
.gkwq_theta_from_fit <- function(object) {
  nm <- object$parts[-1L]
  stats::setNames(lapply(nm, function(p) object$parameter_vectors[[p]]), nm)
}

#' @export
residuals.gkwqregs <- function(object, ...) {
  out <- sapply(object$fits, residuals, ...)
  colnames(out) <- format(object$taus, trim = TRUE)
  out
}
