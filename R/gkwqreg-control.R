#' Control parameters for [gkwqreg()]
#'
#' @param method Optimizer passed to [stats::nlminb()] or [stats::optim()].
#'   `"nlminb"` (the default) is a quasi-Newton method with the AD gradient.
#' @param maxit Maximum optimizer iterations.
#' @param reltol Relative convergence tolerance.
#' @param start_method How starting values are built. `"auto"` runs the full
#'   algorithm with an intercept-only pre-fit for the harder families,
#'   `"ols"` skips the pre-fit, `"intercept"` always runs it, `"gkwdist"` takes
#'   marginal estimates and leaves every slope at zero, `"user"` requires
#'   `start`.
#' @param start Named list of starting coefficients per part, used when
#'   `start_method = "user"`.
#' @param eps_y Response values are clamped into `[eps_y, 1 - eps_y]`. Values at
#'   exactly 0 or 1 lie outside the family's support.
#' @param eps_mu Smooth clamp applied to the fitted conditional quantile inside
#'   the tape. The anchor solves overflow in opposite tails (the beta solve as
#'   `mu -> 0`, the alpha solve as `mu -> 1`; SPEC N6), so this bound is what
#'   keeps either usable near the boundary.
#' @param tiny Floor for log-scale quantities that must stay strictly negative.
#' @param vcov_type Default covariance: `"expected"` from `sdreport`,
#'   `"observed"` from `optimHess` on the AD gradient, or `"sandwich"`.
#' @param hessian Whether to compute the covariance matrix at all.
#' @param warm_start For a vector `tau`, fit the level nearest the median first
#'   and warm-start outward in both directions.
#' @param silent Suppress TMB's tracing output.
#' @param ... Ignored, for forward compatibility.
#'
#' @return An object of class `"gkwq_control"`.
#' @examples
#' gkwq_control(maxit = 500, vcov_type = "sandwich")
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
