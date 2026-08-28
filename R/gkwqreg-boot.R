## ---------------------------------------------------------------------------
## Simulation from a fit, and the bootstrap.
## ---------------------------------------------------------------------------

#' Simulate responses from a fitted quantile regression
#'
#' Draws from each observation's fitted conditional distribution. Because the
#' anchored parameter is reconstructed first, a simulated sample has the fitted
#' conditional `tau`-quantile by construction, which makes this the natural
#' engine for parametric bootstrap and for simulated envelopes.
#'
#' @param object A `"gkwqreg"` fit.
#' @param nsim Number of replicate samples.
#' @param seed Optional seed; handled as in [stats::simulate()].
#' @param ... Unused.
#' @return A data frame with `nsim` columns of length `nobs(object)`.
#' @export
simulate.gkwqreg <- function(object, nsim = 1L, seed = NULL, ...) {
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      stats::runif(1)
    }
    old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  pv <- object$parameter_vectors
  n <- nrow(pv)
  out <- as.data.frame(lapply(seq_len(nsim), function(i) {
    gkwdist::rgkw(n, pv$alpha, pv$beta, pv$gamma, pv$delta, pv$lambda)
  }))
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}

#' Bootstrap a quantile regression fit
#'
#' @param object A `"gkwqreg"` fit.
#' @param R Number of bootstrap replicates.
#' @param type `"pairs"` resamples observations, which is robust to
#'   misspecification of the conditional distribution; `"parametric"` redraws
#'   responses from the fitted distributions, which is more efficient when the
#'   family is right and misleading when it is not.
#' @param seed Optional seed.
#' @param ... Unused.
#'
#' @details
#' `"pairs"` is the default deliberately. The whole trade-off of parametric
#' quantile regression is that it buys efficiency under a distributional
#' assumption; a bootstrap that re-imposes that same assumption cannot tell you
#' anything about whether it holds.
#'
#' @return An object of class `"gkwq_boot"` carrying the replicate coefficient
#'   matrix, its covariance and percentile intervals.
#' @export
gkwq_boot <- function(object, R = 200L, type = c("pairs", "parametric"),
                      seed = NULL, ...) {
  stopifnot(inherits(object, "gkwqreg"))
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(object$model)) {
    stop("the model frame was not retained; refit with model = TRUE to bootstrap.",
         call. = FALSE)
  }

  cl <- object$call
  mf <- object$model
  n <- nrow(mf)
  p <- length(object$coefficients)
  resp <- names(mf)[1L]

  reps <- matrix(NA_real_, R, p,
                 dimnames = list(NULL, names(object$coefficients)))
  sims <- if (type == "parametric") simulate.gkwqreg(object, nsim = R) else NULL

  for (r in seq_len(R)) {
    d <- if (type == "pairs") mf[sample.int(n, n, replace = TRUE), , drop = FALSE]
         else { z <- mf; z[[resp]] <- sims[[r]]; z }
    cl_r <- cl
    cl_r$data <- quote(d)
    cl_r$subset <- NULL
    cl_r$weights <- NULL
    f <- suppressWarnings(try(eval(cl_r, list(d = d), parent.frame()),
                              silent = TRUE))
    if (!inherits(f, "try-error") && isTRUE(f$convergence == 0)) {
      reps[r, ] <- f$coefficients
    }
  }

  ok <- stats::complete.cases(reps)
  if (sum(ok) < 2L) {
    stop("fewer than two bootstrap replicates converged; the fit is too ",
         "unstable to bootstrap.", call. = FALSE)
  }
  V <- stats::cov(reps[ok, , drop = FALSE])

  structure(list(
    replicates = reps, vcov = V, R = R, n_ok = sum(ok), type = type,
    estimate = object$coefficients,
    se = sqrt(diag(V)),
    percentile = t(apply(reps[ok, , drop = FALSE], 2L, stats::quantile,
                         probs = c(0.025, 0.975), names = FALSE)),
    tau = object$tau, family = object$family, anchor = object$anchor
  ), class = "gkwq_boot")
}

#' @export
print.gkwq_boot <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("\nBootstrap for a %s quantile regression (tau = %s, anchor %s)\n",
              x$family, format(x$tau), x$anchor))
  cat(sprintf("%s bootstrap, %d of %d replicates converged\n\n",
              x$type, x$n_ok, x$R))
  tab <- cbind(Estimate = x$estimate, `Boot SE` = x$se,
               `2.5%` = x$percentile[, 1L], `97.5%` = x$percentile[, 2L])
  print(round(tab, digits))
  invisible(x)
}

#' Likelihood-ratio test for nested quantile regression fits
#'
#' A thin wrapper on [anova.gkwqreg()] under the name `lmtest` users expect. The
#' same guards apply: all fits must share the quantile level and the anchor.
#'
#' @param object,... `"gkwqreg"` fits.
#' @return A data frame of class `"anova.gkwqreg"`.
#' @export
lrtest <- function(object, ...) UseMethod("lrtest")

#' @rdname lrtest
#' @export
lrtest.gkwqreg <- function(object, ...) anova(object, ...)
