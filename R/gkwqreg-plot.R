#' Diagnostic plots for a quantile regression fit
#'
#' @param x A `"gkwqreg"` fit.
#' @param which Panels to draw, from 1 to 6.
#' @param nsim Number of simulated envelopes for the Q-Q panel; `0` disables it.
#' @param nbins Bins for the calibration panel.
#' @param ... Passed to the underlying plot calls.
#'
#' @details
#' Panels 4 and 5 have no counterpart in mean regression.
#'
#' Panel 4 plots observed against fitted quantiles and annotates the empirical
#' coverage `mean(y <= fitted)`, which should sit at `tau`.
#'
#' Panel 5 is the calibration plot: `tau`-sign residuals binned by fitted value.
#' Under correct specification `E[1{Y <= Q_tau(X)} | X] = tau` exactly, so any
#' bin whose mean strays outside the band is direct evidence of miscalibration,
#' and the evidence is distribution-free.
#'
#' @return `x`, invisibly.
#' @export
plot.gkwqreg <- function(x, which = 1:6, nsim = 100L, nbins = 10L, ...) {
  which <- intersect(which, 1:6)
  if (!length(which)) return(invisible(x))
  if (is.null(x$y)) {
    stop("the diagnostic plots all need the response; refit with y = TRUE.",
         call. = FALSE)
  }
  show <- rep(FALSE, 6L); show[which] <- TRUE

  n_panel <- length(which)
  if (n_panel > 1L) {
    nc <- min(3L, n_panel); nr <- ceiling(n_panel / nc)
    op <- graphics::par(mfrow = c(nr, nc), mar = c(4.2, 4.2, 3, 1))
    on.exit(graphics::par(op), add = TRUE)
  }

  rq <- residuals(x, type = "quantile")
  rq[!is.finite(rq)] <- NA_real_
  fitv <- x$fitted.values
  y <- x$y
  tau <- x$tau

  if (show[1L]) {
    plot(seq_along(rq), rq, xlab = "observation", ylab = "quantile residual",
         main = "Quantile residuals vs index", pch = 19, cex = 0.5,
         col = "grey30", ...)
    graphics::abline(h = c(-2, 0, 2), lty = c(3, 1, 3), col = c("grey60", "black", "grey60"))
  }

  if (show[2L]) {
    eta <- x$linear.predictors[["mu"]]
    plot(eta, rq, xlab = expression(paste("linear predictor of ", mu[tau])),
         ylab = "quantile residual",
         main = "Quantile residuals vs linear predictor", pch = 19, cex = 0.5,
         col = "grey30", ...)
    graphics::abline(h = 0, lty = 2)
    ok <- is.finite(eta) & is.finite(rq)
    if (sum(ok) > 10L) {
      graphics::lines(stats::lowess(eta[ok], rq[ok]), col = "firebrick", lwd = 2)
    }
  }

  if (show[3L]) {
    stats::qqnorm(rq, main = "Normal Q-Q of quantile residuals", pch = 19,
                  cex = 0.5, col = "grey30")
    stats::qqline(rq, col = "firebrick", lwd = 2)
    if (nsim > 0L) {
      nn <- sum(is.finite(rq))
      env <- replicate(nsim, sort(stats::rnorm(nn)))
      lo <- apply(env, 1L, stats::quantile, 0.025)
      hi <- apply(env, 1L, stats::quantile, 0.975)
      qs <- stats::qnorm(stats::ppoints(nn))
      graphics::lines(qs, lo, lty = 2, col = "steelblue")
      graphics::lines(qs, hi, lty = 2, col = "steelblue")
    }
  }

  if (show[4L]) {
    cov_emp <- mean(y <= fitv)
    plot(fitv, y, xlab = sprintf("fitted %s-quantile", format(tau)),
         ylab = "observed", main = "Observed vs fitted quantile",
         pch = 19, cex = 0.5, col = "grey30", ...)
    graphics::abline(0, 1, col = "firebrick", lwd = 2)
    graphics::legend("topleft", bty = "n", cex = 0.85,
                     legend = sprintf("coverage %.3f (target %.2f)", cov_emp, tau))
  }

  if (show[5L]) {
    rs <- as.numeric(y <= fitv) - tau
    br <- stats::quantile(fitv, probs = seq(0, 1, length.out = nbins + 1L),
                          names = FALSE)
    br <- unique(br)
    if (length(br) > 2L) {
      g <- base::cut(fitv, br, include.lowest = TRUE)
      m <- tapply(rs, g, mean)
      cnt <- tapply(rs, g, length)
      ctr <- tapply(fitv, g, mean)
      ## Binomial band: sd of the bin mean under correct calibration.
      band <- 1.96 * sqrt(tau * (1 - tau) / cnt)
      yl <- range(c(m - band, m + band, 0), na.rm = TRUE)
      plot(ctr, m, ylim = yl, xlab = sprintf("fitted %s-quantile", format(tau)),
           ylab = expression(paste("mean of  ", 1, "{y" <= "Q} - ", tau)),
           main = "Calibration by fitted value", pch = 19, col = "steelblue4",
           ...)
      graphics::arrows(ctr, m - band, ctr, m + band, angle = 90, code = 3,
                       length = 0.03, col = "steelblue")
      graphics::abline(h = 0, col = "firebrick", lwd = 2)
    }
  }

  if (show[6L]) {
    D <- tryCatch({
      V <- vcov(x)
      S <- estfun.gkwqreg(x)
      rowSums((S %*% V) * S) / x$npar
    }, error = function(e) rep(NA_real_, x$nobs))
    plot(seq_along(D), D, type = "h", xlab = "observation",
         ylab = "generalized Cook distance", main = "Influence", col = "grey40",
         ...)
    if (all(is.finite(D))) {
      big <- which(D > 4 / x$nobs)
      if (length(big)) {
        graphics::points(big, D[big], pch = 19, cex = 0.6, col = "firebrick")
        graphics::text(big, D[big], labels = big, pos = 3, cex = 0.6)
      }
    }
  }

  invisible(x)
}
