#' Diagnostic plots for a quantile regression fit
#'
#' Six diagnostic panels for a `"gkwqreg"` fit. Four of them are the familiar
#' regression diagnostics transferred to a quantile-residual scale; two of them,
#' panels 4 and 5, ask a question that mean regression cannot pose at all, namely
#' whether the fitted quantile is a quantile of the right level.
#'
#' @param x A `"gkwqreg"` fit. The response must have been retained, which is the
#'   default (`y = TRUE` in [gkwqreg()]); every panel needs it, and an informative
#'   error is raised if it is absent.
#' @param which Integer vector selecting panels, any subset of `1:6`. Values
#'   outside that range are dropped; an empty selection returns immediately
#'   without drawing anything. When more than one panel is requested the device is
#'   split into a grid at most three columns wide and the previous
#'   [graphics::par()] settings are restored on exit.
#' @param nsim Number of simulated samples used to build the pointwise envelope in
#'   panel 3. `0` draws the Q-Q plot and its reference line without an envelope.
#'   Ignored when panel 3 is not requested.
#' @param nbins Target number of bins for the calibration panel, panel 5. Bins are
#'   the empirical quantiles of the fitted values, so they hold roughly equal
#'   counts. Ties in the fitted values can collapse breakpoints; if fewer than
#'   three distinct breaks survive, panel 5 is skipped silently rather than drawn
#'   from one or two bins.
#' @param ... Passed to the underlying [graphics::plot()] calls in panels 1, 2, 4,
#'   5 and 6. Panel 3 is drawn by [stats::qqnorm()] and does not receive them.
#'   Note that `main`, `xlab`, `ylab`, `pch`, `cex` and `col` are already supplied
#'   by each panel, so passing them again will raise a duplicated-argument error.
#'
#' @details
#' Panels 1 to 3 and 6 are computed from the randomized-quantile residuals of Dunn
#' and Smyth (1996), `residuals(x, type = "quantile")`, which for this family
#' require no randomization because the distribution is continuous:
#'
#' \deqn{r_i = \Phi^{-1}\left\{ F(y_i; \hat{\theta}_i) \right\}.}{
#'       r_i = qnorm(F(y_i; thetahat_i)).}
#'
#' At the true parameter values these are exactly independent standard normal
#' draws, whatever the family, the level or the link, and at the estimates they
#' are so up to the usual estimation error. Non-finite values, which can arise
#' when an observation falls in a tail so extreme that both tails of \eqn{F}
#' underflow, are set to `NA` and omitted from the panels.
#'
#' @section Panel 1 -- quantile residuals against observation index:
#' Plots \eqn{r_i} against \eqn{i}, with solid and dotted horizontal guides at
#' \eqn{0} and \eqn{\pm 2}{+/- 2}. Roughly five percent of the points should lie
#' outside the dotted lines and they should be scattered without pattern. Read it
#' for serial structure and for clusters: if the rows carry an order that means
#' something (time, space, a grouping) a visible run is evidence that the
#' independence assumption behind the likelihood does not hold, and the standard
#' errors from [vcov.gkwqreg()] are then optimistic. Consider a sandwich or
#' clustered covariance through [estfun.gkwqreg()].
#'
#' @section Panel 2 -- quantile residuals against the linear predictor:
#' Plots \eqn{r_i} against \eqn{\eta_i = x_i^{\top}\hat{\beta}}{eta_i}, the linear
#' predictor of the conditional-quantile part `mu`, with a horizontal line at zero
#' and a `lowess` smooth in red when at least eleven finite points are available.
#' The smooth should be flat. Curvature in it is evidence of a misspecified mean
#' structure for the quantile: a missing nonlinear term in a covariate, or a link
#' that does not suit the data. A fan shape is not read the way it would be in a
#' linear model. Quantile residuals are standard normal by construction under a
#' correct model, whatever the conditional dispersion, so spread that changes with
#' \eqn{\eta}{eta} is not heteroscedasticity in the usual sense; it says that the
#' shape parameters, held constant here, should themselves have been regressed on
#' the covariates. Give them a formula part and compare the two fits with
#' [lrtest()].
#'
#' @section Panel 3 -- normal Q-Q plot with a simulated envelope:
#' Ordered quantile residuals against standard normal plotting positions, with the
#' usual Q-Q reference line and, when `nsim > 0`, a dashed pointwise envelope. The
#' envelope is obtained by drawing `nsim` samples of \eqn{n}{n} independent
#' standard normal variates, sorting each and taking the 2.5th and 97.5th
#' percentiles of each order statistic. Two cautions. It is a *pointwise* envelope,
#' so with many observations some points will stray outside it even under a
#' correct model; and it is built from independent standard normal draws, which
#' means it ignores the uncertainty in \eqn{\hat{\theta}}{thetahat} and is
#' consequently a little narrow in small samples. Systematic departure at one end
#' identifies which tail of the response the family is failing to represent.
#'
#' @section Panel 4 -- observed against fitted quantile, with coverage:
#' Plots \eqn{y_i} against \eqn{\hat{Q}_\tau(x_i)}{Qhat_i} with the 45-degree line,
#' and annotates the empirical coverage
#'
#' \deqn{\widehat{\mathrm{cov}} = \frac{1}{n}\sum_i \mathbf{1}\{y_i \le
#'       \hat{Q}_\tau(x_i)\},}{covhat = mean(y <= fitted),}
#'
#' which should sit at \eqn{\tau}{tau}. This panel has no counterpart in mean
#' regression, where the analogous scatter carries no such calibration statement.
#' The points are *not* expected to lie on the 45-degree line -- at
#' \eqn{\tau = 0.5}{tau = 0.5} half of them should fall on each side, and at
#' \eqn{\tau = 0.9}{tau = 0.9} nine tenths should fall below it. The line is a
#' reference for the *proportion* on each side, not for the distances. Coverage
#' is a single global number and a weak check: it can be exactly right while the
#' fit is badly calibrated in every region of the covariate space, with the errors
#' cancelling. Panel 5 is the panel that detects that.
#'
#' @section Panel 5 -- calibration by fitted value:
#' The panel specific to quantile regression, and the reason to look at this plot
#' rather than a mean-regression one. Under correct specification of the
#' conditional \eqn{\tau}{tau}-quantile,
#'
#' \deqn{E\left[ \mathbf{1}\{Y \le Q_\tau(X)\} \mid X \right] = \tau
#'       \quad \text{exactly, for every } X.}{
#'       E[1{Y <= Q_tau(X)} | X] = tau exactly, for every X.}
#'
#' This is an identity, not an approximation and not an asymptotic result: it
#' follows from the definition of a quantile and holds whatever the shape of the
#' conditional distribution. It therefore supplies a **distribution-free**
#' specification check, available in quantile regression and with no analogue in
#' mean regression, where the corresponding statement about
#' \eqn{E[Y - \hat{\mu}(X) \mid X]}{E[Y - muhat(X) | X]} requires the mean model to
#' be correct in the first place.
#'
#' The panel bins the observations into `nbins` groups of roughly equal size by
#' fitted value and plots, at the mean fitted value of each bin, the bin mean of
#' the `tau`-sign residuals
#' \eqn{\mathbf{1}\{y_i \le \hat{Q}_i\} - \tau}{1{y_i <= Qhat_i} - tau}, which is
#' zero in expectation. Each point carries the interval
#'
#' \deqn{\pm 1.96 \sqrt{\tau(1 - \tau) / n_b},}{+/- 1.96 * sqrt(tau (1 - tau) / n_b),}
#'
#' the normal approximation to the standard error of a mean of \eqn{n_b}{n_b}
#' Bernoulli(\eqn{\tau}{tau}) indicators. A bin whose interval excludes the red
#' zero line is direct evidence that the fit is miscalibrated in that region of the
#' covariate space. The intervals are pointwise, so about one bin in twenty strays
#' by chance and an isolated excursion among the default ten bins is unremarkable;
#' what matters is a *pattern* across adjacent bins, such as a systematic tilt from
#' one end of the fitted range to the other. See the examples, where a mixture response passes the global coverage
#' check of panel 4 and fails this one in half its bins.
#'
#' @section Panel 6 -- influence:
#' A generalized Cook distance, plotted as a spike for each observation:
#'
#' \deqn{D_i = \frac{1}{p}\, s_i^{\top} \hat{V} s_i,}{D_i = t(s_i) V s_i / p,}
#'
#' with \eqn{s_i}{s_i} the observation's score contribution from
#' [estfun.gkwqreg()], \eqn{\hat{V}}{V} the estimated covariance matrix from
#' [vcov.gkwqreg()] and \eqn{p}{p} the number of estimated coefficients. This is
#' the usual one-step approximation to the squared change in
#' \eqn{\hat{\beta}}{betahat} induced by deleting observation \eqn{i}, measured in
#' the metric of \eqn{\hat{V}}{V}. Points above the conventional \eqn{4/n}{4/n}
#' threshold are highlighted and labelled with their row number so that they can be
#' looked up in the data. The threshold is a convention with no distributional
#' basis: use it to rank observations, not to test them. This is the one panel
#' that needs the estimated covariance matrix, so it is unavailable for a fit
#' obtained with `gkwq_control(hessian = FALSE)`.
#'
#' @return `x`, invisibly. Called for the plots it draws.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' *Journal of Computational and Graphical Statistics* **5**, 236-244.
#'
#' @seealso [residuals.gkwqreg()] for the residual types used here, including
#'   `type = "tau-sign"` behind panel 5; [pinball()] for an out-of-sample
#'   criterion; [check_crossing()] for monotonicity across levels.
#'
#' @examples
#' ## ---- A correctly specified fit -------------------------------------------
#' set.seed(1)
#' n <- 400
#' x <- runif(n, -2, 2)
#' mu <- plogis(0.4 + 1.1 * x)
#' y <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
#' d <- data.frame(y = y, x = x)
#' fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
#'
#' ## All six panels at once.
#' plot(fit)
#'
#' ## The number panel 4 annotates: it should sit at tau.
#' mean(d$y <= fitted(fit))
#' #> 0.51        (target 0.50)
#'
#' ## Panel 5 alone, with eight bins. Every interval covers zero here.
#' plot(fit, which = 5, nbins = 8)
#'
#' ## Panel 3 with a denser envelope, and panel 1 with none.
#' plot(fit, which = 3, nsim = 200)
#' plot(fit, which = 1, nsim = 0)
#'
#' ## ---- Why panel 5 is worth more than panel 4 ------------------------------
#' ## A two-component beta mixture is bimodal and lies outside the Generalized
#' ## Kumaraswamy family entirely, so the fit below is badly misspecified.
#' set.seed(2)
#' z <- rbinom(n, 1, plogis(0.5 * x))
#' ymix <- ifelse(z == 1, rbeta(n, 8, 2), rbeta(n, 2, 8))
#' dm <- data.frame(y = pmin(pmax(ymix, 1e-8), 1 - 1e-8), x = x)
#' bad <- suppressWarnings(gkwqreg(y ~ x, data = dm, tau = 0.5, family = "kw"))
#'
#' ## Panel 4 sees nothing wrong: global coverage is essentially perfect.
#' mean(dm$y <= fitted(bad))
#' #> 0.505       (target 0.50)
#'
#' ## Panel 5 does. Four of the eight bins have intervals excluding zero. A
#' ## positive bin mean says too many observations fell at or below the fitted
#' ## quantile there, so the fit sits too HIGH where the fitted median is small
#' ## and too LOW where it is large: the fitted curve is flatter than the truth,
#' ## and the two errors cancel exactly in the global average above.
#' plot(bad, which = c(4, 5), nbins = 8)
#'
#' ## The same information as numbers, which is what the panel draws.
#' rs <- residuals(bad, type = "tau-sign")
#' fv <- fitted(bad)
#' br <- unique(quantile(fv, probs = seq(0, 1, length.out = 9), names = FALSE))
#' g <- cut(fv, br, include.lowest = TRUE)
#' round(cbind(bin_mean = tapply(rs, g, mean),
#'             half_width = 1.96 * sqrt(0.5 * 0.5 / tapply(rs, g, length))), 3)
#' #>               bin_mean half_width
#' #> [0.377,0.409]     0.20      0.139     <- excludes zero
#' #> (0.409,0.437]    -0.08      0.139
#' #> (0.437,0.465]     0.14      0.139     <- excludes zero
#' #> ...
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
