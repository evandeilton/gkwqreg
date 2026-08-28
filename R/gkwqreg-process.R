## ---------------------------------------------------------------------------
## Quantile-regression-specific tooling: the process across levels, crossing,
## rearrangement, and check loss.
## ---------------------------------------------------------------------------

#' Pinball (check) loss
#'
#' The loss a quantile estimate actually targets, and the right criterion for
#' choosing a family: AIC compares likelihoods across parametrizations that are
#' not all nested, whereas check loss compares the thing being estimated. Use it
#' out of sample.
#'
#' @param object A `"gkwqreg"` or `"gkwqregs"` fit.
#' @param newdata Optional data to evaluate on; defaults to the fitting data.
#' @param y Optional response for `newdata`; taken from `newdata` if omitted.
#' @return The mean check loss, or one value per quantile level.
#' @export
pinball <- function(object, newdata = NULL, y = NULL) {
  if (inherits(object, "gkwqregs")) {
    return(vapply(object$fits, pinball, numeric(1), newdata = newdata, y = y))
  }
  stopifnot(inherits(object, "gkwqreg"))
  if (is.null(newdata)) {
    q <- object$fitted.values
    yy <- object$y
    w <- object$weights
  } else {
    q <- predict(object, newdata = newdata, type = "quantile")
    yy <- y
    if (is.null(yy)) {
      resp <- all.vars(stats::formula(object$formula, lhs = 1, rhs = 0))[1L]
      if (is.null(resp) || !resp %in% names(newdata)) {
        stop("`newdata` has no response column (", sQuote(resp %||% "?"),
             "); supply `y` instead.", call. = FALSE)
      }
      yy <- newdata[[resp]]
    }
    w <- rep(1, length(yy))
  }
  e <- yy - q
  sum(w * e * (object$tau - (e < 0))) / sum(w)
}

#' Detect quantile crossing
#'
#' @param object A `"gkwqregs"` container, a list of `"gkwqreg"` fits, or a
#'   single `"gkwqreg"` fit.
#' @param newdata Optional data to evaluate the quantiles on.
#' @param taus For a single fit, the levels to read off the fitted distribution.
#' @param tol Tolerance below which a decrease is not counted as a crossing.
#' @param ... Unused.
#'
#' @details
#' There are two genuinely different questions here.
#'
#' Across **independently fitted** levels, nothing forces the fitted quantiles to
#' increase with `tau`, and crossing is the classic failure of per-level quantile
#' regression. That is what this reports for a `"gkwqregs"` container.
#'
#' Within a **single** fit, the quantiles read off one fitted distribution cannot
#' cross: they are quantiles of a proper distribution function. Reporting zero
#' crossings there is not vacuous, it is the argument for parametric quantile
#' regression over the check-function kind, and the printed output says so.
#'
#' @return An object of class `"gkwq_crossing"`.
#' @references
#' Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile and
#' probability curves without crossing. *Econometrica* **78**, 1093-1125.
#' @export
check_crossing <- function(object, newdata = NULL, taus = NULL, tol = 0, ...) {
  if (inherits(object, "gkwqregs")) {
    fits <- object$fits
    tv <- object$taus
    mode <- "separate"
  } else if (is.list(object) && all(vapply(object, inherits, logical(1), "gkwqreg"))) {
    tv <- vapply(object, function(f) f$tau, numeric(1))
    ord <- order(tv)
    fits <- object[ord]
    tv <- tv[ord]
    mode <- "separate"
  } else if (inherits(object, "gkwqreg")) {
    tv <- sort(taus %||% seq(0.05, 0.95, by = 0.05))
    Q <- predict(object, newdata = newdata, type = "quantile", tau = tv)
    mode <- "implied"
    fits <- NULL
  } else {
    stop("`object` must be a \"gkwqregs\", a list of \"gkwqreg\" fits, or one fit.",
         call. = FALSE)
  }

  if (mode == "separate") {
    if (length(tv) < 2L) {
      stop("crossing needs at least two quantile levels.", call. = FALSE)
    }
    Q <- vapply(fits, function(f) predict(f, newdata = newdata,
                                          type = "quantile"),
                numeric(if (is.null(newdata)) fits[[1L]]$nobs else nrow(newdata)))
    if (is.null(dim(Q))) Q <- matrix(Q, ncol = length(tv))
  }
  colnames(Q) <- format(tv, trim = TRUE)

  D <- t(apply(Q, 1L, diff))
  if (is.null(dim(D))) D <- matrix(D, nrow = nrow(Q))
  viol <- D < -tol
  rows <- which(apply(viol, 1L, any))
  pairs <- data.frame(tau_lo = tv[-length(tv)], tau_hi = tv[-1L],
                      n = colSums(viol))

  structure(list(taus = tv, Q = Q, mode = mode,
                 n_crossing = length(rows), frac = length(rows) / nrow(Q),
                 which = rows,
                 worst = if (any(viol)) max(-D[viol]) else 0,
                 pairs = pairs),
            class = "gkwq_crossing")
}

#' @rdname check_crossing
#' @export
qcrossing <- check_crossing

#' @export
print.gkwq_crossing <- function(x, ...) {
  cat("\nQuantile crossing check\n")
  cat(sprintf("levels: %s\n", paste(format(x$taus, trim = TRUE), collapse = ", ")))
  if (x$mode == "implied") {
    cat("source: quantiles implied by ONE fitted distribution\n")
  } else {
    cat("source: independently fitted models, one per level\n")
  }
  cat(sprintf("\n  rows with a crossing: %d of %d (%.2f%%)\n",
              x$n_crossing, nrow(x$Q), 100 * x$frac))
  cat(sprintf("  worst violation     : %.3e\n", x$worst))
  if (x$mode == "implied" && x$n_crossing == 0L) {
    cat("\n  Zero crossings is guaranteed here, not luck: these are quantiles of\n")
    cat("  a single proper distribution function. That guarantee is what a\n")
    cat("  parametric fit buys over independently fitted levels.\n")
  }
  if (x$mode == "separate" && x$n_crossing > 0L) {
    cat("\n  Independently fitted levels carry no monotonicity guarantee.\n")
    cat("  rearrange() applies the Chernozhukov et al. (2010) fix.\n")
  }
  invisible(x)
}

#' Monotone rearrangement of crossing quantiles
#'
#' Sorts each row of the fitted quantile matrix, the rearrangement of
#' Chernozhukov, Fernandez-Val and Galichon (2010). It cannot make a fit worse in
#' the check-loss sense and it restores monotonicity by construction.
#'
#' @param object A `"gkwqregs"` container or a list of `"gkwqreg"` fits.
#' @param newdata Optional data to evaluate on.
#' @param ... Unused.
#' @return A matrix of rearranged quantiles, with a `"crossing"` attribute
#'   recording what was fixed.
#' @export
rearrange <- function(object, newdata = NULL, ...) {
  if (inherits(object, "gkwqreg")) {
    message("a single fit has no crossing to rearrange: its quantiles are those ",
            "of one distribution function. Returning them unchanged.")
    cr <- check_crossing(object, newdata = newdata, ...)
    return(structure(cr$Q, crossing = cr))
  }
  cr <- check_crossing(object, newdata = newdata, ...)
  Q <- cr$Q
  R <- t(apply(Q, 1L, sort))
  dimnames(R) <- dimnames(Q)
  structure(R, crossing = cr)
}

#' The quantile process
#'
#' Fits the model over a grid of quantile levels and collects the coefficient
#' paths. The plot is deliberately in the idiom of `quantreg`'s
#' `plot.summary.rqs`, because that picture is the literature's standard output
#' and readers should recognise it at once.
#'
#' @param object A `"gkwqreg"` fit whose call is reused, or a `"gkwqregs"`
#'   container to harvest.
#' @param taus Grid of quantile levels.
#' @param level Confidence level for the bands.
#' @param ... Passed to [gkwqreg()].
#' @return An object of class `"gkwq_process"`.
#' @export
quantile_process <- function(object, taus = seq(0.05, 0.95, by = 0.05),
                             level = 0.95, ...) {
  if (inherits(object, "gkwqregs")) {
    fits <- object$fits
    taus <- object$taus
  } else {
    stopifnot(inherits(object, "gkwqreg"))
    cl <- object$call
    cl$tau <- taus
    extra <- list(...)
    for (nm in names(extra)) cl[[nm]] <- extra[[nm]]
    many <- eval(cl, parent.frame())
    fits <- many$fits
    taus <- many$taus
  }

  cf <- vapply(fits, function(f) f$coefficients, numeric(length(fits[[1L]]$coefficients)))
  se <- vapply(fits, function(f) f$se, numeric(length(fits[[1L]]$coefficients)))
  if (is.null(dim(cf))) {
    cf <- matrix(cf, nrow = 1L); se <- matrix(se, nrow = 1L)
  }
  rownames(cf) <- rownames(se) <- names(fits[[1L]]$coefficients)
  colnames(cf) <- colnames(se) <- format(taus, trim = TRUE)
  z <- stats::qnorm(1 - (1 - level) / 2)

  structure(list(taus = taus, coef = cf, se = se,
                 lower = cf - z * se, upper = cf + z * se,
                 parts = sub(":.*$", "", rownames(cf)),
                 level = level, family = fits[[1L]]$family,
                 anchor = fits[[1L]]$anchor, fits = fits,
                 call = match.call()),
            class = "gkwq_process")
}

#' @export
print.gkwq_process <- function(x, ...) {
  cat(sprintf("\nQuantile process: family %s, anchor %s, %d levels from %s to %s\n",
              x$family, x$anchor, length(x$taus), format(min(x$taus)),
              format(max(x$taus))))
  cat("\nCoefficients by tau:\n")
  print(round(x$coef, 4))
  invisible(x)
}

#' @export
plot.gkwq_process <- function(x, parm = NULL, parts = "mu", nrow = NULL,
                              ncol = NULL, reference = NULL, ...) {
  keep <- if (is.null(parm)) which(x$parts %in% parts) else match(parm, rownames(x$coef))
  keep <- keep[!is.na(keep)]
  if (!length(keep)) stop("no coefficients selected.", call. = FALSE)

  k <- length(keep)
  if (is.null(nrow) && is.null(ncol)) {
    ncol <- min(3L, k); nrow <- ceiling(k / ncol)
  } else if (is.null(nrow)) nrow <- ceiling(k / ncol) else ncol <- ceiling(k / nrow)

  op <- graphics::par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  for (i in keep) {
    yl <- range(c(x$lower[i, ], x$upper[i, ], reference[i]), na.rm = TRUE)
    plot(x$taus, x$coef[i, ], type = "n", ylim = yl, xlab = expression(tau),
         ylab = "coefficient", main = rownames(x$coef)[i])
    graphics::polygon(c(x$taus, rev(x$taus)),
                      c(x$lower[i, ], rev(x$upper[i, ])),
                      col = grDevices::adjustcolor("steelblue", 0.25),
                      border = NA)
    graphics::lines(x$taus, x$coef[i, ], col = "steelblue4", lwd = 2)
    graphics::points(x$taus, x$coef[i, ], pch = 19, cex = 0.6, col = "steelblue4")
    graphics::abline(h = 0, lty = 3, col = "grey50")
    if (!is.null(reference) && !is.na(reference[i])) {
      graphics::abline(h = reference[i], lty = 2, col = "firebrick")
    }
  }
  invisible(x)
}

#' @export
plot.gkwqregs <- function(x, ...) plot(quantile_process(x), ...)
