## ---------------------------------------------------------------------------
## Model comparison.
##
## One guard runs through all of it: models compared by a likelihood-ratio test
## must share the same quantile level AND the same anchor. Different levels are
## different likelihoods; different anchors are non-nested models of equal
## dimension (docs/adr/0001), so an LR test on either is meaningless even though
## the arithmetic would happily produce a number.
## ---------------------------------------------------------------------------

.gkwq_check_comparable <- function(objects, what = "compared") {
  taus <- vapply(objects, function(o) o$tau, numeric(1))
  if (length(unique(taus)) > 1L) {
    stop("models at different quantile levels cannot be ", what,
         ": each level is a separate likelihood. Levels seen: ",
         paste(format(unique(taus)), collapse = ", "), ".", call. = FALSE)
  }
  anc <- vapply(objects, function(o) o$anchor, character(1))
  if (length(unique(anc)) > 1L) {
    stop("models with different anchors cannot be ", what,
         " by a likelihood-ratio test: they are non-nested models of equal ",
         "dimension. Compare them by AIC, BIC or a Vuong test instead. ",
         "Anchors seen: ", paste(unique(anc), collapse = ", "), ".",
         call. = FALSE)
  }
  n <- vapply(objects, function(o) o$nobs, numeric(1))
  if (length(unique(n)) > 1L) {
    stop("models fitted to different numbers of observations cannot be ", what,
         ".", call. = FALSE)
  }
  invisible(TRUE)
}

#' Likelihood-ratio comparison of nested quantile regression fits
#'
#' The seven families form a genuine nesting (`kw` inside `ekw` inside `gkw`,
#' and so on), so family selection is an ordinary likelihood-ratio test rather
#' than the Vuong test that non-nested collections require.
#'
#' @param object A `"gkwqreg"` fit.
#' @param ... Further `"gkwqreg"` fits, at the same quantile level and anchor.
#' @return A data frame of class `"anova.gkwqreg"`.
#' @export
anova.gkwqreg <- function(object, ...) {
  objects <- c(list(object), list(...))
  if (length(objects) < 2L) {
    stop("anova() needs at least two fits to compare.", call. = FALSE)
  }
  ok <- vapply(objects, inherits, logical(1), "gkwqreg")
  if (!all(ok)) stop("all arguments must be \"gkwqreg\" fits.", call. = FALSE)
  .gkwq_check_comparable(objects)

  ord <- order(vapply(objects, function(o) o$npar, numeric(1)))
  objects <- objects[ord]

  ll <- vapply(objects, function(o) o$loglik, numeric(1))
  df <- vapply(objects, function(o) o$npar, numeric(1))
  fam <- vapply(objects, function(o) o$family, character(1))

  stat <- c(NA_real_, 2 * diff(ll))
  ddf <- c(NA_real_, diff(df))
  p <- ifelse(is.na(stat) | ddf <= 0, NA_real_,
              stats::pchisq(stat, pmax(ddf, 1), lower.tail = FALSE))

  out <- data.frame(family = fam, Df = df, logLik = ll, AIC = -2 * ll + 2 * df,
                    Chisq = stat, `Chi Df` = ddf, `Pr(>Chisq)` = p,
                    check.names = FALSE)
  attr(out, "heading") <- c(
    "Likelihood-ratio test for Generalized Kumaraswamy quantile regression\n",
    sprintf("tau = %s, anchor = %s\n", format(objects[[1L]]$tau),
            objects[[1L]]$anchor)
  )
  if (any(!is.na(ddf) & ddf <= 0)) {
    warning("some compared models have the same dimension; a likelihood-ratio ",
            "test does not apply to them. Use AIC or a Vuong test.",
            call. = FALSE)
  }
  class(out) <- c("anova.gkwqreg", "anova", "data.frame")
  out
}

#' @export
print.anova.gkwqreg <- function(x, ...) {
  cat(attr(x, "heading"), sep = "")
  stats::printCoefmat(as.matrix(x[, -1L, drop = FALSE]), na.print = "",
                      has.Pvalue = TRUE, P.values = TRUE, cs.ind = NULL,
                      tst.ind = 4L, zap.ind = 1:3, digits = 5)
  invisible(x)
}

#' Vuong test for two non-nested quantile regression fits
#'
#' The right tool for comparing two anchors on the same data, or two families
#' that are not nested. Both give models of the same dimension fitted to the
#' same responses, so a likelihood-ratio test does not apply.
#'
#' @param object,object2 Two `"gkwqreg"` fits at the same quantile level, over
#'   the same observations.
#' @param correction Apply the BIC-style dimension correction.
#' @return An object of class `"gkwq_vuong"`.
#' @references
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. *Econometrica* **57**, 307-333.
#' @export
vuong_test <- function(object, object2, correction = TRUE) {
  stopifnot(inherits(object, "gkwqreg"), inherits(object2, "gkwqreg"))
  if (!isTRUE(all.equal(object$tau, object2$tau))) {
    stop("both fits must use the same quantile level.", call. = FALSE)
  }
  if (object$nobs != object2$nobs) {
    stop("both fits must use the same observations.", call. = FALSE)
  }
  n <- object$nobs
  m <- object$loglik_i - object2$loglik_i
  if (correction) {
    m <- m - (object$npar - object2$npar) * log(n) / (2 * n)
  }
  s <- stats::sd(m)
  stat <- sqrt(n) * mean(m) / s
  p <- 2 * stats::pnorm(-abs(stat))
  structure(list(statistic = stat, p.value = p, n = n,
                 model1 = sprintf("%s / anchor %s", object$family, object$anchor),
                 model2 = sprintf("%s / anchor %s", object2$family, object2$anchor),
                 tau = object$tau, correction = correction),
            class = "gkwq_vuong")
}

#' @export
print.gkwq_vuong <- function(x, digits = 4, ...) {
  cat("\nVuong test for non-nested quantile regression fits\n")
  cat(sprintf("tau = %s%s\n", format(x$tau),
              if (x$correction) ", BIC-corrected" else ""))
  cat(sprintf("  model 1: %s\n  model 2: %s\n", x$model1, x$model2))
  cat(sprintf("\n  z = %s, p-value = %s\n", format(x$statistic, digits = digits),
              base::format.pval(x$p.value, digits = digits)))
  cat(sprintf("  %s\n", if (x$p.value > 0.05) "neither model is favoured" else
    if (x$statistic > 0) "model 1 is favoured" else "model 2 is favoured"))
  invisible(x)
}

#' Compare all seven families at one quantile level
#'
#' Refits the same formula under every family and ranks them. Reports check loss
#' alongside AIC, because check loss is the criterion a quantile estimate
#' actually targets, while AIC compares likelihoods across parametrizations that
#' are not all nested.
#'
#' @param object A `"gkwqreg"` fit whose call is reused.
#' @param families Families to try.
#' @param ... Passed to [gkwqreg()].
#' @return A data frame ordered by AIC.
#' @export
compare_families <- function(object,
                             families = c("kw", "ekw", "kkw", "bkw", "gkw",
                                          "mc", "beta"),
                             ...) {
  stopifnot(inherits(object, "gkwqreg"))
  ## Capture the evaluation environment ONCE, here. Reaching for
  ## parent.frame(3L) from inside lapply() depends on the call stack depth and
  ## breaks the moment this is called through another function.
  env <- parent.frame()
  rows <- lapply(families, function(fm) {
    cl <- object$call
    cl$family <- fm
    cl$anchor <- NULL          # each family gets its own default anchor
    f <- suppressWarnings(try(eval(cl, env), silent = TRUE))
    if (inherits(f, "try-error")) {
      return(data.frame(family = fm, anchor = NA_character_, df = NA_integer_,
                        logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                        pinball = NA_real_, converged = FALSE))
    }
    data.frame(family = fm, anchor = f$anchor, df = f$npar, logLik = f$loglik,
               AIC = f$aic, BIC = f$bic, pinball = f$pinball,
               converged = isTRUE(f$convergence == 0))
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$AIC, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
