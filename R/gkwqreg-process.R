## ---------------------------------------------------------------------------
## Quantile-regression-specific tooling: the process across levels, crossing,
## rearrangement, and check loss.
## ---------------------------------------------------------------------------

#' Pinball (check) loss of a fitted quantile regression
#'
#' Evaluates, on the fitting data or on a held-out sample, the loss function that
#' a quantile estimate actually targets. This is the criterion to use when the
#' question is which family to fit and the answer will be acted on as a statement
#' about a quantile.
#'
#' @param object A `"gkwqreg"` fit, or a `"gkwqregs"` container of fits at
#'   several quantile levels.
#' @param newdata Optional data frame on which to evaluate the loss. When
#'   supplied, the fitted quantiles are obtained by
#'   `predict(object, newdata, type = "quantile")`, so `newdata` must contain
#'   every covariate used by any part of the model. When `NULL` (the default) the
#'   loss is computed in sample, from `object$fitted.values` and the stored
#'   response, and is then identical to the value cached in `object$pinball`.
#' @param y Optional numeric response to pair with `newdata`. If omitted, the
#'   response is taken from the column of `newdata` named on the left-hand side
#'   of the model formula; if that column is absent, an error asks for `y`
#'   explicitly. Ignored when `newdata` is `NULL`.
#' @return A single number for a `"gkwqreg"` fit: the weighted mean check loss,
#'   on the scale of the response. For a `"gkwqregs"` container, a numeric vector
#'   with one entry per quantile level, named as the container's elements are
#'   (`"tau=0.1"` and so on).
#'
#' @details
#' Write \eqn{\hat{Q}_\tau(x_i)}{Qhat_i} for the fitted conditional
#' \eqn{\tau}{tau}-quantile of observation \eqn{i} and
#' \eqn{u_i = y_i - \hat{Q}_\tau(x_i)}{u_i = y_i - Qhat_i} for its residual. The
#' pinball, or check, loss is the asymmetrically weighted absolute deviation
#'
#' \deqn{\rho_\tau(u) \;=\; u \, \left( \tau - \mathbf{1}\{u < 0\} \right),}{
#'       rho_tau(u) = u * (tau - 1{u < 0}),}
#'
#' equal to \eqn{\tau u}{tau * u} when the observation lies above the fitted
#' quantile and to \eqn{-(1 - \tau) u}{-(1 - tau) * u} when it lies below.
#' `pinball()` returns the weighted sample mean of that loss,
#'
#' \deqn{L_\tau \;=\; \frac{\sum_i w_i \, \rho_\tau(y_i - \hat{Q}_\tau(x_i))}{\sum_i w_i},}{
#'       L_tau = sum_i w_i rho_tau(y_i - Qhat_i) / sum_i w_i,}
#'
#' with \eqn{w_i}{w_i} the prior weights supplied to [gkwqreg()]. Weights apply
#' only to the in-sample calculation; when `newdata` is given every row counts
#' once.
#'
#' Underprediction is charged \eqn{\tau}{tau} per unit and overprediction
#' \eqn{1 - \tau}{1 - tau} per unit, so at \eqn{\tau = 0.9}{tau = 0.9} falling
#' short of an observation costs nine times as much as overshooting it by the
#' same amount. The loss carries the units of the response and is therefore
#' numerically small for a response confined to \eqn{(0,1)}{(0,1)}; only
#' *differences* between competing fits evaluated on the *same* rows are
#' interpretable.
#'
#' @section Why check loss and not AIC:
#' The population minimiser of \eqn{E[\rho_\tau(Y - q)]}{E[rho_tau(Y - q)]} over
#' \eqn{q} is the \eqn{\tau}{tau}-quantile of \eqn{Y}. This is the defining
#' property of the check function (Koenker and Bassett, 1978) and the reason it is
#' the natural scoring rule for a quantile forecast (Gneiting, 2011). Check loss
#' therefore scores exactly the quantity this package estimates, and nothing else.
#'
#' AIC scores something materially different: the whole conditional density. A
#' family can describe the bulk of the distribution well, be rewarded for it by
#' the likelihood, and still place the quantile of interest worse than a rival
#' that fits the bulk less well. **The two criteria can rank families
#' differently, and when they do it is not a symptom of anything having gone
#' wrong.** In the case study shipped with this package
#' (`vignette("gkwqreg-case-study")`) AIC ranks family `"beta"` first while
#' out-of-sample check loss ranks it last, behind every Kumaraswamy member. When
#' the inferential target is a quantile, act on the check loss; when it is the
#' conditional distribution as a whole, AIC is the relevant criterion. Deciding
#' which of the two applies is a modelling decision, not a computational one.
#'
#' @section Why out of sample:
#' In-sample check loss is not a model-selection criterion. It carries no penalty
#' for dimension, so a family with more free parameters can track the observed
#' rows more closely whether or not the extra flexibility corresponds to anything
#' real, and comparing it with a sub-family in sample mostly rewards that
#' flexibility. Evaluate on rows the model has not seen: pass a held-out
#' `newdata`, or average `pinball()` over cross-validation folds. The in-sample
#' value remains useful as a description of fit -- and it is what
#' [compare_families()] reports, alongside AIC and BIC, which do carry penalties
#' -- but it should not be the basis for choosing between families.
#'
#' @references
#' Koenker, R. and Bassett, G. (1978). Regression quantiles. *Econometrica*.
#'
#' Gneiting, T. (2011). Making and evaluating point forecasts.
#' *Journal of the American Statistical Association*.
#'
#' @seealso [compare_families()], which tabulates AIC, BIC and the in-sample
#'   check loss for every family at once; [residuals.gkwqreg()] with
#'   `type = "check"` for the per-observation contributions to \eqn{L_\tau}{L_tau};
#'   [vuong_test()] for a formal test between two non-nested families.
#'
#' @examples
#' ## A response in (0,1) whose dispersion, not only its location, depends on x.
#' ## Fitting `y ~ x` regresses the quantile alone and leaves the shape constant,
#' ## so every family below is misspecified, to a different degree and in a
#' ## different way -- which is exactly when the choice of criterion matters.
#' set.seed(1)
#' n <- 600
#' x <- runif(n, -2, 2)
#' mu <- plogis(0.3 + 1.0 * x)
#' phi <- exp(1 + 0.8 * x)
#' y <- rbeta(n, phi * mu, phi * (1 - mu))
#' d <- data.frame(y = y, x = x)
#' train <- d[1:300, ]
#' test <- d[301:600, ]
#'
#' fit <- gkwqreg(y ~ x, data = train, tau = 0.5, family = "kw")
#'
#' ## In sample this is the value the likelihood machinery already cached.
#' c(pinball(fit), fit$pinball)
#' #> 0.1116671 0.1116671
#'
#' ## Out of sample, and identical to the definition applied by hand.
#' e <- test$y - predict(fit, newdata = test, type = "quantile")
#' c(pinball(fit, newdata = test), mean(e * (0.5 - (e < 0))))
#' #> 0.1034218 0.1034218
#'
#' ## AIC and out-of-sample check loss need not agree.
#' tab <- t(sapply(c("kw", "ekw", "bkw", "beta", "mc"), function(f) {
#'   m <- suppressWarnings(gkwqreg(y ~ x, data = train, tau = 0.5, family = f))
#'   c(AIC = AIC(m), in_sample = pinball(m), out_of_sample = pinball(m, newdata = test))
#' }))
#' round(tab, 4)
#' #>            AIC in_sample out_of_sample
#' #> kw   -291.5234    0.1117        0.1034
#' #> ekw  -480.8364    0.1003        0.0934
#' #> bkw  -417.9420    0.1010        0.0938
#' #> beta -484.7745    0.1012        0.0940
#' #> mc   -482.7745    0.1012        0.0940
#'
#' rownames(tab)[which.min(tab[, "AIC"])]            # "beta"
#' rownames(tab)[which.min(tab[, "out_of_sample"])]  # "ekw"
#'
#' ## AIC prefers "beta"; the median it delivers is the fourth best of five.
#' ## For a question about the median, "ekw" is the fit to report.
#'
#' ## One value per level from a container of independent fits.
#' fits <- gkwqreg(y ~ x, data = train, tau = c(0.1, 0.5, 0.9), family = "kw")
#' round(pinball(fits, newdata = test), 4)
#' #> tau=0.1 tau=0.5 tau=0.9
#' #>  0.0550  0.1034  0.0389
#' ## The loss shrinks towards the tails because rho_tau is bounded by the
#' ## smaller of tau and 1 - tau; these three numbers are NOT comparable with
#' ## one another, only with a rival fit at the same level.
#'
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
#' Checks whether a set of fitted conditional quantiles is monotone in the
#' quantile level, row by row. A quantile function must be non-decreasing in
#' \eqn{\tau}{tau}; an estimate need not be, and where it is not the fitted object
#' is not a conditional distribution at all. `qcrossing()` is an alias.
#'
#' @param object One of three things, and the answer means something different in
#'   each case (see Details). A `"gkwqregs"` container returned by [gkwqreg()]
#'   with a vector `tau`; a plain list of `"gkwqreg"` fits, which is sorted by
#'   level before checking; or a single `"gkwqreg"` fit, whose implied quantiles
#'   are read off its one fitted conditional distribution.
#' @param newdata Optional data frame on which to evaluate the quantiles. It must
#'   supply every covariate used by any part of the model. Extrapolating beyond
#'   the observed covariate range is legitimate here and often revealing: see the
#'   examples.
#' @param taus Used only when `object` is a single fit: the levels at which to
#'   read the fitted distribution. Defaults to `seq(0.05, 0.95, by = 0.05)`.
#'   Ignored for a container or a list of fits, whose levels are already fixed by
#'   the fits themselves.
#' @param tol Non-negative numeric tolerance. A decrease between adjacent levels
#'   counts as a crossing only if it exceeds `tol` in absolute value. The default
#'   `0` flags any decrease at all; a small positive value (say `1e-8`) is
#'   appropriate if you wish to ignore decreases attributable to floating-point
#'   arithmetic rather than to the fit.
#' @param ... Currently unused, accepted so that the alias and future methods keep
#'   a stable signature.
#'
#' @details
#' Let \eqn{\hat{Q}_{ij}}{Qhat_ij} be the fitted quantile for row \eqn{i} at the
#' \eqn{j}th level of an increasing grid
#' \eqn{\tau_1 < \cdots < \tau_m}{tau_1 < ... < tau_m}. The function forms the
#' adjacent differences
#' \eqn{D_{ij} = \hat{Q}_{i,j+1} - \hat{Q}_{ij}}{D_ij = Qhat_i,j+1 - Qhat_ij} and
#' records a crossing at \eqn{(i,j)}{(i,j)} whenever
#' \eqn{D_{ij} < -\texttt{tol}}{D_ij < -tol}. A row is counted once, however many
#' of its adjacent pairs are violated, so `n_crossing` and `frac` describe
#' *observations affected*, while `pairs` localises the damage in \eqn{\tau}{tau}.
#'
#' @section Two genuinely different questions:
#' The word "crossing" covers two situations that this function deliberately
#' treats together and reports separately.
#'
#' **Across independently fitted levels** (`mode = "separate"`). Each level is a
#' separate optimisation with its own likelihood and its own coefficient vector.
#' Nothing in the estimation ties them together, so nothing forces
#' \eqn{\hat{Q}_\tau(x)}{Qhat_tau(x)} to increase with \eqn{\tau}{tau} at any
#' particular \eqn{x}{x}. This is the classic and well-documented failure of
#' per-level quantile regression, and it is what the check reports for a
#' `"gkwqregs"` container or a list of fits. It typically bites where the fitted
#' lines are least constrained: in sparse regions of the covariate space, near
#' the extreme levels, and above all under extrapolation.
#'
#' **Within a single fit** (`mode = "implied"`). Here the quantiles at every level
#' are read off one fitted conditional distribution function, through
#' \eqn{\hat{Q}(\tau \mid x) = F^{-1}(\tau \mid x)}{Qhat(tau | x) = Finv(tau | x)}
#' with a single parameter vector. Because \eqn{F(\cdot \mid x)}{F(. | x)} is a
#' proper distribution function, its inverse is non-decreasing by construction and
#' crossing is arithmetically impossible.
#'
#' Reporting zero crossings in the second case is therefore not a vacuous
#' tautology dressed up as a result. It is precisely the argument for parametric
#' quantile regression over the check-function kind: the guarantee holds at every
#' covariate value, in sample and out, at every level, without any post-hoc repair.
#' The printed output says so explicitly rather than letting a reader mistake it
#' for good luck. Verified in this package, ten levels on the same data: `0` of
#' the rows cross within a single `gkwqreg` fit, against `4.3` percent of rows for
#' `quantreg::rq` fitted level by level.
#'
#' The price of the guarantee is that it is bought with a distributional
#' assumption. Under misspecification the single fitted distribution can be
#' monotone and wrong, whereas independently fitted levels remain consistent for
#' each quantile separately. Monotonicity is a property worth having, not a
#' certificate of correctness; check it alongside [plot.gkwqreg()] panel 5 and
#' [pinball()], not instead of them.
#'
#' @return An object of class `"gkwq_crossing"`, a list with components
#' \describe{
#'   \item{`taus`}{Numeric vector of the levels checked, in increasing order.}
#'   \item{`Q`}{The \eqn{n \times m}{n x m} matrix of fitted quantiles, rows in
#'     the order of the evaluation data and columns named by level.}
#'   \item{`mode`}{`"separate"` for independently fitted levels, `"implied"` for
#'     the levels of one fitted distribution.}
#'   \item{`n_crossing`}{Number of rows containing at least one crossing.}
#'   \item{`frac`}{`n_crossing` divided by the number of rows.}
#'   \item{`which`}{Integer row indices of the offending rows, so that the
#'     covariate values responsible can be inspected.}
#'   \item{`worst`}{Largest decrease observed, on the scale of the response;
#'     `0` when there is none. Useful for deciding whether a crossing is a
#'     substantive defect or numerical dust.}
#'   \item{`pairs`}{Data frame with one row per adjacent pair of levels, giving
#'     `tau_lo`, `tau_hi` and the count `n` of rows violated by that pair.}
#' }
#' A `print` method summarises the check and states which of the two questions
#' was answered.
#'
#' @references
#' Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile and
#' probability curves without crossing. *Econometrica* **78**, 1093-1125.
#'
#' @seealso [rearrange()] for the post-hoc monotonicity fix; [predict.gkwqreg()]
#'   with a `tau` argument for reading other levels off one fitted distribution.
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' x <- runif(n, -1, 1)
#' mu <- plogis(0.2 + 1.2 * x)
#' y <- rbeta(n, 4 * mu, 4 * (1 - mu))
#' d <- data.frame(y = y, x = x)
#'
#' ## (a) ONE fit, nineteen levels read off its fitted distribution.
#' fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
#' check_crossing(fit, taus = seq(0.05, 0.95, by = 0.05))
#' #>   rows with a crossing: 0 of 120 (0.00%)
#' ## Zero, and it could not have been anything else.
#'
#' ## (b) NINETEEN independent fits, one per level, shape modelled too.
#' fits <- gkwqreg(y ~ x | x, data = d, tau = seq(0.05, 0.95, by = 0.05),
#'                 family = "kw")
#' check_crossing(fits)
#' #>   rows with a crossing: 0 of 120 (0.00%)
#' ## None in sample -- but that is an empirical fact about these 120 rows,
#' ## not a guarantee. Evaluate the same fits where the data do not constrain
#' ## them and the guarantee's absence becomes visible:
#'
#' grid <- data.frame(x = seq(-3, 3, length.out = 201))
#' cr <- check_crossing(fits, newdata = grid)
#' cr
#' #>   rows with a crossing: 44 of 201 (21.89%)
#' #>   worst violation     : 6.146e-03
#'
#' ## Where, in tau: the upper levels, which are the least well determined.
#' cr$pairs[cr$pairs$n > 0, ]
#' #>      tau_lo tau_hi  n
#' #> 0.65   0.60   0.65  2
#' #> 0.70   0.65   0.70 11
#' #> ...
#' #> 0.95   0.90   0.95 44
#'
#' ## Where, in x: entirely outside the observed range of the covariate.
#' range(grid$x[cr$which])   # -3.00 -1.71
#' range(d$x)                # -0.97  0.99
#'
#' ## The same nineteen levels taken from the single fit stay monotone there.
#' check_crossing(fit, newdata = grid, taus = seq(0.05, 0.95, by = 0.05))$n_crossing
#' #> 0
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

#' Monotone rearrangement of crossing quantile curves
#'
#' Restores monotonicity in the quantile level by sorting each row of the fitted
#' quantile matrix. This is the rearrangement operator of Chernozhukov,
#' Fernandez-Val and Galichon (2010), which repairs a non-monotone estimate of a
#' monotone curve without refitting anything and without ever making the estimate
#' worse.
#'
#' @param object A `"gkwqregs"` container, or a plain list of `"gkwqreg"` fits at
#'   different levels. A single `"gkwqreg"` fit is accepted but has nothing to
#'   repair: a message says so and its quantiles are returned unchanged.
#' @param newdata Optional data frame on which to evaluate the quantiles before
#'   rearranging. Passed straight to [check_crossing()], so the same rules apply.
#' @param ... Passed to [check_crossing()]; in particular `taus` and `tol`.
#'
#' @details
#' Fix a row \eqn{i} and let
#' \eqn{\hat{Q}_{i1}, \ldots, \hat{Q}_{im}}{Qhat_i1, ..., Qhat_im} be its fitted
#' quantiles at the increasing levels
#' \eqn{\tau_1 < \cdots < \tau_m}{tau_1 < ... < tau_m}. The rearranged curve is
#' simply the sorted vector,
#'
#' \deqn{\hat{Q}^{*}_{ij} \;=\; \hat{Q}_{i(j)},}{Qhat*_ij = the jth smallest of
#'       Qhat_i1, ..., Qhat_im,}
#'
#' the \eqn{j}th order statistic of the row. Equivalently, and this is the way to
#' see why it is the right operation, \eqn{\hat{Q}^{*}}{Qhat*} is the quantile
#' function of the random variable \eqn{\hat{Q}(U)}{Qhat(U)} obtained by feeding a
#' uniform level \eqn{U} through the unsorted curve: rearrangement replaces a
#' non-monotone curve by the monotone curve with the same distribution of values.
#'
#' @section Why sorting is a legitimate repair, not a cosmetic one:
#' Sorting looks like tidying up an embarrassment, and it is natural to suspect it
#' of hiding a problem rather than solving one. It does not, and the reason is a
#' theorem rather than a convention.
#'
#' The estimand \eqn{Q_0(\cdot \mid x)}{Q0(. | x)} is a genuinely non-decreasing
#' function of \eqn{\tau}{tau}. Chernozhukov, Fernandez-Val and Galichon (2010)
#' show that the rearrangement operator is a contraction towards the set of
#' non-decreasing functions: for every \eqn{p \ge 1}{p >= 1},
#'
#' \deqn{\left\| \hat{Q}^{*} - Q_0 \right\|_p \;\le\; \left\| \hat{Q} - Q_0 \right\|_p,}{
#'       || Qhat* - Q0 ||_p  <=  || Qhat - Q0 ||_p,}
#'
#' and, under their mild regularity conditions, strictly so whenever
#' \eqn{\hat{Q}}{Qhat} actually crosses. Sorting therefore *weakly dominates* the
#' unsorted estimate in every \eqn{L_p}{Lp} distance simultaneously, whatever the
#' true curve happens to be. The intuition is direct: a decreasing stretch of the
#' estimate cannot be tracking a non-decreasing target, so exchanging the two
#' values must move both closer to it. Nothing here is asymptotic, and nothing
#' depends on the estimator that produced \eqn{\hat{Q}}{Qhat}.
#'
#' Two limits are worth stating plainly. First, the guarantee is about the
#' *quantile curve*, not about the coefficients: after rearranging you hold a
#' matrix of fitted values, not a model, and the coefficients of the fits that
#' produced it are unchanged and still non-monotone. Second, monotone is not the
#' same as correct. Rearrangement removes an internal contradiction; it does not
#' repair a misspecified family, and it cannot be read as evidence that the fit is
#' adequate. If a great many rows need repair, that is a diagnostic about the
#' model, and a single parametric fit whose implied quantiles cannot cross at all
#' (see [check_crossing()]) is usually the better answer than a large correction.
#'
#' @return A numeric matrix of rearranged quantiles with the same dimensions and
#'   dimnames as the input quantile matrix: one row per evaluation observation,
#'   one column per level, columns named by level. The matrix carries an attribute
#'   `"crossing"` holding the `"gkwq_crossing"` object computed *before*
#'   rearrangement, so that what was repaired remains recoverable through
#'   `attr(x, "crossing")`.
#'
#' @references
#' Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile and
#' probability curves without crossing. *Econometrica* **78**, 1093-1125.
#'
#' @seealso [check_crossing()] for the diagnostic that motivates this repair.
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' x <- runif(n, -1, 1)
#' mu <- plogis(0.2 + 1.2 * x)
#' y <- rbeta(n, 4 * mu, 4 * (1 - mu))
#' d <- data.frame(y = y, x = x)
#'
#' fits <- gkwqreg(y ~ x | x, data = d, tau = seq(0.05, 0.95, by = 0.05),
#'                 family = "kw")
#'
#' ## Evaluate on a grid that extrapolates well beyond the observed x, which is
#' ## where independently fitted levels are least constrained.
#' grid <- data.frame(x = seq(-3, 3, length.out = 201))
#' R <- rearrange(fits, newdata = grid)
#'
#' ## What was repaired is kept with the result.
#' cr <- attr(R, "crossing")
#' c(rows_repaired = cr$n_crossing, of = nrow(R), worst = cr$worst)
#' #> rows_repaired            of         worst
#' #>  4.400000e+01  2.010000e+02  6.146179e-03
#'
#' ## Monotone in tau afterwards, row by row, by construction.
#' all(apply(R, 1, function(r) all(diff(r) >= 0)))
#' #> TRUE
#'
#' ## The repair is local: only the rows that crossed were touched.
#' Q <- cr$Q
#' identical(which(rowSums(abs(R - Q)) > 0), cr$which)
#' #> TRUE
#'
#' ## And it is small -- the largest fitted value moved by sorting.
#' max(abs(R - Q))
#' #> 0.0100315
#'
#' ## A single fit has nothing to rearrange and says so.
#' fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
#' Q1 <- rearrange(fit, taus = c(0.1, 0.5, 0.9))
#' #> a single fit has no crossing to rearrange: its quantiles are those of one
#' #> distribution function. Returning them unchanged.
#' head(round(Q1, 4), 3)
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

#' The estimated quantile process
#'
#' Fits the model independently at each level of a grid of quantile levels and
#' collects the resulting coefficient paths
#' \eqn{\hat{\beta}_j(\tau)}{betahat_j(tau)} together with pointwise confidence
#' bands. The companion [plot method][plot.gkwq_process] draws them in the visual
#' idiom of `quantreg`'s `plot.summary.rqs`, because that picture is the
#' literature's standard summary of a quantile process and readers should
#' recognise it on sight.
#'
#' @param object Either a `"gkwqreg"` fit, whose stored call is re-evaluated with
#'   `tau` replaced by `taus`, or a `"gkwqregs"` container of fits that have
#'   already been computed and are simply harvested.
#' @param taus Increasing grid of quantile levels in \eqn{(0,1)}{(0,1)}, used only
#'   when `object` is a single fit. **Silently ignored when `object` is a
#'   `"gkwqregs"` container**, whose levels are already fixed by its fits.
#' @param level Confidence level for the pointwise bands, `0.95` by default.
#' @param ... Additional arguments merged into the re-evaluated call, hence passed
#'   to [gkwqreg()]; `family`, `anchor` or `control` may be overridden here.
#'   Ignored when `object` is a container.
#'
#' @details
#' Each level is a separate maximum-likelihood problem with its own likelihood and
#' its own parameter vector, exactly as [gkwqreg()] produces when given a vector
#' `tau`. The process object collects those separate answers side by side; it does
#' not estimate them jointly, and it imposes no smoothness or monotonicity across
#' \eqn{\tau}{tau}.
#'
#' The bands are pointwise Wald intervals computed level by level,
#' \eqn{\hat{\beta}_j(\tau) \pm z_{1 - (1 - \texttt{level})/2}\,
#' \mathrm{se}\{\hat{\beta}_j(\tau)\}}{betahat_j(tau) +/- z * se(betahat_j(tau))},
#' from each fit's own observed information. Three consequences follow and none of
#' them should be glossed over when the picture is presented. The bands are
#' pointwise in \eqn{\tau}{tau} and not simultaneous, so a band that excludes zero
#' at one level out of the nineteen in the default grid is not evidence at the
#' nominal level for the process as a whole. They ignore the dependence between levels, which is substantial
#' because every fit uses the same rows. And they are symmetric on the coefficient
#' scale, which is the scale of the linear predictor, not of the quantile; use
#' [marginal_effects()] for statements about the response.
#'
#' When `object` is a single fit the call is rebuilt and evaluated in the caller's
#' frame, so the data referred to by that call must still be visible from where
#' `quantile_process()` is invoked. Passing a `"gkwqregs"` container avoids the
#' issue entirely and refits nothing.
#'
#' @section A non-flat coefficient path is not evidence of tail heterogeneity:
#' This is the single most consequential misreading of the picture this function
#' produces, and under the present parametrization it is a trap rather than a
#' subtlety.
#'
#' In check-function quantile regression the coefficient path is flat exactly when
#' the covariate shifts the whole conditional distribution without changing its
#' shape, so a sloping path is read, correctly, as tail heterogeneity. **That
#' reading does not carry over here.** The quantity plotted is the slope of
#' \eqn{g(\mu_\tau)}{g(mu_tau)}, the link-transformed conditional
#' \eqn{\tau}{tau}-quantile, and for a fixed conditional distribution the map from
#' a covariate to \eqn{g}{g} of its \eqn{\tau}{tau}-quantile is a different
#' function at each \eqn{\tau}{tau}. The path is therefore *generically*
#' \eqn{\tau}{tau}-varying even when the data-generating process has no tail
#' heterogeneity whatsoever. It is a property of the link and the reparametrization,
#' not a finding about the data.
#'
#' The effect is large, not a second-order curiosity. Simulating from a correctly
#' specified `kw` model with a **constant** \eqn{\alpha}{alpha}, one covariate, a
#' single index driving the entire conditional distribution and a true median slope
#' of `1.1` (\eqn{n = 20000}{n = 20000}), the fitted `mu:x` coefficient runs
#'
#' \tabular{lrrrrr}{
#'   `tau`  \tab 0.10  \tab 0.25  \tab 0.50  \tab 0.75  \tab 0.90 \cr
#'   `mu:x` \tab 0.748 \tab 0.883 \tab 1.100 \tab 1.397 \tab 1.744
#' }
#'
#' The median is recovered exactly. The other four levels are not the same number
#' and were never supposed to be: a path rising by a factor of 2.3 from the first
#' decile to the ninth is what *no* tail heterogeneity looks like on this scale.
#' A reader shown only the picture would report that the covariate "matters more in
#' the upper tail", and would be describing the parametrization.
#'
#' What to do instead. The null hypothesis to compare against is not a horizontal
#' line but the path implied by a homogeneous model, and there are two honest ways
#' to obtain it. Simulate from the fitted homogeneous model with
#' [simulate.gkwqreg()], refit the process on the simulated data, and use that path
#' as the reference; the `reference` argument of [plot.gkwq_process()] takes a
#' value per coefficient for exactly this purpose. Or, better, test the question
#' directly inside the model: give the shape parameter its own regression
#' (`y ~ x | x` rather than `y ~ x`) and compare the two nested fits at a single
#' level with [lrtest()]. Tail heterogeneity in this parametrization means the
#' shape parameters depend on the covariates, and that is a hypothesis with a
#' likelihood-ratio test attached, not something to be read off a slope.
#'
#' The shape coefficients' own paths are the useful diagnostic in the plot: under a
#' homogeneous data-generating process they are flat, because a constant shape
#' parameter is estimated as the same constant at every level.
#'
#' @return An object of class `"gkwq_process"`, a list with components
#' \describe{
#'   \item{`taus`}{Numeric vector of the levels, increasing.}
#'   \item{`coef`}{Coefficient matrix, one row per model coefficient and one
#'     column per level. Rows are named `"part:term"` (for example `"mu:x"`,
#'     `"alpha:(Intercept)"`); columns are named by level.}
#'   \item{`se`}{Standard errors, same shape and dimnames as `coef`.}
#'   \item{`lower`, `upper`}{Pointwise Wald band limits, `coef` plus or minus
#'     `qnorm(1 - (1 - level) / 2)` standard errors.}
#'   \item{`parts`}{Character vector aligned with the rows of `coef`, naming the
#'     formula part each coefficient belongs to; this is what the plot method's
#'     `parts` argument selects on.}
#'   \item{`level`}{The confidence level used for the bands.}
#'   \item{`family`, `anchor`}{Family and anchor, common to every fit.}
#'   \item{`fits`}{The list of underlying `"gkwqreg"` fits, so that any per-level
#'     quantity (log-likelihood, check loss, residuals) remains available.}
#'   \item{`call`}{The matched call.}
#' }
#' A `print` method reports the family, anchor, level range and the rounded
#' coefficient matrix.
#'
#' @seealso [plot.gkwq_process()] for the picture; [check_crossing()] and
#'   [rearrange()] for the monotonicity of independently fitted levels;
#'   [marginal_effects()] for effects on the scale of the response; [lrtest()] for
#'   testing shape heterogeneity properly.
#'
#' @examples
#' ## A correctly specified kw model with a CONSTANT shape parameter: by
#' ## construction there is no tail heterogeneity of any kind here.
#' set.seed(6)
#' n <- 1000
#' x <- runif(n, -2, 2)
#' mu <- plogis(0.4 + 1.1 * x)                    # true median, logit slope 1.1
#' y <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
#' d <- data.frame(y = y, x = x)
#'
#' fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
#' qp <- quantile_process(fit, taus = c(0.1, 0.25, 0.5, 0.75, 0.9))
#' round(qp$coef, 3)
#' #>                     0.10   0.25  0.50  0.75  0.90
#' #> mu:(Intercept)    -1.180 -0.409 0.436 1.342 2.272
#' #> mu:x               0.750  0.891 1.116 1.425 1.789
#' #> alpha:(Intercept)  0.664  0.680 0.688 0.685 0.675
#'
#' ## Read the two informative rows against each other.
#' ##   mu:x rises from 0.75 to 1.79 and recovers 1.1 at the median. That rise
#' ##   is NOT tail heterogeneity: the simulation has none. It is what a single
#' ##   fixed conditional distribution looks like when read at five levels
#' ##   through a logit link.
#' ##   alpha:(Intercept) is flat to within 0.025, correctly reporting that the
#' ##   shape parameter is the same constant at every level. THAT is the row
#' ##   that carries the heterogeneity question.
#'
#' ## Test the question directly instead of eyeballing the slope of mu:x.
#' m0 <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")       # alpha constant
#' m1 <- gkwqreg(y ~ x | x, data = d, tau = 0.5, family = "kw")   # alpha ~ x
#' lrtest(m0, m1)
#' #>      Df logLik     AIC  Chisq Chi Df Pr(>Chisq)
#' #> [1,]  3 521.43 -1036.9
#' #> [2,]  4 521.43 -1034.9 0.0018      1     0.9664
#' ## No evidence at all that the shape depends on x -- which is the right
#' ## answer, and the one the sloping mu:x path would have contradicted.
#'
#' ## Harvesting an existing container refits nothing and ignores `taus`.
#' fits <- gkwqreg(y ~ x, data = d, tau = c(0.25, 0.5, 0.75), family = "kw")
#' quantile_process(fits, taus = seq(0.1, 0.9, by = 0.1))$taus
#' #> 0.25 0.50 0.75
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

#' Plot a quantile process
#'
#' Draws one panel per selected coefficient: the estimated path
#' \eqn{\hat{\beta}_j(\tau)}{betahat_j(tau)} against the quantile level, with a
#' shaded pointwise confidence band, a dotted line at zero and, optionally, a
#' dashed reference line. The layout deliberately reproduces the idiom of
#' `quantreg`'s `plot.summary.rqs` so that the picture is immediately legible to
#' readers of the quantile-regression literature.
#'
#' @param x An object of class `"gkwq_process"` from [quantile_process()]; or, for
#'   the `"gkwqregs"` method, a container of fits, which is passed through
#'   [quantile_process()] first.
#' @param parm Character vector of coefficient names to draw, matched against
#'   `rownames(x$coef)` and therefore written in the `"part:term"` form, for
#'   example `c("mu:x", "mu:(Intercept)")`. Names that match nothing are dropped
#'   silently. When `NULL` (the default) selection falls to `parts`.
#' @param parts Character vector of formula parts to draw when `parm` is `NULL`;
#'   every coefficient belonging to one of these parts is included. Defaults to
#'   `"mu"`, the conditional quantile. Use `parts = unique(x$parts)` for every
#'   coefficient in the model, or `parts = "alpha"` to inspect a shape path.
#' @param nrow,ncol Panel grid. Supply neither and the grid is at most three
#'   columns wide with as many rows as needed; supply one and the other is
#'   derived. The previous [graphics::par()] settings are restored on exit.
#' @param reference Optional numeric vector of horizontal reference values, drawn
#'   as a dashed line and included in the vertical range of the panel. It is
#'   indexed **positionally against all rows of `x$coef`**, not against the
#'   selected subset, so it must have length `nrow(x$coef)` and follow the order
#'   of `rownames(x$coef)`; use `NA` for coefficients that need no line. A shorter
#'   vector is not recycled and simply produces no line.
#' @param ... Currently unused; accepted for compatibility with the `plot`
#'   generic.
#'
#' @details
#' The band drawn is `x$lower` to `x$upper`, the pointwise Wald interval at
#' `x$level` computed separately at each quantile level. It is not simultaneous
#' over \eqn{\tau}{tau} and it ignores the dependence between levels induced by
#' their sharing the same data; see [quantile_process()] for the full statement.
#'
#' **Read the shape of a path with the warning in [quantile_process()] in hand.**
#' Under this parametrization the plotted coefficient is the slope of the
#' link-transformed conditional \eqn{\tau}{tau}-quantile, and that slope is
#' generically \eqn{\tau}{tau}-varying even when the data-generating process has no
#' tail heterogeneity at all. A sloping path is therefore not, by itself, evidence
#' that a covariate matters more in the tail. The `reference` argument exists so
#' that a path can be judged against a value that means something -- a known
#' truth, a single-level estimate, or the path implied by a simulated homogeneous
#' model -- rather than against a horizontal line that carries no null hypothesis.
#'
#' @return `x`, invisibly. Called for the plot it produces.
#'
#' @seealso [quantile_process()], which constructs the object and documents its
#'   components and its interpretation.
#'
#' @examples
#' ## Correctly specified kw data with a CONSTANT shape parameter: the
#' ## data-generating process has no tail heterogeneity at all.
#' set.seed(6)
#' n <- 1000
#' x <- runif(n, -2, 2)
#' mu <- plogis(0.4 + 1.1 * x)                # true median, logit slope 1.1
#' y <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
#' d <- data.frame(y = y, x = x)
#'
#' fits <- gkwqreg(y ~ x, data = d, tau = seq(0.1, 0.9, by = 0.1), family = "kw")
#' qp <- quantile_process(fits)
#'
#' ## Default: every coefficient of the conditional-quantile part.
#' plot(qp)
#'
#' ## One coefficient, with the true median slope drawn in. The path passes
#' ## through the reference near tau = 0.5 and rises steadily away from it in
#' ## both directions -- from 0.750 at tau = 0.1 to 1.789 at tau = 0.9. That is
#' ## correct behaviour and NOT tail heterogeneity: the simulation has none.
#' ref <- rep(NA_real_, nrow(qp$coef))
#' names(ref) <- rownames(qp$coef)
#' ref["mu:x"] <- 1.1
#' plot(qp, parm = "mu:x", reference = ref)
#'
#' ## The shape path is the row that does carry the heterogeneity question, and
#' ## here it is flat to within 0.025, as it should be for a constant alpha.
#' plot(qp, parts = "alpha")
#' round(qp$coef["alpha:(Intercept)", ], 3)
#' #>   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9
#' #> 0.664 0.677 0.683 0.686 0.688 0.688 0.687 0.683 0.675
#'
#' ## Every coefficient of the model, in a single row of panels.
#' plot(qp, parts = unique(qp$parts), nrow = 1)
#'
#' ## A container plots itself by calling quantile_process() first.
#' plot(fits, parm = "mu:x")
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

#' @rdname plot.gkwq_process
#' @export
plot.gkwqregs <- function(x, ...) plot(quantile_process(x), ...)
