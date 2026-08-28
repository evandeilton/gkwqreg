## ---------------------------------------------------------------------------
## Prediction.
## ---------------------------------------------------------------------------

## Rebuild design matrices and offsets for newdata, honouring the factor levels
## and contrasts recorded at fit time.
.gkwq_newdata_matrices <- function(object, newdata, na.action = stats::na.pass) {
  X <- vector("list", length(object$parts))
  O <- vector("list", length(object$parts))
  names(X) <- names(O) <- object$parts
  for (p in object$parts) {
    mt <- stats::delete.response(object$terms[[p]])
    ## Only hand model.frame the levels for variables THIS part actually uses.
    ## Passing the whole xlev list makes it warn "variable 'x' is not a factor"
    ## for every factor that appears in some other part -- noise on essentially
    ## every predict(newdata=) call with a factor in the model.
    xl <- object$levels[intersect(names(object$levels), all.vars(mt))]
    mf <- stats::model.frame(mt, newdata, na.action = na.action, xlev = xl)
    X[[p]] <- stats::model.matrix(mt, mf, contrasts.arg = object$contrasts[[p]])
    o <- stats::model.offset(mf)
    O[[p]] <- if (is.null(o)) rep(0, nrow(X[[p]])) else as.numeric(o)
  }
  list(X = X, offsets = O)
}

.gkwq_eta_mu_theta <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    nd <- .gkwq_newdata_matrices(object, newdata)
    eta <- stats::setNames(lapply(object$parts, function(p) {
      as.numeric(nd$X[[p]] %*% object$coef_list[[p]]) + nd$offsets[[p]]
    }), object$parts)
  }
  mu <- .gkwq_linkinv(eta[["mu"]], object$link[["mu"]], object$link_scale[["mu"]])
  mu <- pmin(pmax(mu, object$control$eps_mu), 1 - object$control$eps_mu)
  theta <- stats::setNames(lapply(object$parts[-1L], function(p) {
    .gkwq_linkinv(eta[[p]], object$link[[p]], object$link_scale[[p]])
  }), object$parts[-1L])
  list(eta = eta, mu = mu, theta = theta)
}

#' Predictions from a fixed-level quantile regression fit
#'
#' Evaluates a fitted [gkwqreg()] model at the estimation data or at `newdata`.
#' The first part of the model *is* the conditional quantile, so the default
#' prediction is the conditional \eqn{\tau}-quantile itself and never a
#' conditional mean. Densities, distribution functions and moments are all
#' derived from the same fitted conditional distribution.
#'
#' @details
#' Write \eqn{g} for the link of the quantile part (`"logit"` by default) and
#' \eqn{o_i} for its offset. The default prediction is
#' \deqn{\widehat{Q}(\tau \mid x_i) \;=\; \hat\mu_\tau(x_i) \;=\;
#'   g^{-1}\!\left(x_i^{\top}\hat\beta_\mu + o_i\right).}
#'
#' Every other `type` is obtained by first reconstructing, row by row, the full
#' parameter vector \eqn{(\alpha, \beta, \gamma, \delta, \lambda)} of the
#' Generalized Kumaraswamy distribution: the non-anchored parameters come from
#' their own linear predictors and links, and the anchored parameter is solved
#' for so that the resulting distribution has \eqn{\hat\mu_\tau(x_i)} as its
#' exact \eqn{\tau}-quantile. A useful consequence is that
#' \deqn{F\!\left(\widehat{Q}(\tau \mid x_i) \mid x_i\right) = \tau}
#' holds to machine precision for every row, which is the cheapest available
#' check that a prediction pipeline has not silently drifted onto some other
#' scale (see the examples).
#'
#' @section Values of `type`:
#'
#' | `type` | What is returned | Typical use |
#' | :----- | :--------------- | :---------- |
#' | `"quantile"` (default), `"response"` | the conditional tau-quantile; with `tau` supplied, several quantiles of the same fitted distribution | the headline prediction |
#' | `"mu"` | the conditional tau-quantile, always ignoring `tau` | when the fitted level is wanted inside code that passes `tau` around |
#' | `"link"` | the linear predictors of every part, each on its own link scale | diagnostics, and plotting on the scale the model is linear in |
#' | `"parameter"` | the reconstructed five distribution parameters per row | handing the fitted conditional law to `gkwdist` functions |
#' | `"mean"` | the conditional mean, by 64-node Gauss-Legendre quadrature of the quantile function | comparison with a mean-parametrized fit |
#' | `"variance"` | the conditional variance, by the same quadrature | assessing conditional dispersion |
#' | `"density"` | the conditional density evaluated at `at` | likelihood displays, simulated envelopes |
#' | `"probability"` | the conditional distribution function evaluated at `at` | exceedance probabilities, probability-integral-transform checks |
#' | `"terms"` | each term's additive contribution to its part's linear predictor | partial-effect plots |
#'
#' `"mean"` and `"variance"` integrate \eqn{\int_0^1 Q(u)^k \, du} rather than
#' the density, because \eqn{Q} is available in closed form and has no boundary
#' singularity to work around. They are provided so that a quantile fit can be
#' put beside a mean-parametrized one; reporting the conditional mean as though
#' it were the estimand of this model defeats the point of the parametrization.
#'
#' @section Reading other levels off the fit versus refitting:
#' Supplying `tau` returns quantiles of the **same fitted conditional
#' distribution** read at other probabilities. This is not the same object as
#' `gkwqreg(..., tau = tau_new)`, which solves a different estimating problem
#' and in general returns different coefficients.
#'
#' * Reading off is a *distributional extrapolation*. It borrows strength from
#'   the parametric family, costs nothing beyond a call to the quantile
#'   function, and is monotone in `tau` by construction, so predictions read off
#'   one fit can never cross.
#' * Refitting is a *level-specific* estimate. It targets the new level
#'   directly, and is the honest thing to report when the family may be wrong
#'   away from the level already fitted.
#'
#' The two coincide only if the family is correctly specified. The gap between
#' them is therefore a specification diagnostic in its own right: a small
#' discrepancy far from the fitted `tau` supports the shape assumption, and a
#' large one says the parametric tail, rather than the systematic part of the
#' model, is doing the work. See [quantile_process()] for the systematic version
#' of this comparison and [check_crossing()] for what can go wrong once several
#' levels are fitted separately.
#'
#' @param object A `"gkwqreg"` fit, as returned by [gkwqreg()].
#' @param newdata Optional data frame in which to evaluate the model. It must
#'   contain every variable used by *every* part of the formula, not only the
#'   quantile part. Factor levels and contrasts recorded at fit time are
#'   reapplied, so a factor whose levels are incompletely represented in
#'   `newdata` still produces columns consistent with the fit. Defaults to the
#'   estimation data.
#' @param type The quantity to predict; see the table below. `"response"` is an
#'   alias for `"quantile"`, retained so that code written against other
#'   regression classes keeps running -- but note that here the "response
#'   scale" prediction is a quantile, not a mean.
#' @param at Numeric vector of points in `(0,1)` at which the density or the
#'   distribution function is evaluated. Required for `type = "density"` and
#'   `type = "probability"`, and ignored otherwise.
#' @param tau Optional numeric vector of quantile levels in `(0,1)`, used only
#'   by `type = "quantile"` and `type = "response"`. **These are further
#'   quantiles of the same fitted conditional distribution, not a refit**; see
#'   the section below. `NULL` (the default) returns the level the model was
#'   fitted at.
#' @param elementwise Logical, governing how `at` is combined with the rows.
#'   `TRUE` pairs the `i`th element of `at` with the `i`th row, recycling `at`
#'   if necessary, and returns a vector; `FALSE` crosses every element of `at`
#'   with every row and returns a matrix. The default is `TRUE` when `at` has
#'   exactly one element per row and there is more than one row, and `FALSE`
#'   otherwise.
#' @param na.action How to treat missing values in `newdata`. The default,
#'   [stats::na.pass()], retains the rows and propagates `NA` into the
#'   predictions, so the result stays aligned with the rows of `newdata`.
#' @param ... Currently unused; present for consistency with the generic.
#'
#' @return
#' The shape depends on `type`, and for `"quantile"` also on `tau`.
#'
#' * `"quantile"`, `"response"`, `"mu"`: a numeric vector with one element per
#'   row of `newdata`, or per observation used in the fit. If `tau` is supplied
#'   to `"quantile"`/`"response"`, an `n` by `length(tau)` matrix whose columns
#'   are named `"tau=<level>"`.
#' * `"link"`: a data frame with `n` rows and one column per model part, named
#'   after the parts (`mu` first), holding the linear predictors.
#' * `"parameter"`: a data frame with `n` rows and the five columns `alpha`,
#'   `beta`, `gamma`, `delta`, `lambda`, in `gkwdist` order.
#' * `"mean"`, `"variance"`: numeric vectors of length `n`.
#' * `"density"`, `"probability"`: an `n` by `length(at)` matrix whose columns
#'   are named after `at`, or a numeric vector of length `n` when
#'   `elementwise = TRUE`.
#' * `"terms"`: a named list with one matrix per model part. Each matrix has `n`
#'   rows and one column per non-intercept term of that part; a part with no
#'   covariates contributes a matrix with zero columns.
#'
#' Applying `predict()` to the `"gkwqregs"` container returned by a vector-valued
#' `tau` predicts from each fit in turn, and binds the results column-wise, named
#' by level, when every level returned a plain vector of the same length.
#'
#' @seealso [gkwqreg()] for the model, [fitted.gkwqreg()] for the in-sample
#'   quantiles, [marginal_effects()] for effects on the quantile scale,
#'   [quantile_process()] and [check_crossing()] for several levels at once.
#'
#' @examples
#' ## A Kumaraswamy sample whose conditional MEDIAN follows a logit model.
#' ## Drawn by inverse transform through gkwq_quantile(), so the data-generating
#' ## process is exactly the one the model assumes.
#' set.seed(2024)
#' n   <- 300
#' x1  <- runif(n, -2, 2)
#' x2  <- rbinom(n, 1, 0.5)
#' mu  <- plogis(0.3 + 0.9 * x1 - 0.5 * x2)      # true conditional median
#' bt  <- log1p(-0.5) / log1p(-mu^2)             # beta anchoring Q(0.5) = mu
#' y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
#' dat <- data.frame(y = y, x1 = x1, x2 = factor(x2))
#'
#' fit <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw")
#'
#' ## -- the default prediction is a conditional quantile --------------------
#' nd <- data.frame(x1 = c(-1, 0, 1), x2 = factor(0, levels = c("0", "1")))
#' predict(fit, nd)
#' ## 0.3651 0.5755 0.7618.  Read this as: among units with x1 = 0 and x2 = 0,
#' ## half are predicted to fall below 0.5755 -- NOT "they average 0.5755".
#'
#' ## -- the arithmetic check ------------------------------------------------
#' predict(fit, nd, type = "probability", at = predict(fit, nd),
#'         elementwise = TRUE)
#' ## 0.5 0.5 0.5 exactly: the fitted distribution puts mass tau below its own
#' ## fitted tau-quantile, by construction of the anchor.
#'
#' ## -- other levels of the SAME fit ---------------------------------------
#' predict(fit, nd, tau = c(0.1, 0.5, 0.9))
#' ## The "tau=0.5" column reproduces the default call above; the other two are
#' ## read off the same conditional law and were never targeted by estimation.
#'
#' ## How much does that extrapolation cost?  Compare with a genuine refit.
#' fit90 <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.9, family = "kw")
#' cbind(read_off = predict(fit, nd, tau = 0.9)[, 1],
#'       refit    = predict(fit90, nd))
#' ##      read_off  refit
#' ## [1,]   0.5945 0.6252
#' ## [2,]   0.8438 0.8533
#' ## [3,]   0.9668 0.9530
#' ## Close, as it should be here: the family is correctly specified. A wide
#' ## gap would indict the shape assumption, not the covariate effects.
#'
#' ## -- the fitted conditional distribution, row by row ---------------------
#' predict(fit, nd, type = "parameter")   # alpha, beta, gamma, delta, lambda
#' predict(fit, nd, type = "link")        # eta for every part
#' predict(fit, nd, type = "mean")        # 0.3720 0.5630 0.7123
#' predict(fit, nd, type = "variance")
#' ## The conditional mean (0.5630) is NOT the conditional median (0.5755):
#' ## the fitted Kumaraswamy law is skewed. This is why fitted() must not be
#' ## read as a mean.
#'
#' ## -- densities and probabilities ----------------------------------------
#' predict(fit, nd, type = "density", at = c(0.25, 0.50, 0.75))
#' ## A VECTOR of length 3, not a matrix: `at` happens to have exactly one
#' ## element per row, so `elementwise` defaults to TRUE and pairs the ith
#' ## point with the ith row. Say so explicitly when that is not the intent:
#' predict(fit, nd, type = "density", at = c(0.25, 0.50, 0.75),
#'         elementwise = FALSE)                                     # 3 x 3
#' predict(fit, nd, type = "probability", at = 0.5)                 # 3 x 1
#'
#' ## -- additive term contributions ----------------------------------------
#' predict(fit, nd, type = "terms")$mu
#' @export
predict.gkwqreg <- function(object, newdata = NULL,
                            type = c("quantile", "response", "link", "parameter",
                                     "mu", "density", "probability", "terms",
                                     "mean", "variance"),
                            at = NULL, tau = NULL, elementwise = NULL,
                            na.action = stats::na.pass, ...) {
  type <- match.arg(type)

  if (type == "terms") {
    return(.gkwq_predict_terms(object, newdata))
  }

  em <- .gkwq_eta_mu_theta(object, newdata)
  if (type == "link") {
    return(as.data.frame(em$eta))
  }

  P <- .gkwq_reconstruct(em$mu, object$tau, object$spec, em$theta)
  pv <- as.data.frame(P)

  switch(type,
    quantile = ,
    response = {
      if (is.null(tau)) {
        em$mu
      } else {
        ## The same fitted distribution, read at other probabilities.
        out <- vapply(tau, function(t) {
          gkwq_quantile(t, pv$alpha, pv$beta, pv$gamma, pv$delta, pv$lambda)
        }, numeric(nrow(pv)))
        if (is.null(dim(out))) out <- matrix(out, nrow = nrow(pv))
        colnames(out) <- paste0("tau=", format(tau, trim = TRUE))
        out
      }
    },
    mu = em$mu,
    parameter = pv,
    mean = .gkwq_moment(pv, 1L),
    variance = .gkwq_variance(pv),
    density = .gkwq_at(object, pv, at, elementwise, "density"),
    probability = .gkwq_at(object, pv, at, elementwise, "probability")
  )
}

.gkwq_at <- function(object, pv, at, elementwise, what) {
  if (is.null(at)) {
    stop("`at` is required for type = ", sQuote(what), ".", call. = FALSE)
  }
  n <- nrow(pv)
  ew <- elementwise %||% (length(at) == n && n > 1L)
  f <- function(a, i) {
    if (what == "density") {
      exp(.gkwq_logdens(a, pv$alpha[i], pv$beta[i], pv$gamma[i], pv$delta[i],
                        pv$lambda[i],
                        delta_is_zero = object$spec$delta_is_zero == 1L))
    } else {
      exp(.gkwq_logcdf(a, pv$alpha[i], pv$beta[i], pv$gamma[i], pv$delta[i],
                       pv$lambda[i]))
    }
  }
  if (ew) {
    at <- rep_len(at, n)
    vapply(seq_len(n), function(i) f(at[i], i), numeric(1))
  } else {
    out <- vapply(at, function(a) vapply(seq_len(n), function(i) f(a, i), numeric(1)),
                  numeric(n))
    if (is.null(dim(out))) out <- matrix(out, nrow = n)
    colnames(out) <- format(at, trim = TRUE)
    out
  }
}

.gkwq_predict_terms <- function(object, newdata) {
  if (is.null(newdata)) {
    X <- lapply(object$parts, function(p) model.matrix(object, part = p))
    names(X) <- object$parts
  } else {
    X <- .gkwq_newdata_matrices(object, newdata)$X
  }
  stats::setNames(lapply(object$parts, function(p) {
    b <- object$coef_list[[p]]
    xm <- X[[p]]
    keep <- setdiff(colnames(xm), "(Intercept)")
    if (!length(keep)) return(matrix(0, nrow(xm), 0))
    sweep(xm[, keep, drop = FALSE], 2, b[keep], `*`)
  }), object$parts)
}

#' @export
predict.gkwqregs <- function(object, newdata = NULL, ...) {
  out <- lapply(object$fits, predict, newdata = newdata, ...)
  ## Bind into a matrix only when every level returned a plain vector of the
  ## same length; types like "parameter" or "link" return data frames, and
  ## forcing those into a matrix would silently mangle them.
  vec <- vapply(out, function(z) is.numeric(z) && is.null(dim(z)), logical(1))
  if (all(vec) && length(unique(lengths(out))) == 1L) {
    out <- do.call(cbind, out)
    colnames(out) <- format(object$taus, trim = TRUE)
  }
  out
}

## ---------------------------------------------------------------------------
## Marginal effects.
## ---------------------------------------------------------------------------

#' Marginal effects on the quantile scale
#'
#' Translates the coefficients of the conditional quantile, which are effects on
#' the *log quantile odds*, into effects on the quantile itself, with delta-method
#' standard errors and confidence intervals.
#'
#' @details
#' Under the default logit link the quantile part of the model is
#' \deqn{\log \frac{\mu_\tau(x)}{1 - \mu_\tau(x)} \;=\; x^{\top}\beta_\mu ,}
#' so \eqn{\exp(\beta_j)} is the multiplicative effect of a one-unit increase in
#' \eqn{x_j} on the **quantile odds** \eqn{\mu_\tau/(1-\mu_\tau)}. It is *not*
#' the odds ratio of an event, and it is *not* an effect on a mean. Reported on
#' its own it is close to uninterpretable, because a fixed multiplicative change
#' in the odds displaces the quantile by an amount that depends entirely on
#' where that quantile already sits.
#'
#' What is interpretable on the scale of the response is the derivative
#' \deqn{\frac{\partial Q(\tau \mid x)}{\partial x_j}
#'   \;=\; \beta_j \, \frac{d\mu}{d\eta}
#'   \;=\; \beta_j \, \mu_\tau (1 - \mu_\tau) \quad \text{(logit)},}
#' which is what this function reports. For links other than the logit the
#' factor \eqn{d\mu/d\eta} changes accordingly and everything else is unaltered.
#'
#' For a factor the derivative is still the quantity reported, evaluated at the
#' dummy column: it is the *marginal* effect of that indicator on the quantile,
#' which for a binary regressor is an accurate local approximation to the
#' discrete contrast rather than the contrast itself.
#'
#' @section Where the effect is evaluated:
#'
#' | `at` | What is reported | Shape |
#' | :--- | :--------------- | :---- |
#' | `"ame"` (default) | average marginal effect: the derivative averaged over the rows | one number per covariate, with a standard error |
#' | `"mem"` | marginal effect at the mean: the derivative evaluated at the average linear predictor | one number per covariate, with a standard error |
#' | `"observed"` | the whole vector of per-row derivatives | an `n` by `p` matrix, no standard errors |
#'
#' `"ame"` is the default because it answers the question actually being asked
#' of a bounded response: how much does the quantile move, on average over the
#' population at hand? `"mem"` is cheaper to explain but describes a unit whose
#' covariates are all at their mean, which need not resemble anybody.
#' `"observed"` is the raw material for a plot of how the effect varies across
#' the covariate space; its column means reproduce the `"ame"` estimates exactly.
#'
#' @section How the standard errors are obtained:
#' For `at = "ame"` the reported effect is
#' \deqn{\widehat{\mathrm{AME}}_j(\beta) \;=\; \beta_j \cdot
#'   \frac{1}{n}\sum_{i=1}^{n} g'\!\left(x_i^{\top}\beta_\mu + o_i\right),}
#' and for `at = "mem"` the same expression with the averaging done inside,
#' \eqn{g'(\bar\eta)}. **Both** factors are functions of the coefficient
#' vector: the slope \eqn{\beta_j} directly, and the averaged derivative through
#' \eqn{\eta = X\beta_\mu}. The standard error is therefore the exact delta
#' method taken over the whole stacked coefficient vector,
#' \deqn{\widehat{\mathrm{Var}}\big(\widehat{\mathrm{AME}}\big) \;=\;
#'   J \, V \, J^{\top}, \qquad
#'   J \;=\; \partial\,\widehat{\mathrm{AME}} / \partial\beta^{\top},}
#' with \eqn{J} evaluated numerically by [numDeriv::jacobian()] and \eqn{V}
#' taken from `vcov.` when supplied and from `vcov(object)` otherwise.
#'
#' Freezing \eqn{g'(\eta)} at its sample average and differentiating only
#' \eqn{\beta_j} -- the obvious shortcut, which yields
#' \eqn{\mathrm{s.e.}(\beta_j) \times \overline{g'}} -- is not the delta method.
#' It discards the second term, which partially cancels the first: raising
#' \eqn{\beta_j} spreads \eqn{\eta}, and a more spread \eqn{\eta} has a smaller
#' \eqn{\overline{g'}} once it is away from the mode of the link derivative. The
#' shortcut consequently *overstates* the standard error, badly for a strong
#' covariate; the examples show the discrepancy on a fit where it is a factor of
#' about 1.6.
#'
#' Whether the interval deserves to be believed is a separate question from
#' whether it is computed correctly. Pass `vcov. = vcov(object, type =
#' "sandwich")` if the conditional distribution may be misspecified, or a
#' clustered covariance matrix if the data are grouped; the delta method is
#' applied to whatever matrix is supplied.
#'
#' @section Covariates that also enter a nuisance part:
#' If \eqn{x_j} also appears in the equation for one of the remaining
#' parameters, the effect on the quantile is still exactly
#' \eqn{\beta_j \, d\mu/d\eta}, with no extra term. This is a genuine
#' convenience of the parametrization rather than an approximation:
#' \eqn{\mu_\tau} *is* the \eqn{\tau}-quantile, so a nuisance part alters the
#' spread and shape of the conditional distribution while leaving the quantile
#' being modelled untouched. The print method says so when it detects the
#' overlap. Note that this is a statement about the \eqn{\tau}-quantile only:
#' the same covariate does move *other* quantiles through the nuisance part.
#'
#' @param object A `"gkwqreg"` fit, as returned by [gkwqreg()].
#' @param variables Character vector naming the covariates to report, using the
#'   coefficient names of the quantile part (`"x1"`, `"x21"` for the second
#'   level of a factor `x2`, and so on). Defaults to every non-intercept
#'   covariate in the quantile part. Names not present are dropped silently;
#'   if nothing remains the function errors.
#' @param at Where to evaluate the derivative: `"ame"` averages it over the
#'   rows, `"mem"` evaluates it at the average linear predictor, and
#'   `"observed"` returns the untouched per-row vector for each covariate. See
#'   the table below.
#' @param newdata Optional data frame over which to average, defaulting to the
#'   estimation data. Useful for reporting the effect for a subpopulation, or
#'   for a counterfactual covariate distribution, without refitting.
#' @param level Confidence level for the reported interval. The interval is
#'   Wald on the effect scale, `estimate` \eqn{\pm} `z * std.error`, and is
#'   ignored when `at = "observed"`.
#' @param vcov. Optional covariance matrix for the full stacked coefficient
#'   vector, dimensioned and named as `vcov(object)`. Typical choices are
#'   `vcov(object, type = "sandwich")` or a clustered matrix such as
#'   `sandwich::vcovCL()`. When `NULL` the model-based `vcov(object)` is used,
#'   and if that is unavailable the standard errors come back as `NA` rather
#'   than the function failing.
#' @param ... Currently unused.
#'
#' @return
#' For `at = "ame"` or `"mem"`, an object of class `"gkwq_meff"`: a list with
#' components `table` (a data frame with one row per covariate and columns
#' `variable`, `effect`, `std.error`, `lower`, `upper`), `at`, `tau`, `level`,
#' `link`, and `overlap` (the covariates that also enter a nuisance part). It
#' has a `print` method; extract `x$table` for further computation.
#'
#' For `at = "observed"`, an object of class `"gkwq_meff_observed"`: a list with
#' components `effects` (an `n` by `length(variables)` numeric matrix of
#' per-row derivatives, columns named after the covariates), `at` and `tau`.
#' There is no print method for this class; use `x$effects` directly.
#'
#' @seealso [gkwqreg()] for the model and the interpretation of its
#'   coefficients, [predict.gkwqreg()] for predictions, [vcov.gkwqreg()] for the
#'   covariance matrices that can be passed to `vcov.`.
#'
#' @examples
#' ## Same simulated design as predict.gkwqreg(): a conditional median that
#' ## follows a logit model, drawn by inverse transform.
#' set.seed(2024)
#' n   <- 300
#' x1  <- runif(n, -2, 2)
#' x2  <- rbinom(n, 1, 0.5)
#' mu  <- plogis(0.3 + 0.9 * x1 - 0.5 * x2)
#' bt  <- log1p(-0.5) / log1p(-mu^2)
#' y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
#' dat <- data.frame(y = y, x1 = x1, x2 = factor(x2))
#'
#' fit <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw")
#'
#' ## -- what the coefficient says, and what it does NOT say -----------------
#' exp(coef(fit, part = "mu"))
#' ## x1 -> 2.3580.  A one-unit rise in x1 multiplies the odds of the median,
#' ## m/(1-m), by 2.36.  That is not an event odds ratio and not a mean effect.
#'
#' me <- marginal_effects(fit)
#' me
#' ## x1  effect  0.16834  (s.e. 0.005758)
#' ## x21 effect -0.10233  (s.e. 0.017803)
#' ## A one-unit rise in x1 lifts the conditional MEDIAN by 0.168 on average,
#' ## on the (0,1) scale of the response.  The factor-2.36 odds statement and
#' ## the 0.168 response-scale statement describe the same coefficient.
#'
#' ## -- the two are linked by the chain rule --------------------------------
#' b   <- coef(fit)[["mu:x1"]]
#' muh <- fitted(fit)                       # fitted conditional medians
#' b * mean(muh * (1 - muh))                # 0.1683372
#' me$table$effect[me$table$variable == "x1"]   # identical
#'
#' ## -- the effect is not constant across the covariate space ---------------
#' obs <- marginal_effects(fit, at = "observed")
#' range(obs$effects[, "x1"])               # 0.0901 to 0.2144
#' colMeans(obs$effects)                    # reproduces the AME exactly
#' ## The AME is an average of effects that differ by more than twofold across
#' ## the sample, which is precisely why the odds-scale coefficient is constant
#' ## and the response-scale effect is not.
#'
#' ## -- at the covariate means instead --------------------------------------
#' marginal_effects(fit, at = "mem")
#' ## Larger in absolute value here (0.2143 for x1): the average linear
#' ## predictor sits nearer the middle of the logistic curve, where mu(1-mu)
#' ## attains its maximum, than the typical observation does.
#'
#' ## -- the delta method versus the shortcut --------------------------------
#' se_shortcut <- sqrt(diag(vcov(fit)))[["mu:x1"]] * mean(muh * (1 - muh))
#' se_exact    <- me$table$std.error[me$table$variable == "x1"]
#' c(shortcut = se_shortcut, exact = se_exact, ratio = se_shortcut / se_exact)
#' ## shortcut 0.009192, exact 0.005758, ratio 1.596.  Freezing mu(1-mu) at its
#' ## average ignores that it, too, depends on the coefficients, and the two
#' ## dependencies partially cancel.
#'
#' ## -- robust and subpopulation variants -----------------------------------
#' marginal_effects(fit, vcov. = vcov(fit, type = "sandwich"))
#' marginal_effects(fit, variables = "x1", newdata = subset(dat, x2 == "1"))
#' @export
marginal_effects <- function(object, variables = NULL,
                             at = c("ame", "mem", "observed"),
                             newdata = NULL, level = 0.95, vcov. = NULL, ...) {
  stopifnot(inherits(object, "gkwqreg"))
  at <- match.arg(at)
  b_mu <- object$coef_list[["mu"]]
  vars <- variables %||% setdiff(names(b_mu), "(Intercept)")
  vars <- intersect(vars, names(b_mu))
  if (!length(vars)) {
    stop("no covariates in the quantile part to report effects for.", call. = FALSE)
  }

  ## Design matrix and offset for the quantile part, so the effect can be
  ## written as an explicit function of the coefficient vector.
  if (is.null(newdata)) {
    Xmu <- stats::model.matrix(object, part = "mu")
    offmu <- object$offsets[["mu"]]
  } else {
    nd <- .gkwq_newdata_matrices(object, newdata)
    Xmu <- nd$X[["mu"]]
    offmu <- nd$offsets[["mu"]]
  }
  lnk <- object$link[["mu"]]
  scl <- object$link_scale[["mu"]]
  idx_mu <- seq_along(b_mu)          # the mu block is first in the stacked vector

  if (at == "observed") {
    eta <- as.numeric(Xmu %*% b_mu) + offmu
    dmu <- .gkwq_mu_eta(eta, lnk, scl)
    out <- vapply(vars, function(v) b_mu[[v]] * dmu, numeric(length(dmu)))
    return(structure(list(effects = out, at = at, tau = object$tau),
                     class = "gkwq_meff_observed"))
  }

  ## The effect as an EXACT function of the whole coefficient vector. Both
  ## factors depend on it: the slope b_j directly, and mean(dmu/deta) through
  ## eta = X b. Holding dmu/deta fixed at its average -- the obvious shortcut --
  ## drops the second term and is NOT the delta method. The two terms partially
  ## cancel (raising b_j spreads eta, which lowers mean(dmu/deta) once eta is off
  ## the mode), so the shortcut overstates the standard error; on a logit fit at
  ## tau = 0.9 it was 2.4 times too large.
  eff_fun <- function(cf) {
    eta <- as.numeric(Xmu %*% cf[idx_mu]) + offmu
    if (at == "mem") eta <- mean(eta)
    dmu <- .gkwq_mu_eta(eta, lnk, scl)
    vapply(vars, function(v) cf[[paste0("mu:", v)]] * mean(dmu), numeric(1))
  }

  cf_full <- object$coefficients
  est <- eff_fun(cf_full)

  V <- if (!is.null(vcov.)) vcov. else tryCatch(vcov(object), error = function(e) NULL)
  se <- rep(NA_real_, length(vars))
  if (!is.null(V)) {
    J <- numDeriv::jacobian(eff_fun, cf_full)
    Vme <- J %*% V %*% t(J)
    se <- sqrt(pmax(diag(Vme), 0))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)

  ## Does any of these also drive a nuisance part?
  overlap <- vars[vapply(vars, function(v) {
    any(vapply(object$parts[-1L], function(p) v %in% names(object$coef_list[[p]]),
               logical(1)))
  }, logical(1))]

  structure(list(
    table = data.frame(variable = vars, effect = est, std.error = se,
                       lower = est - z * se, upper = est + z * se,
                       row.names = NULL),
    at = at, tau = object$tau, level = level, link = object$link[["mu"]],
    overlap = overlap
  ), class = "gkwq_meff")
}

#' @export
print.gkwq_meff <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("\nMarginal effects on the conditional %s-quantile (%s)\n",
              format(x$tau), switch(x$at, ame = "averaged over the sample",
                                    mem = "at covariate means")))
  cat(sprintf("dQ(tau|x)/dx_j, %s link\n\n", x$link))
  print(format(x$table, digits = digits), row.names = FALSE)
  if (length(x$overlap)) {
    cat(sprintf("\nNote: %s also appear(s) in a nuisance part. The effect on the\n",
                paste(x$overlap, collapse = ", ")))
    cat("      quantile is still exactly the above: mu IS the quantile, so the\n")
    cat("      nuisance part changes the spread, not the quantile modelled.\n")
  }
  invisible(x)
}
