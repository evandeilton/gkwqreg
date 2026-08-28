## ---------------------------------------------------------------------------
## S3 surface for class "gkwqreg".
##
## The class deliberately does NOT inherit from "gkwreg". Inheritance would make
## every missing method silently dispatch mean semantics on a quantile object,
## and the most dangerous case is fitted(): here it returns conditional
## quantiles, there conditional means. A missing method must error, not guess.
## ---------------------------------------------------------------------------

#' Extractor methods for a quantile regression fit
#'
#' The standard model accessors, adapted to a model whose coefficients are split
#' across several parts. `coef()`, `terms()` and `model.matrix()` therefore take
#' a `part` argument; the rest behave as their generics lead one to expect.
#'
#' @details
#' A `"gkwqreg"` object carries one linear predictor per modelled parameter. The
#' first part is always `mu`, the conditional \eqn{\tau}-quantile; the remaining
#' parts are the family's free parameters with the anchored one removed, as
#' reported by [gkwq_parts()]. Coefficients are stored in a single stacked
#' vector whose names have the form `"part:term"`, so `coef(object)` returns
#' everything at once and `coef(object, part = "mu")` returns just the quantile
#' block with its bare term names.
#'
#' `logLik()`, `AIC()` and `BIC()` describe the model at **one** quantile level.
#' They are comparable across families, anchors and covariate sets fitted at the
#' same `tau`, and are not comparable across levels: different levels answer
#' different questions of the same data and their likelihoods are not on a
#' common scale.
#'
#' `family()` returns a small description object rather than a
#' [stats::family()] object, because a Generalized Kumaraswamy quantile model has
#' no single variance function or canonical link to report. It records the
#' family, the anchor, the quantile level, the part names and the link used by
#' each part, and prints them compactly.
#'
#' Some accessors depend on what the fit was asked to retain.
#' `model.frame()` requires `model = TRUE` (the default) and errors otherwise;
#' `model.matrix()` uses the stored design matrices when the fit was made with
#' `x = TRUE`, and otherwise rebuilds them from the model frame, so it needs one
#' or the other. Each error message names the argument to set on the refit.
#'
#' @section What `update()` does and does not carry over:
#' `update()` modifies the recorded call and re-evaluates it, so any argument of
#' [gkwqreg()] can be replaced by name -- `tau`, `family`, `anchor`, `data`,
#' `control` and so on. Changing `tau` refits at the new level, which is a
#' different model and not the same thing as reading another level off the
#' present fit; see [predict.gkwqreg()] for that distinction. Changing `family`
#' or `anchor` may change the number of formula parts the model expects, in
#' which case the formula must be updated in the same call.
#'
#' @param object,x,formula A `"gkwqreg"` fit, as returned by [gkwqreg()]. The
#'   three names correspond to the argument names of the respective generics.
#' @param part Which model part to extract, named as in [gkwq_parts()]:
#'   `"mu"` for the conditional quantile, then the family's remaining
#'   parameters. Partial matching applies. `coef()` accepts a vector of parts
#'   and returns a named list in that case, and `part = NULL` (its default)
#'   returns the whole stacked coefficient vector. `terms()` and
#'   `model.matrix()` take exactly one part and default to `"mu"`.
#' @param k The penalty per parameter in `AIC()`. The default `k = 2` gives the
#'   usual Akaike criterion; `k = log(nobs(object))` reproduces `BIC()`.
#' @param level Confidence level recorded by `summary()` and used when its
#'   result is printed or post-processed.
#' @param vcov_type Which covariance estimator `summary()` should use for its
#'   standard errors, one of `"expected"`, `"observed"` or `"sandwich"`; see
#'   [vcov.gkwqreg()]. `NULL`, the default, uses the type recorded in the fit's
#'   control object. If the requested type cannot be computed, `summary()` falls
#'   back to the stored matrix rather than failing.
#' @param formula. A formula update specification for `update()`, in the sense
#'   of [stats::update.formula()], applied to the multi-part model formula.
#' @param evaluate If `TRUE` (the default) `update()` evaluates the modified
#'   call and returns the new fit; if `FALSE` it returns the unevaluated call,
#'   which is useful for inspection or for deferred evaluation.
#' @param ... For `update()`, further named arguments replacing those of the
#'   original call. Unused by the other methods.
#'
#' @return
#' * `coef()`: a named numeric vector. The whole stacked vector, with names
#'   `"part:term"`, when `part` is `NULL`; a single part's coefficients, with
#'   bare term names, when one part is given; a named list of such vectors when
#'   several are given.
#' * `logLik()`: an object of class `"logLik"` with attributes `df`
#'   (the number of estimated coefficients) and `nobs`.
#' * `nobs()`: a single integer, the number of observations used in the fit.
#' * `AIC()`, `BIC()`: a single number.
#' * `family()`: an object of class `"gkwq_family"`, a list with components
#'   `family`, `anchor`, `tau`, `parts` and `link`, with a `print` method.
#' * `formula()`: the multi-part model formula as supplied.
#' * `terms()`: the `"terms"` object of one part.
#' * `model.frame()`: the stored model frame.
#' * `model.matrix()`: the design matrix of one part, with an `"assign"`
#'   attribute and contrasts as recorded at fit time.
#' * `getCall()`: the matched call that produced the fit.
#' * `update()`: a new `"gkwqreg"` fit, or the unevaluated call when
#'   `evaluate = FALSE`.
#' * `summary()`: an object of class `"summary.gkwqreg"` carrying the
#'   coefficient table, the fit statistics, the pinball loss, the pseudo-R1, the
#'   empirical coverage and the information condition number, with a `print`
#'   method.
#'
#' @seealso [gkwqreg()] for the model, [gkwq_parts()] for the part contract,
#'   [vcov.gkwqreg()], [confint.gkwqreg()], [fitted.gkwqreg()].
#'
#' @rdname gkwqreg-extractors
#' @examples
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
#' ## -- coefficients, whole and by part -------------------------------------
#' round(coef(fit), 4)                 # names are "part:term"
#' round(coef(fit, part = "mu"), 4)    # just the conditional-median equation
#' names(coef(fit, part = c("mu", "alpha")))
#'
#' ## -- what the model is -----------------------------------------------------
#' family(fit)
#' gkwq_parts("kw")                    # the part contract this fit obeys
#' formula(fit)
#'
#' ## -- fit statistics, all at this one quantile level ----------------------
#' c(logLik = as.numeric(logLik(fit)), df = attr(logLik(fit), "df"),
#'   n = nobs(fit), AIC = AIC(fit), BIC = BIC(fit))
#' all.equal(AIC(fit, k = log(nobs(fit))), BIC(fit))   # TRUE
#'
#' ## -- design matrices differ by part --------------------------------------
#' fitx <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw",
#'                 x = TRUE)
#' colnames(model.matrix(fitx, part = "mu"))       # intercept, x1, x21
#' colnames(model.matrix(fitx, part = "alpha"))    # intercept, x1
#' ## This is why sandwich::vcovHC() cannot work here: there is no single
#' ## design matrix whose columns match the score matrix.
#'
#' ## -- summary, with a covariance estimator of your choosing ---------------
#' summary(fit, vcov_type = "sandwich")
#'
#' ## -- update() re-evaluates the call --------------------------------------
#' getCall(fit)
#' update(fit, . ~ x1 | x1, evaluate = FALSE)      # inspect before running
#' fit_simpler <- update(fit, . ~ x1 | x1)
#' c(full = AIC(fit), simpler = AIC(fit_simpler))
#' @export
coef.gkwqreg <- function(object, part = NULL, ...) {
  if (is.null(part)) return(object$coefficients)
  part <- match.arg(part, object$parts, several.ok = TRUE)
  if (length(part) == 1L) object$coef_list[[part]] else object$coef_list[part]
}

#' Covariance matrix of a quantile regression fit
#'
#' Returns the estimated covariance matrix of the full stacked coefficient
#' vector -- the quantile part and every nuisance part together -- either from
#' the observed information alone or in sandwich form.
#'
#' @details
#' All three types are built from the same observed information matrix
#' \deqn{\hat{H} \;=\; \left. \frac{\partial^2 \, (-\ell)}
#'   {\partial\theta \, \partial\theta^{\top}} \right|_{\theta = \hat\theta},}
#' evaluated by [stats::optimHess()] applied to the exact automatic-
#' differentiation gradient supplied by `TMB` and stored in the fit at
#' estimation time. Two routes that look shorter are deliberately not taken.
#' The Hessian is never read off `obj$he()`, which would require second-order
#' differentiation through the incomplete-beta atomic; and it is never obtained
#' from the naive \eqn{J^{\top} H J} transformation of the unreparametrized
#' Hessian, which omits the curvature term
#' \eqn{\sum_k g_k \, \nabla^2 \theta_k} and is simply wrong once the parameters
#' vary with covariates.
#'
#' @section The three types:
#'
#' | `type` | Estimator | Appropriate when |
#' | :----- | :-------- | :--------------- |
#' | `"expected"` (default) | \eqn{\hat{H}^{-1}} | the family is taken to be correctly specified |
#' | `"observed"` | \eqn{\hat{H}^{-1}}, identical to the above | provided for interface compatibility |
#' | `"sandwich"` | \eqn{\hat{H}^{-1} \left(\sum_i s_i s_i^{\top}\right) \hat{H}^{-1}} | the conditional distribution may be misspecified |
#'
#' `"expected"` and `"observed"` return **the same matrix**. Nothing here
#' evaluates a Fisher information by expectation: the two names are kept so that
#' code written for other regression classes runs unchanged, and the label is
#' recorded in `summary()` for the record. If the distinction matters to the
#' argument being made, say "observed information" and mean it.
#'
#' The sandwich form replaces the middle of the expression by the empirical
#' outer product of the per-observation scores \eqn{s_i} returned by
#' [estfun.gkwqreg()]. It is consistent for the covariance of the maximum
#' likelihood estimator when the family is wrong but the quantile equation is
#' right, at the cost of a noisier estimate in small samples. A large gap
#' between the two is itself informative: under correct specification the
#' information equality makes them agree up to sampling error, so a systematic
#' discrepancy is evidence against the family.
#'
#' @section Why the whole matrix matters:
#' The conditional quantile is **not** orthogonal to the remaining parameters in
#' the Cox-Reid sense. The information matrix has non-negligible off-diagonal
#' blocks linking \eqn{\beta_\mu} to the nuisance coefficients, so uncertainty in
#' the shape parameters propagates into the quantile coefficients and vice
#' versa. Two consequences follow. Any linear or non-linear function that mixes
#' parts -- a contrast, a marginal effect, a predicted quantile -- must be given
#' the full matrix, not a per-part block. And per-part standard errors reported
#' in isolation understate what is actually known jointly; use
#' [confint.gkwqreg()] with `method = "profile"` or `"boot"` when the Wald
#' approximation is doing too much work.
#'
#' @param object A `"gkwqreg"` fit, as returned by [gkwqreg()].
#' @param type Which covariance estimator to return; see the table below.
#'   Partial matching applies. All three require that the fit carried out the
#'   Hessian computation, that is, `control = gkwq_control(hessian = TRUE)`,
#'   which is the default; otherwise the function errors with instructions to
#'   refit.
#' @param ... Currently unused; present for consistency with the generic.
#'
#' @return A symmetric numeric matrix of dimension `p` by `p`, where
#'   `p = object$npar` is the total number of estimated coefficients across all
#'   parts. Both dimnames are `names(coef(object))`, whose entries have the form
#'   `"part:term"` -- `"mu:(Intercept)"`, `"mu:x1"`, `"alpha:x1"` and so on --
#'   so blocks can be extracted by `grep()` on the part prefix.
#'
#' @seealso [estfun.gkwqreg()] and [bread.gkwqreg()] for the ingredients of the
#'   sandwich and for interoperability with the \pkg{sandwich} package,
#'   [confint.gkwqreg()] for intervals that do not rely on the Wald
#'   approximation, [summary.gkwqreg()] which uses this to build its coefficient
#'   table.
#'
#' @examples
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
#' V <- vcov(fit)
#' dim(V)                                   # 5 x 5: three mu, two alpha
#' rownames(V)                              # "mu:(Intercept)" ... "alpha:x1"
#' identical(V, vcov(fit, type = "observed"))   # TRUE, by construction
#'
#' ## -- the quantile block is NOT orthogonal to the nuisance block ----------
#' round(cov2cor(V), 3)
#' ## cor(mu:(Intercept), alpha:(Intercept)) = 0.319 and
#' ## cor(mu:x1, alpha:x1) = 0.468.  Uncertainty about the shape parameter is
#' ## not separable from uncertainty about the median; this is why anything
#' ## mixing parts must be handed the whole matrix.
#'
#' ## -- model-based versus sandwich standard errors -------------------------
#' Vs <- vcov(fit, type = "sandwich")
#' round(cbind(model = sqrt(diag(V)), sandwich = sqrt(diag(Vs)),
#'             ratio = sqrt(diag(Vs) / diag(V))), 4)
#' ## Ratios between 0.91 and 1.11 here.  The family is correctly specified in
#' ## this simulation, so the information equality holds and the two agree up
#' ## to sampling error, exactly as it should.
#'
#' ## -- a Wald table built by hand from whichever matrix you trust ----------
#' est <- coef(fit); se <- sqrt(diag(Vs))
#' round(cbind(Estimate = est, `Std. Error` = se, z = est / se,
#'             `Pr(>|z|)` = 2 * pnorm(-abs(est / se))), 4)
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
#' Returns the empirical estimating functions of the fit: the gradient of each
#' observation's log-likelihood contribution with respect to the full stacked
#' coefficient vector. This is the "meat" of the sandwich covariance, and
#' supplying it in the conventional shape is what makes `sandwich::sandwich()`,
#' `sandwich::vcovCL()` and `lmtest::coeftest()` work on `"gkwqreg"` fits with no
#' further glue.
#'
#' @details
#' The `i`th row is
#' \deqn{s_i(\hat\theta) \;=\; \left. \frac{\partial \, \log f(y_i \mid x_i;
#'   \theta)}{\partial \theta} \right|_{\theta = \hat\theta},}
#' obtained by applying [numDeriv::jacobian()] to the vector of per-observation
#' log-likelihood contributions that the `TMB` object reports. Numerical
#' differentiation of the reported likelihood is used rather than a second
#' automatic-differentiation pass because the tape returns the summed objective;
#' the per-observation decomposition lives in the report, not on the tape.
#'
#' At the maximum the column sums are zero up to optimizer and differentiation
#' tolerance, which is a useful convergence check in its own right: column sums
#' that are not small say the reported optimum is not a stationary point.
#'
#' The computation goes through the fit's `TMB` object, held in `object$obj`. If
#' that component is absent -- it is dropped by anything that strips the fit down
#' for storage -- the function stops with a message asking for a refit in the
#' current session, rather than returning numbers it cannot stand behind.
#'
#' @section What does not work:
#' `sandwich::vcovHC()` is **not** supported and will error with "cannot match
#' dimension of model.matrix and estfun". It needs working residuals and a
#' single design matrix whose columns line up one-for-one with the scores. This
#' is a multi-part model with one design matrix *per part*, so
#' `model.matrix()` returns the `mu` block while `estfun()` spans every
#' coefficient of every part, and the two cannot be made to line up. There is no
#' sensible thing for a heteroskedasticity-consistent correction to do here.
#' Use `vcov(object, type = "sandwich")`, which is the same estimator computed
#' correctly, or `sandwich::vcovCL()` for clustered data.
#'
#' @param x A `"gkwqreg"` fit, as returned by [gkwqreg()], still carrying its
#'   `TMB` object in `x$obj`.
#' @param ... Currently unused.
#'
#' @return A numeric matrix with `nobs(x)` rows and `length(coef(x))` columns,
#'   the columns named as `names(coef(x))` (`"part:term"`). Row `i` holds the
#'   score contribution of observation `i`.
#'
#' @seealso [bread.gkwqreg()] for the other half of the sandwich,
#'   [vcov.gkwqreg()] for the assembled estimator.
#'
#' @examples
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
#' S <- estfun.gkwqreg(fit)
#' dim(S)                          # 300 x 5: one row per observation
#' round(colSums(S), 6)            # all zero to 1e-4: a stationary point
#'
#' ## -- the sandwich, assembled by hand -------------------------------------
#' B <- bread.gkwqreg(fit)                       # n * H^{-1}
#' Vs <- B %*% crossprod(S) %*% B / nobs(fit)^2
#' all.equal(Vs, vcov(fit, type = "sandwich"), check.attributes = FALSE)
#' ## TRUE.  vcov(type = "sandwich") is exactly H^{-1} (sum s_i s_i') H^{-1},
#' ## and is also what sandwich::sandwich() returns for this class.
#'
#' ## -- the scores carry the clustering information -------------------------
#' grp <- rep(1:30, each = 10)
#' M   <- crossprod(rowsum(S, grp))              # clustered meat
#' Vcl <- (B %*% M %*% B) / nobs(fit)^2
#' round(sqrt(diag(Vcl)), 4)
#' ## The same construction sandwich::vcovCL() performs; shown here so the
#' ## estimator is legible without taking a dependency on that package.
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

#' Bread matrix for the sandwich covariance
#'
#' Returns the inverse observed information scaled by the sample size, in the
#' normalization the \pkg{sandwich} package expects. Together with
#' [estfun.gkwqreg()] this makes `sandwich::sandwich()`, `sandwich::vcovCL()` and
#' `lmtest::coeftest()` work on `"gkwqreg"` fits without any further glue.
#'
#' @details
#' The value is
#' \deqn{\widehat{\mathrm{bread}} \;=\; n \, \hat{H}^{-1}
#'   \;=\; \left(\frac{1}{n}\hat{H}\right)^{-1},}
#' the inverse of the *average* observed information, which is the quantity that
#' converges to a fixed matrix as \eqn{n} grows. That normalization is what lets
#' \pkg{sandwich} assemble
#' \eqn{\widehat{\mathrm{bread}} \, \widehat{\mathrm{meat}} \,
#' \widehat{\mathrm{bread}} / n} and recover the same estimator as
#' `vcov(object, type = "sandwich")`.
#'
#' `sandwich::vcovHC()` is not supported for this class; see
#' [estfun.gkwqreg()] for why.
#'
#' @param x A `"gkwqreg"` fit that carries an observed information matrix, that
#'   is, one fitted with `gkwq_control(hessian = TRUE)` (the default). Otherwise
#'   the function errors.
#' @param ... Currently unused.
#'
#' @return A symmetric numeric matrix of dimension `p` by `p`, where
#'   `p = length(coef(x))`, with dimnames inherited from `vcov(x)`.
#'
#' @seealso [estfun.gkwqreg()], [vcov.gkwqreg()].
#'
#' @examples
#' set.seed(2024)
#' n   <- 300
#' x1  <- runif(n, -2, 2)
#' mu  <- plogis(0.3 + 0.9 * x1)
#' bt  <- log1p(-0.5) / log1p(-mu^2)
#' y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
#' dat <- data.frame(y = y, x1 = x1)
#'
#' fit <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw")
#'
#' B <- bread.gkwqreg(fit)
#' all.equal(B / nobs(fit), vcov(fit), check.attributes = FALSE)   # TRUE
#'
#' ## The bread alone gives the model-based standard errors back:
#' round(sqrt(diag(B / nobs(fit))), 4)
#' round(sqrt(diag(vcov(fit))), 4)
#' @export
bread.gkwqreg <- function(x, ...) {
  if (is.null(x$vcov)) stop("no observed information available.", call. = FALSE)
  x$vcov * x$nobs
}

#' @rdname gkwqreg-extractors
#' @export
logLik.gkwqreg <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$nobs,
            class = "logLik")
}

#' @rdname gkwqreg-extractors
#' @export
nobs.gkwqreg <- function(object, ...) object$nobs

#' @rdname gkwqreg-extractors
#' @export
AIC.gkwqreg <- function(object, ..., k = 2) {
  -2 * object$loglik + k * object$npar
}

#' @rdname gkwqreg-extractors
#' @export
BIC.gkwqreg <- function(object, ...) object$bic

#' Fitted conditional quantiles
#'
#' **Returns the fitted conditional \eqn{\tau}-quantile, not the conditional
#' mean.** This is the single most consequential difference between a
#' `"gkwqreg"` object and the mean-parametrized `gkwreg` fits it otherwise
#' resembles, and the class deliberately does not inherit from that one so that
#' the confusion cannot happen by silent dispatch.
#'
#' @details
#' The default value is
#' \deqn{\hat{y}_i \;=\; \widehat{Q}(\tau \mid x_i) \;=\; \hat\mu_\tau(x_i),}
#' the fitted \eqn{\tau}-quantile for observation \eqn{i}, taken straight from
#' the quantile part of the linear predictor through its inverse link. It is not
#' a prediction of \eqn{y_i} in the least-squares sense and should not be
#' averaged, differenced against \eqn{y_i}, or fed to a residual routine that
#' assumes a mean. Use [residuals.gkwqreg()] for residuals defined on the right
#' scale.
#'
#' The operational meaning is a coverage statement: a fraction \eqn{\tau} of the
#' observations should fall at or below their own fitted value. That check is
#' one line and worth running, and `summary()` reports it as the empirical
#' coverage. A conditional mean satisfies no such property, which is what the
#' examples make concrete.
#'
#' @section Values of `type`:
#'
#' | `type` | What is returned | Shape |
#' | :----- | :--------------- | :---- |
#' | `"quantile"` (default) | the fitted conditional tau-quantile | numeric vector of length `nobs(object)` |
#' | `"mean"` | the conditional mean by 64-node Gauss-Legendre quadrature of the fitted quantile function | numeric vector of length `nobs(object)` |
#' | `"parameter"` | the reconstructed distribution parameters, one row per observation | data frame with columns `alpha`, `beta`, `gamma`, `delta`, `lambda` |
#'
#' `type = "mean"` exists so that a quantile fit can be placed beside a
#' mean-parametrized one, and so that a mean-based diagnostic borrowed from
#' elsewhere can be reproduced deliberately rather than by accident. It is never
#' the default, and it is never what the model was estimated to get right.
#'
#' `type = "parameter"` returns the fitted conditional law itself. For each row
#' the anchored parameter has been solved so that the distribution has
#' \eqn{\hat\mu_\tau(x_i)} as its exact \eqn{\tau}-quantile, so these values can
#' be handed directly to the `gkwdist` density, distribution and random-number
#' functions.
#'
#' @param object A `"gkwqreg"` fit, as returned by [gkwqreg()].
#' @param type Which fitted quantity to return; see the table below. Partial
#'   matching applies.
#' @param ... Currently unused; present for consistency with the generic.
#'
#' @return For `type = "quantile"` and `type = "mean"`, a numeric vector of
#'   length `nobs(object)`, aligned with the rows retained after `na.action`.
#'   For `type = "parameter"`, a data frame with `nobs(object)` rows and the
#'   five columns `alpha`, `beta`, `gamma`, `delta`, `lambda` in `gkwdist`
#'   order.
#'
#' @seealso [predict.gkwqreg()] for the same quantities at new covariate values
#'   and for further quantile levels, [residuals.gkwqreg()] for quantile
#'   residuals, [gkwqreg()] for the model.
#'
#' @examples
#' set.seed(2024)
#' n   <- 300
#' x1  <- runif(n, -2, 2)
#' x2  <- rbinom(n, 1, 0.5)
#' mu  <- plogis(0.3 + 0.9 * x1 - 0.5 * x2)
#' bt  <- log1p(-0.5) / log1p(-mu^2)
#' y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
#' dat <- data.frame(y = y, x1 = x1, x2 = factor(x2))
#'
#' ## Fit the ninth decile, where quantile and mean are far apart.
#' fit90 <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.9, family = "kw")
#'
#' ## -- what fitted() actually promises -------------------------------------
#' mean(dat$y <= fitted(fit90))                  # 0.8900, target tau = 0.9
#' mean(dat$y <= fitted(fit90, type = "mean"))   # 0.5067
#' ## Ninety per cent of the sample lies below its own fitted value, as a 0.9
#' ## quantile must. The conditional mean covers only about half. Reading
#' ## fitted() as a mean would misstate the estimand by four deciles here.
#'
#' mean(fitted(fit90))                # 0.6991
#' mean(fitted(fit90, "mean"))        # 0.4736
#'
#' ## -- the fitted conditional distribution ---------------------------------
#' pv <- fitted(fit90, type = "parameter")
#' head(pv, 3)
#' ## Each row is a Kumaraswamy law (gamma = 1, delta = 0, lambda = 1) whose
#' ## 0.9-quantile is exactly the corresponding fitted value:
#' all.equal(gkwq_quantile(0.9, pv$alpha, pv$beta, pv$gamma, pv$delta,
#'                         pv$lambda),
#'           fitted(fit90))
#' ## TRUE. beta varies across rows because it is the anchor: it is recomputed
#' ## from the fitted quantile rather than estimated.
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

#' @rdname gkwqreg-extractors
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

#' @rdname gkwqreg-extractors
#' @export
formula.gkwqreg <- function(x, ...) x$formula

#' @rdname gkwqreg-extractors
#' @export
terms.gkwqreg <- function(x, part = "mu", ...) {
  part <- match.arg(part, x$parts)
  x$terms[[part]]
}

#' @rdname gkwqreg-extractors
#' @export
model.frame.gkwqreg <- function(formula, ...) {
  if (is.null(formula$model)) {
    stop("the model frame was not retained; refit with model = TRUE.",
         call. = FALSE)
  }
  formula$model
}

#' @rdname gkwqreg-extractors
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

#' @rdname gkwqreg-extractors
#' @export
getCall.gkwqreg <- function(x, ...) x$call

#' @rdname gkwqreg-extractors
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
#' Computes confidence intervals for the coefficients of a `"gkwqreg"` fit by
#' the Wald approximation, by inverting the profile likelihood, or by the
#' bootstrap.
#'
#' @details
#' The intervals are for coefficients on the **link scale**: under the default
#' logit link an interval for a `mu` coefficient is an interval for an effect on
#' the log quantile odds. Exponentiating the limits gives an interval for the
#' multiplicative effect on the quantile odds \eqn{\mu_\tau/(1-\mu_\tau)}, which
#' is legitimate because the transformation is monotone. An interval for the
#' effect on the quantile itself is a different object, since it involves the
#' link derivative as well; use [marginal_effects()] for that.
#'
#' @section The three methods:
#'
#' | `method` | Construction | Cost | Use when |
#' | :------- | :----------- | :--- | :------- |
#' | `"wald"` (default) | estimate plus or minus `z * se`, from the stored standard errors | free | the sample is large and the log-likelihood is close to quadratic |
#' | `"profile"` | inversion of the `TMB` profile likelihood, one coefficient at a time | one profile per coefficient | the Wald approximation is suspect and the sample is small |
#' | `"boot"` | percentile interval from [gkwq_boot()] | `R` refits | the family may be misspecified, or the profile does not resolve |
#'
#' Wald intervals use `object$se`, the square roots of the diagonal of the
#' stored observed-information inverse, and are therefore always symmetric about
#' the estimate. This matters more here than in an ordinary generalized linear
#' model: the conditional quantile is not orthogonal to the remaining parameters
#' in the Cox-Reid sense, so the profile log-likelihood in a `mu` coefficient
#' inherits curvature from the nuisance parameters and can be visibly
#' asymmetric in small samples. Neither `"profile"` nor `"boot"` imposes
#' symmetry.
#'
#' `"boot"` passes `...` through to [gkwq_boot()], whose default resampling
#' scheme is `type = "pairs"`. That default is the right one for interval
#' construction: the appeal of parametric quantile regression is that it buys
#' efficiency under a distributional assumption, and a parametric bootstrap that
#' re-imposes the same assumption cannot say anything about whether it holds.
#' Set `type = "parametric"` only when the family is not in question.
#'
#' `"profile"` is the most accurate of the three when it succeeds. It is
#' not always able to interpolate the likelihood cut-off from the profile it
#' obtains -- the profile may simply not extend far enough on one side. Limits
#' that cannot be obtained come back as `NA` and the function warns, naming how
#' many coefficients were affected, rather than reporting a silently wrong
#' number. Fall back on `"wald"` or `"boot"` when that happens.
#'
#' @param object A `"gkwqreg"` fit, as returned by [gkwqreg()].
#' @param parm Coefficients to report, given either as a character vector of
#'   names from `names(coef(object))` -- which have the form `"part:term"`, as
#'   in `"mu:x1"` or `"alpha:(Intercept)"` -- or as integer positions into that
#'   vector. Defaults to every coefficient of every part.
#' @param level Confidence level, a single number in `(0,1)`. The default,
#'   `0.95`, gives equal-tailed limits at the 2.5th and 97.5th percentiles.
#' @param method How to construct the interval; see the table below. Partial
#'   matching applies.
#' @param R Number of bootstrap replicates, used only when `method = "boot"`.
#'   The default of 200 is adequate for a percentile interval at the usual
#'   levels; raise it for tail levels or for publication.
#' @param ... Passed to [gkwq_boot()] when `method = "boot"`, which is where
#'   `type` and `seed` are set. Ignored by the other methods.
#'
#' @return A numeric matrix with one row per entry of `parm` and two columns,
#'   the lower and upper limits. Row names are the coefficient names; column
#'   names are the percentages of the limits, for example `"2.5 %"` and
#'   `"97.5 %"`. Entries that could not be computed are `NA`.
#'
#' @seealso [vcov.gkwqreg()] for the covariance matrices behind the Wald
#'   interval, [gkwq_boot()] for the bootstrap, [marginal_effects()] for
#'   intervals on the scale of the response, [summary.gkwqreg()] for the
#'   accompanying coefficient table.
#'
#' @examples
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
#' ## -- Wald intervals for every coefficient of every part ------------------
#' round(confint(fit), 4)
#' ## mu:x1 lies in (0.7660, 0.9496); the true value used to simulate was 0.9.
#'
#' ## -- on the quantile-odds scale ------------------------------------------
#' round(exp(confint(fit, parm = "mu:x1")), 4)
#' ## (2.1512, 2.5847): a one-unit rise in x1 multiplies the odds of the median
#' ## by between 2.15 and 2.58. The transformation is monotone, so the limits
#' ## carry over; the effect on the median itself is marginal_effects().
#'
#' ## -- a narrower selection and a different level --------------------------
#' round(confint(fit, parm = c("mu:x1", "mu:x21"), level = 0.9), 4)
#'
#' \donttest{
#' ## -- percentile bootstrap, resampling observations -----------------------
#' round(confint(fit, parm = c("mu:x1", "mu:x21"), method = "boot",
#'               R = 100, seed = 7), 4)
#' ##          2.5 %  97.5 %
#' ## mu:x1   0.7829  0.9529
#' ## mu:x21 -0.7060 -0.2651
#' ## For mu:x1 these sit almost on top of the Wald limits, as they should when
#' ## the family is correctly specified and n = 300. For the factor mu:x21 the
#' ## upper limit is visibly further out; with R = 100 a percentile limit is
#' ## itself an estimate, so raise R before reading much into a single
#' ## discrepancy of this size.
#' }
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
      tryCatch(.gkwq_profile_ci(object$obj, i, level),
               error = function(e) c(NA_real_, NA_real_))
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

#' @rdname gkwqreg-extractors
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

#' Methods for a set of quantile regression fits
#'
#' Accessors for the `"gkwqregs"` container returned by [gkwqreg()] when `tau`
#' has length greater than one. The container holds one entirely independent fit
#' per quantile level in `object$fits`, and these methods arrange their output
#' with one column per level.
#'
#' @details
#' Levels are stored sorted and de-duplicated, and each fit is estimated on its
#' own. Nothing ties the levels together, which is what makes quantile crossing
#' possible and worth testing for; see [check_crossing()] and [rearrange()].
#'
#' `logLik()` deliberately fails. The container holds one log-likelihood per
#' level, and those are neither comparable across levels nor additive: summing
#' them would treat the same observations as independent evidence several times
#' over, and comparing them would compare answers to different questions. Take
#' the likelihood of an individual fit instead, as the error message says.
#'
#' @param object A `"gkwqregs"` container, as returned by [gkwqreg()] with a
#'   vector-valued `tau`.
#' @param ... Passed on to the corresponding method for each individual fit.
#'
#' @return
#' * `coef()`: a numeric matrix with one row per coefficient and one column per
#'   quantile level, the columns named by level.
#' * `fitted()`: a numeric matrix with one row per observation and one column
#'   per level, holding each fit's conditional quantiles.
#' * `summary()`: a list of `"summary.gkwqreg"` objects, one per level.
#' * `logLik()`: never returns; it raises an error explaining why.
#'
#' @seealso [gkwqreg()], [quantile_process()] for the quantile process on a
#'   grid, [check_crossing()] for crossing diagnostics,
#'   [predict.gkwqreg()] for reading several levels off a *single* fit.
#'
#' @rdname gkwqregs-methods
#' @examples
#' set.seed(2024)
#' n   <- 300
#' x1  <- runif(n, -2, 2)
#' mu  <- plogis(0.3 + 0.9 * x1)
#' bt  <- log1p(-0.5) / log1p(-mu^2)
#' y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
#' dat <- data.frame(y = y, x1 = x1)
#'
#' fits <- gkwqreg(y ~ x1, data = dat, tau = c(0.25, 0.5, 0.75),
#'                 family = "kw")
#' fits
#'
#' round(coef(fits), 4)
#' ## One column per level. The mu:x1 row is the covariate effect on the log
#' ## quantile odds at each level; its drift across columns is what a quantile
#' ## regression is for.
#'
#' Q <- fitted(fits)
#' dim(Q)
#' colnames(Q)                   # the levels, formatted: "0.25" "0.50" "0.75"
#' colMeans(dat$y <= Q)          # 0.25, 0.46, 0.73: coverage at each level
#'
#' ## Each level is a separate fit, so nothing forces the columns to be
#' ## ordered; check rather than assume.
#' mean(Q[, 1] <= Q[, 2] & Q[, 2] <= Q[, 3])   # 1: no crossing in this sample
#' @export
coef.gkwqregs <- function(object, ...) {
  cf <- sapply(object$fits, function(f) f$coefficients)
  colnames(cf) <- format(object$taus, trim = TRUE)
  cf
}

#' @rdname gkwqregs-methods
#' @export
fitted.gkwqregs <- function(object, ...) {
  out <- sapply(object$fits, function(f) f$fitted.values)
  colnames(out) <- format(object$taus, trim = TRUE)
  out
}

#' @rdname gkwqregs-methods
#' @export
logLik.gkwqregs <- function(object, ...) {
  stop("a `gkwqregs` container holds one likelihood per tau, and they are not ",
       "comparable or additive. Take logLik() of an individual fit, e.g. ",
       "logLik(object$fits[[1]]).", call. = FALSE)
}

#' @rdname gkwqregs-methods
#' @export
summary.gkwqregs <- function(object, ...) {
  lapply(object$fits, summary, ...)
}

## Profile-likelihood interval for one coefficient.
##
## TMB::tmbprofile() is used for the profile itself, but its confint() method is
## not: for these models it returns a two-column frame whose columns are named
## "value" and NA, with the PARAMETER grid under "value" and the objective under
## the NA name. TMB:::confint.tmbprofile() reads them the other way round and
## every interval comes back NA. The columns are therefore identified by
## content: the parameter grid is the monotone one.
##
## The interval is the usual likelihood-ratio set
##   {theta : 2 (l_max - l_profile(theta)) <= qchisq(level, 1)},
## found by linear interpolation on each side of the minimum of the profiled
## negative log-likelihood.
.gkwq_profile_ci <- function(obj, i, level = 0.95) {
  pr <- TMB::tmbprofile(obj, i, trace = FALSE)
  if (!is.data.frame(pr) || ncol(pr) < 2L) return(c(NA_real_, NA_real_))

  mono <- vapply(pr, function(z) !is.unsorted(z) || !is.unsorted(rev(z)),
                 logical(1))
  par_col <- which(mono)[1L]
  if (is.na(par_col)) return(c(NA_real_, NA_real_))
  obj_col <- setdiff(seq_len(ncol(pr)), par_col)[1L]

  x <- as.numeric(pr[[par_col]])
  f <- as.numeric(pr[[obj_col]])
  ok <- is.finite(x) & is.finite(f)
  x <- x[ok]; f <- f[ok]
  if (length(x) < 3L) return(c(NA_real_, NA_real_))

  o <- order(x); x <- x[o]; f <- f[o]
  k <- which.min(f)
  cut <- f[k] + stats::qchisq(level, df = 1) / 2

  ## Interpolate x as a function of f on each side of the minimum. The
  ## objective DECREASES with x on the left branch, so the pairs must be put in
  ## ascending f order first; `ties = "ordered"` would assert an ordering that
  ## the left branch does not have, and silently return NA.
  branch <- function(ix) {
    if (length(ix) < 2L) return(NA_real_)
    ff <- f[ix]; xx <- x[ix]
    if (max(ff) < cut) return(NA_real_)
    o <- order(ff)
    stats::approx(ff[o], xx[o], xout = cut, ties = mean)$y
  }
  c(branch(seq_len(k)), branch(k:length(x)))
}
