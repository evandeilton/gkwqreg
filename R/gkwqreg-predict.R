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
    mf <- stats::model.frame(mt, newdata, na.action = na.action,
                             xlev = object$levels)
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

#' Predictions from a quantile regression fit
#'
#' @param object A `"gkwqreg"` fit.
#' @param newdata Optional data frame; defaults to the fitting data.
#' @param type One of `"quantile"` (the default, equivalently `"response"`),
#'   `"link"`, `"parameter"`, `"mu"`, `"density"`, `"probability"`, `"terms"`,
#'   `"mean"`, `"variance"`.
#' @param at Evaluation points for `type = "density"` or `"probability"`.
#' @param tau Optional quantile level(s). **This predicts a different quantile of
#'   the same fitted conditional distribution.** Refitting at another level gives
#'   a different model; the two agree only if the family is correctly specified,
#'   and the gap between them is itself a specification diagnostic (see
#'   [check_crossing()]).
#' @param elementwise If `TRUE`, `at` is paired with rows rather than crossed
#'   with them.
#' @param na.action Handling of missing values in `newdata`.
#' @param ... Unused.
#'
#' @return A vector, matrix or data frame depending on `type`.
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
#' The coefficients of the conditional quantile are effects on the log quantile
#' odds, which is the single most likely thing for a reader to misinterpret as a
#' mean effect. This reports the effect on the quantile itself,
#' `dQ(tau|x)/dx_j = beta_j * dmu/deta`, which under the logit link is
#' `beta_j * mu * (1 - mu)`.
#'
#' If a covariate also appears in a nuisance part, the effect on the quantile is
#' still exactly the expression above: `mu` *is* the quantile, so the nuisance
#' part changes the spread of the conditional distribution, not the quantile
#' being modelled. That is a genuine convenience of this parametrization, and
#' the printed output says so when it detects the overlap.
#'
#' @param object A `"gkwqreg"` fit.
#' @param variables Covariates to report; defaults to all in the quantile part.
#' @param at `"ame"` averages the effect over the sample, `"mem"` evaluates it at
#'   the covariate means, `"observed"` returns the full vector per covariate.
#' @param newdata Optional data to evaluate at.
#' @param level Confidence level for the delta-method interval.
#' @param ... Unused.
#'
#' @return An object of class `"gkwq_meff"`.
#' @export
marginal_effects <- function(object, variables = NULL,
                             at = c("ame", "mem", "observed"),
                             newdata = NULL, level = 0.95, ...) {
  stopifnot(inherits(object, "gkwqreg"))
  at <- match.arg(at)
  b_mu <- object$coef_list[["mu"]]
  vars <- variables %||% setdiff(names(b_mu), "(Intercept)")
  vars <- intersect(vars, names(b_mu))
  if (!length(vars)) {
    stop("no covariates in the quantile part to report effects for.", call. = FALSE)
  }

  em <- .gkwq_eta_mu_theta(object, newdata)
  eta <- em$eta[["mu"]]
  if (at == "mem") {
    eta <- mean(eta)
  }
  dmu <- .gkwq_mu_eta(eta, object$link[["mu"]], object$link_scale[["mu"]])

  if (at == "observed") {
    out <- vapply(vars, function(v) b_mu[[v]] * dmu, numeric(length(dmu)))
    return(structure(list(effects = out, at = at, tau = object$tau),
                     class = "gkwq_meff_observed"))
  }

  dbar <- mean(dmu)
  est <- vapply(vars, function(v) b_mu[[v]] * dbar, numeric(1))

  ## Delta method treating dmu/deta as fixed at its average: the leading term,
  ## and the one that matters at the sample sizes this is used at.
  V <- tryCatch(vcov(object), error = function(e) NULL)
  se <- rep(NA_real_, length(vars))
  if (!is.null(V)) {
    idx <- match(paste0("mu:", vars), rownames(V))
    ok <- !is.na(idx)
    se[ok] <- dbar * sqrt(pmax(diag(V)[idx[ok]], 0))
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
