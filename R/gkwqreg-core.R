#' Parametric quantile regression for the Generalized Kumaraswamy family
#'
#' Fits a regression for the conditional quantile of a response in `(0,1)` at a
#' quantile level fixed in advance. For that level the distribution is
#' reparametrized so that one of its parameters *is* the conditional quantile,
#' which is then modelled directly; the remaining parameters may carry their own
#' regressions.
#'
#' @section What `fitted()` returns:
#' **[fitted()] returns the fitted conditional `tau`-quantile, not the
#' conditional mean.** This is the single most important difference from
#' `gkwreg`, which fits the same seven families with mean semantics. A
#' `type = "mean"` escape hatch exists but is never the default.
#'
#' @section The anchor is a modeling choice:
#' `anchor` names the parameter eliminated in favour of the conditional
#' quantile. It looks like a reparametrization and, when every remaining
#' parameter is regressed on at least the covariates used for `mu`, it is one:
#' the likelihood is identical whichever anchor is chosen. It stops being one as
#' soon as `mu` varies with covariates while a nuisance parameter is held
#' constant. Anchoring on `beta` then asserts "`alpha` constant,
#' `beta_i = f(mu_i, alpha)`"; anchoring on `alpha` asserts the reverse, and the
#' two trace different paths through parameter space. In the design study the
#' two differed by 131 in log-likelihood and by 41% in the coefficient of
#' interest.
#'
#' Two anchors on the same data give non-nested models of equal dimension:
#' compare them by AIC, BIC or Vuong's test, **not** by a likelihood-ratio test.
#' Regressing the nuisance parameters on the same covariates as `mu` makes the
#' fit largely anchor-insensitive, and is recommended whenever theory does not
#' dictate an anchor.
#'
#' @section Interpreting the coefficients:
#' Under the default logit link, `exp(beta_j)` is the multiplicative effect of a
#' one-unit increase in `x_j` on the **odds of the conditional `tau`-quantile**,
#' `mu_tau / (1 - mu_tau)`. It is not an effect on the odds of an event and not
#' an effect on a mean. The marginal effect on the quantile itself is
#' `beta_j * mu_tau * (1 - mu_tau)`; see [marginal_effects()].
#'
#' @param formula A multi-part formula, parts separated by `|`. **Part one is
#'   always the conditional quantile.** The remaining parts follow the family's
#'   parameter order with the anchored parameter removed; [gkwq_parts()] reports
#'   the exact contract. Omitted trailing parts default to `~ 1`.
#' @param data,subset,na.action,weights,offset,contrasts Standard model-frame
#'   arguments. Weights and offsets are applied inside the likelihood.
#' @param tau Quantile level in `(0,1)`. A vector returns a `"gkwqregs"`
#'   container of independent fits, one per level. Mandatory in the sense that it is never
#'   estimated: the profile likelihood in `tau` is exactly flat, so `tau` indexes
#'   the question, not the model.
#' @param family One of `"kw"`, `"ekw"`, `"kkw"`, `"bkw"`, `"gkw"`, `"mc"`,
#'   `"beta"`.
#' @param anchor The parameter eliminated in favour of the conditional quantile,
#'   or `NULL` for the family default. See the section above; this is a modeling
#'   argument.
#' @param link,link_scale Link functions and scale factors, either a single
#'   value or a named vector indexed by part. The quantile defaults to `"logit"`
#'   and must map to `(0,1)`; the rest default to `"log"`.
#' @param control An object from [gkwq_control()].
#' @param model,x,y Whether to return the model frame, design matrices and
#'   response in the fitted object.
#' @param ... Passed to [gkwq_control()] when `control` is left at its default.
#'
#' @return An object of class `"gkwqreg"`, or `"gkwqregs"` when `tau` has length
#'   greater than one.
#'
#' @references
#' Mitnik, P. A. and Baek, S. (2013). The Kumaraswamy distribution:
#' median-dispersion re-parameterizations. *Statistical Papers* **54**, 177-192.
#' \doi{10.1007/s00362-011-0417-y}
#'
#' Carrasco, J. M. F., Ferrari, S. L. P. and Cordeiro, G. M. (2010).
#' A generalized Kumaraswamy distribution.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- runif(n, -2, 2)
#' mu <- plogis(0.4 + 1.1 * x)
#' # generate under the beta anchor at the median, with alpha = 2
#' y <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
#'
#' fit <- gkwqreg(y ~ x, tau = 0.5, family = "kw")
#' summary(fit)
#' head(fitted(fit))          # conditional medians, NOT means
#'
#' @seealso [gkwq_parts()] for the formula contract, [marginal_effects()] for
#'   effects on the quantile scale, [quantile_process()] for a grid of levels,
#'   [check_crossing()] for crossing across levels.
#' @export
gkwqreg <- function(formula, data, tau = 0.5,
                    family = c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta"),
                    anchor = NULL, link = NULL, link_scale = NULL,
                    subset = NULL, weights = NULL, offset = NULL,
                    na.action = stats::na.omit, contrasts = NULL,
                    control = gkwq_control(), model = TRUE, x = FALSE,
                    y = TRUE, ...) {
  cl <- match.call()
  family <- match.arg(family)
  if (missing(data)) data <- environment(formula)
  if (length(list(...))) control <- do.call(gkwq_control, list(...))

  if (!is.numeric(tau) || !length(tau) || anyNA(tau) ||
      any(!is.finite(tau)) || any(tau <= 0 | tau >= 1)) {
    stop("`tau` must be one or more finite values strictly inside (0,1).",
         call. = FALSE)
  }

  if (length(tau) > 1L) {
    ev <- .gkwq_eval_mf_args(substitute(subset), substitute(weights),
                             substitute(offset), data, parent.frame())
    return(.gkwq_fit_many(cl, formula, data, sort(unique(tau)), family, anchor,
                          link, link_scale, ev$subset, ev$weights, ev$offset,
                          na.action, contrasts, control, model, x, y))
  }

  spec <- .gkwq_family_info(family, anchor)
  .gkwq_check_identifiability(spec)

  ev <- .gkwq_eval_mf_args(substitute(subset), substitute(weights),
                           substitute(offset), data, parent.frame())
  md <- .gkwq_model_data(formula, data, spec$parts, ev$subset, na.action,
                         ev$weights, ev$offset, contrasts, family, spec$anchor)
  yv <- .gkwq_validate_y(md$y, control$eps_y)

  links <- .gkwq_links(link, spec$parts)
  scales <- .gkwq_link_scales(link_scale, links, spec$parts)

  fit <- .gkwq_fit_one(yv, md, spec, tau, links, scales, control)

  out <- .gkwq_build(fit, cl, formula, md, spec, tau, links, scales, control,
                     model, x, y, yv)
  out
}

## ---------------------------------------------------------------------------
## The TMB data contract. X_k and beta_k are indexed by PART POSITION, so X1 is
## always the conditional quantile's design matrix.
## ---------------------------------------------------------------------------

.gkwq_tmb_data <- function(yv, md, spec, tau, links, scales, control) {
  n <- length(yv)
  parts <- spec$parts
  X <- vector("list", 5L)
  O <- vector("list", 5L)
  link_type <- integer(5L)
  link_scale <- numeric(5L)
  for (k in 1:5) {
    if (k <= length(parts)) {
      p <- parts[k]
      X[[k]] <- md$X[[p]]
      O[[k]] <- md$offsets[[p]]
      link_type[k] <- .GKWQ_LINKS[[links[[p]]]]
      link_scale[k] <- scales[[p]]
    } else {
      X[[k]] <- matrix(0, n, 0)
      O[[k]] <- rep(0, n)
      link_type[k] <- 1L
      link_scale[k] <- 1
    }
  }

  ## Hoist log z out of the observation loop when gamma and delta do not vary.
  ## For the families that need the incomplete beta this turns n atomic
  ## inversions into one, which is the difference between a usable and an
  ## unusable `gkw` fit.
  varying <- vapply(c("gamma", "delta"), function(nm) {
    if (!nm %in% parts) return(FALSE)
    ncol(md$X[[nm]]) > 1L || any(md$offsets[[nm]] != md$offsets[[nm]][1L])
  }, logical(1))
  z_is_scalar <- as.integer(!any(varying))

  list(y = yv, w = md$weights, tau = tau,
       anchor_code = spec$anchor_code, z_mode = spec$z_mode,
       delta_is_zero = spec$delta_is_zero,
       par_id = spec$par_id, fixed_val = spec$fixed_val,
       X1 = X[[1]], X2 = X[[2]], X3 = X[[3]], X4 = X[[4]], X5 = X[[5]],
       O1 = O[[1]], O2 = O[[2]], O3 = O[[3]], O4 = O[[4]], O5 = O[[5]],
       link_type = link_type, link_scale = link_scale,
       z_is_scalar = z_is_scalar, reportScores = 0L,
       eps_mu = control$eps_mu, tiny = control$tiny)
}

.gkwq_tmb_params <- function(start, md, spec) {
  parts <- spec$parts
  out <- vector("list", 5L)
  names(out) <- paste0("beta", 1:5)
  for (k in 1:5) {
    out[[k]] <- if (k <= length(parts)) as.numeric(start[[parts[k]]]) else numeric(0)
  }
  out
}

## ---------------------------------------------------------------------------
## Starting values (SPEC 3.6).
##
## gkwreg starts every coefficient at zero, which puts mu at g^{-1}(0) = 0.5
## regardless of tau -- badly wrong at tau = 0.05. And because each anchor is a
## ratio of logs, a mu that starts near 0 or 1 makes the derived parameter
## explode. Hence the quantile-shifted OLS and the feasibility sweep.
## ---------------------------------------------------------------------------

.gkwq_start <- function(yv, md, spec, tau, links, scales, control) {
  if (control$start_method == "user") {
    st <- control$start
    miss <- setdiff(spec$parts, names(st))
    if (length(miss)) {
      stop("`start` is missing parts: ", paste(sQuote(miss), collapse = ", "),
           call. = FALSE)
    }
    return(list(user = st[spec$parts]))
  }

  parts <- spec$parts
  n <- length(yv)

  ## STEP 0 -- marginal parameter estimates for the family.
  th0 <- try(gkwdist::gkwgetstartvalues(yv, family = spec$family), silent = TRUE)
  neutral <- c(alpha = 1, beta = 1, gamma = 1, delta = 0.5, lambda = 1)
  if (inherits(th0, "try-error") || anyNA(th0) || any(!is.finite(th0))) {
    th0 <- neutral[spec$full]
  }
  th0 <- as.list(th0)
  for (nm in spec$full) {
    if (is.null(th0[[nm]]) || !is.finite(th0[[nm]]) || th0[[nm]] <= 0) {
      th0[[nm]] <- neutral[[nm]]
    }
  }

  ## STEP 2 -- quantile-shifted OLS on the link scale. With an intercept-only
  ## design this returns exactly g(quantile(y, tau)), so STEP 1 and STEP 2 agree
  ## by construction rather than by luck.
  X1 <- md$X[[1L]]
  z <- .gkwq_linkfun(yv, links[["mu"]], scales[["mu"]])
  b1 <- tryCatch(qr.solve(X1, z), error = function(e) {
    c(stats::median(z), rep(0, ncol(X1) - 1L))
  })
  if (anyNA(b1)) b1 <- c(stats::median(z), rep(0, ncol(X1) - 1L))
  r <- as.numeric(z - X1 %*% b1)
  b1[1L] <- b1[1L] + stats::quantile(r, probs = tau, names = FALSE, type = 7)

  start <- vector("list", length(parts))
  names(start) <- parts
  start[["mu"]] <- b1

  ## STEP 3 -- nuisance intercepts from the marginal fit, slopes at zero.
  for (p in parts[-1L]) {
    bk <- numeric(ncol(md$X[[p]]))
    bk[1L] <- .gkwq_linkfun(max(th0[[p]], 1e-6), links[[p]], scales[[p]])
    if (!is.finite(bk[1L])) bk[1L] <- 0
    start[[p]] <- bk
  }

  ## STEP 4 -- feasibility sweep: evaluate the anchor in R at the proposed
  ## start. Cheaper to discover an infeasible start here than as a tape full of
  ## NaNs.
  if (!.gkwq_feasible(start, md, spec, tau, links, scales, control)) {
    mu0 <- stats::quantile(yv, probs = tau, names = FALSE, type = 7)
    mu0 <- min(max(mu0, control$eps_mu), 1 - control$eps_mu)
    start[["mu"]] <- c(.gkwq_linkfun(mu0, links[["mu"]], scales[["mu"]]),
                       rep(0, ncol(X1) - 1L))
    for (p in parts[-1L]) {
      bk <- numeric(ncol(md$X[[p]]))
      bk[1L] <- .gkwq_linkfun(neutral[[p]], links[[p]], scales[[p]])
      start[[p]] <- bk
    }
    if (!.gkwq_feasible(start, md, spec, tau, links, scales, control)) {
      stop("could not find a feasible starting value: the anchor is not finite ",
           "and positive anywhere on the tried starts. This usually means the ",
           "response is concentrated at a boundary; try a different anchor ",
           "(the beta and alpha solves fail in opposite tails) or a simpler ",
           "family.", call. = FALSE)
    }
  }

  ## STEP 5 -- intercept-only pre-fit. Worth ~20 ms for the families with three
  ## or more nuisance parameters, where weak identifiability turns a large share
  ## of bad starts into failures; wasted time for kw and ekw.
  ##
  ## But it is a CANDIDATE, never a replacement. These same families have flat
  ## ridges -- kkw nests ekw nests kw -- and an intercept-only pre-fit can run
  ## off along one, landing at alpha ~ 1e-12 and lambda ~ 1e13. Starting the full
  ## fit there is worse than not pre-fitting at all, so the caller evaluates the
  ## real objective at every candidate and keeps the best.
  candidates <- list(base = start)
  refine <- switch(control$start_method,
    auto = spec$family %in% c("kkw", "bkw", "gkw", "mc", "beta"),
    intercept = TRUE,
    FALSE
  )
  if (refine) {
    pre <- tryCatch(.gkwq_prefit(start, yv, md, spec, tau, links, scales,
                                 control),
                    error = function(e) NULL)
    if (!is.null(pre)) candidates$prefit <- pre
  }
  candidates
}

## Evaluate the anchor at a proposed start using the R implementation.
.gkwq_feasible <- function(start, md, spec, tau, links, scales, control) {
  pars <- .gkwq_eval_parts(start, md, spec, links, scales)
  mu <- pmin(pmax(pars$mu, control$eps_mu), 1 - control$eps_mu)
  P <- try(.gkwq_reconstruct(mu, tau, spec, pars$theta), silent = TRUE)
  if (inherits(P, "try-error")) return(FALSE)
  anc <- P[, match(spec$anchor, c("alpha", "beta", "gamma", "delta", "lambda"))]
  all(is.finite(anc)) && all(anc > 0) && all(anc < 1e8) &&
    all(is.finite(.gkwq_logdens(md$y, P[, 1], P[, 2], P[, 3], P[, 4], P[, 5],
                                delta_is_zero = spec$delta_is_zero == 1L)))
}

## Linear predictors -> conditional quantile and nuisance parameters, in R.
.gkwq_eval_parts <- function(coef_list, md, spec, links, scales) {
  eta <- lapply(spec$parts, function(p) {
    as.numeric(md$X[[p]] %*% coef_list[[p]]) + md$offsets[[p]]
  })
  names(eta) <- spec$parts
  mu <- .gkwq_linkinv(eta[["mu"]], links[["mu"]], scales[["mu"]])
  theta <- lapply(spec$parts[-1L], function(p) {
    .gkwq_linkinv(eta[[p]], links[[p]], scales[[p]])
  })
  names(theta) <- spec$parts[-1L]
  list(eta = eta, mu = mu, theta = theta)
}

## Intercept-only pre-fit: same tape, design matrices collapsed to intercepts.
.gkwq_prefit <- function(start, yv, md, spec, tau, links, scales, control) {
  n <- length(yv)
  md0 <- md
  for (p in spec$parts) {
    md0$X[[p]] <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
    md0$offsets[[p]] <- rep(0, n)
  }
  st0 <- lapply(start, function(b) b[1L])
  obj <- TMB::MakeADFun(
    data = .gkwq_tmb_data(yv, md0, spec, tau, links, scales, control),
    parameters = .gkwq_tmb_params(st0, md0, spec),
    DLL = "gkwqreg", silent = TRUE
  )
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(iter.max = 200L, eval.max = 400L))
  k <- 0L
  for (p in spec$parts) {
    k <- k + 1L
    start[[p]][1L] <- opt$par[[k]]
  }
  start
}

## ---------------------------------------------------------------------------
## Fitting.
## ---------------------------------------------------------------------------

.gkwq_fit_one <- function(yv, md, spec, tau, links, scales, control,
                          start = NULL) {
  candidates <- if (is.null(start)) {
    .gkwq_start(yv, md, spec, tau, links, scales, control)
  } else {
    list(warm = start)
  }

  tmb_data <- .gkwq_tmb_data(yv, md, spec, tau, links, scales, control)
  obj <- TMB::MakeADFun(data = tmb_data,
                        parameters = .gkwq_tmb_params(candidates[[1L]], md, spec),
                        DLL = "gkwqreg", silent = control$silent)

  ## Pick the candidate by the objective it actually produces. The tape is built
  ## once and evaluated at each, so this costs a handful of function calls.
  pars <- lapply(candidates, function(st) {
    unlist(.gkwq_tmb_params(st, md, spec), use.names = FALSE)
  })
  vals <- vapply(pars, function(p) {
    v <- tryCatch(as.numeric(obj$fn(p)), error = function(e) Inf)
    if (!is.finite(v)) Inf else v
  }, numeric(1))
  if (all(!is.finite(vals))) {
    stop("no starting value gives a finite log-likelihood. The response may be ",
         "concentrated at a boundary; try a different anchor (the beta and ",
         "alpha solves fail in opposite tails) or a simpler family.",
         call. = FALSE)
  }
  best <- which.min(vals)
  start <- candidates[[best]]

  opt <- .gkwq_optimize(obj, pars[[best]], control)
  rep <- obj$report(opt$par)

  list(obj = obj, opt = opt, report = rep, start = start, tmb_data = tmb_data,
       start_used = names(candidates)[best])
}

.gkwq_optimize <- function(obj, par, control) {
  if (control$method == "nlminb") {
    opt <- stats::nlminb(par, obj$fn, obj$gr,
                         control = list(iter.max = control$maxit,
                                        eval.max = 2L * control$maxit,
                                        rel.tol = control$reltol))
    opt$value <- opt$objective
  } else {
    opt <- stats::optim(par, obj$fn, obj$gr, method = control$method,
                        control = list(maxit = control$maxit,
                                       reltol = control$reltol))
    opt$objective <- opt$value
  }
  opt
}

## ---------------------------------------------------------------------------
## Assembling the fitted object.
## ---------------------------------------------------------------------------

.gkwq_build <- function(fit, cl, formula, md, spec, tau, links, scales, control,
                        model, x, y, yv) {
  obj <- fit$obj
  opt <- fit$opt
  rep <- fit$report

  cf <- opt$par
  nms <- .gkwq_coef_names(md$X, spec$parts)
  names(cf) <- nms

  coef_list <- list()
  i <- 0L
  for (p in spec$parts) {
    k <- ncol(md$X[[p]])
    coef_list[[p]] <- stats::setNames(cf[i + seq_len(k)], colnames(md$X[[p]]))
    i <- i + k
  }

  ## Observed information from optimHess on the AD gradient. Never obj$he()
  ## (that would need second-order AD through the incomplete-beta atomic) and
  ## never the naive J' H J sandwich of the unreparametrized Hessian, which
  ## omits the curvature term and is simply wrong in a regression (SPEC N5).
  vc <- NULL
  se <- rep(NA_real_, length(cf))
  hess <- NULL
  kappa_I <- NA_real_
  if (control$hessian && length(cf)) {
    hess <- tryCatch(stats::optimHess(opt$par, obj$fn, obj$gr),
                     error = function(e) NULL)
    if (!is.null(hess)) {
      hess <- (hess + t(hess)) / 2
      vc <- tryCatch(solve(hess), error = function(e) NULL)
      kappa_I <- tryCatch(kappa(hess, exact = TRUE), error = function(e) NA_real_)
      if (!is.null(vc)) {
        dimnames(vc) <- list(nms, nms)
        d <- diag(vc)
        se <- ifelse(d > 0, sqrt(d), NA_real_)
        names(se) <- nms
        if (anyNA(se)) {
          warning("the information matrix is not positive definite, so ",
                  sum(is.na(se)), " standard error(s) are unavailable. ",
                  "This is the signature of a parameter drifting along a flat ",
                  "ridge; a smaller sub-family usually fixes it.",
                  call. = FALSE)
        }
      }
    }
  }

  n <- length(yv)
  npar <- length(cf)
  ## Take the log-likelihood from the tape re-evaluated AT opt$par, never from
  ## opt$objective. When nlminb exits without clean convergence its reported
  ## objective can belong to a different evaluation than its reported par -- on
  ## a `gkw` fit here the two differed by 3.28 -- and every downstream quantity
  ## (AIC, BIC, the LR test) would then disagree with the coefficients printed
  ## beside it.
  ll <- -as.numeric(rep$nll)

  pv <- data.frame(alpha = as.numeric(rep$alphaVec), beta = as.numeric(rep$betaVec),
                   gamma = as.numeric(rep$gammaVec), delta = as.numeric(rep$deltaVec),
                   lambda = as.numeric(rep$lambdaVec))

  converged <- isTRUE(opt$convergence == 0)
  if (!converged) {
    warning("the optimizer did not report convergence (code ", opt$convergence,
            "). Treat the estimates as provisional.", call. = FALSE)
  }

  out <- list(
    call = cl, formula = md$formula, family = spec$family, tau = tau,
    anchor = spec$anchor, parts = spec$parts, spec = spec,
    link = links, link_scale = scales,
    coefficients = cf, coef_list = coef_list,
    se = se, vcov = vc, hessian = hess, cond_number = kappa_I,
    loglik = ll, npar = npar, nobs = n,
    aic = -2 * ll + 2 * npar, bic = -2 * ll + log(n) * npar,
    fitted.values = as.numeric(rep$muVec),
    linear.predictors = stats::setNames(
      lapply(spec$parts, function(p) {
        as.numeric(md$X[[p]] %*% coef_list[[p]]) + md$offsets[[p]]
      }), spec$parts),
    parameter_vectors = pv,
    loglik_i = as.numeric(rep$loglik_i),
    pinball = as.numeric(rep$pinball),
    weights = md$weights, offsets = md$offsets,
    terms = md$terms, levels = md$levels, contrasts = md$contrasts,
    convergence = opt$convergence, message = opt$message,
    iterations = opt$iterations %||% NA_integer_,
    control = control, obj = obj, start = fit$start
  )
  if (isTRUE(model)) out$model <- md$mf
  if (isTRUE(x)) out$x <- md$X
  if (isTRUE(y)) out$y <- yv

  class(out) <- "gkwqreg"
  out
}

## ---------------------------------------------------------------------------
## Several quantile levels: a container of independent fits.
##
## Separate levels are separate likelihoods, so keeping per-level objects is
## what makes logLik(), vcov() and AIC() honest. This mirrors the rq/rqs split
## in quantreg.
## ---------------------------------------------------------------------------

.gkwq_fit_many <- function(cl, formula, data, taus, family, anchor, link,
                           link_scale, subset, weights, offset, na.action,
                           contrasts, control, model, x, y) {
  spec <- .gkwq_family_info(family, anchor)
  .gkwq_check_identifiability(spec)
  ev <- .gkwq_eval_mf_args(substitute(subset), substitute(weights),
                           substitute(offset), data, parent.frame())
  md <- .gkwq_model_data(formula, data, spec$parts, ev$subset, na.action,
                         ev$weights, ev$offset, contrasts, family, spec$anchor)
  yv <- .gkwq_validate_y(md$y, control$eps_y)
  links <- .gkwq_links(link, spec$parts)
  scales <- .gkwq_link_scales(link_scale, links, spec$parts)

  ## Warm-start outward from the level nearest the median. Neighbouring levels
  ## have nearly the same coefficients, so each sweep must carry its OWN
  ## neighbour forward: running one chain through the whole grid would hand the
  ## highest level the lowest level's estimates, which is worse than no warm
  ## start at all.
  fits <- vector("list", length(taus))
  k0 <- which.min(abs(taus - 0.5))

  fit_at <- function(idx, st) {
    f <- try(.gkwq_fit_one(yv, md, spec, taus[idx], links, scales, control,
                           start = st), silent = TRUE)
    if (inherits(f, "try-error")) {
      if (is.null(st)) stop(attr(f, "condition"))
      f <- .gkwq_fit_one(yv, md, spec, taus[idx], links, scales, control)
    }
    .gkwq_build(f, cl, formula, md, spec, taus[idx], links, scales, control,
                model, x, y, yv)
  }

  fits[[k0]] <- fit_at(k0, NULL)
  seed <- fits[[k0]]$coef_list
  for (sweep in list(rev(seq_len(k0 - 1L)),
                     seq_along(taus)[-seq_len(k0)])) {
    prev <- seed
    for (idx in sweep) {
      fits[[idx]] <- fit_at(idx, if (control$warm_start) prev else NULL)
      prev <- fits[[idx]]$coef_list
    }
  }
  names(fits) <- paste0("tau=", format(taus, trim = TRUE))

  structure(list(fits = fits, taus = taus, family = spec$family,
                 anchor = spec$anchor, parts = spec$parts, call = cl,
                 nobs = length(yv)),
            class = "gkwqregs")
}

## Evaluate the model-frame arguments here rather than letting model.frame look
## them up by symbol. In the formula's environment `weights` resolves to
## stats::weights -- a closure -- and model.frame then fails with "invalid type
## (closure) for variable '(weights)'". Evaluating against `data` first also
## lets `weights = wcol` name a column, as lm() and betareg() allow.
.gkwq_eval_mf_args <- function(sub_expr, w_expr, o_expr, data, env) {
  ev1 <- function(e) {
    if (is.null(e)) return(NULL)
    v <- try(eval(e, data, env), silent = TRUE)
    if (inherits(v, "try-error") || is.function(v)) NULL else v
  }
  list(subset = ev1(sub_expr), weights = ev1(w_expr), offset = ev1(o_expr))
}
