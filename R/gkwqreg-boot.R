## ---------------------------------------------------------------------------
## Simulation from a fit, and the bootstrap.
## ---------------------------------------------------------------------------

#' Simulate responses from a fitted quantile regression
#'
#' Draws new responses from each observation's fitted conditional distribution.
#' Because the anchored parameter is reconstructed from the fitted conditional
#' quantile before anything is drawn, every simulated sample carries the fitted
#' conditional `tau`-quantile by construction, exactly and not approximately.
#' That property is what makes this function the engine behind the parametric
#' bootstrap of [gkwq_boot()] and the simulated envelopes of [plot.gkwqreg()].
#'
#' @param object A `"gkwqreg"` fit.
#' @param nsim Number of replicate samples to draw, returned as columns.
#'   Defaults to 1.
#' @param seed Optional seed for reproducibility, handled as in
#'   [stats::simulate()]: the current value of `.Random.seed` is saved,
#'   `set.seed(seed)` is called, and the original state is restored when the
#'   function exits. A seeded call is therefore reproducible without disturbing
#'   the surrounding random number stream.
#' @param ... Unused; present for consistency with the generic.
#'
#' @details
#' For observation \eqn{i} the fit stores the reconstructed parameter vector
#' \eqn{(\hat\alpha_i, \hat\beta_i, \hat\gamma_i, \hat\delta_i,
#' \hat\lambda_i)} in `object$parameter_vectors`. It is obtained by evaluating
#' each modelled parameter at its own linear predictor, then solving the anchor
#' equation for the eliminated one so that
#'
#' \deqn{Q(\tau; \hat\alpha_i, \hat\beta_i, \hat\gamma_i, \hat\delta_i,
#'       \hat\lambda_i) \;=\; \hat\mu_i,}
#'
#' \eqn{\hat\mu_i} being the fitted conditional `tau`-quantile returned by
#' `fitted(object)`. Responses are then drawn independently across observations
#' from the corresponding Generalized Kumaraswamy laws,
#'
#' \deqn{y_i^{*} \;\sim\; \mathrm{GKw}(\hat\alpha_i, \hat\beta_i,
#'       \hat\gamma_i, \hat\delta_i, \hat\lambda_i),}
#'
#' using [gkwdist::rgkw()]. Since \eqn{\hat\mu_i} *is* the `tau`-quantile of
#' that law, it follows immediately that
#'
#' \deqn{\Pr\left(y_i^{*} \le \hat\mu_i\right) = \tau,
#'       \qquad i = 1, \ldots, n,}
#'
#' whatever the family, the anchor, the covariates or the level. The example
#' below verifies this numerically, both marginally over the whole sample and
#' observation by observation. It is the simulation counterpart of the fact
#' that a single `gkwqreg` fit cannot produce crossing quantiles: the fitted
#' quantiles are those of one proper distribution function, and this function
#' samples from that distribution.
#'
#' @section What the draws do and do not capture:
#' The parameters are held at their estimates and treated as known. The draws
#' therefore reproduce the sampling variation of the *response* around a fixed
#' fitted model; they say nothing about the estimation uncertainty in
#' \eqn{\hat\theta} itself. Propagating that uncertainty into the coefficients
#' is what [gkwq_boot()] does, and it is why a parametric bootstrap refits the
#' model on each simulated sample instead of stopping here.
#'
#' A corollary worth keeping in view: a sample drawn here is, by construction,
#' precisely what the fitted family says the data should look like. Comparing
#' such a sample with the observed data can reveal a poor *fit* within the
#' family, which is the point of a simulated envelope, but it cannot arbitrate
#' between families, since every family would pass its own test.
#'
#' @return
#' A data frame of `nobs(object)` rows and `nsim` columns, named `sim_1`
#' through `sim_nsim`, the rows in the order of the model frame used for
#' fitting. Each column is one independent replicate sample. Row names are the
#' default integer sequence, not those of the original data.
#'
#' @examples
#' set.seed(11)
#' n  <- 400
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.4 + 1.1 * x)
#' y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.9) / log1p(-mu^2))
#' d  <- data.frame(y = y, x = x)
#'
#' fit <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
#'
#' ## One replicate, of the same length as the data.
#' s1 <- simulate(fit, nsim = 1, seed = 1)
#' dim(s1)
#' ## [1] 400   1
#'
#' ## The defining property: the fitted values are the 0.9-quantiles of the
#' ## distributions being sampled, so about 90 per cent of every simulated
#' ## sample falls below them.
#' s <- as.matrix(simulate(fit, nsim = 200, seed = 1))
#' mean(s < fitted(fit))
#' ## [1] 0.8989625
#'
#' ## Observation by observation, not merely on average. Each row of the matrix
#' ## is one observation across the 200 replicates.
#' summary(rowMeans(s < fitted(fit)))
#' ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' ##   0.835   0.885   0.900   0.899   0.915   0.950
#'
#' ## Seeding is reproducible and leaves the ambient stream untouched.
#' identical(simulate(fit, nsim = 2, seed = 1),
#'           simulate(fit, nsim = 2, seed = 1))
#' ## [1] TRUE
#'
#' @seealso [gkwq_boot()] for the bootstrap built on this; [plot.gkwqreg()] for
#'   simulated envelopes; [residuals.gkwqreg()] for the residuals they envelope;
#'   [stats::simulate()] for the generic and its seed convention.
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
#' Refits the model on `R` resampled data sets and summarizes the resulting
#' spread of the coefficients: a bootstrap standard error and a percentile
#' interval for each one, together with the full replicate matrix and its
#' covariance. It is the tool to reach for when the model-based Wald standard
#' errors are in doubt -- a small sample, a weakly identified family, an
#' information matrix that is barely positive definite, or simply a
#' distributional assumption one would rather not lean on.
#'
#' @param object A `"gkwqreg"` fit. It must have been fitted with
#'   `model = TRUE`, which is the default, because the retained model frame is
#'   what gets resampled; otherwise the function stops with an explanatory
#'   message.
#' @param R Number of bootstrap replicates. The default of 200 is adequate for
#'   a standard error, which is a scale estimate and converges quickly.
#'   Percentile intervals estimate the tails of the replicate distribution and
#'   need considerably more, of the order of 1000 or beyond.
#' @param type Resampling scheme, one of `"pairs"` or `"parametric"`.
#'   `"pairs"`, the default, resamples whole observations -- rows of the model
#'   frame, response and covariates together -- with replacement.
#'   `"parametric"` holds the covariates fixed and redraws only the responses,
#'   from the fitted conditional distributions, via [simulate.gkwqreg()].
#' @param seed Optional seed, passed to [base::set.seed()]. Unlike
#'   [simulate.gkwqreg()], the ambient random number stream is **not** restored
#'   afterwards, so a seeded call does change `.Random.seed`.
#' @param ... Unused; present for future extension.
#'
#' @details
#' Each replicate re-evaluates the original call on the resampled data. Two
#' arguments are stripped from that call because both refer to the original
#' rows and would no longer be aligned: `subset`, which has already been
#' applied in building the model frame, and `weights`. **A weighted fit is
#' therefore bootstrapped unweighted**, and the bootstrap standard errors of
#' such a fit describe a different estimator from the one that was fitted. An
#' offset written into the formula with `offset()` lives in the model frame and
#' does follow the resampling; one supplied as a free-standing vector does not.
#'
#' @section Why pairs is the default:
#' The whole bargain of parametric quantile regression is that it buys
#' efficiency, and interpretable coefficients, by assuming a conditional
#' distribution. A parametric bootstrap re-imposes that very assumption at
#' every replicate: the responses it draws are by construction exactly what the
#' fitted family says they should be, since [simulate.gkwqreg()] samples from
#' the fitted distributions themselves. Its intervals therefore describe the
#' sampling behaviour of the estimator *in a world where the family is
#' correct*, and can say nothing whatever about whether that is the world the
#' data came from. If the distributional assumption is the thing in question,
#' the parametric bootstrap is constitutionally unable to answer.
#'
#' The pairs bootstrap makes no such assumption. It resamples observations and
#' so requires only that they be exchangeable draws from the joint distribution
#' of response and covariates. Its spread reflects whatever generated the data,
#' correctly specified or not, which is why it is the honest default. The price
#' is some efficiency when the family really is right, and a heavier tail in
#' small samples, where a resample can happen to contain too few distinct
#' covariate patterns to identify the model -- a factor variable with a rare
#' level is the usual culprit, and shows up as a low `n_ok`.
#'
#' Use `type = "parametric"` when the family is *not* what is being questioned:
#' to check the accuracy of the Wald approximation for a family one is willing
#' to assume, to study the estimator's behaviour in a simulation, or to build
#' envelopes for a diagnostic plot.
#'
#' @section Failed replicates:
#' A replicate whose refit throws an error, or which the optimizer does not
#' certify as converged, contributes a row of `NA` to `replicates` and is
#' excluded from the covariance, the standard errors and the percentile
#' intervals. The count that survived is reported as `n_ok` and printed
#' alongside `R`. That ratio is itself a diagnostic: a fit that loses a
#' noticeable share of its replicates is telling you that the likelihood
#' surface is fragile under mild perturbation of the data, and the surviving
#' replicates are then a selected sample rather than a random one, so the
#' resulting intervals are optimistic. The function stops outright if fewer
#' than two replicates survive.
#'
#' @section What the output means:
#' With \eqn{\hat\theta^{(r)}} the coefficient vector from replicate \eqn{r} and
#' \eqn{B} the number that converged, the reported covariance is the ordinary
#' sample covariance of the surviving replicates,
#'
#' \deqn{\widehat{V} = \frac{1}{B - 1} \sum_{r} \left(\hat\theta^{(r)} -
#'       \bar{\theta}^{*}\right) \left(\hat\theta^{(r)} -
#'       \bar{\theta}^{*}\right)^{\top},}
#'
#' with \eqn{\bar\theta^{*}} the mean over surviving replicates; `se` is the
#' square root of its diagonal. The percentile interval for coefficient \eqn{j}
#' is the pair of empirical quantiles of \eqn{\hat\theta_j^{(r)}} at 0.025 and
#' 0.975, which is the basic percentile method and is not corrected for bias or
#' acceleration.
#'
#' Note that `estimate` is the coefficient vector of the *original* fit, not
#' the mean of the replicates. The difference between the two estimates the
#' bootstrap bias, and a percentile interval that sits noticeably off-centre
#' from `estimate` is a sign that the sampling distribution is skewed and that
#' the symmetric Wald interval is the wrong shape.
#'
#' @return
#' An object of class `"gkwq_boot"`: a list with components
#' \describe{
#'   \item{`replicates`}{an `R` by \eqn{p} matrix of replicate coefficient
#'     vectors, with column names from the fit and `NA` rows where a replicate
#'     failed.}
#'   \item{`vcov`}{the \eqn{p} by \eqn{p} bootstrap covariance matrix
#'     \eqn{\widehat V}, computed over the surviving replicates only.}
#'   \item{`R`}{the number of replicates requested.}
#'   \item{`n_ok`}{the number that converged and entered the summaries.}
#'   \item{`type`}{the resampling scheme used.}
#'   \item{`estimate`}{the coefficients of the original fit.}
#'   \item{`se`}{bootstrap standard errors, the square root of
#'     `diag(vcov)`.}
#'   \item{`percentile`}{a \eqn{p} by 2 matrix of percentile interval endpoints
#'     at 0.025 and 0.975, one row per coefficient.}
#'   \item{`tau`,`family`,`anchor`}{copied from the fit, for the printed
#'     header.}
#' }
#' The `print()` method shows the original estimates beside the bootstrap
#' standard errors and the percentile interval.
#'
#' @examples
#' set.seed(11)
#' n  <- 400
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.4 + 1.1 * x)
#' y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.9) / log1p(-mu^2))
#' d  <- data.frame(y = y, x = x)
#'
#' fit <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
#'
#' \donttest{
#' ## The default: resample observations. No distributional assumption is
#' ## re-imposed, so this is the version that can disagree with the model-based
#' ## standard errors when the family is wrong.
#' bp <- gkwq_boot(fit, R = 200, seed = 1)
#' bp
#' ##                   Estimate Boot SE   2.5%  97.5%
#' ## mu:(Intercept)      0.3972  0.0586 0.2763 0.5127
#' ## mu:x                1.0582  0.0496 0.9694 1.1681
#' ## alpha:(Intercept)   0.6960  0.0416 0.6192 0.7897
#'
#' ## The parametric alternative, for comparison only. Its replicates are drawn
#' ## from the fitted Kumaraswamy laws themselves, so it assumes precisely what
#' ## a specification check would want to test.
#' bq <- gkwq_boot(fit, R = 200, type = "parametric", seed = 1)
#'
#' ## Here the family is in fact correct -- the data were generated from it --
#' ## so all three sets of standard errors agree closely. That agreement is
#' ## itself the evidence: under misspecification it is the pairs column that
#' ## remains trustworthy.
#' round(cbind(wald       = sqrt(diag(vcov(fit))),
#'             pairs      = bp$se,
#'             parametric = bq$se), 4)
#' ##                     wald  pairs parametric
#' ## mu:(Intercept)    0.0603 0.0586     0.0567
#' ## mu:x              0.0505 0.0496     0.0466
#' ## alpha:(Intercept) 0.0470 0.0416     0.0479
#'
#' ## Every replicate converged, so nothing was silently discarded.
#' c(requested = bp$R, converged = bp$n_ok)
#' ## requested converged
#' ##       200       200
#' }
#'
#' @seealso [simulate.gkwqreg()] for the draws behind `type = "parametric"`;
#'   [vcov.gkwqreg()] and [confint.gkwqreg()] for the model-based alternatives;
#'   [gkwq_control()] for the optimizer settings that govern each refit.
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
  ## Prior weights, if any, live in the model frame under "(weights)".
  wcol <- if ("(weights)" %in% names(mf)) as.numeric(mf[["(weights)"]]) else NULL
  ## An `offset` ARGUMENT is resolved from the caller's environment when the
  ## stored call is re-evaluated, so it comes back at its original length and
  ## is added to resampled rows -- the weights defect below, in another guise.
  ## offset() terms written inside the formula need no such care: they sit in
  ## the model frame and are resampled with it. Recovering the argument on its
  ## own therefore means subtracting those terms from the assembled mu offset.
  ocol <- as.numeric(object$offsets[["mu"]]) -
    .gkwq_part_offset(object$terms[["mu"]], mf, n)
  if (all(ocol == 0)) ocol <- NULL
  p <- length(object$coefficients)
  resp <- names(mf)[1L]

  reps <- matrix(NA_real_, R, p,
                 dimnames = list(NULL, names(object$coefficients)))
  sims <- if (type == "parametric") simulate.gkwqreg(object, nsim = R) else NULL

  for (r in seq_len(R)) {
    idx <- if (type == "pairs") sample.int(n, n, replace = TRUE) else seq_len(n)
    d <- if (type == "pairs") mf[idx, , drop = FALSE]
         else { z <- mf; z[[resp]] <- sims[[r]]; z }
    ## Carry the prior weights into the refit. The model frame stores them in a
    ## column called "(weights)", which is not a usable variable name, so they
    ## are copied to a plain column and named in the call. Dropping them --
    ## which this did before version 0.1.0 -- silently bootstraps a WEIGHTED fit
    ## as though it were unweighted, and the replicates then centre on the wrong
    ## estimate.
    cl_r <- cl
    cl_r$data <- quote(d)
    cl_r$subset <- NULL
    cl_r$weights <- NULL
    if (!is.null(wcol)) {
      d[[".gkwq_w"]] <- wcol[idx]
      cl_r$weights <- as.name(".gkwq_w")
    }
    cl_r$offset <- NULL
    if (!is.null(ocol)) {
      d[[".gkwq_o"]] <- ocol[idx]
      cl_r$offset <- as.name(".gkwq_o")
    }
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
#' A thin wrapper on [anova.gkwqreg()], offered under the name that users of
#' `lmtest::lrtest()` will reach for first. Behaviour, guards and return value
#' are identical: the fits are sorted by dimension, each is tested against the
#' one above it by the statistic
#' \eqn{LR = 2\{\ell_1(\hat\theta_1) - \ell_0(\hat\theta_0)\}} referred to a
#' chi-squared distribution on the difference in dimension, and comparisons
#' across quantile levels or across anchors are refused.
#'
#' @param object A `"gkwqreg"` fit.
#' @param ... One or more further `"gkwqreg"` fits, at the same quantile level,
#'   with the same anchor and the same number of observations.
#'
#' @details
#' `lrtest()` is a generic defined by this package, so that `gkwqreg` need not
#' depend on `lmtest` merely to offer a familiar name. It is a *different*
#' generic from `lmtest::lrtest()`, and the two mask one another when both
#' packages are attached: whichever was attached last wins.
#'
#' The distinction matters, because the two behave differently on these
#' objects. If `lmtest` is attached after `gkwqreg`, a bare call to `lrtest()`
#' reaches `lmtest`'s default method, which knows nothing about quantile levels
#' or anchors. Handed two fits at different levels it will happily compute a
#' difference of log-likelihoods and report a significant result -- the very
#' number that [anova.gkwqreg()] refuses to produce, and for the reasons given
#' there. Call `gkwqreg::lrtest()` or [anova.gkwqreg()] explicitly if there is any doubt
#' about which generic is in scope; both carry the guards.
#'
#' @return
#' A data frame of class `"anova.gkwqreg"`, exactly as returned by
#' [anova.gkwqreg()]; see that help page for the columns and the printed
#' layout.
#'
#' @examples
#' ## Data with a genuine third shape parameter: lambda = 4, so "kw" is too
#' ## small and "ekw" is right.
#' set.seed(2026)
#' n  <- 500
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.3 + 1.2 * x)
#' b  <- log1p(-0.5^(1 / 4)) / log1p(-mu^2)
#' y  <- gkwdist::rekw(n, alpha = 2, beta = b, lambda = 4)
#' d  <- data.frame(y = y, x = x)
#'
#' m_kw  <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
#' m_ekw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "ekw")
#'
#' lrtest(m_kw, m_ekw)
#' ##      Df logLik     AIC   Chisq Chi Df Pr(>Chisq)
#' ## [1,]  3 482.58 -959.16
#' ## [2,]  4 496.90 -985.80  28.641      1  8.713e-08
#'
#' ## Identical to anova(), which is what it calls.
#' identical(lrtest(m_kw, m_ekw), anova(m_kw, m_ekw))
#' ## [1] TRUE
#'
#' ## The guards travel with it: a cross-level comparison is refused here too.
#' m_tau9 <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
#' try(lrtest(m_kw, m_tau9))
#'
#' @seealso [anova.gkwqreg()] for the full account of the test, the nesting of
#'   the families and the two refusals; [vuong_test()] for the non-nested case.
#' @export
lrtest <- function(object, ...) UseMethod("lrtest")

#' @rdname lrtest
#' @export
lrtest.gkwqreg <- function(object, ...) anova(object, ...)
