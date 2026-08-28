## ---------------------------------------------------------------------------
## Model comparison.
##
## One guard runs through all of it: models compared by a likelihood-ratio test
## must share the same quantile level AND the same anchor. Different levels are
## different likelihoods; different anchors are non-nested models of equal
## dimension, so an LR test on either is meaningless even though
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
#' Compares two or more fitted quantile regressions that stand in a nesting
#' relation, by the classical likelihood-ratio test, and tabulates the result in
#' the idiom of [stats::anova()]. Because the seven Generalized Kumaraswamy
#' families are genuine parameter restrictions of one another, choosing among
#' them is here an ordinary hypothesis test. It does not require the Vuong
#' machinery that a collection of mutually non-nested families would force upon
#' it, and which packages built on such a collection have no way to avoid.
#'
#' @param object A `"gkwqreg"` fit.
#' @param ... One or more further `"gkwqreg"` fits. Every fit must use the same
#'   quantile level `tau`, the same anchor, and the same number of observations;
#'   the section *What this method refuses* explains why each requirement is
#'   enforced rather than merely advised. The order of the arguments does not
#'   matter: the fits are sorted internally by number of estimated coefficients.
#'
#' @details
#' Write \eqn{\ell(\hat\theta)} for the maximized log-likelihood of a fit and
#' \eqn{p} for its number of estimated coefficients, reported in the `Df`
#' column. The anchored parameter is not counted in \eqn{p}: it is computed from
#' the conditional quantile by the anchor solve, not estimated.
#'
#' For a model \eqn{M_0} nested inside a model \eqn{M_1}, with
#' \eqn{p_0 < p_1}, the statistic reported in the `Chisq` column is
#'
#' \deqn{LR = 2 \left\{ \ell_1(\hat\theta_1) - \ell_0(\hat\theta_0) \right\},}
#'
#' whose null distribution, under \eqn{M_0} and the usual regularity
#' conditions, is asymptotically chi-squared on \eqn{p_1 - p_0} degrees of
#' freedom. The columns `Chi Df` and `Pr(>Chisq)` hold \eqn{p_1 - p_0} and the
#' upper-tail probability \eqn{\Pr(\chi^2_{p_1 - p_0} > LR)}. The `AIC` column
#' reports \eqn{-2\ell(\hat\theta) + 2p} for each row; unlike the test itself,
#' it is comparable down the whole table and requires no nesting.
#'
#' With more than two fits the table is **sequential**: after sorting by
#' dimension, each row is tested against the row immediately above it, never
#' against the first row. Nothing verifies that consecutive rows really are
#' nested; that remains the analyst's judgement, and the map in the next section
#' is the guide.
#'
#' @section The families are genuinely nested:
#' Every family is the five-parameter `gkw` distribution with some of its
#' parameters held at neutral values, so each is a restriction of the ones
#' above it in the lattice.
#'
#' | family | free parameters | held fixed |
#' |:--|:--|:--|
#' | `"gkw"` | `alpha`, `beta`, `gamma`, `delta`, `lambda` | none |
#' | `"bkw"` | `alpha`, `beta`, `gamma`, `delta` | \eqn{\lambda = 1} |
#' | `"kkw"` | `alpha`, `beta`, `delta`, `lambda` | \eqn{\gamma = 1} |
#' | `"ekw"` | `alpha`, `beta`, `lambda` | \eqn{\gamma = 1}, \eqn{\delta = 0} |
#' | `"kw"` | `alpha`, `beta` | \eqn{\gamma = 1}, \eqn{\delta = 0}, \eqn{\lambda = 1} |
#' | `"mc"` | `gamma`, `delta`, `lambda` | \eqn{\alpha = \beta = 1} |
#' | `"beta"` | `gamma`, `delta` | \eqn{\alpha = \beta = \lambda = 1} |
#'
#' The restrictions compose into chains, and it is along a chain that this
#' method is valid:
#'
#' 1. `"kw"` is `"ekw"` with \eqn{\lambda = 1}; `"ekw"` is `"kkw"` with
#'    \eqn{\delta = 0}; `"kkw"` is `"gkw"` with \eqn{\gamma = 1}.
#' 2. `"kw"` is also `"bkw"` with \eqn{\gamma = 1} and \eqn{\delta = 0}, and
#'    `"bkw"` is `"gkw"` with \eqn{\lambda = 1}.
#' 3. `"beta"` is `"mc"` with \eqn{\lambda = 1}, and `"mc"` is `"gkw"` with
#'    \eqn{\alpha = \beta = 1}; `"beta"` is likewise `"bkw"` with
#'    \eqn{\alpha = \beta = 1}.
#'
#' The reparametrization does not disturb any of this. The anchor solve
#' eliminates the same parameter from both members of a pair, so a restriction
#' on the original family is a restriction on the anchored one and \eqn{LR} is
#' the same quantity it always was.
#'
#' Two consequences are worth stating. First, a larger family cannot attain a
#' smaller maximum log-likelihood than a family nested inside it, so the
#' `logLik` column is non-decreasing down the table as a matter of mathematics,
#' not of luck; the package's test suite checks that invariant. A `logLik` that
#' *falls* as `Df` rises is therefore not evidence about the data but evidence
#' that one of the fits has not reached its maximum, and the offending fit
#' should be refitted or discarded. Second, not every pair of families is
#' nested: `"kw"` and `"mc"`, for instance, restrict disjoint sets of
#' parameters and neither contains the other. Such a pair belongs to
#' [vuong_test()] or [compare_families()], not here.
#'
#' A practical corollary of the anchor guard is that `"mc"` admits only the
#' `lambda` anchor and `"beta"` only the `gamma` anchor, whereas every
#' Kumaraswamy member defaults to `beta`. Those two families therefore cannot
#' reach this method in company with a Kumaraswamy member at all, even where
#' the nesting exists on paper.
#'
#' @section A caution at the boundary:
#' The restriction \eqn{\delta = 0} -- the step from `"kkw"` down to `"ekw"`, or
#' from `"bkw"` down to `"kw"` -- lies on the boundary of the parameter space:
#' under the default log link \eqn{\delta} is strictly positive and reaches zero
#' only in the limit. For such a restriction the null distribution of \eqn{LR}
#' is not chi-squared on one degree of freedom but a mixture with an atom at
#' zero (Chernoff 1954; Self and Liang 1987), so the p-value reported here is
#' conservative and errs towards keeping the smaller model. Its signature in
#' practice is an \eqn{LR} of essentially zero with a p-value close to one, as
#' in the third row of the example below. Restrictions setting \eqn{\gamma} or
#' \eqn{\lambda} to 1 are interior and need no such caveat.
#'
#' @section What this method refuses, and why:
#' Two comparisons raise an error instead of returning a number. Both refusals
#' are deliberate, and in both cases the arithmetic would have produced a
#' perfectly presentable statistic that means nothing.
#'
#' **Fits at different quantile levels.** The level `tau` indexes the question
#' being asked, not the model: it never enters the parameter vector, and the
#' profile likelihood in `tau` is exactly flat. Each level therefore defines its
#' own likelihood, and two likelihoods for two different questions are not
#' comparable at all -- neither nested, nor non-nested, but incommensurable. A
#' difference of log-likelihoods between a `tau = 0.5` fit and a `tau = 0.9` fit
#' is a difference between the fit of one model to one target and another model
#' to another target. Note that a generic likelihood-ratio routine written for
#' other model classes has no way to know this and will happily report a
#' significant result; see [lrtest()].
#'
#' **Fits with different anchors.** Two anchors on the same data give
#' *non-nested models of equal dimension*. Neither is a restriction of the
#' other, so no null hypothesis relates them; and \eqn{p_1 - p_0 = 0} leaves no
#' degrees of freedom against which to refer a statistic. The difference in
#' log-likelihood between them can nonetheless be substantial, because the
#' anchor is a modeling choice and not a relabelling whenever the nuisance
#' parameters are unsaturated (see [gkwqreg()]). The correct instrument is
#' [vuong_test()], which is built for models of equal dimension that do not
#' nest; AIC and BIC are also legitimate.
#'
#' A third case is warned about rather than refused. If two fits of *equal*
#' dimension reach the test by some other route -- same family, same anchor,
#' same level, but different covariates of the same count -- a warning is issued
#' and the corresponding `Chisq` and `Pr(>Chisq)` entries are `NA`.
#'
#' @return
#' A data frame of class `"anova.gkwqreg"`, inheriting from `"anova"` and
#' `"data.frame"`, with one row per fit sorted by increasing `Df` and a
#' `"heading"` attribute recording the common `tau` and anchor. Its columns are
#' \describe{
#'   \item{`family`}{the family of the fit in that row.}
#'   \item{`Df`}{\eqn{p}, the number of estimated coefficients.}
#'   \item{`logLik`}{\eqn{\ell(\hat\theta)}, the maximized log-likelihood.}
#'   \item{`AIC`}{\eqn{-2\ell(\hat\theta) + 2p}.}
#'   \item{`Chisq`}{the likelihood-ratio statistic against the row above; `NA`
#'     in the first row.}
#'   \item{`Chi Df`}{the difference in `Df` against the row above; `NA` in the
#'     first row.}
#'   \item{`Pr(>Chisq)`}{the upper-tail chi-squared probability; `NA` in the
#'     first row and wherever `Chi Df` is not positive.}
#' }
#' The `print()` method formats the table through [stats::printCoefmat()] and
#' prepends the heading.
#'
#' @examples
#' ## A response generated from an EKW distribution: the true lambda is 4, so
#' ## the two-parameter Kumaraswamy family is genuinely wrong, and its larger
#' ## relatives are genuinely right.
#' set.seed(2026)
#' n  <- 500
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.3 + 1.2 * x)          # the true conditional median
#' b  <- log1p(-0.5^(1 / 4)) / log1p(-mu^2)   # the beta anchor solve at tau = 0.5
#' y  <- gkwdist::rekw(n, alpha = 2, beta = b, lambda = 4)
#' d  <- data.frame(y = y, x = x)
#'
#' m_kw  <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
#' m_ekw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "ekw")
#' m_kkw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kkw")
#'
#' ## The chain kw inside ekw inside kkw, tested one step at a time.
#' anova(m_kw, m_ekw, m_kkw)
#' ##      Df logLik     AIC   Chisq Chi Df Pr(>Chisq)
#' ## [1,]  3 482.58 -959.16
#' ## [2,]  4 496.90 -985.80 28.6406      1  8.713e-08
#' ## [3,]  5 496.91 -983.81  0.0093      1      0.923
#' ##
#' ## Row 2 tests lambda = 1, that is kw inside ekw, and rejects it decisively:
#' ## the extra shape parameter is real, exactly as the simulation was built.
#' ## Row 3 tests delta = 0, that is ekw inside kkw, and finds nothing -- note
#' ## the statistic of essentially zero, the boundary signature discussed above.
#' ## Note also that logLik never falls as Df rises, as the nesting requires.
#'
#' ## The refusals. Neither of these returns a number.
#' m_tau9 <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
#' try(anova(m_kw, m_tau9))       # different levels: separate likelihoods
#'
#' m_alpha <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw",
#'                    anchor = "alpha")
#' try(anova(m_kw, m_alpha))      # different anchors: non-nested, equal dimension
#'
#' ## The comparison anova() refused is exactly the one vuong_test() makes.
#' vuong_test(m_kw, m_alpha)
#' ##   z = 7.642, p-value = 2.137e-14 -- model 1, the beta anchor, is favoured.
#'
#' @seealso [lrtest()] for the same test under the name `lmtest` users expect;
#'   [vuong_test()] for non-nested comparisons, including two anchors;
#'   [compare_families()] to sweep every family at once; [gkwqreg()] for what
#'   the anchor is and why it is a modeling argument.
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
#' Tests whether two models fitted to the same responses, at the same quantile
#' level, are equally close to the true data-generating distribution. This is
#' the right instrument for the two comparisons a likelihood-ratio test cannot
#' make: two anchors of the same family, and two families that do not nest.
#' Both produce models of equal dimension fitted to the same data, for which no
#' null hypothesis of restriction exists and no degrees of freedom remain.
#'
#' @param object,object2 Two `"gkwqreg"` fits at the same quantile level and
#'   over the same observations. The order matters only through the sign of the
#'   statistic: `object` is model 1 and `object2` is model 2. The two may differ
#'   in family, in anchor, in the covariates entering any part, or in all three.
#' @param correction Logical. If `TRUE`, the default, apply the BIC-style
#'   dimension correction described below, which penalizes whichever model
#'   estimates more coefficients. If `FALSE`, compare the raw per-observation
#'   log-likelihood contributions. The correction is exactly zero when the two
#'   fits have the same number of coefficients, which is the case for a
#'   comparison of two anchors.
#'
#' @details
#' Let \eqn{f_1} and \eqn{f_2} denote the two fitted conditional densities and
#' \eqn{h} the unknown true one. Vuong's null hypothesis is that the two models
#' are *equally close to the truth* in Kullback-Leibler divergence,
#'
#' \deqn{H_0 : \quad E\left[ \log \frac{f_1(Y_i \mid x_i; \theta_1^{*})}
#'                                {f_2(Y_i \mid x_i; \theta_2^{*})} \right] = 0,}
#'
#' the expectation taken under \eqn{h} and \eqn{\theta_k^{*}} being the
#' pseudo-true parameter of model \eqn{k}. The alternatives are one-sided in
#' either direction: \eqn{H_1}, that model 1 is closer to the truth, and
#' \eqn{H_2}, that model 2 is. Crucially, **neither model is assumed correct**.
#' The test is about relative distance to the truth, so rejecting \eqn{H_0} says
#' that one model is nearer, never that it is right.
#'
#' @section The statistic:
#' Write \eqn{\ell_{ki}(\hat\theta_k)} for the log-likelihood contribution of
#' observation \eqn{i} under model \eqn{k}, stored in `object$loglik_i`, and let
#'
#' \deqn{m_i = \ell_{1i}(\hat\theta_1) - \ell_{2i}(\hat\theta_2), \qquad
#'       i = 1, \ldots, n.}
#'
#' With `correction = TRUE` these differences are shifted by the per-observation
#' share of the dimension penalty,
#'
#' \deqn{\tilde m_i = m_i - \frac{(p_1 - p_2)\log n}{2n},}
#'
#' where \eqn{p_k} is the number of coefficients estimated by model \eqn{k}; with
#' `correction = FALSE`, \eqn{\tilde m_i = m_i}. The reported statistic is the
#' studentized mean
#'
#' \deqn{Z = \frac{\sqrt{n}\, \bar{\tilde m}}{s}, \qquad
#'       s^2 = \frac{1}{n - 1} \sum_{i=1}^{n} (\tilde m_i - \bar{\tilde m})^2,}
#'
#' with \eqn{s} the ordinary sample standard deviation, as computed by
#' [stats::sd()]. Under \eqn{H_0}, and for strictly non-nested models,
#' \eqn{Z} converges in distribution to a standard normal, so the two-sided
#' p-value is \eqn{2\,\Phi(-|Z|)}.
#'
#' The denominator is what distinguishes this from a naive comparison of
#' log-likelihoods. A difference in total log-likelihood says nothing on its own
#' about sampling variability; \eqn{s} measures how consistently one model beats
#' the other across observations, so a small total advantage accumulated
#' steadily over every observation can be decisive, while a large advantage
#' driven by a handful of points is not.
#'
#' @section Reading the sign:
#' The sign of \eqn{Z} points at the favoured model, and the p-value says
#' whether the data support the preference at all.
#'
#' * \eqn{Z > 0} with a small p-value: model 1, that is `object`, is favoured.
#' * \eqn{Z < 0} with a small p-value: model 2, that is `object2`, is favoured.
#' * a large p-value: the data do not discriminate between the two. This is a
#'   genuine third conclusion and not a failure to establish a preference: it
#'   states that both models sit at the same Kullback-Leibler distance from the
#'   truth as far as these data can tell. For an anchor comparison it is the
#'   expected outcome once the nuisance parameters are saturated, and the
#'   example below shows exactly that.
#'
#' The `print()` method states the verdict at the conventional five per cent
#' level.
#'
#' @section The dimension correction:
#' The uncorrected statistic compares maximized log-likelihoods, and therefore
#' rewards the more flexible model for its flexibility alone. The correction
#' subtracts half the difference of BIC penalties, spread evenly over the
#' observations. Summing the corrected differences shows what it buys:
#'
#' \deqn{\sum_{i=1}^{n} \tilde m_i \;=\; \ell_1 - \ell_2
#'       - \frac{(p_1 - p_2)\log n}{2}
#'       \;=\; \frac{\mathrm{BIC}_2 - \mathrm{BIC}_1}{2},}
#'
#' since \eqn{\mathrm{BIC}_k = -2\ell_k + p_k \log n}. The corrected numerator is
#' thus one half of the BIC difference, oriented so that a positive value favours
#' model 1, and the test is a formal significance statement about a BIC
#' comparison rather than an informal reading of one. When \eqn{p_1 = p_2} the
#' correction term vanishes identically and the two settings of `correction`
#' agree to the last digit.
#'
#' @section When this test does not apply:
#' The standard normal reference requires the two models to be *strictly*
#' non-nested. For a genuinely nested pair, one family restricting another, the
#' variance \eqn{s^2} tends to zero under the null and the limit of \eqn{Z} is
#' not normal; the correct reference there is chi-squared, which is what
#' [anova.gkwqreg()] uses. Use `anova()` for a nested pair and this function for
#' everything else.
#'
#' The two fits must also share the quantile level and the number of
#' observations, both of which are checked. Fits at different levels answer
#' different questions and their likelihoods are not commensurable; fits on
#' different samples give differences \eqn{m_i} that are not paired.
#'
#' @return
#' An object of class `"gkwq_vuong"`: a list with components
#' \describe{
#'   \item{`statistic`}{the value of \eqn{Z}.}
#'   \item{`p.value`}{the two-sided p-value \eqn{2\,\Phi(-|Z|)}.}
#'   \item{`n`}{the number of observations, common to both fits.}
#'   \item{`model1`,`model2`}{character labels of the form
#'     `"family / anchor a"` identifying the two fits.}
#'   \item{`tau`}{the common quantile level.}
#'   \item{`correction`}{the value of `correction` that was used.}
#' }
#' The `print()` method reports the two labels, the statistic, the p-value, and
#' the verdict at the five per cent level.
#'
#' @examples
#' ## Data generated under the "beta" anchor of the Kumaraswamy family: alpha is
#' ## constant at 2 and beta is the quantity that varies with the covariate.
#' set.seed(11)
#' n  <- 400
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.4 + 1.1 * x)                  # true conditional 0.9-quantile
#' y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.9) / log1p(-mu^2))
#' d  <- data.frame(y = y, x = x)
#'
#' ## Two anchors, same family, same level, same data: three coefficients each,
#' ## and neither nested in the other. anova() refuses them; this does not.
#' f_beta  <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
#' f_alpha <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw",
#'                    anchor = "alpha")
#' vuong_test(f_beta, f_alpha)
#' ##   z = 6.316, p-value = 2.679e-10
#' ##   model 1 is favoured
#' ## z is large and positive, so model 1 -- the beta anchor, the one the data
#' ## were generated under -- is favoured decisively.
#'
#' ## Both fits estimate three coefficients, so the dimension correction is
#' ## identically zero and switching it off changes nothing.
#' identical(vuong_test(f_beta, f_alpha)$statistic,
#'           vuong_test(f_beta, f_alpha, correction = FALSE)$statistic)
#' ## [1] TRUE
#'
#' ## Saturating the nuisance parameter -- regressing alpha, respectively beta,
#' ## on the same covariate as the quantile -- makes the anchor choice almost
#' ## immaterial, which is the recommended practice when theory does not
#' ## dictate an anchor.
#' s_beta  <- gkwqreg(y ~ x | x, data = d, tau = 0.9, family = "kw")
#' s_alpha <- gkwqreg(y ~ x | x, data = d, tau = 0.9, family = "kw",
#'                    anchor = "alpha")
#' vuong_test(s_beta, s_alpha)
#' ##   z = 1.62, p-value = 0.1052
#' ##   neither model is favoured
#' ## The evidence has evaporated. Saturation has turned two different models
#' ## into two parametrizations of very nearly the same one.
#'
#' ## Non-nested families work the same way. "kw" restricts alpha and beta while
#' ## "beta" restricts them away entirely, so neither contains the other.
#' f_bet <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "beta")
#' vuong_test(f_beta, f_bet)
#'
#' @seealso [anova.gkwqreg()] for the nested case; [compare_families()] for a
#'   sweep over all seven families by AIC, BIC and check loss; [gkwqreg()] for
#'   why two anchors are two models.
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

#' Compare the seven families at one quantile level
#'
#' Refits the same model under each family in turn and tabulates the results
#' side by side, reporting AIC, BIC and check loss together. The three criteria
#' answer different questions and can genuinely disagree; the sections below
#' explain when that happens and which number to act on.
#'
#' @param object A fitted `"gkwqreg"` model whose call is reused. Everything is
#'   held fixed except the family: the same formula, data, quantile level,
#'   links, weights, offsets and control settings are used for every refit.
#' @param families Character vector of families to try, any subset of
#'   `c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")`. Defaults to all seven.
#'   Restricting it is often sensible: `"gkw"` is weakly identified in every
#'   parametrization and warns accordingly, and a family that cannot represent
#'   the data at all merely costs time.
#' @param ... Currently unused; present for future extension. Arguments given
#'   here do not reach [gkwqreg()]. To vary anything other than the family,
#'   change it in `object` and call this function again.
#'
#' @details
#' The stored call of `object` is re-evaluated once per family, with `family`
#' replaced and `anchor` deleted, in the environment from which
#' `compare_families()` was itself called. Two practical consequences follow.
#' The data and any other objects named in the call must be visible from that
#' environment, exactly as they were when `object` was fitted. And because the
#' anchor is dropped rather than carried over, **each family is fitted under its
#' own default anchor**: `beta` for the five Kumaraswamy members, `lambda` for
#' `"mc"`, `gamma` for `"beta"`. That is not a lapse. Not every anchor is
#' admissible in every family -- `"mc"` and `"beta"` each admit exactly one --
#' so insisting on a common anchor would simply fail. It does mean the rows can
#' differ in more than the family, which is why `anchor` is a column of the
#' result and worth reading.
#'
#' A family that fails to fit is not silently dropped. It contributes a row of
#' `NA` values with `converged = FALSE`, so an absent family is visible as an
#' absence rather than as nothing at all. The `converged` column carries the
#' optimizer's own verdict, and a row reporting `FALSE` should be excluded from
#' the ranking rather than interpreted: an unconverged fit's log-likelihood is
#' a lower bound on its maximum, so its AIC is not comparable with the rest.
#'
#' @section The three criteria:
#' With \eqn{\ell(\hat\theta)} the maximized log-likelihood, \eqn{p} the number
#' of estimated coefficients and \eqn{n} the number of observations,
#'
#' \deqn{\mathrm{AIC} = -2\ell(\hat\theta) + 2p, \qquad
#'       \mathrm{BIC} = -2\ell(\hat\theta) + p \log n.}
#'
#' Both are computed from the whole conditional density and differ only in how
#' hard they penalize dimension: \eqn{\log n} exceeds 2 for any \eqn{n \ge 8},
#' so BIC prefers smaller models and is consistent for the true model when it
#' is among the candidates, whereas AIC targets predictive accuracy of the
#' density and tends to over-select.
#'
#' The `pinball` column is the mean check loss at the fitted quantiles
#' \eqn{\hat Q_i(\tau)},
#'
#' \deqn{L_\tau = \frac{\sum_{i=1}^n w_i \,(y_i - \hat Q_i(\tau))
#'        \left\{ \tau - \mathbf{1}(y_i < \hat Q_i(\tau)) \right\}}
#'        {\sum_{i=1}^n w_i},}
#'
#' with \eqn{w_i} the prior weights. The summand is
#' \eqn{\tau\,(y_i - \hat Q_i)} when the observation lies above the fitted
#' quantile and \eqn{(1 - \tau)(\hat Q_i - y_i)} when it lies below, so the
#' loss is asymmetric in exactly the proportion \eqn{\tau : 1 - \tau} that
#' defines the quantile. Its population minimizer is the true conditional
#' \eqn{\tau}-quantile and nothing else; see [pinball()].
#'
#' @section Why AIC and check loss can disagree:
#' They measure different things, and the disagreement is not a paradox but the
#' point. AIC and BIC score the entire conditional density: how well the model
#' describes the response at every point of its support, the bulk of the
#' distribution included. Check loss scores one feature of that density, the
#' conditional \eqn{\tau}-quantile, and is indifferent to everything else. A
#' family can capture the bulk of the distribution well, which is what earns it
#' the likelihood, and still place the \eqn{\tau}-quantile worse than a rival
#' that fits the bulk less well -- especially at an extreme level, where the
#' quantile is determined by a region of the support that contributes little to
#' the likelihood.
#'
#' This is not hypothetical. In the case study shipped with the package,
#' `vignette("gkwqreg-case-study")`, AIC ranks the `"beta"` family first while
#' out-of-sample check loss ranks it last, behind every Kumaraswamy member. The
#' example below reproduces a milder version of the same reversal on simulated
#' data.
#'
#' The resolution is to let the question decide. If the object of inference is
#' the conditional quantile, the check loss is the criterion that matches the
#' goal, and AIC is answering a question that was not asked.
#'
#' @section Compute the check loss out of sample:
#' One qualification is essential. The `pinball` column returned here is
#' computed **in sample**, on the same observations that were used to fit, and
#' an in-sample loss rewards flexibility for its own sake: a larger family can
#' always drive it down a little by tracking noise. It is reported because it is
#' free and because it exposes gross failures, not because it settles anything.
#'
#' A check loss used to *choose* a family must be evaluated on data the fit has
#' not seen -- a holdout sample or a cross-validation, both shown in the example
#' below -- via [pinball()] with its `newdata` argument. The out-of-sample check
#' loss is the number to act on for a quantile question.
#'
#' @return
#' A data frame with one row per entry of `families`, ordered by increasing
#' `AIC` with failures last, and with row names reset. Its columns are
#' \describe{
#'   \item{`family`}{the family fitted in that row.}
#'   \item{`anchor`}{the anchor that family defaulted to; `NA` on failure.}
#'   \item{`df`}{the number of estimated coefficients, \eqn{p}.}
#'   \item{`logLik`}{the maximized log-likelihood, \eqn{\ell(\hat\theta)}.}
#'   \item{`AIC`}{\eqn{-2\ell(\hat\theta) + 2p}.}
#'   \item{`BIC`}{\eqn{-2\ell(\hat\theta) + p\log n}.}
#'   \item{`pinball`}{the in-sample mean check loss \eqn{L_\tau}.}
#'   \item{`converged`}{`TRUE` when the optimizer reported convergence.}
#' }
#'
#' @examples
#' ## Data from an EKW distribution, split into a training and a holdout part.
#' set.seed(2026)
#' n  <- 500
#' x  <- runif(n, -2, 2)
#' mu <- plogis(0.3 + 1.2 * x)
#' b  <- log1p(-0.9^(1 / 4)) / log1p(-mu^2)      # beta anchor solve at tau = 0.9
#' y  <- gkwdist::rekw(n, alpha = 2, beta = b, lambda = 4)
#' d  <- data.frame(y = y, x = x)
#' train <- d[1:350, ]
#' test  <- d[351:500, ]
#'
#' fit <- gkwqreg(y ~ x, data = train, tau = 0.9, family = "kw")
#' cmp <- compare_families(fit, families = c("kw", "ekw", "kkw", "beta"))
#' cmp
#' ##   family anchor df   logLik       AIC       BIC    pinball converged
#' ## 1    ekw   beta  4 391.9588 -775.9176 -760.4858 0.01504029      TRUE
#' ## 2    kkw   beta  5 392.7234 -775.4469 -756.1572 0.01500626      TRUE
#' ## 3     kw   beta  3 383.0529 -760.1058 -748.5320 0.01507971      TRUE
#' ## 4   beta  gamma  3 297.2013 -588.4027 -576.8289 0.01994690      TRUE
#'
#' ## Now recompute the check loss on the holdout, which is the criterion that
#' ## matches a quantile question.
#' holdout <- vapply(cmp$family, function(fm)
#'   pinball(gkwqreg(y ~ x, data = train, tau = 0.9, family = fm),
#'           newdata = test),
#'   numeric(1))
#' cbind(cmp[, c("family", "AIC", "pinball")], holdout = round(holdout, 5))
#' ##      family       AIC    pinball holdout
#' ## ekw     ekw -775.9176 0.01504029 0.01859
#' ## kkw     kkw -775.4469 0.01500626 0.01871
#' ## kw       kw -760.1058 0.01507971 0.01854
#' ## beta   beta -588.4027 0.01994690 0.02362
#' ##
#' ## The orderings differ. AIC puts "kw" third and the in-sample check loss puts
#' ## it third as well, both rewarding the larger families -- which is what extra
#' ## parameters are for. On data the fits have not seen, the smallest family
#' ## predicts the 0.9-quantile best. Only "beta" is last by every criterion,
#' ## and it is the one family here that is genuinely excluded.
#'
#' @seealso [pinball()] for the check loss, in or out of sample;
#'   [anova.gkwqreg()] for a formal test along a nesting chain;
#'   [vuong_test()] for a formal test between two non-nested fits;
#'   [gkwqreg()] for the family and anchor arguments themselves.
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
  extra <- list(...)
  rows <- lapply(families, function(fm) {
    cl <- object$call
    cl$family <- fm
    cl$anchor <- NULL          # each family gets its own default anchor
    ## Anything in `...` overrides the original call, so a single
    ## compare_families(fit, tau = 0.9) refits every family at that level.
    for (nm in names(extra)) cl[[nm]] <- extra[[nm]]
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
