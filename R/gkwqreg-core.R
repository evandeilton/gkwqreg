#' Parametric quantile regression for the Generalized Kumaraswamy family
#'
#' @description
#' Fits a regression model for the conditional \eqn{\tau}-quantile of a response
#' confined to the open unit interval \eqn{(0,1)}, at a level \eqn{\tau} that the
#' analyst fixes in advance. The Generalized Kumaraswamy (GKw) distribution is
#' reparametrized so that one of its five parameters *is* the conditional
#' \eqn{\tau}-quantile \eqn{\mu_\tau}; that parameter is then modelled directly
#' through a link function, while the parameters left free may carry regressions
#' of their own. Estimation is by maximum likelihood with exact gradients from
#' automatic differentiation.
#'
#' Seven nested families are available (`kw`, `ekw`, `kkw`, `bkw`, `gkw`, `mc`,
#' `beta`), so the shape of the conditional distribution is chosen by model
#' selection rather than assumed. Unlike check-function quantile regression, a
#' single fit carries a complete conditional distribution: every other quantile,
#' the density, the mean and the variance are all available from it.
#'
#' @section The distribution and its quantile function:
#' For \eqn{0<y<1}, with \eqn{\alpha,\beta,\gamma,\lambda>0} and
#' \eqn{\delta\ge 0}, the GKw density (Carrasco, Ferrari and Cordeiro, 2010) is
#'
#' \deqn{f(y;\alpha,\beta,\gamma,\delta,\lambda)=
#'   \frac{\lambda\alpha\beta}{B(\gamma,\delta+1)}\,
#'   y^{\alpha-1}\,(1-y^{\alpha})^{\beta-1}\,
#'   \bigl[1-(1-y^{\alpha})^{\beta}\bigr]^{\gamma\lambda-1}\,
#'   \bigl\{1-\bigl[1-(1-y^{\alpha})^{\beta}\bigr]^{\lambda}\bigr\}^{\delta}.}{
#'   f(y) = (lambda*alpha*beta / B(gamma, delta+1)) * y^(alpha-1) *
#'          (1 - y^alpha)^(beta-1) * [1 - (1 - y^alpha)^beta]^(gamma*lambda-1) *
#'          [1 - (1 - (1 - y^alpha)^beta)^lambda]^delta.}
#'
#' Writing \eqn{v=1-(1-y^{\alpha})^{\beta}}{v = 1 - (1 - y^alpha)^beta}, the distribution function is the
#' regularized incomplete beta function \eqn{F(y)=I_{v^{\lambda}}(\gamma,\delta+1)}{F(y) = I[v^lambda](gamma, delta+1)},
#' which inverts in closed form. With
#' \eqn{z_\tau=I^{-1}_\tau(\gamma,\delta+1)}{z_tau = I^-1_tau(gamma, delta+1)}, that is `qbeta(tau, gamma, delta + 1)`,
#'
#' \deqn{Q(\tau)=\Bigl\{1-\bigl[1-z_\tau^{1/\lambda}\bigr]^{1/\beta}\Bigr\}^{1/\alpha}.}{
#'   Q(tau) = ( 1 - [1 - z_tau^(1/lambda)]^(1/beta) )^(1/alpha).}
#'
#' Two special cases keep the incomplete beta out of the computation entirely:
#' \eqn{z_\tau=\tau} when \eqn{\gamma=1,\delta=0} (families `kw`, `ekw`), and
#' \eqn{z_\tau=1-(1-\tau)^{1/(\delta+1)}}{z_tau = 1 - (1-tau)^(1/(delta+1))} when \eqn{\gamma=1} (family `kkw`).
#'
#' @section The quantile reparametrization (the anchor):
#' Fix \eqn{\tau} and set \eqn{\mu=Q(\tau)}. Rearranging the quantile function
#' gives the master identity
#'
#' \deqn{\mu^{\alpha}=1-\bigl(1-z_\tau^{1/\lambda}\bigr)^{1/\beta},}{
#'   mu^alpha = 1 - (1 - z_tau^(1/lambda))^(1/beta),}
#'
#' from which any one of \eqn{\alpha}, \eqn{\beta} or \eqn{\lambda} can be solved
#' for in closed form:
#'
#' \deqn{\alpha=\frac{\log\bigl(1-(1-z_\tau^{1/\lambda})^{1/\beta}\bigr)}{\log\mu},
#'   \qquad
#'   \beta=\frac{\log\bigl(1-z_\tau^{1/\lambda}\bigr)}{\log(1-\mu^{\alpha})},
#'   \qquad
#'   \lambda=\frac{\log z_\tau}{\log\bigl(1-(1-\mu^{\alpha})^{\beta}\bigr)}.}{
#'   alpha  = log(1 - (1 - z_tau^(1/lambda))^(1/beta)) / log(mu),
#'   beta   = log(1 - z_tau^(1/lambda)) / log(1 - mu^alpha),
#'   lambda = log(z_tau) / log(1 - (1 - mu^alpha)^beta).}
#'
#' The parameter eliminated this way is the **anchor**. Each solve is a ratio of
#' two logarithms of numbers in \eqn{(0,1)}, hence strictly positive and finite
#' for every \eqn{\mu,\tau\in(0,1)}{mu, tau in (0,1)} and every admissible value of the remaining
#' parameters. No box constraint is ever needed, which is what makes an
#' unconstrained optimizer safe here.
#'
#' Setting \eqn{\gamma=1}, \eqn{\delta=0}, \eqn{\lambda=1} gives \eqn{z_\tau=\tau}
#' and reduces the \eqn{\beta}-solve to
#' \eqn{\beta=\log(1-\tau)/\log(1-\mu^{\alpha})}{beta = log(1-tau) / log(1-mu^alpha)}, which is exactly the
#' Kumaraswamy quantile (median-dispersion) reparametrization of Mitnik and Baek
#' (2013). The general construction therefore contains the established
#' two-parameter case as a slice.
#'
#' \eqn{\gamma} and \eqn{\delta} admit no closed-form solve: they enter only
#' through \eqn{z_\tau}. Family `beta` is anchored on \eqn{\gamma} by a bracketed
#' root of \eqn{I_\mu(\gamma,\delta+1)=\tau}{I_mu(gamma, delta+1) = tau}, which always exists because the
#' left-hand side sweeps all of \eqn{(0,1)} as \eqn{\gamma} ranges over
#' \eqn{(0,\infty)}. The \eqn{\delta}-solve is inadmissible in every family: it
#' would require \eqn{\tau>\mu^{\gamma}}{tau > mu^gamma}, which nothing guarantees.
#'
#' @section The regression and its estimating equations:
#' For observation \eqn{i}, with \eqn{\tau} fixed,
#'
#' \deqn{g_\mu(\mu_{\tau i}) = \mathbf{x}_i^{\top}\boldsymbol\beta_\mu,
#'   \qquad
#'   g_k(\theta_{ki}) = \mathbf{x}_{ki}^{\top}\boldsymbol\beta_k,}{
#'   g_mu(mu_i) = x_i' b_mu,      g_k(theta_ki) = x_ki' b_k,}
#'
#' where \eqn{g_\mu} maps \eqn{(0,1)\to\mathbb{R}}{(0,1) -> R} (logit by default) and each
#' \eqn{g_k} maps \eqn{(0,\infty)\to\mathbb{R}}{(0,Inf) -> R} (log by default), one equation per
#' non-anchored parameter \eqn{k}. The anchored parameter is *computed* from
#' \eqn{(\mu_{\tau i},\tau,\theta_i)}{(mu_i, tau, theta_i)} by the closed form above; it is never a free
#' coefficient. The log-likelihood is the ordinary GKw log-likelihood evaluated at
#' the reconstructed parameter vector,
#'
#' \deqn{\ell(\boldsymbol\beta)=\sum_{i=1}^{n} w_i \,
#'   \log f\bigl(y_i;\alpha_i,\beta_i,\gamma_i,\delta_i,\lambda_i\bigr),}{
#'   l(b) = sum_i w_i * log f(y_i; alpha_i, beta_i, gamma_i, delta_i, lambda_i),}
#'
#' and the score \eqn{\partial\ell/\partial\boldsymbol\beta}{dl/db} is obtained by
#' automatic differentiation through the anchor solve, so the chain rule through
#' the reparametrization is exact rather than approximated. Maximization uses
#' [stats::nlminb()] or [stats::optim()]; see [gkwq_control()].
#'
#' Standard errors come from the observed information, [stats::optimHess()]
#' applied to the automatic-differentiation gradient. The naive \eqn{J^\top H J}{J' H J}
#' transformation of the unreparametrized Hessian is **wrong** in a regression
#' because it omits the curvature term \eqn{\sum_k g_k \nabla^2\theta_k}{sum_k g_k grad^2 theta_k}, and is
#' never used. See [vcov.gkwqreg()].
#'
#' `tau` is fixed, never estimated. The profile likelihood in \eqn{\tau} is
#' exactly flat -- for any other level there is a parameter value reproducing the
#' same distribution -- so \eqn{\tau} indexes the *question asked of the data*,
#' not the model fitted to it.
#'
#' @section The multi-part formula:
#' The right-hand side is split into parts by `|`, as in the \pkg{Formula}
#' package. The contract is fixed and is reported by [gkwq_parts()]:
#'
#' 1. **Part one is always the conditional \eqn{\tau}-quantile**, `mu`.
#' 2. The remaining parts follow the family's own parameter order
#'    (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta}, \eqn{\lambda})
#'    **with the anchored parameter removed**, because that one is computed
#'    rather than estimated.
#' 3. Omitted trailing parts default to `~ 1`, so `y ~ x` and `y ~ x | 1` are the
#'    same model.
#' 4. Supplying more parts than the family has is an error, not a silent
#'    truncation; the message names the parts the model expects.
#'
#' Because the anchor changes which parameter disappears, it also changes what
#' the later parts mean. Always confirm with `gkwq_parts(family, anchor)`:
#'
#' ```
#' gkwq_parts("kw")            # "mu" "alpha"
#' gkwq_parts("kw", "alpha")   # "mu" "beta"
#' gkwq_parts("gkw")           # "mu" "alpha" "gamma" "delta" "lambda"
#' ```
#'
#' So for `family = "gkw"`, `y ~ x1 | x2 | 1 | x3` puts `x1` in the quantile,
#' `x2` in \eqn{\alpha}, an intercept in \eqn{\gamma}, `x3` in \eqn{\delta}, and
#' an intercept in \eqn{\lambda}. Offsets for a specific part are written as an
#' `offset()` term inside that part.
#'
#' @section The seven families:
#' Every family is a constrained GKw, so they nest: `kw` inside `ekw` inside
#' `kkw` inside `gkw`, and `bkw` and `beta` inside `gkw` as well. Nested pairs
#' can be compared by a likelihood-ratio test ([anova.gkwqreg()]); non-nested
#' pairs by AIC, BIC or [vuong_test()].
#'
#' | family | free parameters | constraints | default anchor | other anchors |
#' |:---|:---|:---|:---|:---|
#' | `kw`   | \eqn{\alpha,\beta}{alpha, beta} | \eqn{\gamma=1,\ \delta=0,\ \lambda=1}{gamma=1, delta=0, lambda=1} | \eqn{\beta}{beta} | \eqn{\alpha}{alpha} |
#' | `ekw`  | \eqn{\alpha,\beta,\lambda}{alpha, beta, lambda} | \eqn{\gamma=1,\ \delta=0}{gamma=1, delta=0} | \eqn{\beta}{beta} | \eqn{\alpha}{alpha}, \eqn{\lambda}{lambda} |
#' | `kkw`  | \eqn{\alpha,\beta,\delta,\lambda}{alpha, beta, delta, lambda} | \eqn{\gamma=1}{gamma=1} | \eqn{\beta}{beta} | \eqn{\alpha}{alpha}, \eqn{\lambda}{lambda} |
#' | `bkw`  | \eqn{\alpha,\beta,\gamma,\delta}{alpha, beta, gamma, delta} | \eqn{\lambda=1}{lambda=1} | \eqn{\beta}{beta} | \eqn{\alpha}{alpha} |
#' | `gkw`  | \eqn{\alpha,\beta,\gamma,\delta,\lambda}{all five} | none | \eqn{\beta}{beta} | \eqn{\alpha}{alpha}, \eqn{\lambda}{lambda} |
#' | `mc`   | \eqn{\gamma,\delta,\lambda}{gamma, delta, lambda} | \eqn{\alpha=\beta=1}{alpha=beta=1} | \eqn{\lambda}{lambda} | none |
#' | `beta` | \eqn{\gamma,\delta}{gamma, delta} | \eqn{\alpha=\beta=\lambda=1}{alpha=beta=lambda=1} | \eqn{\gamma}{gamma} | none |
#'
#' `mc` fixes \eqn{\alpha=\beta=1}{alpha=beta=1}, so neither is available and the
#' \eqn{\lambda}-solve is forced; it collapses to the exact form
#' \eqn{\lambda=\log z_\tau/\log\mu}{lambda = log(z_tau) / log(mu)}. `beta` fixes \eqn{\alpha=\beta=\lambda=1}{alpha=beta=lambda=1},
#' leaving \eqn{\gamma} as the only admissible anchor and the only one solved
#' numerically. [gkwq_anchors()] lists the admissible anchors for a family, the
#' default first.
#'
#' Two identifiability facts are enforced rather than left to the user. When
#' \eqn{\delta=0} with \eqn{\gamma} free, \eqn{I_u(\gamma,1)=u^{\gamma}}{I_u(gamma,1) = u^gamma} and only
#' the product \eqn{\gamma\lambda}{gamma*lambda} is identified, so a \eqn{\lambda}-anchor there
#' is refused with an error. And `family = "gkw"` is weakly identified in *any*
#' parametrization, with information-matrix condition numbers of order
#' \eqn{10^{8}}{1e8} to \eqn{10^{11}}{1e11}; it warns, and `summary()` prints the condition
#' number so the warning can be checked.
#'
#' @section The anchor is a modeling choice:
#' `anchor` names the parameter eliminated in favour of the conditional quantile.
#' It looks like a reparametrization and, when every remaining parameter is
#' regressed on at least the covariates used for `mu`, it is one: the likelihood
#' is identical whichever anchor is chosen. It stops being one as soon as `mu`
#' varies with covariates while a nuisance parameter is held constant. Anchoring
#' on `beta` then asserts "`alpha` constant, `beta_i = f(mu_i, alpha)`";
#' anchoring on `alpha` asserts the reverse, and the two trace different paths
#' through parameter space. In the design study the two differed by 131 in
#' log-likelihood and by 41% in the coefficient of interest.
#'
#' Two anchors on the same data give non-nested models of equal dimension:
#' compare them by AIC, BIC or Vuong's test, **not** by a likelihood-ratio test.
#' [anova.gkwqreg()] refuses the comparison rather than printing a number that
#' would mean nothing. Regressing the nuisance parameters on the same covariates
#' as `mu` makes the fit largely anchor-insensitive, and is recommended whenever
#' theory does not dictate an anchor.
#'
#' The two closed-form solves also fail in **opposite tails**: the
#' \eqn{\beta}-solve overflows as \eqn{\mu\to0}{mu -> 0} and the
#' \eqn{\alpha}-solve as \eqn{\mu\to1}{mu -> 1}. Neither dominates, so where the response mass lies is a
#' legitimate reason to prefer one; see the `eps_mu` entry of [gkwq_control()]
#' for the arithmetic.
#'
#' @section Interpreting the coefficients:
#' Under the default logit link for the quantile part,
#' \eqn{\exp(\beta_{\mu j})}{exp(b_j)} is the multiplicative effect of a one-unit increase
#' in \eqn{x_j} on the **odds of the conditional \eqn{\tau}-quantile**,
#' \eqn{\mu_\tau/(1-\mu_\tau)}{mu/(1 - mu)}. It is *not* an effect on the odds of an event, and
#' *not* an effect on a mean. `summary()` labels the block accordingly.
#'
#' On the scale of the quantile itself the effect is
#'
#' \deqn{\frac{\partial Q_\tau(x)}{\partial x_j}
#'   = \beta_{\mu j}\,\frac{d\mu}{d\eta}
#'   = \beta_{\mu j}\,\mu_\tau(1-\mu_\tau) \quad\text{(logit link)},}{
#'   dQ(tau | x) / dx_j = b_j * (dmu/deta) = b_j * mu_tau * (1 - mu_tau)
#'   under the logit link,}
#'
#' reported with delta-method standard errors by [marginal_effects()]. A
#' covariate that also appears in a nuisance part does not change this
#' expression: `mu` *is* the quantile, so a nuisance part alters the spread of
#' the conditional distribution, not the quantile being modelled. That separation
#' is a genuine convenience of this parametrization.
#'
#' Nuisance coefficients are on the log scale by default, so `exp()` of one is a
#' multiplicative effect on the corresponding shape parameter.
#'
#' @section What `fitted()` returns:
#' **[fitted()] returns the fitted conditional `tau`-quantile, not the
#' conditional mean.** This is the single most important difference from
#' `gkwreg`, which fits the same seven families with mean semantics. A
#' `type = "mean"` escape hatch exists but is never the default. For the same
#' reason the class deliberately does not inherit from `"gkwreg"`: a missing
#' method must error rather than silently dispatch mean semantics on a quantile
#' object.
#'
#' @param formula A multi-part formula, parts separated by `|`. **Part one is
#'   always the conditional quantile.** The remaining parts follow the family's
#'   parameter order with the anchored parameter removed; [gkwq_parts()] reports
#'   the exact contract. Omitted trailing parts default to `~ 1`.
#' @param data An optional data frame, list or environment holding the variables
#'   in `formula`. If missing, they are looked up in the environment of
#'   `formula`.
#' @param tau Quantile level in `(0,1)`. A vector returns a `"gkwqregs"`
#'   container of independent fits, one per level, fitted from the level nearest
#'   the median outward with warm starts. Mandatory in the sense that it is never
#'   estimated: the profile likelihood in `tau` is exactly flat, so `tau` indexes
#'   the question, not the model.
#' @param family One of `"kw"`, `"ekw"`, `"kkw"`, `"bkw"`, `"gkw"`, `"mc"`,
#'   `"beta"`. See the table above for the constraints each imposes.
#' @param anchor The parameter eliminated in favour of the conditional quantile,
#'   or `NULL` for the family default. See the section above; this is a modeling
#'   argument, not an internal detail. [gkwq_anchors()] lists what a family
#'   allows.
#' @param link Link functions, as a single name or a named vector indexed by
#'   part. The quantile part defaults to `"logit"` and must map to `(0,1)`, so it
#'   is one of `"logit"`, `"probit"`, `"cauchy"`, `"cloglog"`. The remaining
#'   parts default to `"log"` and may also be `"identity"`, `"sqrt"`,
#'   `"inverse"` or `"inverse-square"`. A single *unnamed* link applies to the
#'   nuisance parts only, never silently to the quantile; an unnamed vector as
#'   long as `gkwq_parts()` applies part by part.
#' @param link_scale Multiplicative scale factors for the inverse links, as a
#'   single value or a named vector indexed by part. Only the bounded links use
#'   them, rescaling the image from `(0,1)` to `(0, scale)`; the conditional
#'   quantile is a probability and is never rescaled, so its scale is forced
#'   to 1.
#' @param subset An optional expression selecting a subset of rows, evaluated in
#'   `data` first and then in the calling frame.
#' @param weights Optional prior weights \eqn{w_i}, entering the likelihood as
#'   multipliers of the per-observation log-density; the reported check loss is
#'   weighted to match. May name a column of `data`. Negative weights are an
#'   error.
#' @param offset An optional numeric offset added to the linear predictor of the
#'   **conditional quantile**, which is the part an offset is normally meant for.
#'   Offsets for another part are written as an `offset()` term inside that part
#'   of `formula`, and the two add.
#' @param na.action How to treat missing values in the model frame; defaults to
#'   [stats::na.omit()]. Responses at exactly 0 or 1 lie outside the support and
#'   are an error, never silently clamped.
#' @param contrasts An optional list of contrasts for factor covariates, as in
#'   [stats::model.matrix()].
#' @param control An object from [gkwq_control()] holding the optimizer,
#'   starting-value and covariance settings.
#' @param model,x,y Whether to keep the model frame, the per-part design matrices
#'   and the response in the fitted object. `y = TRUE` (the default) is required
#'   by [residuals.gkwqreg()], [plot.gkwqreg()] and the goodness-of-fit lines of
#'   `summary()`; `x = TRUE` avoids rebuilding design matrices in
#'   [marginal_effects()].
#' @param ... Passed to [gkwq_control()]. If any are supplied, the resulting
#'   control object **replaces** `control` entirely, so use one route or the
#'   other rather than both.
#'
#' @return
#' For a single `tau`, an object of class `"gkwqreg"`: a list whose components
#' include
#'
#' \describe{
#'   \item{`coefficients`}{the stacked coefficient vector, named `"part:term"`.}
#'   \item{`coef_list`}{the same coefficients split by part; `coef(fit, part =)`
#'     reads from it.}
#'   \item{`se`, `vcov`, `hessian`, `cond_number`}{standard errors, the inverse
#'     observed information, the observed information itself, and its exact
#'     condition number. `cond_number` above `1e8` is the signature of a weakly
#'     identified fit.}
#'   \item{`fitted.values`}{the fitted conditional \eqn{\tau}-quantiles
#'     \eqn{\mu_{\tau i}} -- **quantiles, not means**.}
#'   \item{`linear.predictors`}{a named list of linear predictors, one per part,
#'     offsets included.}
#'   \item{`parameter_vectors`}{an `n` by 5 data frame of the reconstructed
#'     `alpha`, `beta`, `gamma`, `delta`, `lambda` per observation, the anchored
#'     column included.}
#'   \item{`loglik`, `loglik_i`, `npar`, `nobs`, `aic`, `bic`}{the maximized
#'     log-likelihood (re-evaluated at the reported coefficients, never taken
#'     from the optimizer's own record), its per-observation contributions, the
#'     number of estimated coefficients, the sample size, and the two
#'     information criteria.}
#'   \item{`pinball`}{the in-sample check loss
#'     \eqn{n^{-1}\sum_i (y_i-\mu_{\tau i})\bigl(\tau-\mathbf{1}(y_i<\mu_{\tau i})\bigr)}{mean_i (y_i - mu_i) (tau - 1(y_i < mu_i))}, the
#'     criterion a quantile estimate actually targets.}
#'   \item{`family`, `tau`, `anchor`, `parts`, `spec`}{the family name, the
#'     quantile level, the anchored parameter, the part names in order, and the
#'     full internal family specification.}
#'   \item{`link`, `link_scale`}{the link name and scale used for each part.}
#'   \item{`call`, `formula`, `terms`, `levels`, `contrasts`}{the matched call,
#'     the multi-part formula, per-part terms objects, the factor levels and the
#'     contrasts recorded at fit time, all used by `predict(newdata =)`.}
#'   \item{`weights`, `offsets`}{the prior weights and the per-part offsets.}
#'   \item{`convergence`, `message`, `iterations`}{the optimizer's exit code,
#'     message and iteration count. A non-zero code triggers a warning at fit
#'     time and a note in `summary()`.}
#'   \item{`control`, `start`, `obj`}{the control object, the starting values
#'     actually used, and the TMB object. `obj` is what makes
#'     [estfun.gkwqreg()] and profile-likelihood intervals possible; it holds
#'     external pointers, so it is valid only in the session that built it.}
#'   \item{`model`, `x`, `y`}{present only when the matching argument was `TRUE`:
#'     the model frame, the list of design matrices, and the (clamped) response.}
#' }
#'
#' For a vector `tau`, an object of class `"gkwqregs"` with components `fits` (a
#' named list of `"gkwqreg"` objects, one per level), `taus`, `family`, `anchor`,
#' `parts`, `call` and `nobs`. Each level is a separate likelihood, so the
#' container has no single `logLik()`; take it from an individual fit.
#'
#' @references
#' Carrasco, J. M. F., Ferrari, S. L. P. and Cordeiro, G. M. (2010).
#' A generalized Kumaraswamy distribution.
#'
#' Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile and
#' probability curves without crossing. *Econometrica* **78**, 1093-1125.
#'
#' Cox, D. R. and Reid, N. (1987). Parameter orthogonality and approximate
#' conditional inference. *Journal of the Royal Statistical Society B* **49**,
#' 1-39.
#'
#' Koenker, R. and Bassett, G. (1978). Regression quantiles. *Econometrica*
#' **46**, 33-50.
#'
#' Mitnik, P. A. and Baek, S. (2013). The Kumaraswamy distribution:
#' median-dispersion re-parameterizations. *Statistical Papers* **54**, 177-192.
#' \doi{10.1007/s00362-011-0417-y}
#'
#' @examples
#' ## -------------------------------------------------------------------------
#' ## 1. A Kumaraswamy median regression, end to end
#' ## -------------------------------------------------------------------------
#' ## Simulate a response whose conditional MEDIAN is logit-linear in x1 and x2.
#' ## The generator is the beta anchor itself: a Kumaraswamy with
#' ## beta = log(1 - tau) / log(1 - mu^alpha) has its tau-quantile exactly at mu.
#' set.seed(2024)
#' n  <- 300
#' x1 <- runif(n, -2, 2)
#' x2 <- rbinom(n, 1, 0.5)
#' mu <- plogis(0.4 + 1.1 * x1 - 0.6 * x2)      # the true conditional median
#' y  <- gkwdist::rkw(n, alpha = 2, beta = log(1 - 0.5) / log(1 - mu^2))
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' fit <- gkwqreg(y ~ x1 + x2, data = dat, tau = 0.5, family = "kw")
#' summary(fit)
#' ## The "mu" block should recover (0.4, 1.1, -0.6) up to sampling error, and
#' ## the "alpha" block log(2) = 0.693 -- alpha is modelled on the log scale.
#' ## The footer reports "Anchored parameter: beta", the reminder that beta was
#' ## computed from the fitted median rather than estimated, plus the pinball
#' ## loss and the empirical coverage, which should sit near tau.
#'
#' ## -------------------------------------------------------------------------
#' ## 2. Reading the coefficients: quantile odds, not event odds
#' ## -------------------------------------------------------------------------
#' exp(coef(fit, part = "mu"))
#' ## exp(b_x1) should be close to exp(1.1) = 3.0: a one-unit rise in x1
#' ## multiplies the ODDS OF THE CONDITIONAL MEDIAN, mu/(1 - mu), by about 3.
#' ## It is not the odds of an event, and not an effect on a mean.
#'
#' ## On the quantile scale itself the effect is b_j * mu * (1 - mu), i.e. in
#' ## probability units, averaged over the sample:
#' marginal_effects(fit)
#'
#' ## -------------------------------------------------------------------------
#' ## 3. Fitted values are quantiles
#' ## -------------------------------------------------------------------------
#' head(fitted(fit))                     # conditional medians
#' head(fitted(fit, type = "mean"))      # conditional means, by quadrature
#' mean(dat$y <= fitted(fit))            # empirical coverage; target is tau
#'
#' ## -------------------------------------------------------------------------
#' ## 4. Prediction: one fit, a whole conditional distribution
#' ## -------------------------------------------------------------------------
#' nd <- data.frame(x1 = c(-1, 0, 1), x2 = 0)
#' predict(fit, newdata = nd)                        # conditional medians
#' ## Other quantiles read off the SAME fitted distribution -- these cannot
#' ## cross, because they are quantiles of one proper distribution function:
#' predict(fit, newdata = nd, tau = c(0.1, 0.5, 0.9))
#'
#' \donttest{
#' ## -------------------------------------------------------------------------
#' ## 5. The formula grammar
#' ## -------------------------------------------------------------------------
#' gkwq_parts("kw")             # "mu" "alpha" -- part 1 is always the quantile
#' gkwq_parts("kw", "alpha")    # "mu" "beta"  -- the anchor is what disappears
#' gkwq_parts("gkw")            # five parts, in family order
#'
#' ## Give the shape parameter a regression of its own. Omitted trailing parts
#' ## default to ~ 1, so `y ~ x1 + x2` is `y ~ x1 + x2 | 1`.
#' fit_sat <- gkwqreg(y ~ x1 + x2 | x1 + x2, data = dat, tau = 0.5,
#'                    family = "kw")
#' anova(fit, fit_sat)          # nested, same anchor: an ordinary LR test
#'
#' ## -------------------------------------------------------------------------
#' ## 6. The anchor is a modeling choice
#' ## -------------------------------------------------------------------------
#' fb <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw", anchor = "beta")
#' fa <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw", anchor = "alpha")
#' c(beta = as.numeric(logLik(fb)), alpha = as.numeric(logLik(fa)))
#' ## Same number of parameters, different likelihood: these are NON-NESTED
#' ## models, so compare them by AIC/BIC or a Vuong test, never by an LR test.
#' vuong_test(fb, fa)
#'
#' ## Saturating the nuisance part makes the choice nearly immaterial, which is
#' ## the recommended default when theory does not dictate an anchor:
#' fbs <- gkwqreg(y ~ x1 | x1, data = dat, tau = 0.5, family = "kw",
#'                anchor = "beta")
#' fas <- gkwqreg(y ~ x1 | x1, data = dat, tau = 0.5, family = "kw",
#'                anchor = "alpha")
#' c(beta = coef(fbs)[["mu:x1"]], alpha = coef(fas)[["mu:x1"]])
#'
#' ## -------------------------------------------------------------------------
#' ## 7. A grid of quantile levels
#' ## -------------------------------------------------------------------------
#' fits <- gkwqreg(y ~ x1, data = dat, tau = c(0.1, 0.25, 0.5, 0.75, 0.9),
#'                 family = "kw")
#' round(coef(fits), 3)   # one column per level: the estimated quantile process
#' ## The slope on x1 growing with tau is heteroskedasticity that a mean
#' ## regression would average away.
#' ## Independently fitted levels are not ordered by construction; here none
#' ## cross, but nothing guaranteed that, which is exactly what this reports:
#' check_crossing(fits)
#'
#' ## -------------------------------------------------------------------------
#' ## 8. Choosing a family
#' ## -------------------------------------------------------------------------
#' fit_ekw <- gkwqreg(y ~ x1 + x2, data = dat, tau = 0.5, family = "ekw")
#' anova(fit, fit_ekw)    # kw is nested in ekw: one degree of freedom
#' ## Check loss is the criterion a quantile fit targets, so read it beside
#' ## AIC -- and read the `converged` column before either of them, since a
#' ## richer family that failed to converge can still post the best logLik:
#' compare_families(fit, families = c("kw", "ekw", "kkw", "beta"))
#' }
#'
#' @seealso
#' [gkwq_parts()] and [gkwq_anchors()] for the formula and anchor contracts;
#' [gkwq_control()] for optimizer, starting-value and covariance settings;
#' [predict.gkwqreg()], [fitted.gkwqreg()] and [residuals.gkwqreg()] for what
#' comes out of a fit; [marginal_effects()] for effects on the quantile scale;
#' [vcov.gkwqreg()] and [confint.gkwqreg()] for inference;
#' [anova.gkwqreg()], [vuong_test()] and [compare_families()] for model choice;
#' [quantile_process()], [check_crossing()], [rearrange()] and [pinball()] for
#' quantile-specific tooling; [plot.gkwqreg()] for diagnostics;
#' [gkwq_quantile()] for the closed-form quantile function on its own.
#'
#' @family model fitting
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
  ## The part is constant across observations only when every column of its
  ## design matrix is constant. `ncol(X) > 1` is a different test: a part
  ## written without an intercept -- `~ z - 1`, `~ 0 + z` -- has one column
  ## that varies with i, and hoisting log z there evaluates a likelihood that
  ## is not the model's.
  varying <- vapply(c("gamma", "delta"), function(nm) {
    if (!nm %in% parts) return(FALSE)
    X <- md$X[[nm]]
    any(vapply(seq_len(ncol(X)), function(j) any(X[, j] != X[1L, j]), logical(1))) ||
      any(md$offsets[[nm]] != md$offsets[[nm]][1L])
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
